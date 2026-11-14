"""Test wheel selection and runtime pinning without native builds or downloads."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


RUNNER = Path(__file__).with_name('test_in_container.sh')

# Run the real shell wrapper, replacing only pip/native modules with fixtures.
# The wrapper's inline Python checks are executed, including version validation.
PYTHON_FIXTURE = r'''
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import types

state_path = Path(os.environ['RUNNER_STATE'])
state = json.loads(state_path.read_text()) if state_path.exists() else {'version': '2025.9.1'}
args = sys.argv[1:]
event = {'args': args, 'cwd': os.getcwd()}
if args[:3] == ['-m', 'pip', 'install']:
    if '--constraint' in args:
        constraint = Path(args[args.index('--constraint') + 1]).read_text().strip()
        event['constraint'] = constraint
        state['version'] = constraint.removeprefix('rdkit==')
    if args[-1].endswith('[pdbqt]') and os.environ.get('RUNNER_SIMULATE_DRIFT'):
        state['version'] = '2099.9.1'
    state_path.write_text(json.dumps(state))
elif args[:2] == ['-m', 'unittest']:
    copied_tests = Path(args[args.index('-s') + 1])
    assert (copied_tests / 'test_fixture.py').is_file()
    assert Path.cwd() != Path(os.environ['RUNNER_SOURCE'])

with Path(os.environ['RUNNER_LOG']).open('a') as log:
    log.write(json.dumps(event) + '\n')

if not args or args == ['-']:
    for name in ('rdkit', 'msani', 'msani_confgen_cpp', 'msani_stereoisomers', 'amsolcpp', 'amsolcpp._amsolcpp'):
        module = types.ModuleType(name)
        module.__file__ = '/fixture/site-packages/' + name.replace('.', '/') + '.py'
        sys.modules[name] = module
    rd_base = types.ModuleType('rdkit.rdBase')
    parts = state['version'].split('.')
    rd_base.rdkitVersion = '.'.join([parts[0], parts[1].zfill(2), parts[2]])
    sys.modules['rdkit'].rdBase = rd_base
    # The wrapper also checks the RDKit the extension was compiled against.
    sys.modules['msani_confgen_cpp'].RDKIT_CPP_VERSION = '2026.09.2'
    sys.modules['msani_confgen_cpp'].EMBEDDING_PRESETS = ('sretkdgv3', 'etkdgv4', 'etkdgv4-ic')
    sys.modules['amsolcpp']._amsolcpp =sys.modules['amsolcpp._amsolcpp']
    importlib.metadata.version = lambda name: state['version']
    exec(compile(sys.stdin.read(), '<wheel-runtime-check>', 'exec'))
'''


class WheelTestRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root / 'bin').mkdir()
        python = self.root / 'bin/python'
        python.write_text(f'#!{sys.executable}\n' + PYTHON_FIXTURE)
        python.chmod(0o755)
        (self.root / 'tests').mkdir()
        (self.root / 'tests/test_fixture.py').write_text('# Copied test input\n')
        (self.root / 'wheels').mkdir()
        self.log = self.root / 'calls.jsonl'
        self.environment = {
            **os.environ,
            'PATH': str(self.root / 'bin') + os.pathsep + os.environ['PATH'],
            'MSANI_TAG': 'cp312',
            'MSANI_WHEELS': 'wheels',
            'MSANI_TESTS': 'tests',
            'RUNNER_STATE': str(self.root / 'state.json'),
            'RUNNER_LOG': str(self.log),
            'RUNNER_SOURCE': str(self.root),
        }
        self.environment.pop('MSANI_TEST_RDKIT', None)

    def wheel(self, name='molsanitizer-0.7.5-cp312-cp312-manylinux_2_28_x86_64.whl'):
        (self.root / 'wheels' / name).touch()

    def run_wrapper(self, version=None, drift=False):
        env = self.environment.copy()
        if version:
            env['MSANI_TEST_RDKIT'] = version
        if drift:
            env['RUNNER_SIMULATE_DRIFT'] = '1'
        return subprocess.run(['bash', str(RUNNER)], cwd=self.root, env=env,
                              text=True, capture_output=True)

    def calls(self):
        return [json.loads(line) for line in self.log.read_text().splitlines()]

    def test_pins_all_selected_versions_through_optional_dependencies(self):
        self.wheel()
        for version in ('2025.3.1', '2025.9.1', '2026.9.1', '2027.3.1', '2027.3.1.post1'):
            with self.subTest(version=version):
                self.log.unlink(missing_ok=True)
                result = self.run_wrapper(version)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                installs = [call for call in self.calls() if call['args'][:3] == ['-m', 'pip', 'install']]
                self.assertEqual(len(installs), 2)
                self.assertTrue(all(call['constraint'] == f'rdkit=={version}' for call in installs))
                self.assertTrue(all('--only-binary=rdkit' in call['args'] for call in installs))
                self.assertTrue(installs[1]['args'][-1].endswith('.whl[pdbqt]'))
                suite = self.calls()[-1]
                self.assertEqual(suite['args'][:3], ['-m', 'unittest', 'discover'])
                self.assertFalse(Path(suite['args'][suite['args'].index('-s') + 1]).exists())

    def test_unpinned_release_test_remains_supported(self):
        self.wheel()
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        installs = [call for call in self.calls() if call['args'][:3] == ['-m', 'pip', 'install']]
        self.assertTrue(all('--constraint' not in call['args'] for call in installs))

    def test_runtime_drift_fails_before_running_the_suite(self):
        self.wheel()
        result = self.run_wrapper('2025.9.1', drift=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Expected Python RDKit 2025.9.1, found 2099.9.1', result.stderr)
        self.assertFalse(any(call['args'][:2] == ['-m', 'unittest'] for call in self.calls()))

    def test_missing_wheel_fails_without_installing_from_source(self):
        result = self.run_wrapper('2025.9.1')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Expected exactly one cp312 wheel', result.stderr)
        self.assertFalse(self.log.exists())

    def test_multiple_matching_wheels_are_rejected(self):
        self.wheel()
        self.wheel('molsanitizer-0.7.6-cp312-cp312-manylinux_2_28_x86_64.whl')
        result = self.run_wrapper('2025.9.1')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Expected exactly one cp312 wheel', result.stderr)
        self.assertFalse(self.log.exists())


if __name__ == '__main__':
    unittest.main()
