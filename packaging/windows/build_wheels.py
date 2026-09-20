"""Local Windows x64 wheel matrix; reuse one SDK across all Python versions."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from windows_toolchain import select_msvc, wheel_compiler_settings


def run(args, *, env, cwd, log):
    args = [str(a) for a in args]
    log.write('\n> ' + subprocess.list2cmdline(args) + '\n')
    log.flush()
    with subprocess.Popen(args, env=env, cwd=cwd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True,
                          encoding='utf-8', errors='replace') as process:
        for line in process.stdout:
            print(line, end='')
            log.write(line)
        code = process.wait()
    log.flush()
    if code:
        raise subprocess.CalledProcessError(code, args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sdk', required=True)
    parser.add_argument('--rdkit-version', required=True,
                        help='Python RDKit version matching the Windows C++ SDK')
    parser.add_argument('--extra-prefixes', default='')
    parser.add_argument('--versions', nargs='+', default=['3.10', '3.11', '3.12', '3.13', '3.14'])
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--python-provider', choices=['conda', 'launcher', 'venv'],
                        default='conda',
                        help="'venv' builds from --base-python; use it in CI, where a\n"
                             'single interpreter is provisioned per job.')
    parser.add_argument('--base-python', type=Path,
                        help='Interpreter for --python-provider venv (one version per run).')
    args = parser.parse_args()
    if sys.platform != 'win32':
        parser.error('Run on native Windows, not WSL.')
    if args.jobs < 1 or any(v not in ('3.10', '3.11', '3.12', '3.13', '3.14') for v in args.versions):
        parser.error('Use positive jobs and Python versions 3.10 through 3.14.')
    if not shutil.which('cl'):
        parser.error('Start from the VS x64 developer shell (cl.exe must be available).')
    launcher = shutil.which('py') if args.python_provider == 'launcher' else None
    conda = os.environ.get('CONDA_EXE') or shutil.which('conda.exe')
    if args.python_provider == 'launcher' and not launcher:
        parser.error('Install the Windows Python launcher and standalone CPython versions.')
    if args.python_provider == 'conda' and (not conda or not Path(conda).is_file()):
        parser.error('Activate conda in this Command Prompt; CONDA_EXE must identify conda.exe.')
    if args.python_provider == 'venv':
        if not args.base_python or not args.base_python.is_file():
            parser.error('--python-provider venv requires --base-python PATH_TO_PYTHON_EXE.')
        # One interpreter cannot seed environments for several versions.
        if len(set(args.versions)) != 1:
            parser.error('--python-provider venv builds exactly one --versions entry.')
    repo = Path(__file__).resolve().parents[2]
    sdk = Path(args.sdk).resolve()
    # Accept either a conda environment or an independent SDK prefix. Using
    # Library directly avoids discovering the conda environment's Python.
    if (sdk / 'Library').is_dir():
        sdk = sdk / 'Library'
    if not (sdk / 'include').is_dir() or not any((sdk / d).is_dir() for d in ('bin', 'lib')):
        parser.error('SDK must contain include/ and bin/ or lib/ (or a conda Library/ directory).')
    sdk_extras = []
    manifest = sdk / 'msani-sdk.json'
    if manifest.exists():
        metadata = json.loads(manifest.read_text(encoding='utf-8'))
        if metadata.get('rdkit_version') != args.rdkit_version:
            parser.error('The requested Python RDKit version does not match msani-sdk.json.')
        sdk_extras = [Path(p) for p in metadata.get('extra_prefixes', [])]
        if not all(p.is_dir() for p in sdk_extras):
            parser.error('An SDK dependency directory recorded in msani-sdk.json is missing.')
    prefixes = list(dict.fromkeys([sdk, *sdk_extras,
                    *[Path(p).resolve() for p in args.extra_prefixes.split(';') if p]]))
    dll_dirs = [p / d for p in prefixes for d in ('bin', 'lib') if (p / d).is_dir()]
    (repo / 'dist').mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='windows-wheels-', dir=repo / 'dist'))
    wheels = output / 'wheels'
    wheels.mkdir()
    logs = output / 'logs'
    logs.mkdir()
    print(f'Output: {output}')
    base_env = os.environ.copy()
    for key in list(base_env):
        if key.upper().startswith(('CONDA', 'PYTHON', 'CMAKE', 'PIP_')) or key.upper() in ('RDBASE', 'VIRTUAL_ENV'):
            base_env.pop(key, None)
    base_env['PYTHONUTF8'] = '1'
    base_env['PYTHONNOUSERSITE'] = '1'
    compiler, compiler_options = select_msvc(base_env)
    print(f'Pinned C/C++ compiler: {compiler}')
    # Keep the SDK outside interpreter environments. Conda caches package
    # downloads automatically, while each interpreter environment stays fresh.
    def create_environment(prefix, version, env, cwd, log):
        if args.python_provider == 'venv':
            # --seed is unavailable here; ensurepip provides pip in the venv.
            run([args.base_python, '-m', 'venv', prefix], env=env, cwd=cwd, log=log)
            python = prefix / 'Scripts/python.exe'
            run([python, '-m', 'pip', 'install', '--upgrade', 'pip'], env=env, cwd=cwd, log=log)
            return python
        if args.python_provider == 'conda':
            conda_env = os.environ.copy()
            conda_env['CONDA_SUBDIR'] = 'win-64'
            conda_env.pop('PYTHONPATH', None)
            conda_env.pop('PYTHONHOME', None)
            run([conda, 'create', '--yes', '--prefix', prefix,
                 f'python={version}', 'pip'], env=conda_env, cwd=cwd, log=log)
            return prefix / 'python.exe'
        run([launcher, f'-{version}', '-m', 'venv', prefix], env=env, cwd=cwd, log=log)
        return prefix / 'Scripts/python.exe'

    def interpreter_path(prefix):
        paths = [prefix, prefix / 'Scripts']
        if args.python_provider == 'venv':
            return os.pathsep.join(map(str, paths))
        if args.python_provider == 'conda':
            paths.append(prefix / 'Library/bin')
        return os.pathsep.join(map(str, paths))
    results = []
    (output / 'build-settings.json').write_text(json.dumps({
        'sdk': str(sdk), 'rdkit_version': args.rdkit_version,
        'versions': args.versions, 'extra_prefixes': args.extra_prefixes,
        'python_provider': args.python_provider,
        'base_python': str(args.base_python) if args.base_python else None,
        'compiler_options': compiler_options,
    }, indent=2), encoding='utf-8')
    for version in dict.fromkeys(args.versions):
        tag = 'cp' + version.replace('.', '')
        phase = 'BUILD'
        try:
            with (logs / f'{tag}.log').open('w', encoding='utf-8') as log:
                with tempfile.TemporaryDirectory(prefix=f'msani-{tag}-') as tmp:
                    work = Path(tmp)
                    check = (
                        "import sys,struct,platform,pathlib; "
                        f"assert sys.version_info[:2] == {tuple(map(int, version.split('.')))!r}; "
                        "assert platform.python_implementation() == 'CPython'; "
                        "assert struct.calcsize('P') == 8 and platform.machine().lower() in ('amd64','x86_64'); "
                        "assert not __import__('sysconfig').get_config_var('Py_GIL_DISABLED'); "
                        "print(sys.executable, sys.version)"
                    )
                    build = work / 'build-env'
                    python = create_environment(build, version, base_env, work, log)
                    env = base_env.copy()
                    env.update(RDBASE=sdk.as_posix(), CMAKE_GENERATOR='Ninja',
                               CMAKE_BUILD_PARALLEL_LEVEL=str(args.jobs))
                    env['PATH'] = os.pathsep.join([interpreter_path(build), *map(str, dll_dirs), env['PATH']])
                    env['PATH'] = str(Path(compiler).parent) + os.pathsep + env['PATH']
                    run([python, '-c', check], env=env, cwd=work, log=log)
                    run([python, '-m', 'pip', 'install', 'scikit-build-core', 'pybind11', 'ninja',
                         'cmake', 'delvewheel', 'twine', f'rdkit=={args.rdkit_version}'], env=env, cwd=work, log=log)
                    run([python, '-m', 'pip', 'freeze'], env=env, cwd=work, log=log)
                    raw = work / 'raw'
                    run([python, '-m', 'pip', 'wheel', repo, '--no-build-isolation', '--no-deps',
                         f'-Cbuild-dir={work / "native"}', '-Ccmake.define.MSANI_RDKIT_LINKAGE=shared',
                         *wheel_compiler_settings(compiler_options),
                         '-Ccmake.define.CMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL',
                         '-Ccmake.define.CMAKE_PREFIX_PATH=' + ';'.join(p.as_posix() for p in prefixes),
                         f'-Ccmake.define.Python3_EXECUTABLE={python.as_posix()}',
                         f'-Ccmake.define.Python_EXECUTABLE={python.as_posix()}',
                         '--wheel-dir', raw], env=env, cwd=work, log=log)
                    candidates = list(raw.glob(f'*-{tag}-{tag}-win_amd64.whl'))
                    if len(candidates) != 1:
                        raise RuntimeError(f'Expected one {tag} wheel, found {candidates}')
                    wheel = candidates[0]
                    for command in ('show', 'repair'):
                        cmd = [python, '-m', 'delvewheel', command, wheel,
                               '--add-path', os.pathsep.join(map(str, dll_dirs))]
                        if command == 'repair':
                            cmd += ['-w', wheels]
                        run(cmd, env=env, cwd=work, log=log)
                    repaired = wheels / wheel.name
                    run([python, '-m', 'twine', 'check', '--strict', repaired], env=env, cwd=work, log=log)
                    phase = 'TEST'
                    test_env = base_env.copy()
                    # Do not expose build tools, conda, or SDK DLLs during tests.
                    system = Path(os.environ['SystemRoot'])
                    test_env['PATH'] = os.pathsep.join(map(str, [system / 'System32', system]))
                    clean = work / 'test-env'
                    test_python = create_environment(clean, version, test_env, work, log)
                    test_env['PATH'] = interpreter_path(clean) + os.pathsep + test_env['PATH']
                    run([test_python, '-c', check], env=test_env, cwd=work, log=log)
                    run([test_python, '-m', 'pip', 'install', repaired], env=test_env, cwd=work, log=log)
                    smoke = (
                        'from pathlib import Path; import msani, msani_confgen_cpp, msani_stereoisomers, amsolcpp; '
                        'from amsolcpp import _amsolcpp; '
                        'modules=(msani,msani_confgen_cpp,msani_stereoisomers,amsolcpp,_amsolcpp); '
                        'assert all("site-packages" in Path(m.__file__).parts for m in modules); '
                        'from rdkit import rdBase; print("Runtime RDKit:",rdBase.rdkitVersion); '
                        'print("Native imports passed")'
                    )
                    run([test_python, '-c', smoke], env=test_env, cwd=work, log=log)
                    # Install the wheel's declared optional dependencies, not an
                    # unrelated MolSanitizer release from the package index.
                    run([test_python, '-m', 'pip', 'install', str(repaired) + '[pdbqt]'], env=test_env, cwd=work, log=log)
                    run([test_python, '-m', 'pip', 'check'], env=test_env, cwd=work, log=log)
                    run([test_python, '-m', 'pip', 'freeze'], env=test_env, cwd=work, log=log)
                    tests = work / 'tests'
                    shutil.copytree(repo / 'test', tests, ignore=shutil.ignore_patterns('__pycache__'))
                    run([test_python, '-m', 'unittest', 'discover', '-s', tests, '-p', 'test*.py', '-v'],
                        env=test_env, cwd=work, log=log)
            results.append(f'{version} PASSED')
        except Exception as error:
            results.append(f'{version} {phase}_FAILED: {error}')
            print(results[-1], file=sys.stderr)
        (output / 'status.txt').write_text('\n'.join(results) + '\n', encoding='utf-8')
    hashes = [f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}' for p in sorted(wheels.glob('*.whl'))]
    (output / 'SHA256SUMS').write_text('\n'.join(hashes) + '\n', encoding='ascii')
    print('\n'.join(results))
    print(f'Results: {output}')
    return 0 if all(r.endswith(' PASSED') for r in results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
