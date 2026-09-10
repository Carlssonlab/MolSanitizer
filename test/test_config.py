"""Configuration discovery and CLI defaults without modifying user settings."""
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import yaml

from msani import config
from msani.io import parsers


class TestConfig(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name)
        directory_patch = patch.object(config, 'config_directory', return_value=self.path / 'settings')
        directory_patch.start()
        self.addCleanup(directory_patch.stop)
        environment_patch = patch.dict(os.environ)
        environment_patch.start()
        self.addCleanup(environment_patch.stop)
        os.environ.pop('MSANI_CONFIG', None)
        self.stdout = StringIO()
        self.stderr = StringIO()
        for redirect in (redirect_stdout(self.stdout), redirect_stderr(self.stderr)):
            redirect.__enter__()
            self.addCleanup(redirect.__exit__, None, None, None)

    def test_missing_uses_defaults_without_writes(self):
        values, sources = config.load_defaults()
        self.assertEqual(values['NUMCONFS'], 2000)
        self.assertEqual(sources['NUMCONFS'], 'bundled defaults')
        self.assertEqual(self.stdout.getvalue(), '')
        self.assertEqual(self.stderr.getvalue(), '')
        self.assertFalse((self.path / 'settings').exists())

    def test_init_and_remember(self):
        path = self.path / 'custom directory' / 'defaults.yaml'
        config.main(['init', '--path', str(path)])
        pointer = json.loads((self.path / 'settings/settings.json').read_text())
        self.assertEqual(pointer['config_path'], str(path))
        path.write_text('NUMCONFS: 123\n')
        values, sources = config.load_defaults()
        self.assertEqual(values['NUMCONFS'], 123)
        self.assertEqual(values['TIMEOUT'], 2)
        self.assertEqual(sources['NUMCONFS'], str(path))
        with self.assertRaises(SystemExit):
            config.main(['init', '--path', str(path)])
        self.assertEqual(config.load_defaults()[0]['NUMCONFS'], 123)

    def test_env_override_does_not_change_pointer(self):
        config.main(['init'])
        pointer = (self.path / 'settings/settings.json').read_text()
        path = self.path / 'env.yaml'
        path.write_text('NUMCONFS: 77\n')
        os.environ['MSANI_CONFIG'] = str(path)
        self.assertEqual(config.load_defaults()[0]['NUMCONFS'], 77)
        self.assertEqual((self.path / 'settings/settings.json').read_text(), pointer)
        path.unlink()
        with self.assertRaisesRegex(config.ConfigError, 'env.yaml'):
            config.load_defaults()

    def test_validation(self):
        for text in ['NUMCONFS: wrong', 'UNKNOWN: 3', '- list',
                     'WHOLE_NODE: 1', 'TIMEOUT: -1', 'PH: .nan']:
            with self.subTest(text=text):
                path = self.path / 'bad.yaml'
                path.write_text(text)
                with self.assertRaisesRegex(config.ConfigError, 'bad.yaml'):
                    config.read_config(path)

    def test_precedence_and_reload(self):
        config.main(['init'])
        path = config.selected_path()[0]
        path.write_text('NUMCONFS: 123\nWHOLE_NODE: true\nEMBED_METHOD: obabel\n')
        args, _ = parsers.parseArguments([], batch_mode=True)
        self.assertEqual((args.numconfs, args.whole_node, args.method), (123, True, 'obabel'))
        job = self.path / 'job.yaml'
        job.write_text('numconfs: 456\n')
        self.assertEqual(parsers.parseArguments(['--config', str(job)]).numconfs, 456)
        self.assertEqual(
            parsers.parseArguments(['--config', str(job), '--numconfs', '789']).numconfs,
            789,
        )
        path.write_text('NUMCONFS: 321\n')
        self.assertEqual(parsers.parseArguments([]).numconfs, 321)

    def test_use_reset_and_missing_remembered(self):
        path = self.path / 'other.yaml'
        path.write_text('NUMCONFS: 55\n')
        config.main(['use', str(path)])
        path.unlink()
        with self.assertRaises(config.ConfigError):
            config.load_defaults()
        config.main(['reset-path'])
        self.assertEqual(config.load_defaults()[0]['NUMCONFS'], 2000)

    def test_snapshot(self):
        from msani.batchmode import write_single_job_script

        config.main(['init'])
        config.selected_path()[0].write_text('NUMCONFS: 123\n')
        previous_directory = Path.cwd()
        self.addCleanup(os.chdir, previous_directory)
        os.chdir(self.path)
        with patch('msani.batchmode.subprocess.run', return_value=SimpleNamespace(stdout='/bin/msani')):
            write_single_job_script('#!/bin/bash\n', 'MSANI_PATH -i input.smi\n')
        self.assertEqual(yaml.safe_load(Path('msani_defaults.yaml').read_text())['NUMCONFS'], 123)
        self.assertIn('SLURM_SUBMIT_DIR', Path('submit_msani.sh').read_text())

    def test_repeated_parsing_is_quiet(self):
        for _ in range(3):
            self.assertEqual(parsers.parseArguments([]).numconfs, 2000)
        self.assertEqual(self.stdout.getvalue(), '')
        self.assertEqual(self.stderr.getvalue(), '')

    def test_help_explains_personal_defaults(self):
        for batch_mode in [False, True]:
            with self.subTest(batch_mode=batch_mode):
                self.stdout.seek(0)
                self.stdout.truncate(0)
                with self.assertRaises(SystemExit) as error:
                    parsers.parseArguments(['--help'], batch_mode=batch_mode)
                self.assertEqual(error.exception.code, 0)
                self.assertIn('msani config init', self.stdout.getvalue())
                self.assertIn('msani config show', self.stdout.getvalue())
                self.assertEqual(self.stderr.getvalue(), '')


if __name__ == '__main__':
    unittest.main()
