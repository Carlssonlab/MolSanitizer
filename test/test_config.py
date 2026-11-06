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
        os.environ.pop('MSANI_BATCH_WORKER', None)
        self.stdout = StringIO()
        self.stderr = StringIO()
        for redirect in (redirect_stdout(self.stdout), redirect_stderr(self.stderr)):
            redirect.__enter__()
            self.addCleanup(redirect.__exit__, None, None, None)

    def test_missing_uses_defaults_without_writes(self):
        values, sources = config.load_defaults()
        self.assertEqual(values['NUMCONFS'], 2000)
        self.assertEqual(values['PH_RANGE'], 1)
        self.assertEqual(sources['NUMCONFS'], 'bundled defaults')
        self.assertEqual(self.stdout.getvalue(), '')
        self.assertEqual(self.stderr.getvalue(), '')
        self.assertFalse((self.path / 'settings').exists())

    def test_init_and_remember(self):
        path = self.path / 'custom directory' / 'defaults.yaml'
        config.main(['init', '--path', str(path)])
        pointer = json.loads((self.path / 'settings/settings.json').read_text())
        # Configuration stores canonical paths; macOS /var is a symlink to
        # /private/var, including the default temporary directory.
        self.assertEqual(pointer['config_path'], str(path.resolve()))
        path.write_text('NUMCONFS: 123\n')
        values, sources = config.load_defaults()
        self.assertEqual(values['NUMCONFS'], 123)
        self.assertEqual(values['TIMEOUT'], 2)
        self.assertEqual(sources['NUMCONFS'], str(path.resolve()))
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
        self.assertIn('export MSANI_BATCH_WORKER=1', Path('submit_msani.sh').read_text())

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
                self.assertIn('bundled default changed from 0 to 1 in 0.8.0', self.stdout.getvalue())
                self.assertEqual(self.stderr.getvalue(), '')

    def test_personal_ph_settings_do_not_activate_protonation(self):
        path = self.path / 'personal.yaml'
        path.write_text('PH: 5\nPH_RANGE: 0\n')
        os.environ['MSANI_CONFIG'] = str(path)
        for batch_mode in [False, True]:
            with self.subTest(batch_mode=batch_mode):
                parsed = parsers.parseArguments([], batch_mode=batch_mode)
                args = parsed[0] if batch_mode else parsed
                parsers.Sanitycheck(args)
                self.assertEqual((args.pH, args.pH_range), (5, 0))
                self.assertFalse(args.protonation)
                self.assertFalse(args.neutralize)

    def test_missing_personal_ph_range_inherits_new_default(self):
        path = self.path / 'personal.yaml'
        path.write_text('PH: 5\n')
        os.environ['MSANI_CONFIG'] = str(path)
        args = parsers.parseArguments(['--protonation'])
        self.assertEqual(args.pH_range, 1)
        self.assertEqual(args._pH_sources['pH_range'], 'bundled defaults')

    def test_explicit_ph_options_activate_protonation_even_at_default(self):
        for flags in [['-p', '7'], ['--pH=7'], ['-r', '0'], ['-r', '1'],
                      ['--pH_range=3'], ['-r2']]:
            with self.subTest(flags=flags):
                args = parsers.Sanitycheck(parsers.parseArguments(flags))
                self.assertTrue(args.protonation)
                self.assertTrue(args.neutralize)

    def test_job_yaml_and_cli_ph_sources(self):
        personal = self.path / 'personal.yaml'
        personal.write_text('PH_RANGE: 0\n')
        os.environ['MSANI_CONFIG'] = str(personal)
        job = self.path / 'job.yaml'
        for settings in ['pH_range: 1\n', 'pH: 7\n']:
            with self.subTest(settings=settings):
                job.write_text(settings)
                args = parsers.Sanitycheck(parsers.parseArguments(['--config', str(job)]))
                self.assertTrue(args.protonation)
        job.write_text('pH_range: 2\n')
        args = parsers.parseArguments(['--config', str(job)])
        self.assertEqual(args.pH_range, 2)
        self.assertEqual(args._pH_sources['pH_range'], f'job configuration {job.resolve()}')
        args = parsers.parseArguments(['--config', str(job), '-r', '0'])
        self.assertEqual(args.pH_range, 0)
        self.assertEqual(args._pH_sources['pH_range'], 'command line (-r)')

    def test_notice_describes_bundled_and_preserved_defaults(self):
        from msani.io.loggers import protonation_notice

        args = parsers.parseArguments(['--protonation'])
        protonation_notice(args)
        self.assertIn('pH_range=1 from bundled defaults', self.stderr.getvalue())
        self.assertIn('This can change protonation results', self.stderr.getvalue())
        self.assertIn('--pH_range 0', self.stderr.getvalue())
        self.stderr.seek(0)
        self.stderr.truncate(0)
        path = self.path / 'personal.yaml'
        path.write_text('PH_RANGE: 0\n')
        os.environ['MSANI_CONFIG'] = str(path)
        protonation_notice(parsers.parseArguments(['--protonation']))
        self.assertIn(f'pH_range=0 from {path.resolve()}', self.stderr.getvalue())
        self.assertIn('Your configured value is preserved', self.stderr.getvalue())
        self.assertEqual(path.read_text(), 'PH_RANGE: 0\n')

    def test_cli_notice_is_printed_and_logged_before_processing(self):
        from msani import cli

        def check_notice(args):
            self.assertIn('pH_range=1 from bundled defaults', self.stderr.getvalue())
            self.assertIn('NOTICE:', (self.path / 'run.log').read_text())

        # Keep the real logging setup and clean up its handlers afterwards.
        logger = cli.logger
        original_handlers = logger.handlers[:]
        original_level = logger.level
        def restore_logger():
            for handler in logger.handlers[:]:
                if handler not in original_handlers:
                    logger.removeHandler(handler)
                    handler.close()
            logger.setLevel(original_level)
        self.addCleanup(restore_logger)
        input_file = self.path / 'input.smi'
        input_file.write_text('CCO ethanol\n')
        argv = ['msani', '-i', str(input_file), '--protonation', '--prefix', str(self.path / 'run')]
        with patch('sys.argv', argv), patch.object(cli, 'clean_data', side_effect=check_notice) as run:
            cli.main()
        run.assert_called_once()

    def test_batch_notice_precedes_submission(self):
        from msani import batchmode

        input_file = self.path / 'input.smi'
        input_file.write_text('CCO ethanol\n')
        def check_notice(args, parser):
            self.assertIn('pH_range=1 from bundled defaults', self.stderr.getvalue())

        # --lazy activates protonation in workers without needing Sanitycheck
        # to normalize the submitter's other options.
        argv = ['msani_batch', '-i', str(input_file), '--lazy']
        with patch('sys.argv', argv), patch.object(batchmode, 'Split_Submit_jobs', side_effect=check_notice) as submit:
            batchmode.main()
        submit.assert_called_once()

    def test_worker_notice_only_goes_to_log(self):
        from msani.io import loggers

        os.environ['MSANI_BATCH_WORKER'] = '1'
        args = parsers.parseArguments(['--protonation'])
        with patch.object(loggers.logger, 'info') as log:
            loggers.protonation_notice(args)
        self.assertEqual(self.stderr.getvalue(), '')
        self.assertIn('pH_range=1', log.call_args.args[0])
        loggers.protonation_notice(args, batch_submission=True)
        self.assertIn('NOTICE:', self.stderr.getvalue())

    def test_runs_without_protonation_and_standardization_are_quiet(self):
        from msani.io.loggers import protonation_notice

        for flags in [[], ['--standardize', '--protonation']]:
            with self.subTest(flags=flags):
                protonation_notice(parsers.parseArguments(flags))
        self.assertEqual(self.stderr.getvalue(), '')

    def test_version_help_and_templates_have_no_migration_notice(self):
        from msani import cli, batchmode

        for main, flags in [(cli.main, ['--version', '-r', '0']),
                            (batchmode.main, ['--version', '-r', '0'])]:
            with self.subTest(main=main):
                with patch('sys.argv', ['msani', *flags]):
                    main()
        with patch('sys.argv', ['msani', '--create_config', '-r', '0']), \
                patch.object(cli, 'generateCustomTemplate') as generate:
            cli.main()
        generate.assert_called_once()
        self.assertEqual(self.stderr.getvalue(), '')

    def test_batch_workers_preserve_explicit_protonation_intent(self):
        import shlex
        from msani.batchmode import parse_flags_single_job

        for flags in [[], ['-r', '1'], ['-p', '7'], ['-r', '0'], ['--lazy']]:
            with self.subTest(flags=flags):
                args, parser = parsers.parseArguments(flags, batch_mode=True)
                worker_flags = parse_flags_single_job(args, parser)
                worker = parsers.Sanitycheck(parsers.parseArguments(shlex.split(worker_flags)))
                self.assertEqual(worker.pH_range, args.pH_range)
                self.assertEqual(worker.protonation, bool(flags))

    def test_python_and_ionizer_command_defaults(self):
        from msani.api import Msani
        from msani.moltransform import ionizer

        # Personal defaults are CLI settings, not Python constructor defaults.
        personal = self.path / 'personal.yaml'
        personal.write_text('PH_RANGE: 0\n')
        os.environ['MSANI_CONFIG'] = str(personal)
        self.assertEqual(Msani().pH_range, 1)
        self.assertEqual(Msani(pH_range=0).pH_range, 0)
        # This compound acquires another protonation state with the new range.
        smiles = 'NCC(F)(F)F'
        self.assertEqual(set(ionizer.Ionizer().ionize(smiles=smiles)),
                         {smiles, '[NH3+]CC(F)(F)F'})
        self.assertEqual(ionizer.Ionizer(pH_range=0).ionize(smiles=smiles), [smiles])
        for flags, expected in [([], 1), (['-r', '0'], 0)]:
            with self.subTest(flags=flags):
                with patch('sys.argv', ['m_ionizer', '-s', 'CCO', *flags]), \
                        patch.object(ionizer, 'Ionizer') as constructor:
                    ionizer.main()
                self.assertEqual(constructor.call_args.kwargs['pH_range'], expected)


if __name__ == '__main__':
    unittest.main()
