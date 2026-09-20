"""Persistent user defaults, independent of the installed package directory."""
import argparse
import json
import math
import os
from importlib.resources import files
from pathlib import Path
import sys

from platformdirs import user_config_path
import yaml


class ConfigError(ValueError):
    """An invalid or inaccessible configuration."""


def config_directory():
    return user_config_path('msani', appauthor=False, roaming=False)


def absolute_path(value):
    return Path(value).expanduser().resolve()


def bundled_defaults():
    return yaml.safe_load(files('msani').joinpath('msani_configurations.yaml').read_text())


def selected_path():
    override = os.environ.get('MSANI_CONFIG')
    if override is not None:
        if not override.strip():
            raise ConfigError('MSANI_CONFIG must contain a configuration filename.')
        return absolute_path(override), True
    pointer = config_directory() / 'settings.json'
    if pointer.exists():
        try:
            value = json.loads(pointer.read_text(encoding='utf-8'))['config_path']
            if not isinstance(value, str) or not value.strip():
                raise ValueError('config_path must be a nonempty string')
            path = Path(value).expanduser()
            if not path.is_absolute():
                raise ValueError('config_path must be absolute')
            return path, True
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise ConfigError(f'{pointer}: {exc}') from exc
    return config_directory() / 'msani_configurations.yaml', False


def read_config(path):
    try:
        data = yaml.safe_load(path.read_text(encoding='utf-8'))
        if data is None:
            data = {}
        if not isinstance(data, dict):
            raise ValueError('expected a YAML mapping of setting names to values')
        defaults = bundled_defaults()
        positive = {'NUMCONFS', 'MAX_STEREOISOMERS', 'TIMEOUT', 'LINES_PER_JOB',
                    'TIME_LIMIT', 'MAX_ARRAY_SIZE', 'MAX_JOBS', 'MAX_LIMIT_PROJECT',
                    'WHOLE_NODE_CORES'}
        for key, value in data.items():
            if key not in defaults:
                raise ValueError(f'unknown setting {key!r}')
            default = defaults[key]
            valid = True
            if key in {'CORINA', 'SLURM_ACCOUNT', 'SLURM_PARTITION'}:
                valid = value is None or (isinstance(value, str) and bool(value.strip()))
            elif isinstance(default, bool):
                valid = isinstance(value, bool)
            elif key in {'PH', 'PH_RANGE', 'TIMEOUT', 'ENERGY_WINDOW'}:
                valid = type(value) in (int, float) and math.isfinite(value)
            elif isinstance(default, int):
                valid = type(value) is int
            else:
                valid = isinstance(value, str)
            if not valid:
                raise ValueError(f'{key}: invalid value {value!r}')
            if key in positive and value <= 0:
                raise ValueError(f'{key}: must be greater than zero')
            if key in {'PH_RANGE', 'ENERGY_WINDOW'} and value < 0:
                raise ValueError(f'{key}: must be nonnegative')
            if key == 'EMBED_METHOD' and value not in {'rdkit', 'obabel', 'corina'}:
                raise ValueError(f'{key}: expected rdkit, obabel, or corina')
        return data
    except (OSError, ValueError, yaml.YAMLError) as exc:
        raise ConfigError(f'{path}: {exc}') from exc


def load_defaults():
    defaults = bundled_defaults()
    path, explicit = selected_path()
    sources = dict.fromkeys(defaults, 'bundled defaults')
    if explicit or path.exists():
        overrides = read_config(path)
        defaults.update(overrides)
        sources.update(dict.fromkeys(overrides, str(path)))
    return defaults, sources


def remember(path):
    directory = config_directory()
    directory.mkdir(parents=True, exist_ok=True)
    # Replace atomically so another process cannot read a partial pointer.
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=directory,
                                     delete=False) as handle:
        temporary = Path(handle.name)
        json.dump({'config_path': str(path)}, handle, indent=2)
        handle.write('\n')
    try:
        temporary.replace(directory / 'settings.json')
    finally:
        temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(prog='msani config', description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    init = commands.add_parser('init', help='Create personal defaults')
    init.add_argument('--path', type=absolute_path)
    init.add_argument('--force', action='store_true', help='Overwrite an existing file')
    use = commands.add_parser('use', help='Remember an existing configuration')
    use.add_argument('path', type=absolute_path)
    commands.add_parser('path', help='Print the active configuration location')
    commands.add_parser('show', help='Show defaults and their sources')
    commands.add_parser('reset-path', help='Return to the standard configuration location')
    args = parser.parse_args(argv)
    try:
        if args.command == 'init':
            path = args.path or config_directory() / 'msani_configurations.yaml'
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open('w' if args.force else 'x', encoding='utf-8') as handle:
                handle.write(files('msani').joinpath('msani_configurations.yaml').read_text())
            remember(path)
            print(f'Created configuration: {path}')
        elif args.command == 'use':
            read_config(args.path)
            remember(args.path)
            print(f'Using configuration: {args.path}')
        elif args.command == 'reset-path':
            (config_directory() / 'settings.json').unlink(missing_ok=True)
            print(f'Default configuration location: {config_directory() / "msani_configurations.yaml"}')
        elif args.command == 'path':
            print(selected_path()[0])
        else:
            values, sources = load_defaults()
            print(yaml.safe_dump({'settings': values, 'sources': sources}, sort_keys=False), end='')
        if args.command in {'init', 'use', 'reset-path'} and 'MSANI_CONFIG' in os.environ:
            print('MSANI_CONFIG is set and still overrides the remembered location.', file=sys.stderr)
    except (ConfigError, OSError) as exc:
        parser.exit(2, f'Configuration error: {exc}\n')


def dispatch(argv):
    """Handle configuration commands before normal argument parsing."""
    if argv and argv[0] == 'config':
        main(argv[1:])
        return True
    return False
