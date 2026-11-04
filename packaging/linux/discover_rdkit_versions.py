"""Discover stable RDKit releases with wheels compatible with this interpreter."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import sys
import time
from urllib.error import URLError
from urllib.request import Request, urlopen

from packaging.specifiers import InvalidSpecifier, SpecifierSet
from packaging.tags import Tag, sys_tags
from packaging.utils import InvalidWheelFilename, canonicalize_name, parse_wheel_filename
from packaging.version import Version


INDEX_URL = 'https://pypi.org/simple/rdkit/'


def select_versions(index: dict, minimum: str, supported_tags: set[Tag],
                    python_version: str) -> list[str]:
    lower_bound = Version(minimum)
    selected = set()
    for file in index['files']:
        # An empty string is also a yanked reason, not an unyanked file.
        yanked = file.get('yanked', False)
        if yanked is True or isinstance(yanked, str):
            continue
        try:
            name, version, _, tags = parse_wheel_filename(file['filename'])
        except InvalidWheelFilename:
            continue
        if canonicalize_name(name) != 'rdkit':
            continue
        if version < lower_bound or version.is_prerelease or version.is_devrelease or version.local:
            continue
        if not tags.intersection(supported_tags):
            continue
        try:
            requires_python = SpecifierSet(file.get('requires-python') or '')
        except InvalidSpecifier as error:
            raise ValueError(f"Invalid Requires-Python for {file['filename']}") from error
        if requires_python.contains(python_version):
            selected.add(version)

    if not selected:
        raise ValueError(f'No compatible stable RDKit wheels found from {minimum} onward')
    if len(selected) > 256:
        raise ValueError('More than 256 RDKit versions found; split the GitHub Actions matrix')
    return [str(version) for version in sorted(selected)]


def fetch_index() -> dict:
    request = Request(INDEX_URL, headers={'Accept': 'application/vnd.pypi.simple.v1+json'})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except (URLError, TimeoutError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError('Unreachable')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum', default='2025.3.1')
    parser.add_argument('--index-json', type=Path, help='Read a saved index response instead of PyPI')
    args = parser.parse_args()
    if (sys.version_info[:2] != (3, 12) or platform.python_implementation() != 'CPython'
            or platform.system() != 'Linux' or platform.machine() != 'x86_64'):
        parser.error('Run discovery on Linux x86_64 with CPython 3.12, matching the test container')
    index = json.loads(args.index_json.read_text()) if args.index_json else fetch_index()
    versions = select_versions(index, args.minimum, set(sys_tags()), platform.python_version())
    encoded = json.dumps(versions)
    print(encoded)
    if output := os.environ.get('GITHUB_OUTPUT'):
        with Path(output).open('a') as stream:
            stream.write(f'rdkit={encoded}\n')
    if summary := os.environ.get('GITHUB_STEP_SUMMARY'):
        with Path(summary).open('a') as stream:
            stream.write(f'### RDKit compatibility: {len(versions)} releases\n\n')
            stream.write(f'Linux x86_64 / CPython {platform.python_version()}, minimum {args.minimum}.\n\n')
            stream.write('\n'.join(f'- `{version}`' for version in versions) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
