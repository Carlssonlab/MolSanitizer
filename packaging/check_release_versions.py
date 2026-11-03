"""Check the versions recorded inside a wheel and source distribution bundle."""

from __future__ import annotations

import argparse
from email.parser import BytesParser
from pathlib import Path
import re
import tarfile
import zipfile


def distribution_version(path: Path) -> str:
    if path.suffix == '.whl':
        with zipfile.ZipFile(path) as archive:
            names = [name for name in archive.namelist()
                     if name.count('/') == 1 and name.endswith('.dist-info/METADATA')]
            if len(names) != 1:
                raise ValueError(f'{path.name}: expected one wheel METADATA file')
            content = archive.read(names[0])
    else:
        with tarfile.open(path, 'r:gz') as archive:
            members = [member for member in archive.getmembers()
                       if member.name.count('/') == 1 and member.name.endswith('/PKG-INFO')
                       and member.isfile()]
            if len(members) != 1:
                raise ValueError(f'{path.name}: expected one top-level PKG-INFO file')
            with archive.extractfile(members[0]) as metadata_file:
                content = metadata_file.read()

    metadata = BytesParser().parsebytes(content)
    if metadata.get('Name', '').lower() != 'molsanitizer':
        raise ValueError(f'{path.name}: expected MolSanitizer package metadata')
    version = metadata.get('Version')
    if not version:
        raise ValueError(f'{path.name}: missing Version metadata')
    return version


def verify_distributions(directory: Path, expected_version: str | None = None) -> str:
    wheels = sorted(directory.glob('*.whl'))
    sdists = sorted(directory.glob('*.tar.gz'))
    if not wheels or len(sdists) != 1:
        raise ValueError('The bundle must contain wheels and exactly one source distribution')
    if expected_version is not None and not re.fullmatch(
        r'[0-9]+\.[0-9]+\.[0-9]+(?:(?:a|b|rc)[0-9]+)?', expected_version
    ):
        raise ValueError(f'Invalid bare release tag: {expected_version!r}')

    expected = expected_version or distribution_version(sdists[0])
    for path in [*sdists, *wheels]:
        actual = distribution_version(path)
        if actual != expected:
            raise ValueError(f'{path.name}: version {actual!r} does not match {expected!r}')
    return expected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--expected-version', help='Release tag that every artifact must match')
    args = parser.parse_args()
    try:
        version = verify_distributions(args.directory, args.expected_version)
    except (ValueError, OSError, tarfile.TarError, zipfile.BadZipFile) as error:
        parser.exit(1, f'{error}\n')
    print(f'All distributions match MolSanitizer {version}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
