"""Reject repaired wheels with external native dependencies or incorrect architectures."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import zipfile


# Architectures `lipo -archs` should report for each macOS wheel platform tag.
TAG_ARCHS = {
    'x86_64': {'x86_64'},
    'arm64': {'arm64'},
    'universal2': {'x86_64', 'arm64'},
    'intel': {'i386', 'x86_64'},
}


def expected_archs(wheel):
    """Derive the architecture set a wheel must contain from its platform tag."""
    plat_tags = Path(wheel).stem.split('-')[-1].split('.')
    archs = set()
    for tag in plat_tags:
        match = re.fullmatch(r'macosx_\d+_\d+_(.+)', tag)
        if not match:
            raise RuntimeError(f'Not a macOS wheel platform tag: {tag}')
        if match.group(1) not in TAG_ARCHS:
            raise RuntimeError(f'Unrecognised macOS architecture tag: {tag}')
        archs |= TAG_ARCHS[match.group(1)]
    return archs


def check_wheel(wheel):
    wanted = expected_archs(wheel)
    with tempfile.TemporaryDirectory(prefix='msani-wheel-audit-') as tmp:
        root = Path(tmp).resolve()
        with zipfile.ZipFile(wheel) as archive:
            archive.extractall(root)
        binaries = [p for p in root.rglob('*') if p.suffix in ('.so', '.dylib')]
        if len([p for p in binaries if p.suffix == '.so']) != 3:
            raise RuntimeError('Expected the three native MolSanitizer extensions.')
        for binary in binaries:
            archs = set(subprocess.check_output(['lipo', '-archs', binary], text=True).split())
            if archs != wanted:
                raise RuntimeError(
                    f'{binary.name}: expected {sorted(wanted)}, found {sorted(archs)}')
            identity = (subprocess.check_output(['otool', '-D', binary], text=True).splitlines()[1:]
                        if binary.suffix == '.dylib' else [])
            lines = subprocess.check_output(['otool', '-L', binary], text=True).splitlines()[1:]
            for line in lines:
                dependency = line.strip().split(' (compatibility version')[0]
                if dependency.startswith(('/usr/lib/', '/System/Library/')):
                    continue
                # A dylib's first entry is its install ID, not a dependency.
                if dependency in identity:
                    continue
                if not dependency.startswith('@loader_path/'):
                    raise RuntimeError(f'{binary.name}: non-local dependency {dependency}')
                resolved = (binary.parent / dependency[len('@loader_path/'):]).resolve()
                if root not in resolved.parents or not resolved.is_file():
                    raise RuntimeError(f'{binary.name}: missing or escaping dependency {dependency}')
        print(f'Native wheel audit passed: {len(binaries)} '
              f'{"/".join(sorted(wanted))} binaries; '
              f'all dependencies wheel-local or system.')


if __name__ == '__main__':
    check_wheel(sys.argv[1])