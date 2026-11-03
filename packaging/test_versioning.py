"""Exercise version inference and release artifacts without an RDKit build."""

from __future__ import annotations

from contextlib import contextmanager
from email.parser import BytesParser
import io
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from scikit_build_core.build import build_sdist, prepare_metadata_for_build_wheel
from scikit_build_core.metadata.setuptools_scm import dynamic_metadata

from check_release_versions import verify_distributions


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@contextmanager
def in_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


class GitVersionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.workspace = Path(self.temporary.name)
        self.root = self.workspace / 'repo'
        self.root.mkdir()
        # Use the production tag and version schemes with a minimal project.
        scm_config = (PROJECT_ROOT / 'pyproject.toml').read_text().split(
            '[tool.setuptools_scm]\n', 1
        )[1]
        (self.root / 'pyproject.toml').write_text('''
[build-system]
requires = ["scikit-build-core>=0.9", "setuptools-scm>=8,<9"]
build-backend = "scikit_build_core.build"
[project]
name = "MolSanitizer"
dynamic = ["version"]
[tool.scikit-build]
metadata.version.provider = "scikit_build_core.metadata.setuptools_scm"
wheel.packages = ["msani"]
sdist.exclude = [".git_archival.txt"]
[tool.setuptools_scm]
''' + scm_config)
        (self.root / 'msani').mkdir()
        (self.root / 'msani/__init__.py').write_text('')
        (self.root / 'CMakeLists.txt').write_text(
            'cmake_minimum_required(VERSION 3.20)\nproject(version_fixture LANGUAGES NONE)\n'
        )
        for name in ('.git_archival.txt', '.gitattributes'):
            (self.root / name).write_text((PROJECT_ROOT / name).read_text())
        self.git('init', '--quiet')
        self.git('config', 'user.email', 'version-test@example.invalid')
        self.git('config', 'user.name', 'Version Test')
        self.git('add', '.')
        self.git('commit', '--quiet', '-m', 'Initial release')
        self.git('tag', '0.7.5')

    def git(self, *arguments):
        return subprocess.check_output(
            ['git', '-C', str(self.root), *arguments], text=True,
            stderr=subprocess.STDOUT,
        ).strip()

    def inferred_version(self, root=None):
        with in_directory(root or self.root):
            return dynamic_metadata('version')

    def next_commit(self):
        (self.root / 'feature.txt').write_text('A feature after the release\n')
        self.git('add', '.')
        self.git('commit', '--quiet', '-m', 'Add a feature')

    def test_clean_release_and_detached_tag(self):
        self.assertEqual(self.inferred_version(), '0.7.5')
        self.git('checkout', '--quiet', '--detach', '0.7.5')
        self.assertEqual(self.inferred_version(), '0.7.5')

    def test_later_commit_ignores_sdk_tags(self):
        self.next_commit()
        self.git('tag', 'sdk-2025.09.1')
        self.assertRegex(self.inferred_version(), r'^0\.7\.6\.dev1\+g[0-9a-f]+$')

    def test_dirty_release_is_identified(self):
        (self.root / 'msani/__init__.py').write_text('# Local changes\n')
        self.assertRegex(self.inferred_version(), r'^0\.7\.6\.dev0\+g[0-9a-f]+\.d[0-9]{8}$')

    def test_release_candidate(self):
        self.next_commit()
        self.git('tag', '0.7.6rc1')
        self.assertEqual(self.inferred_version(), '0.7.6rc1')

    def test_sdist_preserves_version_without_git(self):
        self.next_commit()
        expected = self.inferred_version()
        out = self.workspace / 'dist'
        out.mkdir()
        with in_directory(self.root):
            filename = build_sdist(str(out))
        extracted = self.workspace / 'extracted'
        extracted.mkdir()
        with tarfile.open(out / filename) as archive:
            self.assertFalse(any(name.endswith('.git_archival.txt')
                                 for name in archive.getnames()))
            archive.extractall(extracted, filter='data')
        source = next(extracted.iterdir())
        self.assertFalse((source / '.git').exists())
        self.assertEqual(self.inferred_version(source), expected)
        metadata_out = self.workspace / 'metadata'
        metadata_out.mkdir()
        with in_directory(source):
            metadata_directory = prepare_metadata_for_build_wheel(str(metadata_out))
        metadata = BytesParser().parsebytes(
            (metadata_out / metadata_directory / 'METADATA').read_bytes()
        )
        self.assertEqual(metadata['Version'], expected)

    def test_git_archive_preserves_release(self):
        archive_path = self.workspace / 'source.tar'
        with archive_path.open('wb') as output:
            subprocess.run(['git', '-C', str(self.root), 'archive', '0.7.5'],
                           stdout=output, check=True)
        extracted = self.workspace / 'archive'
        extracted.mkdir()
        with tarfile.open(archive_path) as archive:
            archive.extractall(extracted, filter='data')
        self.assertEqual(self.inferred_version(extracted), '0.7.5')

    def test_linux_copy_preserves_original_version(self):
        self.next_commit()
        expected = self.inferred_version()
        copied = self.workspace / 'copied-source'
        shutil.copytree(self.root, copied, ignore=shutil.ignore_patterns('.git', '.git_archival.txt'))
        with self.assertRaisesRegex(ValueError, 'unable to detect version'):
            self.inferred_version(copied)
        with patch.dict(os.environ, {'SETUPTOOLS_SCM_PRETEND_VERSION_FOR_MOLSANITIZER': expected}):
            self.assertEqual(self.inferred_version(copied), expected)


class ReleaseBundleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def bundle(self, wheel_version='0.7.5', sdist_version='0.7.5', name='MolSanitizer'):
        metadata = lambda version: f'Metadata-Version: 2.2\nName: {name}\nVersion: {version}\n'
        with zipfile.ZipFile(self.root / 'molsanitizer-0.7.5-py3-none-any.whl', 'w') as wheel:
            wheel.writestr('molsanitizer-0.7.5.dist-info/METADATA', metadata(wheel_version))
        content = metadata(sdist_version).encode()
        with tarfile.open(self.root / 'molsanitizer-0.7.5.tar.gz', 'w:gz') as sdist:
            member = tarfile.TarInfo('molsanitizer-0.7.5/PKG-INFO')
            member.size = len(content)
            sdist.addfile(member, io.BytesIO(content))

    def test_matching_release(self):
        self.bundle()
        self.assertEqual(verify_distributions(self.root, '0.7.5'), '0.7.5')

    def test_manual_development_bundle(self):
        version = '0.7.6.dev12+gabc1234'
        self.bundle(version, version)
        self.assertEqual(verify_distributions(self.root), version)

    def test_checks_metadata_instead_of_filenames(self):
        self.bundle(wheel_version='0.7.4')
        with self.assertRaisesRegex(ValueError, 'does not match'):
            verify_distributions(self.root, '0.7.5')

    def test_sdist_must_match_release_tag(self):
        self.bundle('0.7.6.dev1+gabc1234', '0.7.6.dev1+gabc1234')
        with self.assertRaisesRegex(ValueError, 'does not match'):
            verify_distributions(self.root, '0.7.5')

    def test_rejects_another_package(self):
        self.bundle(name='AnotherPackage')
        with self.assertRaisesRegex(ValueError, 'expected MolSanitizer'):
            verify_distributions(self.root, '0.7.5')

    def test_rejects_prefixed_release_tag(self):
        self.bundle()
        with self.assertRaisesRegex(ValueError, 'Invalid bare release tag'):
            verify_distributions(self.root, 'v0.7.5')

    def test_requires_complete_bundle(self):
        self.bundle()
        next(self.root.glob('*.tar.gz')).unlink()
        with self.assertRaisesRegex(ValueError, 'exactly one source distribution'):
            verify_distributions(self.root, '0.7.5')


if __name__ == '__main__':
    unittest.main()
