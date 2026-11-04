"""Check release filtering without contacting PyPI."""

import unittest

from packaging.tags import Tag

from discover_rdkit_versions import select_versions


class RDKitDiscoveryTests(unittest.TestCase):
    tags = {Tag('cp312', 'cp312', 'manylinux_2_28_x86_64'),
            Tag('cp311', 'abi3', 'manylinux_2_28_x86_64')}

    def file(self, version, tag='cp312-cp312-manylinux_2_28_x86_64', **fields):
        return {'filename': f'rdkit-{version}-{tag}.whl', 'yanked': False, **fields}

    def select(self, *files):
        return select_versions({'files': list(files)}, '2025.3.1', self.tags, '3.12.14')

    def test_includes_every_patch_in_numeric_order_and_future_release(self):
        self.assertEqual(self.select(self.file('2025.3.10'), self.file('2025.3.2'),
                                     self.file('2025.3.1'), self.file('2027.9.1')),
                         ['2025.3.1', '2025.3.2', '2025.3.10', '2027.9.1'])

    def test_excludes_older_prerelease_dev_local_and_yanked_files(self):
        self.assertEqual(self.select(
            self.file('2025.3.1'), self.file('2024.9.6'), self.file('2025.3.2rc1'),
            self.file('2025.3.2.dev1'), self.file('2025.3.2+local'),
            self.file('2025.3.2', yanked=True), self.file('2025.3.3', yanked=''),
            self.file('2025.3.4', yanked='broken release')),
            ['2025.3.1'])

    def test_requires_a_compatible_wheel_and_python_requirement(self):
        self.assertEqual(self.select(
            self.file('2025.3.1'), self.file('2025.3.2', 'cp313-cp313-manylinux_2_28_x86_64'),
            self.file('2025.3.3', 'cp312-cp312-win_amd64'),
            self.file('2025.3.4', 'cp312-cp312-manylinux_2_28_aarch64'),
            self.file('2025.3.5', **{'requires-python': '>=3.13'}),
            {'filename': 'rdkit-2025.3.6.tar.gz'}), ['2025.3.1'])

    def test_handles_duplicate_and_abi3_wheels_and_postreleases(self):
        self.assertEqual(self.select(self.file('2025.3.1'), self.file('2025.3.1'),
                                     self.file('2025.3.2', 'cp311-abi3-manylinux_2_28_x86_64'),
                                     self.file('2025.3.2.post1')),
                         ['2025.3.1', '2025.3.2', '2025.3.2.post1'])

    def test_does_not_select_a_release_when_only_its_foreign_wheel_is_unyanked(self):
        self.assertEqual(self.select(self.file('2025.3.1'), self.file('2025.3.2', yanked=True),
                                     self.file('2025.3.2', 'cp312-cp312-win_amd64')),
                         ['2025.3.1'])

    def test_empty_selection_fails_instead_of_skipping_all_tests(self):
        with self.assertRaisesRegex(ValueError, 'No compatible stable RDKit wheels'):
            self.select(self.file('2024.9.6'))


if __name__ == '__main__':
    unittest.main()
