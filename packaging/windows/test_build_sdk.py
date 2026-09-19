"""Host-independent tests for SDK safeguards; no compiler/network required."""
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import json

import build_sdk


class TestSDKHelpers(unittest.TestCase):
    def test_msvc_selected_despite_conda_gcc_environment(self):
        env = {'PATH': 'conda-mingw', 'CC': 'gcc', 'CXX': 'g++',
               'CXXFLAGS': '-std=gnu++11', 'INCLUDE': 'VS headers'}
        with patch.object(build_sdk.shutil, 'which', return_value='/VS Tools/bin/cl.exe'):
            compiler, options = build_sdk.select_msvc(env)
        self.assertEqual(options, {'CMAKE_C_COMPILER': compiler,
                                   'CMAKE_CXX_COMPILER': compiler})
        self.assertTrue(env['PATH'].startswith('/VS Tools/bin'))
        self.assertEqual(env['INCLUDE'], 'VS headers')
        for key in ('CC', 'CXX', 'CXXFLAGS'):
            self.assertNotIn(key, env)

    def test_missing_msvc_rejected(self):
        with patch.object(build_sdk.shutil, 'which', return_value=None):
            with self.assertRaisesRegex(RuntimeError, 'MSVC'):
                build_sdk.select_msvc({'PATH': 'conda-mingw'})

    def test_pinned_boost_accepted(self):
        build_sdk.validate_packages([{'name': name, 'version': build_sdk.BOOST_VERSION}
                                     for name in ('libboost', 'libboost-devel')])

    def test_missing_or_wrong_boost_rejected(self):
        for packages in ([], [{'name': 'libboost', 'version': '1.86.0'}]):
            with self.assertRaisesRegex(RuntimeError, 'must be'):
                build_sdk.validate_packages(packages)

    def test_python_packages_rejected(self):
        packages = [{'name': name, 'version': build_sdk.BOOST_VERSION}
                    for name in ('libboost', 'libboost-devel')]
        for name in ('rdkit', 'libboost-python', 'libboost-python-devel'):
            with self.subTest(name=name), self.assertRaisesRegex(RuntimeError, 'must not include'):
                build_sdk.validate_packages(packages + [{'name': name, 'version': '1'}])

    def test_python_dependencies_rejected(self):
        for name in ('python312.dll', 'python3.dll', 'boost_python312-vc143-mt-x64-1_85.dll',
                     'boost_numpy314.dll'):
            with self.subTest(name=name), patch.object(build_sdk, 'capture', return_value=name):
                with self.assertRaisesRegex(RuntimeError, 'Python-dependent'):
                    build_sdk.reject_python_dlls([Path('RDKitGraphMol.dll')], {}, io.StringIO())

    def test_native_dependencies_accepted(self):
        with patch.object(build_sdk, 'capture', return_value='KERNEL32.dll\nboost_serialization.dll'):
            build_sdk.reject_python_dlls([Path('RDKitGraphMol.dll')], {}, io.StringIO())

    def test_copy_license_from_exact_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            meta = root / 'env/conda-meta'
            meta.mkdir(parents=True)
            licenses = root / 'cache/libboost-exact/info/licenses'
            licenses.mkdir(parents=True)
            (licenses / 'LICENSE_1_0.txt').write_text('fixture')
            (meta / 'libboost-exact.json').write_text(json.dumps({
                'name': 'libboost', 'link': {'source': str(licenses.parent.parent)}}))
            build_sdk.copy_boost_licenses(root / 'env', root / 'notices')
            self.assertEqual((root / 'notices/libboost/LICENSE_1_0.txt').read_text(), 'fixture')

    def test_missing_license_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RuntimeError, 'license texts missing'):
                build_sdk.copy_boost_licenses(Path(tmp), Path(tmp) / 'notices')


if __name__ == '__main__':
    unittest.main()
