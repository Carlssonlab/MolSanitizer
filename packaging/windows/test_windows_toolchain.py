"""Compiler settings passed to pip must pin both languages, including spaced paths."""
import unittest
from unittest.mock import patch

import build_sdk
import build_wheels
import windows_toolchain


class TestWheelToolchain(unittest.TestCase):
    def test_builders_share_compiler_selection(self):
        self.assertIs(build_sdk.select_msvc, build_wheels.select_msvc)

    def test_pip_settings_pin_both_languages(self):
        env = {'PATH': 'anaconda/mingw/bin', 'CC': 'gcc', 'CXX': 'g++'}
        with patch.object(windows_toolchain.shutil, 'which', return_value='/VS Tools/bin/cl.exe'):
            _, options = windows_toolchain.select_msvc(env)
        self.assertEqual(windows_toolchain.wheel_compiler_settings(options), [
            '-Ccmake.define.CMAKE_C_COMPILER=/VS Tools/bin/cl.exe',
            '-Ccmake.define.CMAKE_CXX_COMPILER=/VS Tools/bin/cl.exe',
        ])
        self.assertNotIn('CC', env)
        self.assertNotIn('CXX', env)


if __name__ == '__main__':
    unittest.main()
