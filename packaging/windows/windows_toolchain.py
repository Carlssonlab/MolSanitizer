"""Compiler selection shared by the Windows SDK and wheel builders."""
import os
from pathlib import Path
import shutil


def select_msvc(env):
    compiler = shutil.which('cl', path=env['PATH'])
    if not compiler or Path(compiler).name.lower() != 'cl.exe':
        raise RuntimeError('MSVC cl.exe not found; run the build .cmd wrapper in an x64 VS shell.')
    # Ninja does not imply MSVC. Remove inherited GCC compiler selections and
    # flags, and prefer VS tools (including link.exe) over Anaconda's MinGW.
    for key in list(env):
        if key.upper() in ('CC', 'CXX', 'CFLAGS', 'CXXFLAGS', 'CPPFLAGS', 'LDFLAGS'):
            env.pop(key)
    env['PATH'] = str(Path(compiler).parent) + os.pathsep + env['PATH']
    return compiler, {
        'CMAKE_C_COMPILER': Path(compiler).as_posix(),
        'CMAKE_CXX_COMPILER': Path(compiler).as_posix(),
    }


def wheel_compiler_settings(options):
    return [f'-Ccmake.define.{key}={value}' for key, value in options.items()]
