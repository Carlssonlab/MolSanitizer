"""Bootstrap msani_build and a reusable Windows RDKit 2025.9.1 shared SDK."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from build_wheels import run
from windows_toolchain import select_msvc

RDKIT_REF = 'Release_2025_09_1'
BOOST_VERSION = '1.85.0'
ENV_PACKAGES = ['python=3.12', 'cmake>=3.28,<4', 'ninja', 'git', 'eigen=3.4',
                f'libboost={BOOST_VERSION}', f'libboost-devel={BOOST_VERSION}']


def capture(command, env=None):
    return subprocess.check_output(command, env=env, text=True, encoding='utf-8').strip()


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def validate_packages(packages):
    versions = {p['name']: p['version'] for p in packages}
    for name in ('libboost', 'libboost-devel'):
        if versions.get(name) != BOOST_VERSION:
            raise RuntimeError(f'{name} must be {BOOST_VERSION}; use a new SDK root.')
    if any(name in versions for name in ('libboost-python', 'libboost-python-devel',
                                        'libboost-numpy', 'rdkit')):
        raise RuntimeError('SDK tools must not include Python RDKit or Boost.Python packages.')


def copy_boost_licenses(tools, destination):
    # Conda retains license texts under info/licenses in its extracted cache.
    # Locate the exact installed package via its conda-meta record, not a glob
    # that might select another Boost version from the cache.
    copied = False
    for meta in (tools / 'conda-meta').glob('libboost*.json'):
        record = json.loads(meta.read_text(encoding='utf-8'))
        source = record.get('link', {}).get('source')
        if not source:
            continue
        licenses = Path(source) / 'info/licenses'
        if licenses.is_dir():
            shutil.copytree(licenses, destination / record['name'], dirs_exist_ok=True)
            copied = True
    if not copied:
        raise RuntimeError('Boost license texts missing from conda package cache. '
                           'Preserve the cache until SDK creation completes.')


def reject_python_dlls(dlls, env, log):
    for dll in dlls:
        report = capture(['dumpbin', '/nologo', '/dependents', str(dll)], env)
        log.write(f'\n{dll}\n{report}\n')
        if re.search(r'(?i)\b(?:python\d*|(?:lib)?boost_(?:python|numpy)[\w.-]*)\.dll\b', report):
            raise RuntimeError(f'Python-dependent SDK DLL: {dll}; see dependency report')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path,
                        default=Path.home() / 'msani-sdk' / 'windows-2025.9.1-conda')
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--conda-lock', type=Path, help='Exact win-64 conda export from a previous build')
    parser.add_argument('--rdkit-commit', help='Require the recorded RDKit commit when rebuilding elsewhere')
    args = parser.parse_args()
    if args.conda_lock and not args.conda_lock.is_file():
        parser.error('--conda-lock file does not exist')
    if sys.platform != 'win32':
        parser.error('Run this script in native Windows Command Prompt, not WSL.')
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    conda = os.environ.get('CONDA_EXE') or shutil.which('conda.exe')
    if not conda or not Path(conda).is_file():
        parser.error('Activate conda first; CONDA_EXE must identify conda.exe.')
    if not shutil.which('cl') or not shutil.which('dumpbin'):
        parser.error('Use build_sdk.cmd or an x64 Native Tools Command Prompt.')
    arch = os.environ.get('VSCMD_ARG_TGT_ARCH', '')
    if arch.lower() not in ('x64', 'amd64'):
        parser.error('Initialize VS tools for x64 (VSCMD_ARG_TGT_ARCH must be x64).')

    root = args.root.resolve()
    owner = root / 'msani-sdk-build.json'
    if root.exists() and any(root.iterdir()) and not owner.exists():
        parser.error(f'{root} is not an SDK workspace created by this script; choose an empty --root.')
    root.mkdir(parents=True, exist_ok=True)
    identity = {'schema': 2, 'rdkit_ref': RDKIT_REF, 'boost_version': BOOST_VERSION,
                'architecture': 'x64', 'runtime': '/MD', 'rdkit_linkage': 'shared',
                'boost_linkage': 'shared', 'boost_provider': 'conda-forge'}
    if owner.exists() and json.loads(owner.read_text()) != identity:
        parser.error('SDK settings changed. Use a new --root to avoid mixing SDK artifacts.')
    json_write(owner, identity)
    with (root / 'build-sdk.log').open('a', encoding='utf-8') as log:
        env = os.environ.copy()
        env['CONDA_SUBDIR'] = 'win-64'
        for key in ('PYTHONPATH', 'PYTHONHOME'):
            env.pop(key, None)
        tools = root / 'msani_build'
        if not (tools / 'conda-meta/history').exists():
            if tools.exists() and any(tools.iterdir()):
                raise RuntimeError('Incomplete SDK conda environment; choose a new --root.')
            dependencies = (['--file', str(args.conda_lock.resolve())] if args.conda_lock
                            else ['--override-channels', '--channel', 'conda-forge',
                                  '--strict-channel-priority', *ENV_PACKAGES])
            run([conda, 'create', '--yes', '--prefix', tools, *dependencies],
                env=env, cwd=root, log=log)
        packages = json.loads(capture([conda, 'list', '--prefix', str(tools), '--json'], env))
        validate_packages(packages)
        explicit = capture([conda, 'list', '--prefix', str(tools), '--explicit', '--md5'], env)
        if args.conda_lock:
            # Ignore machine-specific comments in explicit exports.
            entries = lambda text: {line.strip() for line in text.splitlines()
                                    if line.strip() and not line.startswith('#')}
            if entries(explicit) != entries(args.conda_lock.read_text(encoding='utf-8')):
                raise RuntimeError('Installed SDK environment differs from --conda-lock.')
        (root / 'conda-win-64.lock.txt').write_text(explicit + '\n', encoding='utf-8')
        # Never silently change an existing environment's dependencies.
        cmake = tools / 'Library/bin/cmake.exe'
        ninja = tools / 'Library/bin/ninja.exe'
        git = tools / 'Library/bin/git.exe'
        if not git.exists():
            git = tools / 'Library/bin/git'  # diagnostic below handles absent tools
        if not all(p.exists() for p in (cmake, ninja, git, tools / 'python.exe', tools / 'Library/include/eigen3')):
            raise RuntimeError('Existing msani_build lacks required tools. Install python=3.12, '
                               'cmake>=3.28,<4, ninja, git and eigen=3.4 there, then retry.')
        for key in list(env):
            if key.upper().startswith(('PYTHON', 'CMAKE', 'CONDA')) or key.upper() in ('RDBASE', 'VIRTUAL_ENV'):
                env.pop(key, None)
        env['PATH'] = os.pathsep.join([str(tools), str(tools / 'Scripts'),
                                      str(tools / 'Library/bin'), env['PATH']])
        env['PYTHONNOUSERSITE'] = '1'
        env['PYTHONUTF8'] = '1'
        compiler, compiler_options = select_msvc(env)
        print(f'Pinned C/C++ compiler: {compiler}')
        toolchain = {'compiler': compiler, 'vc_tools_version': os.environ.get('VCToolsVersion'),
                     'windows_sdk': os.environ.get('WindowsSDKVersion'),
                     'cmake': capture([str(cmake), '--version'], env),
                     'conda_packages': packages}
        record = root / 'toolchain.json'
        if record.exists() and json.loads(record.read_text()) != toolchain:
            raise RuntimeError('Toolchain/environment changed. Use a new --root for a clean SDK build.')
        json_write(record, toolchain)
        boost = tools / 'Library'
        rdkit_src = root / 'rdkit-source'
        if not rdkit_src.exists():
            # Clone atomically so a failed download does not poison retries.
            with tempfile.TemporaryDirectory(dir=root) as tmp:
                clone = Path(tmp) / 'source'
                run([git, 'clone', '--depth', '1', '--branch', RDKIT_REF,
                     'https://github.com/rdkit/rdkit.git', clone], env=env, cwd=root, log=log)
                clone.rename(rdkit_src)
        commit = capture([str(git), '-C', str(rdkit_src), 'rev-parse', 'HEAD'], env)
        if args.rdkit_commit and commit.lower() != args.rdkit_commit.lower():
            raise RuntimeError('RDKit commit differs from --rdkit-commit; refusing to build.')
        expected = capture([str(git), '-C', str(rdkit_src), 'rev-parse', RDKIT_REF + '^{commit}'], env)
        if commit != expected or capture([str(git), '-C', str(rdkit_src), 'diff', 'HEAD', '--'], env):
            raise RuntimeError('Cached RDKit checkout differs from the pinned tag; use a fresh --root.')
        sdk = root / 'rdkit'
        # Do not reuse the original directory: older bootstraps could cache
        # Anaconda's MinGW compiler there. Keep it intact for diagnostics.
        build = root / 'rdkit-build-msvc'
        prefix = boost.as_posix()
        options = {
            **compiler_options,
            'CMAKE_BUILD_TYPE': 'Release', 'CMAKE_INSTALL_PREFIX': sdk.as_posix(),
            'CMAKE_PREFIX_PATH': prefix, 'CMAKE_MAKE_PROGRAM': ninja.as_posix(),
            'CMAKE_MSVC_RUNTIME_LIBRARY': 'MultiThreadedDLL',
            'CMAKE_POLICY_DEFAULT_CMP0091': 'NEW',
            'Python3_EXECUTABLE': (tools / 'python.exe').as_posix(),
            'RDK_INSTALL_INTREE': 'OFF', 'RDK_INSTALL_DLLS_MSVC': 'ON',
            'RDK_INSTALL_STATIC_LIBS': 'OFF', 'RDK_BUILD_STATIC_LIBS_ONLY': 'OFF',
            'RDK_INSTALL_DEV_COMPONENT': 'ON', 'RDK_BUILD_PYTHON_WRAPPERS': 'OFF',
            'RDK_BUILD_SWIG_WRAPPERS': 'OFF', 'RDK_BUILD_CPP_TESTS': 'OFF',
            'BUILD_TESTING': 'OFF', 'RDK_BUILD_THREADSAFE_SSS': 'ON',
            'RDK_BUILD_INCHI_SUPPORT': 'OFF', 'RDK_BUILD_FREETYPE_SUPPORT': 'OFF',
            'RDK_INSTALL_COMIC_FONTS': 'OFF', 'RDK_BUILD_MAEPARSER_SUPPORT': 'OFF',
            'RDK_BUILD_PUBCHEMSHAPE_SUPPORT': 'OFF', 'RDK_BUILD_COMPRESSED_SUPPLIERS': 'OFF',
            'Boost_USE_STATIC_LIBS': 'OFF', 'Boost_USE_STATIC_RUNTIME': 'OFF',
        }
        run([cmake, '-S', rdkit_src, '-B', build, '-G', 'Ninja',
             *[f'-D{k}={v}' for k, v in options.items()]], env=env, cwd=root, log=log)
        run([cmake, '--build', build, '--parallel', args.jobs], env=env, cwd=root, log=log)
        run([cmake, '--install', build], env=env, cwd=root, log=log)
        # RDKit/Boost may install DLLs under lib rather than bin on Windows.
        dll_dirs = [p for prefix_dir in (sdk, boost) for p in (prefix_dir / 'bin', prefix_dir / 'lib')
                    if p.exists()]
        rdkit_dlls = list(sdk.rglob('RDKit*.dll'))
        boost_dlls = [p for p in boost.rglob('*.dll') if 'boost_' in p.name.lower()]
        if not rdkit_dlls or not boost_dlls:
            raise RuntimeError('Expected installed shared RDKit and Boost DLLs')
        runtime_env = env.copy()
        runtime_env['PATH'] = os.pathsep.join([*map(str, dll_dirs), env['PATH']])
        probe = Path(__file__).resolve().parent / 'sdk_probe'
        probe_build = root / 'probe-build-msvc'
        run([cmake, '-S', probe, '-B', probe_build, '-G', 'Ninja',
             *[f'-D{k}={v}' for k, v in compiler_options.items()],
             f'-DCMAKE_MAKE_PROGRAM={ninja.as_posix()}',
             f'-DCMAKE_PREFIX_PATH={sdk.as_posix()};{prefix}',
             '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL'],
            env=runtime_env, cwd=root, log=log)
        run([cmake, '--build', probe_build], env=runtime_env, cwd=root, log=log)
        run([probe_build / 'msani_sdk_probe.exe'], env=runtime_env, cwd=root, log=log)
        reject_python_dlls(rdkit_dlls + boost_dlls, runtime_env, log)
        notices = sdk / 'share/msani-sdk-licenses'
        notices.mkdir(parents=True, exist_ok=True)
        shutil.copy2(rdkit_src / 'license.txt', notices / 'RDKit.txt')
        copy_boost_licenses(tools, notices)
        json_write(sdk / 'msani-sdk.json', {
            **identity, 'rdkit_version': '2025.9.1', 'rdkit_commit': commit,
            'extra_prefixes': [str(boost)],
            'cmake_options': options, 'validation': 'C++ probe and DLL Python-dependency checks passed',
        })
        print(f'\nSDK validated: {sdk}')
        print('Next, from the MolSanitizer repository:')
        print(f'packaging\\windows\\build_wheels.cmd --sdk "{sdk}" --rdkit-version 2025.9.1 --versions 3.12')
        print('Omit --versions 3.12 after that passes to build the full matrix.')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f'SDK build failed: {error}', file=sys.stderr)
        raise SystemExit(1)
