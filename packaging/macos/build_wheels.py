"""Build, delocate, and clean-install test a native macOS wheel matrix.

Supports arm64 and x86_64. Builds are always native: the target architecture
must match the host, because every wheel is import-tested and run against the
repository suite before it is accepted.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile


VERSIONS = ['3.10', '3.11', '3.12', '3.13', '3.14']
# conda subdir and delocate/CMake spellings for each supported host.
ARCHITECTURES = {'arm64': 'osx-arm64', 'x86_64': 'osx-64'}


def run(command, *, env, cwd, log):
    command = list(map(str, command))
    log.write('\n> ' + shlex.join(command) + '\n')
    log.flush()
    with subprocess.Popen(command, env=env, cwd=cwd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True,
                          encoding='utf-8', errors='replace') as process:
        for line in process.stdout:
            print(line, end='', flush=True)
            log.write(line)
        code = process.wait()
    log.flush()
    if code:
        raise subprocess.CalledProcessError(code, command)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sdk', required=True, type=Path, help='Native arm64 RDKit C++ SDK prefix')
    parser.add_argument('--rdkit-version', required=True, help='Matching PyPI RDKit version')
    parser.add_argument('--versions', nargs='+', choices=VERSIONS, default=VERSIONS)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--arch', choices=sorted(ARCHITECTURES), default=platform.machine(),
                        help='Target architecture; must match the host.')
    parser.add_argument('--deployment-target', default='11.0')
    parser.add_argument('--uv', default='uv', help='uv executable (manages standalone CPython)')
    args = parser.parse_args()
    if sys.platform != 'darwin':
        parser.error('Run on macOS.')
    if platform.machine() != args.arch:
        parser.error(f'--arch {args.arch} needs a native {args.arch} host and Python, '
                     f'not {platform.machine()} (check for Rosetta).')
    if args.jobs < 1:
        parser.error('--jobs must be positive.')
    try:
        target = tuple(map(int, args.deployment_target.split('.')))
        if len(target) != 2 or target < (11, 0):
            raise ValueError
    except ValueError:
        parser.error('--deployment-target must be MAJOR.MINOR, at least 11.0.')
    uv = shutil.which(args.uv)
    if not uv:
        parser.error('Install uv, or pass --uv /path/to/uv.')
    sdk = args.sdk.resolve()
    if not (sdk / 'lib/cmake/rdkit/rdkit-config.cmake').is_file():
        parser.error('SDK is missing lib/cmake/rdkit/rdkit-config.cmake.')
    for metadata in (sdk / 'conda-meta').glob('librdkit-*.json'):
        record = json.loads(metadata.read_text())
        if record['name'] == 'librdkit':
            if (tuple(map(int, record['version'].split('.'))) !=
                    tuple(map(int, args.rdkit_version.split('.')))):
                parser.error('--rdkit-version does not match the SDK.')
            if record.get('subdir') != ARCHITECTURES[args.arch]:
                parser.error(f'SDK must be built for {ARCHITECTURES[args.arch]}.')
    subprocess.run(['xcrun', '--find', 'clang++'], check=True)
    repo = Path(__file__).resolve().parents[2]
    cache = repo / f'build/macos-{args.arch}'
    cache.mkdir(parents=True, exist_ok=True)
    (repo / 'dist').mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix=f'macos-{args.arch}-wheels-', dir=repo / 'dist'))
    wheels, logs = output / 'wheels', output / 'logs'
    wheels.mkdir()
    logs.mkdir()
    print(f'Output: {output}', flush=True)
    base = os.environ.copy()
    for key in list(base):
        if key.startswith(('CONDA', 'PYTHON', 'CMAKE', 'PIP_', 'UV_', 'DYLD_')) or key in (
                'RDBASE', 'VIRTUAL_ENV', 'ARCHFLAGS', 'MACOSX_DEPLOYMENT_TARGET',
                'CC', 'CXX', 'CFLAGS', 'CXXFLAGS', 'CPPFLAGS', 'LDFLAGS',
                'CPATH', 'LIBRARY_PATH', 'SDKROOT'):
            base.pop(key, None)
    base.update(PATH='/usr/bin:/bin:/usr/sbin:/sbin', PYTHONNOUSERSITE='1',
                UV_CACHE_DIR=str(cache / 'uv-cache'),
                UV_PYTHON_INSTALL_DIR=str(cache / 'python'),
                UV_PYTHON_PREFERENCE='only-managed',
                PIP_CACHE_DIR=str(cache / 'pip-cache'))
    (output / 'build-settings.json').write_text(json.dumps({
        'sdk': str(sdk), 'rdkit_version': args.rdkit_version, 'arch': args.arch,
        'versions': args.versions, 'deployment_target': args.deployment_target,
        'jobs': args.jobs,
    }, indent=2) + '\n')
    results = []
    for version in dict.fromkeys(args.versions):
        tag = 'cp' + version.replace('.', '')
        phase = 'BUILD'
        try:
            with (logs / f'{tag}.log').open('w') as log, tempfile.TemporaryDirectory(prefix=f'msani-{tag}-') as tmp:
                work = Path(tmp)
                def create_env(name):
                    prefix = work / name
                    run([uv, 'venv', '--seed', '--python', version, prefix], env=base, cwd=work, log=log)
                    python = prefix / 'bin/python'
                    run([python, '-c', 'import platform,sys,sysconfig; '
                         f"assert platform.machine() == '{args.arch}'; "
                         "assert platform.python_implementation() == 'CPython'; "
                         f'assert sys.version_info[:2] == {tuple(map(int, version.split(".")))}; '
                         "assert not sysconfig.get_config_var('Py_GIL_DISABLED'); print(sys.version)"],
                        env=base, cwd=work, log=log)
                    return python
                python = create_env('build-env')
                env = base.copy()
                env.update(PATH=str(python.parent) + ':' + base['PATH'], RDBASE=str(sdk),
                           CMAKE_GENERATOR='Ninja', CMAKE_BUILD_PARALLEL_LEVEL=str(args.jobs),
                           MACOSX_DEPLOYMENT_TARGET=args.deployment_target,
                           ARCHFLAGS=f'-arch {args.arch}',
                           CC='/usr/bin/clang', CXX='/usr/bin/clang++')
                run([python, '-m', 'pip', 'install', 'scikit-build-core', 'pybind11', 'cmake',
                     'ninja', 'delocate', 'twine', f'rdkit=={args.rdkit_version}'], env=env, cwd=work, log=log)
                run([python, '-m', 'pip', 'freeze'], env=env, cwd=work, log=log)
                raw = work / 'raw'
                run([python, '-m', 'pip', 'wheel', repo, '--no-build-isolation', '--no-deps',
                     f'-Cbuild-dir={work / "native"}', '-Ccmake.define.MSANI_RDKIT_LINKAGE=shared',
                     f'-Ccmake.define.CMAKE_PREFIX_PATH={sdk}',
                     f'-Ccmake.define.CMAKE_OSX_ARCHITECTURES={args.arch}',
                     f'-Ccmake.define.CMAKE_OSX_DEPLOYMENT_TARGET={args.deployment_target}',
                     f'-Ccmake.define.Python3_EXECUTABLE={python}',
                     f'-Ccmake.define.Python_EXECUTABLE={python}', '--wheel-dir', raw],
                    env=env, cwd=work, log=log)
                candidates = list(raw.glob(f'*-{tag}-{tag}-macosx_*_{args.arch}.whl'))
                if len(candidates) != 1:
                    raise RuntimeError(f'Expected one {args.arch} {tag} wheel, found {candidates}')
                # CMake deliberately removes build-machine RPATHs. Give only the
                # repair step a library search path, then remove it for tests.
                repair_env = env.copy()
                repair_env['DYLD_LIBRARY_PATH'] = str(sdk / 'lib')
                run([python.parent / 'delocate-wheel', '--require-archs', args.arch,
                     '--sanitize-rpaths', '-v', '-w', wheels, candidates[0]],
                    env=repair_env, cwd=work, log=log)
                repaired = list(wheels.glob(f'*-{tag}-{tag}-macosx_*_{args.arch}.whl'))
                if len(repaired) != 1:
                    raise RuntimeError(f'Expected one repaired wheel, found {repaired}')
                wheel = repaired[0]  # delocate may raise the minimum macOS tag.
                run([python.parent / 'delocate-listdeps', wheel], env=base, cwd=work, log=log)
                run([python, repo / 'packaging/macos/check_wheel.py', wheel], env=base, cwd=work, log=log)
                run([python, '-m', 'twine', 'check', '--strict', wheel], env=base, cwd=work, log=log)
                phase = 'TEST'
                clean_python = create_env('test-env')
                run([clean_python, '-m', 'pip', 'install', wheel], env=base, cwd=work, log=log)
                smoke = (
                    'from pathlib import Path; import msani, msani_confgen_cpp, msani_stereoisomers, amsolcpp; '
                    'from amsolcpp import _amsolcpp; '
                    'modules=(msani,msani_confgen_cpp,msani_stereoisomers,amsolcpp,_amsolcpp); '
                    'assert all("site-packages" in Path(m.__file__).parts for m in modules); '
                    'from rdkit import rdBase; print("Runtime RDKit:",rdBase.rdkitVersion); '
                    'print("Native imports passed")'
                )
                run([clean_python, '-c', smoke], env=base, cwd=work, log=log)
                run([clean_python, '-m', 'pip', 'install', str(wheel) + '[pdbqt]'], env=base, cwd=work, log=log)
                run([clean_python, '-m', 'pip', 'check'], env=base, cwd=work, log=log)
                run([clean_python, '-m', 'pip', 'freeze'], env=base, cwd=work, log=log)
                tests = work / 'tests'
                shutil.copytree(repo / 'test', tests, ignore=shutil.ignore_patterns('__pycache__'))
                run([clean_python, '-m', 'unittest', 'discover', '-s', tests, '-p', 'test*.py', '-v'],
                    env=base, cwd=work, log=log)
            results.append(f'{version} PASSED')
        except Exception as error:
            results.append(f'{version} {phase}_FAILED: {error}')
            print(results[-1], file=sys.stderr)
        (output / 'status.txt').write_text('\n'.join(results) + '\n')
    hashes = [f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}' for p in sorted(wheels.glob('*.whl'))]
    (output / 'SHA256SUMS').write_text('\n'.join(hashes) + '\n')
    print('\n'.join(results))
    print(f'Results: {output}')
    return 0 if all(r.endswith(' PASSED') for r in results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
