"""Shared multiprocessing context for MolSanitizer worker processes."""

import multiprocessing as mp
import platform

# forkserver avoids inheriting parent locks (RDKit allocator, logging) on Linux.
# Windows only supports 'spawn'; macOS prefers 'spawn' too (fork is deprecated).
START_METHOD = 'forkserver' if platform.system() == 'Linux' else 'spawn'

# Without a preload, every forkserver child re-imports MolSanitizer, RDKit,
# SciPy and pandas and rebuilds the torsion library (~0.4 s per worker).
# Importing them once in the clean, single-purpose server process makes each
# child a cheap fork of that server instead.
_FORKSERVER_PRELOAD = ['msani.cli', 'msani.conformers.conformers']


def get_context():
    """Return the context for worker processes, preloading the forkserver."""
    if START_METHOD == 'forkserver':
        # The server is shared and started once; this only takes effect if it
        # is set before the first forkserver worker starts.
        mp.set_forkserver_preload(_FORKSERVER_PRELOAD)
    return mp.get_context(START_METHOD)
