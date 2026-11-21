"""Timeout handling for the RDKit embedding worker and the Open Babel fallback."""
from dataclasses import replace
from multiprocessing import Pipe
from queue import Queue
from types import SimpleNamespace
from unittest.mock import patch
import os
import time
import unittest

from msani.conformers import conformers
from msani.io import parsers


def _config(**overrides):
    args = parsers.parseArguments(['-s', 'CCO', '--gen3d', '--format', 'sdf'])
    return replace(conformers.ConformerRunConfig.from_args(args), **overrides)


class TestEmbeddingTimeouts(unittest.TestCase):
    def test_reader_forwards_start_message_then_result(self):
        receiver, sender = Pipe(duplex=False)
        sender.send(conformers._WORKER_STARTED)
        sender.send(('result',))
        sender.close()
        messages = Queue()
        conformers._receive_embedding_result(receiver, messages)
        self.assertEqual(messages.get_nowait(), (True, conformers._WORKER_STARTED))
        self.assertEqual(messages.get_nowait(), (True, ('result',)))
        self.assertTrue(messages.empty())

    def test_reader_reports_worker_exit_without_result(self):
        receiver, sender = Pipe(duplex=False)
        sender.send(conformers._WORKER_STARTED)
        sender.close()
        messages = Queue()
        conformers._receive_embedding_result(receiver, messages)
        messages.get_nowait()
        received, _ = messages.get_nowait()
        self.assertFalse(received)

    def test_obabel_timeout_option_reaches_run_config(self):
        self.assertEqual(_config().obabel_timeout, 30)
        args = parsers.parseArguments(['-s', 'CCO', '--gen3d', '-tob', '12.5'])
        self.assertEqual(conformers.ConformerRunConfig.from_args(args).obabel_timeout, 12.5)

    def test_embedder_is_resolved_once_into_run_config(self):
        # The worker uses the embedder its parent resolved, not its own environment.
        self.assertEqual(_config().embedder, conformers._DEFAULT_EMBEDDER)
        module = SimpleNamespace(EMBEDDING_PRESETS=('sretkdgv3', 'etkdgv4-ic'))
        with patch.object(conformers, 'cpp_sampler', module, create=True), \
                patch.object(conformers, 'CPP_AVAILABLE', True):
            with patch.dict(os.environ, {'MSANI_EMBEDDER': 'ETKDGv4-IC'}):
                config = _config()
            self.assertEqual(config.embedder, 'etkdgv4-ic')
            # An explicit name wins over the environment of the calling process.
            self.assertEqual(conformers._modern_embedder(config.embedder), 'etkdgv4-ic')
            with patch.dict(os.environ, {'MSANI_EMBEDDER': 'etkdgv4-ic'}):
                self.assertIsNone(conformers._modern_embedder('sretkdgv3'))

    def test_embedding_worker_returns_result(self):
        amsol, ring_confs, _, error = conformers._run_embedding_worker('CCO', 'ethanol', _config())
        self.assertIsNone(error)
        self.assertIsNotNone(amsol)
        self.assertEqual(len(ring_confs), 1)

    def test_embedding_timeout_raises_and_stops_worker(self):
        # A near-zero timeout expires during embedding, not during start-up.
        start = time.monotonic()
        with self.assertRaisesRegex(TimeoutError, 'exceeded its deadline'):
            conformers._run_embedding_worker(
                'CC(=O)N1CC(C(=O)N2CCC[C@](C)(O)C2)C1', 'slow', _config(embedding_timeout=1e-6))
        self.assertLess(time.monotonic() - start, conformers._EMBED_WORKER_STARTUP_TIMEOUT)


if __name__ == '__main__':
    unittest.main()
