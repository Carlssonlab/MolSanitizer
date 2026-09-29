import bz2
import gzip
import io
import lzma
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from rdkit import Chem

from msani import cli
from msani.io.readers import detect_input_separator, validate_input_chunk


class TestInputReaders(unittest.TestCase):
    CXSMILES = 'C[C@H](O)F |&1:1|'

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.path = Path(directory.name) / 'input.smi'

    def read(self, contents, extended=False, synthon=False):
        self.path.write_text(contents)
        with cli.read_input_file(self.path, extended, synthon) as chunks:
            self.assertEqual(chunks.engine, 'c')
            self.assertEqual(chunks.chunksize, 100_000)
            return next(chunks)

    def test_plain_smiles_and_extra_columns(self):
        for separator in (' ', '   ', '\t', ' \t'):
            with self.subTest(separator=separator):
                frame = self.read(f'\nCCO{separator}ethanol\textra\nCC{separator}ethane\textra\n')
                # Whitespace within a tab-delimited field is preserved, as
                # with explicit -e; RDKit accepts the trailing space.
                trailing = ' ' if separator == ' \t' else ''
                self.assertEqual(frame.to_dict('records'), [
                    {'smiles': 'CCO' + trailing, 'ids': 'ethanol'},
                    {'smiles': 'CC' + trailing, 'ids': 'ethane'},
                ])

    def test_cxsmiles_compressed_and_plain_match_explicit_flag(self):
        contents = f'\nCCO\tethanol\n{self.CXSMILES}\tstereo\n'.encode()
        for opener in (open, gzip.open, bz2.open, lzma.open):
            with self.subTest(opener=opener):
                with opener(self.path, 'wb') as stream:
                    stream.write(contents)
                with cli.read_input_file(self.path, False, False) as chunks:
                    auto = next(chunks)
                with cli.read_input_file(self.path, True, False) as chunks:
                    explicit = next(chunks)
                self.assertTrue(auto.equals(explicit))
                self.assertEqual(auto.iloc[1].to_dict(), {'smiles': self.CXSMILES, 'ids': 'stereo'})
                mol = Chem.MolFromSmiles(auto.iloc[1]['smiles'])
                groups = mol.GetStereoGroups()
                self.assertEqual(len(groups), 1)
                self.assertEqual(groups[0].GetGroupType(), Chem.StereoGroupType.STEREO_AND)
                self.assertEqual([atom.GetIdx() for atom in groups[0].GetAtoms()], [1])

    def test_cxsmiles_after_sample_and_chunk_boundary(self):
        self.path.write_text('CCO\tethanol\n' * 100_000 + f'{self.CXSMILES}\tstereo\n')
        with cli.read_input_file(self.path, False, False) as chunks:
            first = next(chunks)
            second = next(chunks)
        self.assertEqual(len(first), 100_000)
        self.assertEqual(second.to_dict('records'), [{'smiles': self.CXSMILES, 'ids': 'stereo'}])

    def test_quoted_cxsmiles_preserves_stereo_and_ids(self):
        for separator in (' ', '   ', '\t'):
            for prefix in ('', f'CCO{separator}ethanol\n'):
                for opener in (open, gzip.open, bz2.open, lzma.open):
                    with self.subTest(separator=separator, prefix=prefix, opener=opener):
                        contents = prefix + f'"{self.CXSMILES}"{separator}stereo\n'
                        with opener(self.path, 'wb') as stream:
                            stream.write(contents.encode())
                        with cli.read_input_file(self.path, False, False) as chunks:
                            frame = next(chunks)
                        validate_input_chunk(frame, self.path)
                        self.assertEqual(frame.iloc[-1].to_dict(),
                                         {'smiles': self.CXSMILES, 'ids': 'stereo'})
                        mol = Chem.MolFromSmiles(frame.iloc[-1]['smiles'])
                        groups = mol.GetStereoGroups()
                        self.assertEqual(len(groups), 1)
                        self.assertEqual(groups[0].GetGroupType(), Chem.StereoGroupType.STEREO_AND)
                        self.assertEqual([atom.GetIdx() for atom in groups[0].GetAtoms()], [1])

    def test_synthon_keeps_three_columns(self):
        frame = self.read('CCO ethanol cap\nCC\tethane\tother\n', synthon=True)
        self.assertEqual(frame.to_dict('records'), [
            {'smiles': 'CCO', 'ids': 'ethanol', 'longname': 'cap'},
            {'smiles': 'CC', 'ids': 'ethane', 'longname': 'other'},
        ])

    def test_explicit_flag_skips_detection(self):
        with patch.object(cli, 'detect_input_separator', side_effect=AssertionError):
            frame = self.read(f'{self.CXSMILES}\tstereo\n', extended=True)
        self.assertEqual(frame.iloc[0]['smiles'], self.CXSMILES)

    def test_detection_read_is_bounded(self):
        stream = io.BytesIO(b'CCO\tethanol\n' * 100_000)
        with patch('msani.io.readers.open') as opener:
            opener.return_value.__enter__.return_value = stream
            self.assertEqual(detect_input_separator(self.path, None), '\t')
        self.assertEqual(stream.tell(), 64 * 1024)

    def test_no_complete_record_in_sample(self):
        self.path.write_text(' ' * (64 * 1024) + 'CCO\tethanol\n')
        with self.assertRaisesRegex(ValueError, 'first 64 KiB'):
            detect_input_separator(self.path, None)

    def test_single_record_without_final_newline(self):
        frame = self.read(f'{self.CXSMILES}\tstereo')
        self.assertEqual(frame.iloc[0]['ids'], 'stereo')

    def test_invalid_delimiters_are_caught_before_processing(self):
        cases = (
            f'CCO ethanol\n{self.CXSMILES} stereo\n',
            f'CCO\tethanol\n{self.CXSMILES} stereo\n',
            'CCO\tethanol\nCC\t\n',
            'CCO\tethanol\nCC\t   \n',
        )
        args = SimpleNamespace(
            input_files=[self.path], standardize=False, extended=False,
            synthon=False, prefix=str(self.path.parent / 'out'), metal_error_file=None,
            gen3d=False, test=True,
        )
        for contents in cases:
            with self.subTest(contents=contents):
                self.path.write_text(contents)
                processor = Mock()
                with self.assertRaisesRegex(ValueError, 'input record 2:'):
                    cli.process_files(processor, args)
                processor.run.assert_not_called()

    def test_validation_reports_record_offset(self):
        frame = self.read(f'CCO ethanol\n{self.CXSMILES} stereo\n')
        with self.assertRaisesRegex(ValueError, 'input record 100002: CXSMILES'):
            validate_input_chunk(frame, self.path, 100_000)


if __name__ == '__main__':
    unittest.main()
