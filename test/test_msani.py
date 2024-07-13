import unittest
import tempfile
from pathlib import Path
from types import SimpleNamespace
from os import system

from pandas.testing import assert_frame_equal
import MolSanitizer.parsers as parsers
import MolSanitizer.filters as filters
import MolSanitizer.molSanitizer as molSanitizer

class TestMolSanitizer(unittest.TestCase):
    def test_single_input(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_d2_2col_100.txt'], [], temp_dir)
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_single_input.txt')
    
    def test_multiple_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_d2_2col_100.txt', 'goldenData/in_stereo.txt'], [], temp_dir)
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_multiple_inputs.txt')
    
    def test_removesalts(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_salt.txt'], ['removesalts'], temp_dir)
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_salt_clean.txt')

    def test_tautomers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_tautomers.txt'], ['tautomers'], temp_dir)
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_tautomers.txt')

    def test_painsfilter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_pains.txt'], ['pains'], temp_dir)
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_pains.txt')

    def test_unwanted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_unwanted.txt'], [], temp_dir)
            # Test all filters work together
            args.unwanted = ['Regular','Special','Optional']
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_unwanted_all_clean.txt')
            self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', 'goldenData/out_unwanted_all_rejected.txt')
         
            args.unwanted = ['Regular','Optional']
            system(f'rm {temp_dir}/*.txt')
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_unwanted_regular_optional_clean.txt')
            self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', 'goldenData/out_unwanted_regular_optional_rejected.txt')

            # Test if filters works together
            args.unwanted = ['Regular']
            system(f'rm {temp_dir}/*.txt')
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_unwanted_regular_clean.txt')
            self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', 'goldenData/out_unwanted_regular_rejected.txt')

            args.unwanted = ['Special']
            system(f'rm {temp_dir}/*.txt')
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_unwanted_special_clean.txt')
            self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', 'goldenData/out_unwanted_special_rejected.txt')

            args.unwanted = ['Optional']
            system(f'rm {temp_dir}/*.txt')
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_unwanted_optional_clean.txt')
            self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', 'goldenData/out_unwanted_optional_rejected.txt')
 
    def test_create_customfile(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_unwanted.txt'], ['create_custom'], temp_dir)
            molSanitizer.generateCustomTemplate(args)
            self.assertTrue(Path(f"{args.prefix}.tsv").exists(), "Output file was not created.")

    def test_stereoisomers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_stereo.txt'], ['stereoisomers'], temp_dir)
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_stereo.txt')

    def test_protonation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_protonation.txt'], ['protonation'], temp_dir)
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_protonation.txt')

    def test_integrity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments(['goldenData/in_protonation.txt'], ['protonation'], temp_dir)
            molSanitizer.cleanData(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', 'goldenData/out_protonation.txt')

    def compareFiles(self, newfile: str, goldenfile: str):
    
        with open(goldenfile, 'rb') as goldenFile, open(newfile, 'rb') as newFile:
            
            goldenfile_content = goldenFile.read()
            newFile_content = newFile.read()
            # Assert that the contents are the same
            self.assertEqual(goldenfile_content, newFile_content, "Files' contents differ")

    def generate_mock_arguments(self, in_files: list, modes: list, temp_dir: tempfile.TemporaryDirectory):
    
        output_prefix = Path(temp_dir) / 'dummy_output'

        args = {
            'input_files': in_files,
            'removesalts' : False, 
            'tautomers' : False, 
            'pains': False,
            'unwanted': None,
            'create_custom': False,
            'stereoisomers': False, 
            'protonation': False, 
            "neutralize": False, 
            "debug": False, 
            "custom":None, 
            "prefix":output_prefix, 
            "max_isomers":0
         } 
        for mode in modes: 
            if (mode not in ['unwanted','custom']): args[mode] = True
        return SimpleNamespace(**args)
    
if __name__ == '__main__':
    unittest.main()
