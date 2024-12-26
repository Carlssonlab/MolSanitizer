import unittest
import tempfile
import os

from pathlib import Path
from types import SimpleNamespace
from os import system
import shutil

import pandas as pd
import MolSanitizer.molSanitizer as molSanitizer

class TestMolSanitizer(unittest.TestCase):
    if os.path.basename(os.getcwd()) == 'test':
        path = os.path.join(os.getcwd(), 'goldenData')
    else:
        path = os.path.join(os.getcwd(), 'test', 'goldenData')
    
    def test_single_input(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'], ['test'], temp_dir)
            molSanitizer.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_single_input.txt')
    
    def test_multiple_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt', f'{self.path}/in_stereo.txt'], ['test'], temp_dir)
            molSanitizer.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_multiple_inputs.txt')
    
    def test_removesalts(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_salt.txt'], ['removesalts', 'test'], temp_dir)
            molSanitizer.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_salt_clean.txt')

    def test_tautomers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_tautomers.txt'], ['tautomers', 'test'], temp_dir)
            molSanitizer.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_tautomers.txt')

    def test_painsfilter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_pains.txt'], ['pains', 'test'], temp_dir)
            molSanitizer.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_pains.txt')

    def test_unwanted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_unwanted.txt'], ['test'], temp_dir)
            # Test all filters work together
            args.unwanted = ['Regular','Special','Optional']
            molSanitizer.clean_data(args)
            os.system(f'ls {temp_dir}')
            #self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_all_clean.txt')
            self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_all_rejected.txt')
         
            args.unwanted = ['Regular','Optional']
            system(f'rm {temp_dir}/*.txt')
            molSanitizer.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_regular_optional_clean.txt')
            self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_regular_optional_rejected.txt')

            # Test if filters works together
            args.unwanted = ['Regular']
            system(f'rm {temp_dir}/*.txt')
            molSanitizer.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_regular_clean.txt')
            self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_regular_rejected.txt')

            args.unwanted = ['Special']
            system(f'rm {temp_dir}/*.txt')
            molSanitizer.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_special_clean.txt')
            self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_special_rejected.txt')

            args.unwanted = ['Optional']
            system(f'rm {temp_dir}/*.txt')
            molSanitizer.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_optional_clean.txt')
            self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_optional_rejected.txt')
 
    def test_create_customfile(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_unwanted.txt'], ['create_custom', 'test'], temp_dir)
            molSanitizer.generateCustomTemplate(args)
            self.assertTrue(Path(f"{args.prefix}.tsv").exists(), "Output file was not created.")

    def test_stereoisomers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_stereo.txt'], ['stereoisomers', 'test', ], temp_dir)
            args.max_stereoisomers = 128
            molSanitizer.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_stereo.txt')

    def test_protonation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_protonation.txt'], ['protonation', 'test'], temp_dir)
            args.pH = 7
            molSanitizer.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/ph7_clean.txt')
            system(f'rm {temp_dir}/*.txt')
            
            args.pH = 5
            molSanitizer.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/ph5_clean.txt')
            system(f'rm {temp_dir}/*.txt')

            args.pH = 9
            molSanitizer.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/ph9_clean.txt')
            system(f'rm {temp_dir}/*.txt')
    
    def test_integrity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_enamine.txt'], ['enamine','lazy', 'test'], temp_dir)
            args = molSanitizer.Sanitycheck(args)
            molSanitizer.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_enamine_clean.txt')
            self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_enamine_rejected.txt')
    
    def test_db2_generation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_db2.smi'], ['db2', 'test', 'enrichment'], temp_dir)
            args.prefix = Path(temp_dir)

            # If AMSOL is installed
            self.assertTrue((Path(__file__).parent.parent / 'MolSanitizer' / 'amsol' / 'amsol7.1').exists(), "amsol7.1 file does not exist. Check MolSanitizer/amsol directory for AMSOL installation")
            molSanitizer.clean_data(args)
            # If the file was produced
            self.assertTrue(Path(f"{temp_dir}/db2/3,4-diclorophenol.db2").exists(), "DB2 file was not created.")

            # If produce 2 conformers
            with open(f"{temp_dir}/db2/3,4-diclorophenol.db2") as db2_file:
                first_line = db2_file.readline()
                self.assertEqual(first_line.split()[7],'2', "DB2 file was not created correctly.")   


    def compare_relative(self, newfile: str, goldenfile: str):
        # Read the files into dataframes
        df1 = pd.read_csv(newfile, header=None, sep=r'\s+')
        df2 = pd.read_csv(goldenfile, header=None, sep=r'\s+')

        # Extract the first column from both dataframes
        column1_df1 = df1.iloc[:, 0]
        column1_df2 = df2.iloc[:, 0]

        # Convert the first column of df1 and df2 into sets
        set1 = set(column1_df1)
        set2 = set(column1_df2)

        # Check if the sets are equal
        self.assertEqual(set1, set2, "Files' contents differ")

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
            'enamine': False,
            'lazy': False,
            'removesalts' : False, 
            'tautomers' : False, 
            'pains': False,
            'unwanted': None,
            'create_custom': False,
            'stereoisomers': False, 
            'protonation': False,
            'pH': 7,
            'pH_range': 0,
            "neutralize": False, 
            "debug": False, 
            "custom":None, 
            "prefix":output_prefix, 
            "max_stereoisomers": 8,
            "numcores": 4,
            "test": False,
            "smiles": None,
            "db2": False,
            "numconfs": 2000,
            "cleanup": True,
            "randomSeed": 42,
            "enrichment": False,
            "energywindow": 25,
            "timeout": 2,
            "ignoretorlib":False,
            "timing":False,
            'corina': False,
            "synthon": False,
            "taurdkit": True
         } 
        for mode in modes: 
            if (mode not in ['unwanted','custom']): args[mode] = True
        return SimpleNamespace(**args)
 

    
       
if __name__ == '__main__':
        unittest.main()
