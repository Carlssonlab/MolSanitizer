import unittest
import tempfile
import os
import shutil

from pathlib import Path
from types import SimpleNamespace
from os import system
import platform


from pandas import read_csv
from msani import cli
from msani.batchmode import Split_Submit_jobs
from msani.io import parsers

OS = platform.system()

'''
This is only a test for reproducibility of MolSanitizer.
It should not be regarded that the expected output on DrugBank is correct.
'''
class Test_DrugBank(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up class-level paths before running tests."""
        cls.path = Path(__file__).parent / "goldenData"
        try:
            os.chdir(cls.path)  # Ensure test runs in the correct directory
        except FileNotFoundError:
            print(f"Warning: Directory {cls.path} not found, using default")

    
    def test_tautomer(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_DB.txt'],
                                                ['test', 'tautomers'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_DB_tauto.txt')
            
    def test_tautomer_protonation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_DB.txt'],
                                                ['test', 'tautomers', 'protonation'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_DB_tauto_prot.txt')
    
    def compare_relative(self, newfile: str, goldenfile: str):
        # Read the files into dataframes
        df1 = read_csv(newfile, header=None, sep=r'\s+')
        df2 = read_csv(goldenfile, header=None, sep=r'\s+')

        # Extract the first column from both dataframes
        column1_df1 = df1.iloc[:, 0]
        column1_df2 = df2.iloc[:, 0]

        # Convert the first column of df1 and df2 into sets
        set1 = set(column1_df1)
        set2 = set(column1_df2)

        # Check if the sets are equal
        self.assertEqual(set1, set2, "Files' contents differ")

    def generate_mock_arguments(self, in_files: list, modes: list, temp_dir: tempfile.TemporaryDirectory):
    
        output_prefix = Path(temp_dir) / 'dummy_output'

        args = {
            'input_files': in_files,
            'format': None,
            'extended': False,
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
            "neutralize": True, 
            "debug": False, 
            "custom":None, 
            "prefix":output_prefix, 
            "max_stereoisomers": 8,
            "numcores": 4,
            "test": False,
            "smiles": None,
            "gen3d": False,
            "format": None,
            "method": "rdkit",
            "mode": "fixed",
            "numconfs": 2000,
            "cleanup": True,
            "randomSeed": 42,
            "energywindow": 25,
            "timeout": 2,
            "tolerance": 30,
            "ignoretorlib":False,
            "timing":False,
            "synthon": False,
            "taurdkit": True,
            "standardize": False,
            "ha": None,
            "logp": None,
            "hba": None,
            "hbd": None,
            "mw": None,
            "chiral": None,
            "nringconfs": 1,
            "allowNonring": False,
            "eps": 1,
            "rigid": None,
            "protlib": None,
            "taulib": None,
            "create_protlib": False,
            "create_taulib": False,
         } 
        for mode in modes: 
            if (mode not in ['unwanted','custom']): args[mode] = True
        return SimpleNamespace(**args)
 
if __name__ == '__main__':
        unittest.main()