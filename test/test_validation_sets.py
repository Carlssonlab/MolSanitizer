import unittest
import tempfile
import os

from pathlib import Path
from types import SimpleNamespace
import platform


from pandas import read_csv
from msani import cli

OS = platform.system()

'''
This is only a test for reproducibility of MolSanitizer.
It should not be regarded that the expected output on DrugBank/TautoBase/Drug-like set is correct.
'''
class Test_ValidationSets(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up class-level paths before running tests."""
        cls.path = Path(__file__).parent / "validationSets"
        try:
            os.chdir(cls.path)  # Ensure test runs in the correct directory
        except FileNotFoundError:
            print(f"Warning: Directory {cls.path} not found, using default")

    
    def test_tau_drugbank(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_DB.txt'],
                                                ['test', 'tautomers'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_DB_tauto.txt',
                              f'{self.path}/in_DB.txt')
            
    def test_tau_prot_drugbank(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_DB.txt'],
                                                ['test', 'tautomers', 'protonation'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_DB_tauto_prot7.txt',
                              f'{self.path}/in_DB.txt')
    
    def test_tau_prot_REAL(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_sample_REAL1M.txt'],
                                                ['test', 'tautomers', 'protonation'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_sample_REAL1M.txt',
                              f'{self.path}/in_sample_REAL1M.txt')
    
    def test_TautoBase(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_TautoBase.txt'],
                                                ['test', 'tautomers'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_TautoBase.txt',
                              f'{self.path}/in_TautoBase.txt')
            
    def test_protonation_monoprotic(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_monoprotic.txt'],
                                                ['test', 'protonation'], temp_dir)
            args.pH_range = 1
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_monoprotic.txt',
                              f'{self.path}/in_monoprotic.txt')
            
    def test_protonation_druglikesets(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_druglikesets.txt'],
                                                ['test', 'protonation'], temp_dir)
            args.pH_range = 1
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_druglikesets.txt', 
                              f'{self.path}/in_druglikesets.txt')    
            
    def compare_relative(self, newfile: str, goldenfile: str, inputfile:str):
        # Read the files into dataframes
        df1 = read_csv(newfile, header=None, sep=r'\s+', names=['smiles', 'ids'], dtype={'smiles': str, 'ids': str})
        df2 = read_csv(goldenfile, header=None, sep=r'\s+', names=['smiles', 'ids'], dtype={'smiles': str, 'ids': str})
        inputdf = read_csv(inputfile, header=None, sep=r'\s+', names=['smiles', 'ids'], dtype={'smiles': str, 'ids': str})
        # Extract the first column from both dataframes
        column1_df1 = df1.loc[:, 'smiles']
        column1_df2 = df2.loc[:, 'smiles']

        # Convert the first column of df1 and df2 into sets
        set1 = set(column1_df1)
        set2 = set(column1_df2)

        diff = set1.difference(set2).union(set2.difference(set1))
        if set1 != set2:
            print(str(len(diff)) + " differences found between the files")
            # Get rows where column 0 is in diff, extract column 1 values
            diff_set1 = set()
            diff_set2 = set()
            
            if not df1.empty and not df1[df1.loc[:, 'smiles'].isin(diff)].empty:
                # Extract column 1 values, split at "_" and take first part
                diff_set1 = set(df1[df1.loc[:, 'smiles'].isin(diff)].loc[:, 'ids'].apply(lambda x: str(x)[:-2] if '_' in x else str(x)))
            
            if not df2.empty and not df2[df2.loc[:, 'smiles'].isin(diff)].empty:
                # Do the same for df2
                diff_set2 = set(df2[df2.loc[:, 'smiles'].isin(diff)].loc[:,  'ids'].apply(lambda x: str(x)[:-2] if '_' in x else str(x)))
            total_diff = diff_set1.union(diff_set2)
            print(f"Difference in files:")
            inputdf['mismatch'] = inputdf.loc[:,  'ids'].apply(lambda x: any(str(x) == y for y in total_diff))
            mismatch_df_input = inputdf[inputdf['mismatch'] == True]
            for _, row in mismatch_df_input.iterrows():
                print(row['smiles'] + " " + row['ids'])
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