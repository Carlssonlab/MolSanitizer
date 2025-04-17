import unittest
import tempfile
import os
import shutil

from pathlib import Path
from types import SimpleNamespace
from os import system
import platform


import pandas as pd
import EirVS.cli as cli
from EirVS.batchmode import Split_Submit_jobs
from EirVS.io import parsers

OS = platform.system()

class TestEirVS(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up class-level paths before running tests."""
        cls.path = Path(__file__).parent / "goldenData"
        try:
            os.chdir(cls.path)  # Ensure test runs in the correct directory
        except FileNotFoundError:
            print(f"Warning: Directory {cls.path} not found, using default")

    
    def test_single_input(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'], ['test'], temp_dir)
            cli.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_single_input.txt')
    
    def test_multiple_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt', f'{self.path}/in_stereo.txt'], ['test'], temp_dir)
            cli.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_multiple_inputs.txt')
    
    def test_removesalts(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_salt.txt'], ['removesalts', 'test'], temp_dir)
            cli.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_salt_clean.txt')

    def test_tautomers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_tautomers.txt'], ['tautomers', 'test'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_tautomers.txt')

    def test_painsfilter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_pains.txt'], ['pains', 'test'], temp_dir)
            cli.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_pains.txt')

    def test_unwanted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_unwanted.txt'], ['test'], temp_dir)
            # Test all filters work together
            args.unwanted = ['Regular','Special','Optional']
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted all:"):
            #self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_all_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_all_rejected.txt')
         
            args.unwanted = ['Regular','Optional']
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Regular and Optional:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_regular_optional_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_regular_optional_rejected.txt')

            # Test if filters works together
            args.unwanted = ['Regular']
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Regular:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_regular_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_regular_rejected.txt')

            args.unwanted = ['Special']
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Special:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_special_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_special_rejected.txt')

            args.unwanted = ['Optional']
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Optional:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_optional_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_optional_rejected.txt')

    def test_descriptor_filter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'], ['test'], temp_dir)
            # Test all filters work together
            args.ha = '17-25'
            cli.clean_data(args)
            with self.subTest(msg="Checking HA 17-25:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_ha1725_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_ha1725_rejected.txt')
            
            args.ha = '>24'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking HA >24:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_ha_over24_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_ha_over24_rejected.txt')

            args.ha = None
            system(f'rm {temp_dir}/*.txt')
            args.logp = '100-200'
            cli.clean_data(args)
            with self.subTest(msg="Checking logP 100-200:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_logp_100200_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_logp_100200_rejected.txt')

            args.logp = '<=350'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking logP <=350:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_logp_350_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_logp_350_rejected.txt')
            
            args.logp = None
            args.mw = '>=300'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking MW >=300:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_mw_300_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_mw_300_rejected.txt')

            args.mw = None
            args.hba = '<=4'
            args.hbd = '<=2'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking HBA <=4 HBD <=2:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_hba4_hbd2_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_hba4_hbd2_rejected.txt')

            args = self.generate_mock_arguments([f'{self.path}/in_chiral.txt'], ['test'], temp_dir)
            args.hba = None
            args.hbd = None
            args.chiral = '<=2'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking chiral <=2:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_chiral_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_chiral_rejected.txt')


    def test_create_filter_customfile(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_unwanted.txt'], ['create_custom', 'test'], temp_dir)
            cli.generateCustomTemplate(args)
            with self.subTest(msg="Checking creation of customized filter file:"):
                self.assertTrue(Path(f"{args.prefix}.txt").exists(), "Output file was not created.")
            args.create_custom = False
            args.custom = f"{args.prefix}.txt"
            cli.clean_data(args)
            with self.subTest(msg="Checking if customized filter file is applied:"):
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_unwanted_all_rejected.txt')


    def test_stereoisomers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_stereo.txt'], ['stereoisomers', 'test', ], temp_dir)
            args.max_stereoisomers = 128
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_stereo.txt')


    def test_protonation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_protonation.txt'], ['protonation', 'test'], temp_dir)
            args.pH = 7
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 7:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/ph7_clean.txt')

            system(f'rm {temp_dir}/*.txt')
            
            args.pH = 5
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 5:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/ph5_clean.txt')
            system(f'rm {temp_dir}/*.txt')

            args.pH = 9
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 9:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/ph9_clean.txt')
            system(f'rm {temp_dir}/*.txt')

    def test_integrity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            from EirVS.io import parsers
            args = self.generate_mock_arguments([f'{self.path}/in_enamine.txt'], ['enamine','lazy', 'test'], temp_dir)
            args = parsers.Sanitycheck(args)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_enamine_clean.txt')
            self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt', f'{self.path}/out_enamine_rejected.txt')
        os.chdir(self.path)

    def test_db2_generation(self):
        tmp_obj = tempfile.TemporaryDirectory()
        temp_dir = tmp_obj.name
        args = self.generate_mock_arguments([f'{self.path}/in_db2.smi'], ['db2', 'test', 'enrichment'], temp_dir)
        args.prefix = Path(temp_dir)

        cli.clean_data(args)
        # If the file was produced
        self.assertTrue(Path(f"{temp_dir}/db2/3,4-diclorophenol.db2").exists(), "DB2 file was not created.")

        # If produce 2 conformers
        with open(f"{temp_dir}/db2/3,4-diclorophenol.db2") as db2_file:
            first_line = db2_file.readline()
            self.assertEqual(first_line.split()[7],'2', "DB2 file was not created correctly.")   
            del db2_file
        
        os.chdir(self.path)
        shutil.rmtree(f"{temp_dir}")
        tmp_obj.cleanup()

    def test_pdbqt_generation(self):
        #with tempfile.TemporaryDirectory() as temp_dir:
            try:
                import meeko
            except ImportError:
                print("""The Meeko program is not installed. PDBQT options are not tested""")
                return
            tmp_obj = tempfile.TemporaryDirectory()
            temp_dir = tmp_obj.name
            args = self.generate_mock_arguments([f'{self.path}/in_pdbqt.smi'], ['protonation', 'pdbqt', 'test'], temp_dir)
            args.prefix = Path(temp_dir)
            cli.clean_data(args)
            self.assertTrue(Path(f"{temp_dir}/pdbqt/salicylic_acid.pdbqt").exists(), "PDBQT file was not created.")

            with open(f"{temp_dir}/pdbqt/salicylic_acid.pdbqt") as pdbqt_file:
                lines = pdbqt_file.readlines()
                del pdbqt_file
            if lines[-1] == '\n': lines.pop(-1)
            #shutil.rmtree(f"{temp_dir}/pdbqt")
            
            
            self.assertEqual(lines[0],"REMARK SMILES O=C([O-])c1ccccc1O\n", "PDBQT file was not created correctly.")
            self.assertEqual(lines[-1],"TORSDOF 2\n", "PDBQT file was not created correctly.")
            
            os.chdir(self.path)
            shutil.rmtree(temp_dir)
            tmp_obj.cleanup()
    
    @unittest.skipIf(OS in ["Windows","Darwin"], "Skipping test on Windows due to incompatible `split` command.")
    def test_batch(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args, parser = parsers.parseArguments([], batch_mode=True)
            applied_flags = ['protonation', 'tautomers', 'db2', 'test']

            # Generate arguments
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'], applied_flags, temp_dir)
            args.prefix = Path(temp_dir)
            args.proj_name = 'dummy_output'
            args.max_jobs = 2
            args.lines = 50
            args.timelimit = 96

            # Run batch submission
            Split_Submit_jobs(args, parser)

            batch_dir = Path(temp_dir) / 'in_data100'
            os.chdir(batch_dir)

            # Check for required files
            required_files = ['submit_eirvs.sh', 'in0000.smi', 'in0001.smi']
            for file in required_files:
                with self.subTest(file=file):
                    self.assertTrue(os.path.exists(file), f"{file} was not created.")

            # Validate submit_eirvs.sh content
            expected_template = [
                '#!/bin/bash\n',
                '#SBATCH -A dummy_output\n',
                '#SBATCH -n 1\n',
                '#SBATCH -J eirvs_3d\n',
                '#SBATCH -t 96:00:00\n',
                '#SBATCH --mail-type=FAIL\n'
            ]

            with open('submit_eirvs.sh', 'r') as f:
                file_contents = f.readlines()
                self.assertEqual(file_contents[:6], expected_template, "submit_eirvs.sh header is incorrect.")

                # Extract and verify flags
                command_line = file_contents[16].strip()
                extracted_flags = command_line.split('/eirvs -i $smiles_file ')[-1].split(' --')

                # Ensure applied_flags match extracted_flags
                with self.subTest(msg="Checking applied flags"):
                    self.assertTrue(set(applied_flags).issubset(set(extracted_flags)), "Flags were not passed correctly.")            
            os.chdir(self.path)

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
            "neutralize": True, 
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
            "taurdkit": True,
            "standardize": False,
            "ha": None,
            "logp": None,
            "hba": None,
            "hbd": None,
            "mw": None,
            "chiral": None,
            "pdbqt": False,
            "nringconfs": 1,
            "rigid": None
         } 
        for mode in modes: 
            if (mode not in ['unwanted','custom']): args[mode] = True
        return SimpleNamespace(**args)
 

    
       
if __name__ == '__main__':
        unittest.main()
