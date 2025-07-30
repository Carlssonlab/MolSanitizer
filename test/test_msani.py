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

class Test_MolSanitizer(unittest.TestCase):
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
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'],
                                                ['test'], temp_dir)
            cli.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_single_input.txt')
    
    def test_multiple_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt',
                                                 f'{self.path}/in_stereo.txt'],
                                                 ['test'], temp_dir)
            cli.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_multiple_inputs.txt')
    
    def test_removesalts(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_salt.txt'],
                                                ['removesalts', 'test'], temp_dir)
            cli.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_salt_clean.txt')

    def test_tautomers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_tautomers.txt'],
                                                ['tautomers', 'test'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_tautomers.txt')
            
    def test_tautomers_pseudo_chiralities(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_pseudochiral.txt'],
                                                ['tautomers', 'test'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_pseudochiral.txt')

    def test_painsfilter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_pains.txt'],
                                                ['pains', 'test'], temp_dir)
            cli.clean_data(args)
            self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                              f'{self.path}/out_pains.txt')

    def test_unwanted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_unwanted.txt'],
                                                ['test'], temp_dir)
            # Test all filters work together
            args.unwanted = ['Regular','Special','Optional']
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted all:"):
            #self.compareFiles(f'{temp_dir}/dummy_output_clean.txt', f'{self.path}/out_unwanted_all_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_unwanted_all_rejected.txt')
         
            args.unwanted = ['Regular','Optional']
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Regular and Optional:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_unwanted_regular_optional_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_unwanted_regular_optional_rejected.txt')

            # Test if filters works together
            args.unwanted = ['Regular']
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Regular:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_unwanted_regular_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt',
                                  f'{self.path}/out_unwanted_regular_rejected.txt')

            args.unwanted = ['Special']
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Special:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_unwanted_special_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt',
                                  f'{self.path}/out_unwanted_special_rejected.txt')

            args.unwanted = ['Optional']
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Optional:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_unwanted_optional_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt',
                                  f'{self.path}/out_unwanted_optional_rejected.txt')

    def test_descriptor_filter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'],
                                                ['test'], temp_dir)
            # Test all filters work together
            args.ha = '17-25'
            cli.clean_data(args)
            with self.subTest(msg="Checking HA 17-25:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_ha1725_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_ha1725_rejected.txt')
            
            args.ha = '>24'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking HA >24:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_ha_over24_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_ha_over24_rejected.txt')

            args.ha = None
            system(f'rm {temp_dir}/*.txt')
            args.logp = '100-200'
            cli.clean_data(args)
            with self.subTest(msg="Checking logP 100-200:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_logp_100200_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_logp_100200_rejected.txt')

            args.logp = '<=350'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking logP <=350:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_logp_350_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_logp_350_rejected.txt')
            
            args.logp = None
            args.mw = '>=300'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking MW >=300:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_mw_300_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_mw_300_rejected.txt')

            args.mw = None
            args.hba = '<=4'
            args.hbd = '<=2'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking HBA <=4 HBD <=2:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_hba4_hbd2_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_hba4_hbd2_rejected.txt')

            args = self.generate_mock_arguments([f'{self.path}/in_chiral.txt'], ['test'], temp_dir)
            args.hba = None
            args.hbd = None
            args.chiral = '<=2'
            system(f'rm {temp_dir}/*.txt')
            cli.clean_data(args)
            with self.subTest(msg="Checking chiral <=2:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_chiral_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_chiral_rejected.txt')

    def test_create_customfiles(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_unwanted.txt'], 
                                                ['test'], temp_dir)
            
            with self.subTest(msg="Checking creation of customized filter file:"):
                cli.generateCustomTemplate(args, 'filter_out.txt')
                self.assertTrue(Path(f"{args.prefix}.txt").exists(),
                                "Custom unwanted filter file was not created.")
                
            with self.subTest(msg="Checking if customized filter file is applied:"):
                args.custom = f"{args.prefix}.txt"
                cli.clean_data(args)
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_unwanted_all_rejected.txt')
                os.remove(f"{args.prefix}.txt")

            with self.subTest(msg="Checking if customized protonation file is generated:"):
                args.custom = None
                cli.generateCustomTemplate(args, 'ionizations_v3.txt')
                self.assertTrue(Path(f"{args.prefix}.txt").exists(),
                                "Custom protonation file was not created.")
            
            with self.subTest(msg="Checking if customized protonation file is applied:"):
                args.protonation_library = f"{args.prefix}.txt"
                args.protonation = True
                args.pH = 7
                args.input_files = [f'{self.path}/in_protonation.txt']
                cli.clean_data(args)
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/ph7_clean.txt')
                os.remove(f"{args.prefix}.txt")

            with self.subTest(msg="Checking if customized tautomer file is generated:"):
                args.protonation = False
                cli.generateCustomTemplate(args, 'tautomers_v3.txt')
                self.assertTrue(Path(f"{args.prefix}.txt").exists(),
                                "Custom tautomer file was not created.")
                
            with self.subTest(msg="Checking if customized tautomer file is applied:"):
                args.tautomer_library = f"{args.prefix}.txt"
                args.tautomers = True
                args.input_files = [f'{self.path}/in_tautomers.txt']
                cli.clean_data(args)
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_tautomers.txt')
                os.remove(f"{args.prefix}.txt")
            
    def test_stereoisomers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_stereo.txt'],
                                                ['stereoisomers', 'test', ], temp_dir)
            args.max_stereoisomers = 128
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_stereo.txt')


    def test_protonation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_protonation.txt'],
                                                ['protonation', 'test'], temp_dir)
            args.pH = 7
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 7:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/ph7_clean.txt')

            system(f'rm {temp_dir}/*.txt')
            
            args.pH = 5
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 5:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/ph5_clean.txt')
            system(f'rm {temp_dir}/*.txt')

            args.pH = 9
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 9:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/ph9_clean.txt')
            system(f'rm {temp_dir}/*.txt')

    def test_integrity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            from msani.io import parsers
            args = self.generate_mock_arguments([f'{self.path}/in_enamine.txt'],
                                                ['extended','lazy', 'test'], temp_dir)
            args = parsers.Sanitycheck(args)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_enamine_clean.txt')
            self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                  f'{self.path}/out_enamine_rejected.txt')
        os.chdir(self.path)

    def test_db2_generation(self):
        tmp_obj = tempfile.TemporaryDirectory()
        temp_dir = tmp_obj.name
        args = self.generate_mock_arguments([f'{self.path}/in_db2.smi'],
                                            ['protonation', 'gen3d', 'test'], temp_dir)
        args.format = ['db2']        
        args.prefix = Path(temp_dir)

        cli.clean_data(args)
        # If the file was produced
        self.assertTrue(Path(f"{temp_dir}/db2/3,4-diclorophenol.db2").exists(),
                        "DB2 file was not created.")
        
        # If produce 2 conformers
        with open(f"{temp_dir}/db2/3,4-diclorophenol.db2") as db2_file:
            first_line = db2_file.readline()
            self.assertEqual(first_line.split()[7],'2',
                             "DB2 file was not created correctly.")   
            del db2_file
        
        os.chdir(self.path)
        shutil.rmtree(f"{temp_dir}")
        tmp_obj.cleanup()

    def test_db2_sulfonamide(self):
        tmp_obj = tempfile.TemporaryDirectory()
        temp_dir = tmp_obj.name
        args = self.generate_mock_arguments([f'{self.path}/in_sulfonamide.smi'],
                                            ['gen3d', 'test'], temp_dir)
        args.format = ['db2']        
        args.prefix = Path(temp_dir)

        cli.clean_data(args)
        # If the file was produced
        self.assertTrue(Path(f"{temp_dir}/db2/N-Methylbenzenesulfonamide.db2").exists(),
                        "DB2 file was not created.")
        
        # If produce 1 conformer
        n_rigid = 0
        with open(f"{temp_dir}/db2/N-Methylbenzenesulfonamide.db2") as db2_file:
            for line in db2_file:
                if line.startswith('M '): n_rigid += 1
            del db2_file
        n_rigid /= 5 # 5 lines per rigid scaffold
        # Check if the two regioisomers were generated
        self.assertEqual(n_rigid, 2, 'DB2 file was not created correctly.')
        os.chdir(self.path)
        shutil.rmtree(f"{temp_dir}")
        tmp_obj.cleanup()
        
    def test_pdbqt_generation(self):
        #with tempfile.TemporaryDirectory() as temp_dir:
            try:
                import meeko
            except ImportError:
                print("""The Meeko program is not installed.
                      PDBQT options are not tested""")
                return
            tmp_obj = tempfile.TemporaryDirectory()
            temp_dir = tmp_obj.name
            args = self.generate_mock_arguments([f'{self.path}/in_pdbqt.smi'],
                                                ['protonation', 'gen3d', 'test'], temp_dir)
            args.format = ['pdbqt']
            args.prefix = Path(temp_dir)
            cli.clean_data(args)
            self.assertTrue(Path(f"{temp_dir}/pdbqt/salicylic_acid.pdbqt").exists(),
                            "PDBQT file was not created.")

            with open(f"{temp_dir}/pdbqt/salicylic_acid.pdbqt") as pdbqt_file:
                lines = pdbqt_file.readlines()
                del pdbqt_file
            if lines[-1] == '\n': lines.pop(-1)   
            
            self.assertEqual(lines[0],"REMARK SMILES O=C([O-])c1ccccc1O\n",
                             "PDBQT file was not created correctly.")
            
            self.assertEqual(lines[-1],"TORSDOF 2\n",
                             "PDBQT file was not created correctly.")
            
            os.chdir(self.path)
            shutil.rmtree(temp_dir)
            tmp_obj.cleanup()
    
    @unittest.skipIf(OS in ["Windows","Darwin"],
                     "Skipping test on Windows due to incompatible `split` command.")
    def test_batch(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args, parser = parsers.parseArguments([], batch_mode=True)
            applied_flags = ['protonation', 'tautomers', 'gen3d', 'test']

            # Generate arguments
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'],
                                                applied_flags, temp_dir)
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
            required_files = ['submit_msani.sh', 'in0000.smi', 'in0001.smi']
            for file in required_files:
                with self.subTest(file=file):
                    self.assertTrue(os.path.exists(file),
                                    f"{file} was not created.")

            # Validate submit_msani.sh content
            expected_template = [
                '#!/bin/bash\n',
                '#SBATCH -A dummy_output\n',
                '#SBATCH -n 1\n',
                '#SBATCH -J msani_3d\n',
                '#SBATCH -t 96:00:00\n',
                '#SBATCH --mail-type=FAIL\n'
            ]

            with open('submit_msani.sh', 'r') as f:
                file_contents = f.readlines()
                self.assertEqual(file_contents[:6], expected_template,
                                 "submit_msani.sh header is incorrect.")

                # Extract and verify flags
                command_line = file_contents[16].strip()
                extracted_flags = command_line.split('/msani -i $smiles_file ')[-1].split(' --')

                # Ensure applied_flags match extracted_flags
                with self.subTest(msg="Checking applied flags"):
                    self.assertTrue(set(applied_flags).issubset(set(extracted_flags)),
                                    "Flags were not passed correctly.")            
            os.chdir(self.path)

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
            "torsion": None
         } 
        for mode in modes: 
            if (mode not in ['unwanted','custom']): args[mode] = True
        return SimpleNamespace(**args)
 
if __name__ == '__main__':
        unittest.main()
