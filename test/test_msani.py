import unittest
import tempfile
import os
import shutil
import tarfile

from pathlib import Path
import platform
from types import SimpleNamespace


from pandas import read_csv
from rdkit import Chem

from msani import cli
from msani.batchmode import Split_Submit_jobs
from msani.conformers import mol2writer
from msani.db2.db2conv import db2converter
from msani.io import parsers

OS = platform.system()
machine = platform.machine().lower()
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
    def test_tautomers_extended(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_extended_tautomers.txt'],
                                                ['tautomers','extended-tautomers', 'test'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_extended_tautomers.txt')

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
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Regular and Optional:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_unwanted_regular_optional_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_unwanted_regular_optional_rejected.txt')

            # Test if filters works together
            args.unwanted = ['Regular']
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Regular:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_unwanted_regular_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt',
                                  f'{self.path}/out_unwanted_regular_rejected.txt')

            args.unwanted = ['Special']
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Special:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_unwanted_special_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt',
                                  f'{self.path}/out_unwanted_special_rejected.txt')

            args.unwanted = ['Optional']
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking unwanted Optional:"):
                self.compareFiles(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_unwanted_optional_clean.txt')
                self.compareFiles(f'{temp_dir}/dummy_output_rejected.txt',
                                  f'{self.path}/out_unwanted_optional_rejected.txt')

    def test_descriptor_filter(self):
        os.chdir(self.path)
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
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking HA >24:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_ha_over24_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_ha_over24_rejected.txt')

            args.ha = None
            self.remove_temp_text_files(temp_dir)
            args.logp = '100-200'
            cli.clean_data(args)
            with self.subTest(msg="Checking logP 100-200:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_logp_100200_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_logp_100200_rejected.txt')

            args.logp = '<=350'
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking logP <=350:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_logp_350_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_logp_350_rejected.txt')
            
            args.logp = None
            args.mw = '>=300'
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking MW >=300:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_mw_300_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_mw_300_rejected.txt')

            args.mw = None
            args.hba = '<=4'
            args.hbd = '<=2'
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking HBA <=4 HBD <=2:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_hba4_hbd2_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_hba4_hbd2_rejected.txt')

            args.hba = None
            args.hbd = None
            args.fsp3 = '0.2-0.7'
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking FSP3 0.2-0.7:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_fsp3_0.2-0.7_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_fsp3_0.2-0.7_rejected.txt')

            args.fsp3 = None
            args.tpsa = '50-100'
            self.remove_temp_text_files(temp_dir)
            cli.clean_data(args)
            with self.subTest(msg="Checking TPSA 50-100:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/out_tpsa50-100_clean.txt')
                self.compare_relative(f'{temp_dir}/dummy_output_rejected.txt',
                                      f'{self.path}/out_tpsa50-100_rejected.txt')

            args = self.generate_mock_arguments([f'{self.path}/in_chiral.txt'], ['test'], temp_dir)
            args.hba = None
            args.hbd = None
            args.chiral = '<=2'
            self.remove_temp_text_files(temp_dir)
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
            args.max_isomers = 128
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_stereo.txt')
    
    def test_stereoisomers_complex(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_stereo_complex.txt'],
                                                ['stereoisomers', 'test'], temp_dir)
            args.stereo_timeout = 1
            args.max_isomers = 128
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_stereo_complex.txt')

    def test_neutralize(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_neutralize.txt'],
                                                ['neutralize', 'test', ], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_neutralize.txt')

    def test_protonation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_protonation.txt'],
                                                ['protonation', 'test'], temp_dir)
            args.pH = 7
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 7:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/ph7_clean.txt')

            self.remove_temp_text_files(temp_dir)
            
            args.pH = 5
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 5:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/ph5_clean.txt')
            self.remove_temp_text_files(temp_dir)

            args.pH = 9
            cli.clean_data(args)
            with self.subTest(msg="Checking pH 9:"):
                self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                      f'{self.path}/ph9_clean.txt')
            self.remove_temp_text_files(temp_dir)

    def test_standardization(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = self.generate_mock_arguments([f'{self.path}/in_standardize.txt'],
                                                ['standardize', 'test'], temp_dir)
            cli.clean_data(args)
            self.compare_relative(f'{temp_dir}/dummy_output_clean.txt',
                                  f'{self.path}/out_standardize.txt')

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

    def test_sdf_generation(self):
        tmp_obj = tempfile.TemporaryDirectory()
        temp_dir = tmp_obj.name
        with self.subTest(msg="Generating SDF file with 1 ring conformation:"):
            args = self.generate_mock_arguments([f'{self.path}/in_confgen.smi'],
                                                ['gen3d', 'test'], temp_dir)
            args.prefix = Path(temp_dir)
            args.format = ['sdf']
            cli.clean_data(args)
            self.assertTrue(Path(f"{temp_dir}/sdf/3,4-diclorophenol.sdf").exists(),
                        "SDF file was not created.")
            # If produce 2 conformers
            with open(f"{temp_dir}/sdf/3,4-diclorophenol.sdf") as sdf_file:
                conf = 0
                for line in sdf_file:
                    if line.strip().endswith('M  END'):
                        conf += 1
                del sdf_file
            self.assertEqual(conf, 2, "SDF file was not created correctly.")

        with self.subTest(msg="Generating SDF file with 2 ring conformations (sulfonamide/cyclohexane):"):
            args = self.generate_mock_arguments([f'{self.path}/in_sulfonamide.smi'],
                                                ['gen3d', 'test'], temp_dir)
            args.prefix = Path(temp_dir)
            args.format = ['sdf']
            cli.clean_data(args)
            self.assertTrue(Path(f"{temp_dir}/sdf/N-Methylbenzenesulfonamide.nr0.sdf").exists(),
                        "SDF file was not created.")
            self.assertTrue(Path(f"{temp_dir}/sdf/N-Methylbenzenesulfonamide.nr1.sdf").exists(),
                        "SDF file was not created.")
            # If sdf file can be read back
            with open(f"{temp_dir}/sdf/N-Methylbenzenesulfonamide.nr0.sdf") as sdf_file:
                lines = ''
                for line in sdf_file: 
                    if line.strip().startswith('M  END'):
                        lines += line
                        break
                    else: lines += line
                mol = Chem.MolFromMolBlock(lines)
                del sdf_file
                if mol is None:
                    raise ValueError("The molecule was not written correctly")
                
            with open(f"{temp_dir}/sdf/N-Methylbenzenesulfonamide.nr1.sdf") as sdf_file:
                lines = ''
                for line in sdf_file: 
                    if line.strip().endswith('M  END'):
                        lines += line
                        break
                    else: lines += line
                mol = Chem.MolFromMolBlock(lines)
                del sdf_file
                if mol is None:
                    raise ValueError("The molecule was not written correctly")
        os.chdir(self.path)
        shutil.rmtree(f"{temp_dir}")
        tmp_obj.cleanup()

    def test_mol2_generation(self):
        tmp_obj = tempfile.TemporaryDirectory()
        temp_dir = tmp_obj.name
        with self.subTest(msg="Generating mol2 file with 1 ring conformation:"):
            args = self.generate_mock_arguments([f'{self.path}/in_confgen.smi'],
                                                ['gen3d', 'test'], temp_dir)
            args.prefix = Path(temp_dir)
            args.format = ['mol2']
            cli.clean_data(args)
            self.assertTrue(Path(f"{temp_dir}/mol2/3,4-diclorophenol.mol2").exists(),
                        "mol2 file was not created.")
            # If the mol2 file can be read
            with open(f"{temp_dir}/mol2/3,4-diclorophenol.mol2") as mol2_file:
                lines = ''
                for line in mol2_file: 
                    if line.strip().startswith('@<TRIPOS>MOLECULE'):
                        if lines == '': lines += line
                        else: break
                    else: lines += line
                mol = Chem.MolFromMol2Block(lines)
                del mol2_file
                if mol is None:
                    raise ValueError("The molecule was not written correctly")
            

        with self.subTest(msg="Generating mol2 file with 2 ring conformations (sulfonamide/cyclohexane):"):
            args = self.generate_mock_arguments([f'{self.path}/in_sulfonamide.smi'],
                                                ['gen3d', 'test'], temp_dir)
            args.prefix = Path(temp_dir)
            args.format = ['mol2']
            cli.clean_data(args)
            self.assertTrue(Path(f"{temp_dir}/mol2/N-Methylbenzenesulfonamide.nr0.mol2").exists(),
                        "mol2 file was not created.")
            self.assertTrue(Path(f"{temp_dir}/mol2/N-Methylbenzenesulfonamide.nr1.mol2").exists(),
                        "mol2 file was not created.")
            # If the mol2 file can be read
            with open(f"{temp_dir}/mol2/N-Methylbenzenesulfonamide.nr0.mol2") as mol2_file:
                lines = ''
                for line in mol2_file: 
                    if line.strip().startswith('@<TRIPOS>MOLECULE'):
                        if lines == '': lines += line
                        else: break
                    else: lines += line
                mol = Chem.MolFromMol2Block(lines)
                if mol is None:
                    raise ValueError("The molecule was not written correctly")
                del mol2_file
                
            with open(f"{temp_dir}/mol2/N-Methylbenzenesulfonamide.nr1.mol2") as mol2_file:
                lines = ''
                for line in mol2_file: 
                    if line.strip().startswith('@<TRIPOS>MOLECULE'):
                        if lines == '': lines += line
                        else: break
                    else: lines += line
                mol = Chem.MolFromMol2Block(lines)
                if mol is None:
                    raise ValueError("The molecule was not written correctly")
                del mol2_file
        os.chdir(self.path)
        shutil.rmtree(f"{temp_dir}")
        tmp_obj.cleanup()

    
    def test_db2_generation(self):
        tmp_obj = tempfile.TemporaryDirectory()
        temp_dir = tmp_obj.name
        args = self.generate_mock_arguments([f'{self.path}/in_confgen.smi'],
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

    def test_db2converter_from_pregenerated_mol2(self):
        """Convert the fixed ZINC conformer archive without embedding."""

        mol2_path = self.path / "ZINCpg000027Y0hd.mol2.gz"
        with tarfile.open(mol2_path, "r:gz") as archive:
            mol2_file = archive.extractfile("ZINCpg000027Y0hd.mol2")
            self.assertIsNotNone(mol2_file)
            mol2_text = mol2_file.read().decode("ascii")

        mol2_blocks = [
            "@<TRIPOS>MOLECULE" + block
            for block in mol2_text.split("@<TRIPOS>MOLECULE")[1:]
        ]
        self.assertEqual(len(mol2_blocks), 1196)

        molecule = Chem.MolFromMol2Block(
            mol2_blocks[0], sanitize=True, removeHs=False
        )
        self.assertIsNotNone(molecule, "Could not read MOL2 topology")
        atom_count = molecule.GetNumAtoms()
        for mol2_block in mol2_blocks[1:]:
            lines = mol2_block.splitlines()
            atom_start = lines.index("@<TRIPOS>ATOM") + 1
            atom_end = lines.index("@<TRIPOS>BOND")
            coordinates = [
                tuple(map(float, line.split()[2:5]))
                for line in lines[atom_start:atom_end]
            ]
            self.assertEqual(len(coordinates), atom_count)
            conformer = Chem.Conformer(atom_count)
            for atom_index, xyz in enumerate(coordinates):
                conformer.SetAtomPosition(atom_index, xyz)
            molecule.AddConformer(conformer, assignId=True)

        topology = mol2writer.Mol2Writer(
            molecule, mol2_template=mol2_blocks[0]
        ).to_db2_topology(
            name="ZINCpg000027Y0hd",
            smiles="C#CCN(CCF)C(=O)[C@@H]1C[C@H](OC)CN1C(=O)[C@@H](C)OCC(C)C",
            longname="fake",
        )

        expected = (self.path / "ZINCpg000027Y0hd.db2").read_text()
        expected_lines = expected.splitlines()
        atom_records = [
            line.split() for line in expected_lines if line.startswith("A ")
        ]
        total_values = list(map(float, expected_lines[1].split()[1:]))
        solvation = SimpleNamespace(
            charge=[float(record[6]) for record in atom_records],
            polarSolv=[float(record[7]) for record in atom_records],
            apolarSolv=[float(record[8]) for record in atom_records],
            solv=[float(record[9]) for record in atom_records],
            surface=[float(record[10]) for record in atom_records],
            totalCharge=total_values[0],
            totalPolarSolv=total_values[1],
            totalApolarSolv=total_values[2],
            totalSolv=total_values[3],
            totalSurface=total_values[4],
        )
        molecule_data = mol2writer.Mol2Writer.with_db2_conformers(
            topology, molecule
        )
        observed_lines = db2converter(molecule_data, solvation).splitlines()

        self.assertEqual(len(observed_lines), len(expected_lines))
        for line_number, (observed, reference) in enumerate(
            zip(observed_lines, expected_lines), 1
        ):
            if observed == reference:
                continue
            # The MOL2 archive stores four decimal places, while the golden
            # DB2 was written from the original higher-precision coordinates.
            # Every non-coordinate record must therefore remain byte-exact;
            # only the last printed decimal of an X coordinate may differ.
            observed_fields = observed.split()
            reference_fields = reference.split()
            self.assertEqual(observed_fields[0], "X", f"DB2 line {line_number}")
            self.assertEqual(
                observed_fields[:4], reference_fields[:4], f"DB2 line {line_number}"
            )
            for observed_xyz, reference_xyz in zip(
                observed_fields[4:], reference_fields[4:]
            ):
                self.assertAlmostEqual(
                    float(observed_xyz),
                    float(reference_xyz),
                    delta=0.00011,
                    msg=f"DB2 line {line_number}",
                )

    def test_amsol(self):
        os.chdir(self.path)
        reference_dir = self.path / 'db2_zinc'

        def read_atom_charges(db2_path):
            charges = []
            with open(db2_path) as db2_file:
                for line in db2_file:
                    columns = line.split()
                    if columns and columns[0] == 'A':
                        charges.append((columns[1], columns[2], float(columns[6])))
            return charges

        with tempfile.TemporaryDirectory() as temp_dir:
            try:
                args = self.generate_mock_arguments(
                    [f'{self.path}/in_db2_zinc.txt'],
                    ['gen3d', 'test', 'no-stereoisomers'],
                    temp_dir,
                )
                args.format = ['db2']
                args.numconfs = 1
                args.prefix = Path(temp_dir)

                cli.clean_data(args)
                self.assertFalse(
                    (Path(temp_dir) / 'solv').exists(),
                    "AMSOLcpp should not create a legacy solv directory.",
                )

                reference_files = sorted(reference_dir.glob('*.db2'))
                self.assertTrue(reference_files, "No AMSOLcpp reference DB2 files found.")

                for reference_file in reference_files:
                    with self.subTest(file=reference_file.name):
                        generated_file = Path(temp_dir) / 'db2' / reference_file.name
                        self.assertTrue(
                            generated_file.exists(),
                            f"DB2 file was not created: {reference_file.name}",
                        )

                        expected = read_atom_charges(reference_file)
                        observed = read_atom_charges(generated_file)
                        self.assertEqual(
                            [(atom_id, atom_name) for atom_id, atom_name, _ in observed],
                            [(atom_id, atom_name) for atom_id, atom_name, _ in expected],
                            f"Atoms differ in {reference_file.name}",
                        )
                        charge_pairs = zip(expected, observed)
                        for expected_atom, observed_atom in charge_pairs:
                            atom_id, atom_name, expected_charge = expected_atom
                            observed_charge = observed_atom[2]
                            self.assertAlmostEqual(
                                observed_charge,
                                expected_charge,
                                delta=0.15,
                                msg=(f"Partial charge differs for atom {atom_id} "
                                     f"({atom_name}) in {reference_file.name}"),
                            )
            finally:
                os.chdir(self.path)
    
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
                from meeko.preparation import MoleculePreparation
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
            applied_flags = ['protonation', 'tautomers', 'gen3d', 'test', 'no-neutralize', 'no-stereoisomers']

            # Generate arguments
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'],
                                                applied_flags, temp_dir)
            args.prefix = Path(temp_dir)
            args.proj_name = 'dummy_output'
            args.partition = None
            args.whole_node = False
            args.whole_node_cores = 1
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
                '\n',
                '#SBATCH -n 1\n',
                '#SBATCH -J msani_3d\n',
                '#SBATCH -t 96:00:00\n',
                '#SBATCH --mail-type=FAIL\n'
            ]

            with open('submit_msani.sh', 'r') as f:
                file_contents = f.readlines()
                self.assertEqual(file_contents[:7], expected_template,
                                 "submit_msani.sh header is incorrect.")

                # Extract and verify flags
                command_line = file_contents[17].strip()
                extracted_flags = command_line.split('/msani -i $smiles_file ')[-1].split(' --')
                # Ensure applied_flags match extracted_flags
                with self.subTest(msg="Checking applied flags"):
                    self.assertTrue(set(applied_flags).issubset(set(extracted_flags)),
                                    "Flags were not passed correctly.")            
            os.chdir(self.path)

    @unittest.skipIf(OS in ["Windows","Darwin"],
                     "Skipping test on Windows due to incompatible `split` command.")
    def test_batch_whole_node(self):
        """Test whole-node mode: header uses tetralith partition + nodes/ntasks,
        and the script body contains the chunked for-loop with & and wait."""
        with tempfile.TemporaryDirectory() as temp_dir:
            args, parser = parsers.parseArguments([], batch_mode=True)
            applied_flags = ['protonation', 'tautomers', 'gen3d', 'test', 'no-neutralize', 'no-stereoisomers']

            # Generate arguments — 100 lines / 50 per job = 2 jobs
            # whole_node_cores=2 → ceil(2/2) = 1 array task
            args = self.generate_mock_arguments([f'{self.path}/in_data100.txt'],
                                                applied_flags, temp_dir)
            args.prefix = Path(temp_dir)
            args.proj_name = 'dummy_output'
            args.partition = 'tetralith'
            args.whole_node = True
            args.whole_node_cores = 2
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
            with open('submit_msani.sh', 'r') as f:
                contents = f.read()
                lines = contents.splitlines(keepends=True)

            # --- Header checks ---
            expected_header_lines = [
                '#!/bin/bash\n',
                '#SBATCH -A dummy_output\n',
                '#SBATCH --partition=tetralith\n',
                '#SBATCH --nodes=1\n',
                '#SBATCH --ntasks=2\n',
                '#SBATCH -J msani_3d\n',
                '#SBATCH -t 96:00:00\n',
                '#SBATCH --mail-type=FAIL\n',
            ]
            with self.subTest(msg="Checking whole-node header"):
                self.assertEqual(lines[:8], expected_header_lines,
                                 "submit_msani.sh whole-node header is incorrect.")

            # --- Script body checks ---
            with self.subTest(msg="Checking chunked for-loop"):
                self.assertIn('START_IDX=$(( TASK_ID * 2 ))', contents,
                              "START_IDX calculation missing or incorrect.")
                self.assertIn('END_IDX=$(( START_IDX + 2 - 1 ))', contents,
                              "END_IDX calculation missing or incorrect.")
                self.assertIn('for i in $(seq $START_IDX $END_IDX); do', contents,
                              "Chunked for-loop missing.")

            with self.subTest(msg="Checking background & operator"):
                self.assertIn(' &', contents,
                              "Background '&' operator missing from script.")

            with self.subTest(msg="Checking wait command"):
                self.assertIn('\nwait\n', contents,
                              "'wait' command missing from script.")

            os.chdir(self.path)

    def clear_temp_txt(self, temp_dir: str):
        for path in Path(temp_dir).glob("*.txt"):
            path.unlink()

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

    @staticmethod
    def remove_temp_text_files(temp_dir):
        for text_file in Path(temp_dir).glob('*.txt'):
            text_file.unlink()

    def compareFiles(self, newfile: str, goldenfile: str):
        with open(goldenfile, 'rb') as goldenFile, open(newfile, 'rb') as newFile:
            
            goldenfile_content = goldenFile.read()
            newFile_content = newFile.read()
            # Assert that the contents are the same
            self.assertEqual(goldenfile_content, newFile_content, "Files' contents differ")

    def generate_mock_arguments(self, in_files: list, modes: list, temp_dir: tempfile.TemporaryDirectory):
        """Build test arguments from the production parser's current defaults."""
        args = parsers.parseArguments([])
        args.input_files = in_files
        args.prefix = Path(temp_dir) / 'dummy_output'
        args.max_isomers = 16

        for mode in modes:
            if mode in {'unwanted', 'custom'}:
                continue
            if mode.startswith('no-'):
                setattr(args, mode.removeprefix('no-').replace('-', '_'), False)
            else:
                setattr(args, mode.replace('-', '_'), True)
        return args
        

if __name__ == '__main__':
        unittest.main()
