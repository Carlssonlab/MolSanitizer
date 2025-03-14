from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import  SaltRemover, rdMolDescriptors
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

import pandas as pd
import logging

logger = logging.getLogger('eirvs')
RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit

class Filters():
    def __init__(self, 
                 removesalts = False, 
                 ha = None,
                 logp = None,
                 hba = None,
                 hbd = None,
                 mw = None,
                 custom = None,
                 unwanted = None,
                 pains = None,
                 rejectedFile = 'rejected_entries.txt'):
        self.removesalts = removesalts
        self.ha = ha
        self.logp = logp
        self.hba = hba
        self.hbd = hbd
        self.mw = mw
        self.custom = custom
        self.unwanted = unwanted
        self.pains = pains
        self.rejectedFile = rejectedFile

    def __repr__(self):
        cls_name = self.__class__.__name__
        attrs = ', '.join(f'{k}={v!r}' for k, v in self.__dict__.items())
        return f'{cls_name}({attrs})'
    
    @staticmethod
    def remove_invalid_SMILES(df:pd.DataFrame) -> pd.DataFrame:
        """Remove rows with invalid SMILES from the input DataFrame.

        Args:
            df (pd.DataFrame): Input DataFrame with 'smiles' column containing SMILES strings.

        Returns:
            pd.DataFrame: A new DataFrame chunk with valid SMILES strings.
        """
        # Log rows where 'mol' is None before dropping
        invalid_rows = df[df['mol'].isna()]
        for index, row in invalid_rows.iterrows():
            logger.warning(f"INVALID SMILES: Removing row at index {index}: SMILES: {row['smiles']} - ID: {row['ids']}")

        # Remove rows where 'mol' is None
        df_cleaned = df.dropna(subset=['mol'])

        return df_cleaned
    
    @staticmethod
    def remove_exotic_chem_to_db2(df:pd.DataFrame) -> pd.DataFrame:
        '''
        These are substructures that are not supported by either MMFF94(s)
        and Mol2 format so could not be DB2-compatible.
        It mainly includes: hypervalent atoms, halogens with more than 1 bond, atoms with more than 5 bonds,
        and atom types not supported by the Mol2 format.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules that are DB2-compatible.

        '''
        exotic_chems = Chem.MolFromSmarts('[$([*X{5-}]),$([#35!X1]),$([#53!X1]),$([*;!$([#1,#6,#7,#8,#9,#14,#15,#16,#17,#35,#53,#26,#3,#11,#19,#30,#20,#29,#12])])]')
        df_cleaned = df.copy()
        df_cleaned['mol'] = df['mol'].apply(lambda x: x if x.HasSubstructMatch(exotic_chems) == False else None)
        exotic_chem_df = df_cleaned[df_cleaned['mol'].isna()]
        for _, row in exotic_chem_df.iterrows():
            logger.warning(f"Removed entries that are not DB2-compatible: SMILES: {row['smiles']} - ID: {row['ids']}")
        df_cleaned = df_cleaned.dropna(subset=['mol'])
        return df_cleaned

    @staticmethod
    def stripSMILESsalt(mol, molRemover, debug = False):
        """Strip salts from the input molecule using the RDKit SaltRemover class.

        Args:
            mol (rdkit mol object): The input molecule.
            molRemover: The RDKit SaltRemover object.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            rdkit mol object: The stripped molecule.
        """
        res, deleted = molRemover.StripMolWithDeleted(mol)

        if debug and len(deleted) > 0: logger.info(f"Stripped salt {Chem.MolToSmiles(mol)}:  Salts:{' '.join([Chem.MolToSmiles(m) for m in deleted])}")
        
        if len(Chem.rdmolops.GetMolFrags(res)) > 1:
            # If still contains more than one fragment, retains the largest one
            rdMolStandardize.FragmentParentInPlace(res)
        return res

    @staticmethod
    def saltstripping(df: pd.DataFrame, debug = False) -> pd.DataFrame:
        """Remove salts from the input molecules using the RDKit SaltRemover class.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with salt-stripped molecules.
        """
        # Get the absolute path to the template SMARTS file using pathlib
        smartsFile = Path(__file__).parent.parent / 'Data' / 'salt_stripping.txt'
        if debug: logger.info(f'Parsing salts SMARTS file: {smartsFile.resolve()}')

        remover = SaltRemover.SaltRemover(defnFilename=smartsFile)
        filtered_df = df.copy()
        # Remove salts in the list
        filtered_df['mol'] = df['mol'].apply(lambda x: Filters.stripSMILESsalt(x, remover))

        filtered_df['smiles'] = filtered_df['mol'].apply(lambda x: Chem.MolToSmiles(x))
        filtered_df=filtered_df[filtered_df['smiles']!=''] #Remove purely salt molecules
        filtered_df['mol'].apply(lambda x: Chem.SanitizeMol(x))
        return filtered_df
    
    @staticmethod
    def convert_to_query(condition: str, column: str) -> str:
        if "-" in condition:  # Range condition
            lower, upper = map(int, condition.split('-'))
            return f"{column} >= {lower} and {column} <= {upper}"
        elif '>' in condition or '<' in condition:  # Single value condition
            return column + condition
        elif condition.startswith('='):  # Exact value condition
            return f"{column} == {condition.split('=')[1]}"
        elif column != 'logp':  # Exact value condition for heavy atoms
            return f"{column} == {condition}"
        elif column == 'logp':  # Filter out molecules with logP value upto the defined value
            return f"{column} <= {condition}"
        else:
            raise ValueError(f"Invalid condition: {condition}, supported formats: range (e.g. 1-5), greater than (or equal to) (e.g. >= 5), less than or equal to (e.g. <= 5), equal to (e.g. 5).")

    @staticmethod    
    def filter_by_ha(df, filter_query, rejectedFile, debug = False) -> pd.DataFrame:
        """Filter out molecules with heavy atoms only using the RDKit Mol.GetNumHeavyAtoms() function.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        df['ha'] = df['mol'].apply(lambda x: x.GetNumHeavyAtoms())
        query = Filters.convert_to_query(filter_query, 'ha')
        rejected_df = df.query(f'not ({query})').copy()
        rejected_df['ha'] = rejected_df['ha'].apply(lambda x: f'ha{x}')
        rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','ha'], sep = ' ', header=False)
        if debug: logger.info(f"Removed {len(rejected_df)} molecules with heavy atoms requirements: {rejected_df['ha'].values}")
        df.query(query, inplace=True)
        return df

    @staticmethod
    def filter_by_logp(df, filter_query, rejectedFile, debug = False) -> pd.DataFrame:
        """Filter out molecules with a logP value upto the defined value using the RDKit Crippen logP calculation.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (_type_): Path to the file to save rejected molecules.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules containing the specified logP value.
        """
        df['logp'] = df['mol'].apply(lambda x: (MolLogP(x))*100)
        query = Filters.convert_to_query(filter_query, 'logp')
        rejected_df = df.query(f'not ({query})').copy()
        rejected_df['logp'] = rejected_df['logp'].apply(lambda x: f'logp {x/100:.2f}')
        rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','logp'], sep = ' ', header=False)
        if debug: logger.info(f"Removed {len(rejected_df)} molecules with logP requirements: {rejected_df['logp'].values}")
        df.query(query, inplace=True)

        return df
    
    def filter_by_hba(df, filter_query, rejectedFile, debug = False) -> pd.DataFrame:
        """Filter out molecules with required number of H-bond acceptors using the RDKit CalcNumHBA().
        NOTE: It is by intention that the function uses CalcNumHBA() instead of CalcNumLipinskiHBA() was used as we believe that it better represents the chemistry.
        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        df['hba'] = df['mol'].apply(lambda x: rdMolDescriptors.CalcNumHBA(x))
        query = Filters.convert_to_query(filter_query, 'hba')
        rejected_df = df.query(f'not ({query})').copy()
        rejected_df['hba'] = rejected_df['hba'].apply(lambda x: f'hba {x}')
        rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','hba'], sep = ' ', header=False)
        if debug: logger.info(f"Removed {len(rejected_df)} molecules with number of H-bond acceptors requirements: {rejected_df['hba'].values}")
        df.query(query, inplace=True)
        return df
    
    def filter_by_hbd(df, filter_query, rejectedFile, debug = False) -> pd.DataFrame:
        """Filter out molecules with required number of H-bond donors using the RDKit CalcNumLipinskiHBD().
        NOTE: It is by intention that the function uses CalcNumLipinskiHBD() instead of CalcNumHBD() was used as we believe that it better represents the chemistry.
        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        df['hbd'] = df['mol'].apply(lambda x: rdMolDescriptors.CalcNumLipinskiHBD(x))
        query = Filters.convert_to_query(filter_query, 'hbd')
        rejected_df = df.query(f'not ({query})').copy()
        rejected_df['hbd'] = rejected_df['hbd'].apply(lambda x: f'hbd {x}')
        rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','hbd'], sep = ' ', header=False)
        if debug: logger.info(f"Removed {len(rejected_df)} molecules with number of H-bond donors requirements: {rejected_df['hbd'].values}")
        df.query(query, inplace=True)
        return df

    def filter_by_mw(df, filter_query, rejectedFile, debug = False) -> pd.DataFrame:
        """Filter out molecules with required molecular weight using the RDKit GetMolWt().
        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        df['mw'] = df['mol'].apply(lambda x: rdMolDescriptors.CalcExactMolWt(x))
        query = Filters.convert_to_query(filter_query, 'mw')
        rejected_df = df.query(f'not ({query})').copy()
        rejected_df['mw'] = rejected_df['mw'].apply(lambda x: f'mw {x:.2f}')
        rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','mw'], sep = ' ', header=False)
        if debug: logger.info(f"Removed {len(rejected_df)} molecules with molecular weight requirements: {rejected_df['mw'].values}")
        df.query(query, inplace=True)
        return df

    @staticmethod
    def applyStandarizeFilters(mol, params):

        taut_uncharged_parent_clean_mol = None

        try:

            clean_mol = rdMolStandardize.Cleanup(mol, params) 

            # if many fragments, get the "parent"
            parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)

            # try to neutralize molecule
            uncharger = rdMolStandardize.Uncharger()
            uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
            
            # tautomer enumerator
            te = rdMolStandardize.TautomerEnumerator(params) 
            taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)

        except:

            logger.info(f'Molecule NOT processed: {Chem.MolToSmiles(mol)}')

        return taut_uncharged_parent_clean_mol

    @staticmethod
    def standarizeFilters(df: pd.DataFrame) -> pd.DataFrame:


        # follows the steps in
        # https://github.com/greglandrum/RSC_OpenScience_Standardization_202104/blob/main/MolStandardize%20pieces.ipynb
        # https://www.youtube.com/watch?v=eWTApNX8dJQ
        # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule

        #Sc, Y, In, Sn, W, Ac are not included for now (need to check how they can hit other structures).
        organometallics={'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Ga', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',\
        'Cd', 'La', 'Hf ', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Ce', 'Pr', 'Nd',\
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf ', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Ge', 'Sb'}


        # Filter out rows where 'smiles' contains any organometallics
        filtered_df = df[~df['smiles'].apply(lambda x: any(om in x for om in organometallics))]
        if filtered_df.empty:
            return filtered_df
        params = rdMolStandardize.CleanupParameters()
        params.tautomerRemoveSp3Stereo = False
        params.tautomerRemoveBondStereo = False
        params.tautomerRemoveIsotopicHs = False

        filtered_df['mol'] = filtered_df['mol'].apply(lambda x:  Filters.applyStandarizeFilters(x, params))
        filtered_df['smiles'] = filtered_df['mol'].apply(lambda x: Chem.MolToSmiles(x))

        return filtered_df
    @staticmethod
    def detect_and_label_pains(mol, catalog):
        """Detect PAINS functional groups in a molecule.

        Args:
            mol (Rdkit Mol object): The molecule to be checked.
            catalog (FilterCatalog): The PAINS catalog.

        Returns:
            str: Return 'OK' if no PAINS functional groups are detected, otherwise return the description of the first PAINS functional group detected.
        """
        entry = catalog.GetFirstMatch(mol)  # Get the first matching PAINS
        if entry is not None:
            return ('PAINS violation: ' + str(entry.GetDescription().capitalize()))  # Indicate a PAINS violation
        else:
            return 'OK'
        
    @staticmethod
    def painsFilter(df: pd.DataFrame, rejectedFile: str, debug: bool = False) -> pd.DataFrame:
        """
        Detect and filter out molecules with PAINS functional groups using the RDKit PAINS catalog.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (str): Path to the file to save rejected molecules.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules that passed the PAINS filter.
        """
        # Set up the PAINS catalog
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)


        # Create a copy of the DataFrame to avoid modifying the input
        df = df.copy()

        # Detect PAINS violations
        df['reason'] = df['mol'].apply(lambda x: Filters.detect_and_label_pains(x, catalog))

        # Separate rejected molecules
        rejected_df = df[df['reason'] != 'OK']
        if not rejected_df.empty:
            rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles', 'ids', 'reason'], sep=' ', header=False)

        # Filter passed molecules
        passed_df = df[df['reason'] == 'OK'].drop(columns=['reason'])

        if debug:
            print(f"Rejected {len(rejected_df)} molecules due to PAINS violations.")

        return passed_df
    
    @staticmethod
    def check_header(file_path):
        """
        Check if the file has a header containing 'SMARTS' in the first column.
        Args:
            file_path (Path): Path to the CSV file.
            
        Returns:
            bool: True if the file has a header, False otherwise.
        """
        with open(file_path, 'r') as file:
            first_line = file.readline().strip().upper()
            return 'SMARTS' in first_line
        
    @staticmethod
    def loadSMARTSdata(smartsFile: str, unwanted_option=None) -> pd.DataFrame:
        """Load SMARTS patterns from a file and convert them to RDKit molecule objects.

        Args:
            smartsFile (str): Path to the file containing the SMARTS patterns.
            unwanted_option (list): The mode input by the user. Defaults to None.

        Returns:
            pd.DataFrame: A DataFrame containing the SMARTS patterns and their corresponding RDKit molecule objects.
        """
        if unwanted_option is not None: 
            # Using default substructure file 
            smarts_df = pd.read_csv(smartsFile, sep=r'\s+', header=0, names=['smarts','label', 'reason', 'mode', 'ref'])
            smarts_df = smarts_df[smarts_df["mode"].isin(unwanted_option)]
        else:
            # Using customized substructure file
            has_header = Filters.check_header(smartsFile)
            if has_header:
                logger.info(f'Found header in {smartsFile}')
                smarts_df = pd.read_csv(smartsFile, sep=r'\s+', header=0, usecols=[0,1], names=['smarts','label'])
            else:
                smarts_df = pd.read_csv(smartsFile, sep=r'\s+', header=None, usecols=[0,1], names=['smarts','label'])

        smarts_df['mol'] = smarts_df['smarts'].apply(lambda x: Chem.MolFromSmarts(x)) #do we need mergeHs here?
        
        return smarts_df

    @staticmethod
    def filterbysmarts(mol, smarts_df: pd.DataFrame) -> str:
        for _, substructure in smarts_df.iterrows():
            if mol.HasSubstructMatch(substructure.mol):
                return substructure.label
        return 'OK'
    
    @staticmethod
    def unwantedFilter(df: pd.DataFrame, rejectedFile, unwanted_option, debug = False) -> pd.DataFrame:
        """Filter out unwanted substructures using a default list of SMARTS patterns.

            Args:
                df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
                rejectedFile (str): Path to the file to save rejected molecules.
                unwanted_option (list): The mode input by thle user.
                debug (bool, optional): Debug mode. Defaults to False.

            Returns:
                pd.DataFrame: A new DataFrame chunk with molecules that passed the filter.
        """
        # Get the absolute path to the template SMARTS file using pathlib
        smartsFile = Path(__file__).parent.parent / 'Data' / 'filter_out.csv'

        # Load smarts to clean  from file
        unwanted_df = Filters.loadSMARTSdata(smartsFile.resolve(), unwanted_option)
        logger.info(f'Parsed {len(unwanted_df)} substructures from: {smartsFile}')
        # Apply reactions to each SMILES in the DataFrame
        df_clean = df.copy()
        df_clean['reason'] = df['mol'].apply(lambda x: Filters.filterbysmarts(x, unwanted_df))
        rejected_df=df_clean[df_clean['reason']!='OK']
        rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','reason'], sep = ' ', header=False)
        return df_clean[df_clean['reason']=='OK']


    @staticmethod
    def customFilter(df: pd.DataFrame, rejectedFile, smartsFile, debug = False) -> pd.DataFrame:
        """Filter out unwanted substructures using a customized list of SMARTS patterns.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (_type_): Path to the file to save rejected molecules.
            smartsFile (_type_): Path to the file containing the customized list of SMARTS patterns.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules that passed the filter.
        """
        # Load smarts to clean  from file
        unwanted_df = Filters.loadSMARTSdata(smartsFile)

        # Apply reactions to each SMILES in the DataFrame
        df_clean = df.copy()
        df_clean['reason'] = df['mol'].apply(lambda x: Filters.filterbysmarts(x, unwanted_df))
        rejected_df=df_clean[df_clean['reason']!='OK']
        rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','reason'], sep = ' ', header=False)
        return df_clean[df_clean['reason']=='OK']
    
    def filter_df(self, df: pd.DataFrame, rejectedFile: str, debug: bool = False) -> pd.DataFrame:
        """Apply the filters to the input DataFrame.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (str): Path to the file to save rejected molecules.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules that passed all the filters.
        """
        if self.removesalts:
            df = Filters.saltstripping(df, debug)
        if self.ha:
            df = Filters.filter_by_ha(df, self.ha, rejectedFile, debug)
        if self.logp:
            df = Filters.filter_by_logp(df, self.logp, rejectedFile, debug)
        if self.hba:
            df = Filters.filter_by_hba(df, self.hba, rejectedFile, debug)
        if self.hbd:
            df = Filters.filter_by_hbd(df, self.hbd, rejectedFile, debug)
        if self.mw:
            df = Filters.filter_by_mw(df, self.mw, rejectedFile, debug)
        if self.custom:
            df = Filters.customFilter(df, rejectedFile, self.custom, debug)
        if self.unwanted:
            df = Filters.unwantedFilter(df, rejectedFile, self.unwanted, debug)
        if self.pains:
            df = Filters.painsFilter(df, rejectedFile, debug)
        return df