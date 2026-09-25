from pathlib import Path
import re

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import FilterCatalog as RDFilterCatalog
from rdkit.Chem import SaltRemover, rdMolDescriptors
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.MolStandardize import rdMolStandardize

from pandas import DataFrame, Series, read_csv, concat
import logging

logger = logging.getLogger('msani')
uncharger = rdMolStandardize.Uncharger()
saltfile = Path(__file__).parent.parent / 'Data' / 'salt_stripping.txt'
salt_remover = SaltRemover.SaltRemover(defnFilename=saltfile)        
_pains_catalog = None

RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit

_NUMBER_PATTERN = r'-?(?:\d+(?:\.\d*)?|\.\d+)'


def _format_scaled_logp(value: str) -> str:
    """Return a logP value on the legacy cLogP*100 comparison scale."""
    numeric_value = float(value)
    scaled_value = numeric_value if abs(numeric_value) >= 100 else numeric_value * 100
    return f'{scaled_value:g}'


def _normalize_logp_condition(condition) -> str:
    """Accept cLogP units (3.5) as well as UCSF cLogP*100 units (350)."""
    condition = str(condition).strip().replace(' ', '')
    operator = re.fullmatch(f'(<=|>=|<|>|=)({_NUMBER_PATTERN})', condition)
    if operator:
        return operator.group(1) + _format_scaled_logp(operator.group(2))

    number = re.fullmatch(_NUMBER_PATTERN, condition)
    if number:
        return _format_scaled_logp(number.group(0))

    interval = re.fullmatch(
        f'({_NUMBER_PATTERN})-({_NUMBER_PATTERN})', condition
    )
    if interval:
        return '-'.join(_format_scaled_logp(value) for value in interval.groups())

    raise ValueError(
        f'Invalid logP condition {condition!r}; use 3.5, 350, <=3.5, '
        'or a range such as 1.0-3.5'
    )


class Filters():
    """A class to apply various filtering operations on molecular data.
        
    Parameters
    ---------
        
        removesalts: bool, default False.
            Remove salts from the input molecules. Only the largest fragment will be retained if the salt is not in the database.
        
        ha: str, default None.
            Filter molecules by the number of heavy atoms. Accepts a string in the format '1-5', '>=5', '<=5', or '5'.
        
        logp: str, default None.
            Filter molecules by the logP value. Accepts a string in the format '1-5', '>=5', '<=5', or '5'.
        
        hba: str, default None.
            Filter molecules by the number of H-bond acceptors. Accepts a string in the format '1-5', '>=5', '<=5', or '5'.

        hbd: str, default None.
            Filter molecules by the number of H-bond donors. Accepts a string in the format '1-5', '>=5', '<=5', or '5'.

        mw: str, default None.
            Filter molecules by the molecular weight. Accepts a string in the format '1-5', '>=5', '<=5', or '5'.

        tpsa: str, default None.
            Filter molecules by topological polar surface area (TPSA). Accepts a string in the format '20-80', '>=20', '<=80', or '50'.

        fsp3: str, default None.
            Filter molecules by the fraction of sp3 carbon atoms (FSP3). Accepts a string in the format '0.2-0.8', '>=0.2', '<=0.8', or '0.5'.

        chiral: str, default None.
            Filter molecules by the number of unspecified chiral centers. Accepts a string in the format '1-5', '>=5', '<=5', or '5'.

        custom: str, default None.
            Path to a custom SMARTS file for filtering molecules. The file should contain SMARTS patterns, one per line.

        unwanted: list of str. Default None.
            Filter molecules by unwanted substructures. Acceptes in ['regular', 'optional', 'special']

        pains: bool, default False.
            Apply the PAINS filter to the molecules. If True, it will check for PAINS functional groups.

        numcores: int, default 1.
            Number of native threads used by the PAINS catalog search.

        rejectedFile: str, default 'rejected_entries.txt'.
            Path to the file where rejected molecules will be saved. The file will be created if it does not exist.

        debug: bool, default False.
            Print debug message

    Examples
    ----------

    There are two ways to use the Filters class. The first is to create an instance of the class with the desired parameters and then call the filter_df method on a DataFrame:    
    
    >>> from msani.filtering.filters import Filters
    >>> filters = Filters(removesalts=True,
                            ha='>=5',
                            logp='<=3.5',
                            hba='1-3',
                            hbd='1-2',
                            mw='200-500',
                            tpsa='20-100',
                            fsp3='>=0.25',
                            chiral='0-2',
                            custom='path/to/custom.smarts',
                            unwanted=['regular'],
                            pains=True,
                            rejectedFile='rejected.txt',
                            debug=True)
    >>> filtered_df = filters.filter_df(df)
    
    the second way is to use the static methods of the class directly on a DataFrame:

    >>> from msani.filtering.filters import Filters
    >>> rejectedFile = 'rejected.txt'
    >>> df = Filters.filter_by_ha(df, '>=5', rejectedFile, debug)
    >>> df = Filters.filter_by_logp(df, '<=3.5', rejectedFile, debug)
    >>> df = Filters.filter_by_hba(df, '1-3', rejectedFile, debug)
    >>> df = Filters.filter_by_hbd(df, '1-2', rejectedFile, debug)
    >>> df = Filters.filter_by_mw(df, '200-500', rejectedFile, debug)
    >>> df = Filters.filter_by_tpsa(df, '20-100', rejectedFile, debug)
    >>> df = Filters.filter_by_fsp3(df, '>=0.25', rejectedFile, debug)
    >>> df = Filters.filter_by_chiralcenters(df, '0-2', rejectedFile, debug)
    >>> df = Filters.customFilter(df, rejectedFile, 'custom_smarts.txt', debug)
    >>> df = Filters.unwantedFilter(df, rejectedFile, 'all', debug)
    >>> df = Filters.painsFilter(df, rejectedFile, debug)
    
    """
    def __init__(self, 
                 removesalts = False, 
                 ha = None,
                 logp = None,
                 hba = None,
                 hbd = None,
                 mw = None,
                 chiral = None,
                 custom = None,
                 unwanted = None,
                 pains = None,
                 rejectedFile = 'rejected_entries.txt',
                 tpsa = None,
                 fsp3 = None,
                 numcores = 1):
        self.removesalts = removesalts
        self.ha = ha
        self.logp = logp
        self.hba = hba
        self.hbd = hbd
        self.mw = mw
        self.tpsa = tpsa
        self.fsp3 = fsp3
        self.chiral = chiral

        self.custom = custom
        self.unwanted = unwanted

        if self.custom is not None or self.unwanted is not None:
            if self.custom: temp_df_custom = loadSMARTSdata(self.custom)
            if self.unwanted:
                smartsFile = Path(__file__).parent.parent / 'Data' / 'filter_out.txt'
                temp_df_unwanted = loadSMARTSdata(smartsFile.resolve(), self.unwanted)
            self.unwanted_df = concat([temp_df_custom, temp_df_unwanted]) if self.custom and self.unwanted \
                else temp_df_custom if self.custom else temp_df_unwanted
            self.unwanted_catalog = Filters.buildSMARTScatalog(self.unwanted_df)
            logger.info(f"Loaded {len(self.unwanted_df)} SMARTS patterns for unwanted filtering.")
        else:
            self.unwanted_df = None
            self.unwanted_catalog = None

        self.pains = pains
        self.rejectedFile = rejectedFile
        self.numcores = numcores

    def __repr__(self):
        cls_name = self.__class__.__name__
        attrs = ', '.join(f'{k}={v!r}' for k, v in self.__dict__.items())
        return f'{cls_name}({attrs})'
    
    @staticmethod
    def remove_invalid_SMILES(df:DataFrame) -> DataFrame:
        """Remove rows with invalid SMILES from the input DataFrame.

        Args:
            df (DataFrame): Input DataFrame with 'smiles' column containing SMILES strings.

        Returns:
            DataFrame: A new DataFrame chunk with valid SMILES strings.
        """
        # Log rows where 'mol' is None before dropping
        invalid_rows = df[df['mol'].isna()]
        for index, row in invalid_rows.iterrows():
            logger.warning(f"INVALID SMILES: Removing row at index {index}: SMILES: {row['smiles']} - ID: {row['ids']}")

        # Remove rows where 'mol' is None
        df_cleaned = df.dropna(subset=['mol'])

        return df_cleaned
    
    @staticmethod
    def remove_exotic_chem_to_db2(df:DataFrame) -> DataFrame:
        '''
        These are substructures that are not supported by either MMFF94(s)
        and Mol2 format so could not be DB2-compatible.
        It mainly includes: hypervalent atoms, halogens with more than 1 bond, atoms with more than 5 bonds,
        and atom types not supported by the Mol2 format.

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
        Returns:
            DataFrame: A new DataFrame chunk with molecules that are DB2-compatible.

        '''
        exotic_chems = Chem.MolFromSmarts('[$([*X{5-}]),$([#35!X1]),$([#53!X1]),$([*;!$([#1,#6,#7,#8,#9,#14,#15,#16,#17,#35,#53,#26,#3,#11,#19,#30,#20,#29,#12])])]')
        df_cleaned = df.copy()
        df_cleaned['mol'] = df['mol'].apply(lambda x: x if x.HasSubstructMatch(exotic_chems) == False else None)
        exotic_chem_df = df_cleaned[df_cleaned['mol'].isna()]
        for _, row in exotic_chem_df.iterrows():
            logger.warning(f"Removed entries that are not DB2/Mol2/PDBQT-compatible: SMILES: {row['smiles']} - ID: {row['ids']}")
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

        if debug and len(deleted) > 0: 
            logger.info(f"Stripped salt {Chem.MolToSmiles(mol)}:  Salts:{' '.join([Chem.MolToSmiles(m) for m in deleted])}")
        if len(Chem.rdmolops.GetMolFrags(res)) > 1:
            # If still contains more than one fragment, retains the largest one
            rdMolStandardize.FragmentParentInPlace(res)

        return res

    @staticmethod
    def stripalkali(mol, debug = False):
        """Strip alkali metals (Na, K) from the input molecule.
        Args:
            mol (rdkit mol object): The input molecule.
            debug (bool, optional): Debug mode. Defaults to False.
        Returns:
            rdkit mol object: The stripped molecule.
        """
        res = Chem.Mol(mol)
        rdMolStandardize.FragmentParentInPlace(res)
        if debug: 
            logger.info(f"Stripped alkali metals: {Chem.MolToSmiles(mol)} -> {Chem.MolToSmiles(res)}")
        return res
    
    @staticmethod
    def saltstripping(df: DataFrame, debug = False) -> DataFrame:
        """Remove salts from the input molecules using the RDKit SaltRemover class.

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with salt-stripped molecules.
        """

        global salt_remover

        # We will track molecules and smiles directly
        # To avoid a full df copy, work directly on the series. 
        # But since we'll drop rows, it's safer to work on copies of the columns.
        new_mols = df['mol'].copy()
        new_smiles = df['smiles'].copy()

        # Only process the entries with '.' in the SMILES (multiple )
        has_dot = new_smiles.str.contains('.', regex=False, na=False)
        idx = df.index[has_dot]
        
        updated_mols = new_mols.loc[idx].apply(lambda m: Filters.stripSMILESsalt(m, salt_remover, debug))
        new_mols.loc[idx] = updated_mols
        new_smiles.loc[idx] = updated_mols.map(Chem.MolToSmiles)

        # Disconnect alkali metals (Na, K) from the molecules Issue #33
        has_alkali = new_smiles.str.contains(r'Na|K', regex=True, na=False)
        idx_alkali = df.index[has_alkali]
        
        updated_mols_alkali = new_mols.loc[idx_alkali].apply(lambda m: Filters.stripalkali(m, debug))
        new_mols.loc[idx_alkali] = updated_mols_alkali
        new_smiles.loc[idx_alkali] = updated_mols_alkali.map(Chem.MolToSmiles)

        # Filter out empty SMILES
        valid_mask = (new_smiles != '')
        
        if debug: 
            dropped_count = len(df) - valid_mask.sum()
            logger.info(f"Removed pure {dropped_count} salts from the input molecules.")
            print(f"Removed pure {dropped_count} salts from the input molecules.")

        # Sanitize and keep valid ones
        errors = new_mols[valid_mask].apply(lambda x: Chem.SanitizeMol(x, catchErrors=True))
        
        # Log errors
        error_mask = (errors != 0)
        if error_mask.any():
            for idx_err in errors[error_mask].index:
                logger.warning(f'''SANITIZE FAILED: SMILES: {new_smiles[idx_err]} - ID: {df.loc[idx_err, 'ids']}''')
        
        # Final passing mask
        final_mask = valid_mask.copy()
        final_mask.loc[valid_mask] = ~error_mask
        
        # Create the resulting dataframe only once
        filtered_df = df.loc[final_mask].copy()
        filtered_df['mol'] = new_mols.loc[final_mask]
        filtered_df['smiles'] = new_smiles.loc[final_mask]
        
        return filtered_df
    
    @staticmethod
    def convert_to_query(condition: str, column: str) -> str:
        if "-" in condition:  # Range condition
            lower, upper = map(float, condition.split('-'))
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
    def filter_mask(values: Series, condition: str, column: str) -> Series:
        """Evaluate a descriptor condition directly against a Series."""
        condition = str(condition).strip().replace(' ', '')
        operator = re.fullmatch(f'(<=|>=|<|>|=)({_NUMBER_PATTERN})', condition)
        if operator:
            value = float(operator.group(2))
            comparisons = {
                '>=': values.ge,
                '<=': values.le,
                '>': values.gt,
                '<': values.lt,
                '=': values.eq,
            }
            return comparisons[operator.group(1)](value)

        interval = re.fullmatch(
            f'({_NUMBER_PATTERN})-({_NUMBER_PATTERN})', condition
        )
        if interval:
            lower, upper = map(float, interval.groups())
            return (values >= lower) & (values <= upper)

        number = re.fullmatch(_NUMBER_PATTERN, condition)
        if number:
            value = float(number.group(0))
            return values.le(value) if column == 'logp' else values.eq(value)

        raise ValueError(f'Invalid condition: {condition}')

    @staticmethod    
    def filter_by_ha(df, filter_query, rejectedFile, debug = False) -> DataFrame:
        """Filter out molecules with heavy atoms only using the RDKit Mol.GetNumHeavyAtoms() function.

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        ha_series = df['mol'].apply(lambda x: x.GetNumHeavyAtoms())
        mask = Filters.filter_mask(ha_series, filter_query, 'ha')
        rejected_mask = ~mask
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['ha'] = ha_series.loc[rejected_mask]
            rejected_df['ha'] = rejected_df['ha'].apply(lambda x: f'ha{x}')
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug: 
                logger.info(f"Removed {len(rejected_df)} molecules with heavy atoms requirements: {rejected_df['ha'].values}")
                print(f"Removed {len(rejected_df)} molecules with heavy atoms requirements: {rejected_df['ha'].values}")
        return df.loc[mask].copy()

    @staticmethod
    def filter_by_logp(df, filter_query, rejectedFile, debug = False) -> DataFrame:
        """Filter out molecules with a logP value upto the defined value using the RDKit Crippen logP calculation.

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (_type_): Path to the file to save rejected molecules.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with molecules containing the specified logP value.
        """
        filter_query = _normalize_logp_condition(filter_query)
        logp_series = df['mol'].apply(lambda mol: MolLogP(mol) * 100)
        mask = Filters.filter_mask(logp_series, filter_query, 'logp')
        rejected_mask = ~mask
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['logp'] = logp_series.loc[rejected_mask]
            rejected_df['logp'] = rejected_df['logp'].apply(lambda x: f'logp {x/100:.2f}')
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug: 
                logger.info(f"Removed {len(rejected_df)} molecules with logP requirements: {rejected_df['logp'].values}")
                print(f"Removed {len(rejected_df)} molecules with logP requirements: {rejected_df['logp'].values}")
        
        return df.loc[mask].copy()
    
    def filter_by_hba(df, filter_query, rejectedFile, debug = False) -> DataFrame:
        """Filter out molecules with required number of H-bond acceptors using the RDKit CalcNumHBA().

        NOTE: It is recently updated to use the Python function CalcNumHBA() within this file to remain the consistency between the RDKit versions.
        The previous version of RDKit CalcNumHBA() was updated in 2025.09.3 so that we may have different results between the different versions of RDKit.
        SMARTS pattern: [$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]
        
        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        hba_series = df['mol'].apply(CalcNumHBA)
        mask = Filters.filter_mask(hba_series, filter_query, 'hba')
        rejected_mask = ~mask
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['hba'] = hba_series.loc[rejected_mask]
            rejected_df['hba'] = rejected_df['hba'].apply(lambda x: f'hba {x}')
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug: 
                logger.info(f"Removed {len(rejected_df)} molecules with number of H-bond acceptors requirements: {rejected_df['hba'].values}")
                print(f"Removed {len(rejected_df)} molecules with number of H-bond acceptors requirements: {rejected_df['hba'].values}")
        
        return df.loc[mask].copy()
    
    def filter_by_hbd(df, filter_query, rejectedFile, debug = False) -> DataFrame:
        """Filter out molecules with required number of H-bond donors using the RDKit CalcNumLipinskiHBD().

        NOTE: It is by intention that the function uses CalcNumLipinskiHBD() instead of CalcNumHBD() was used as we believe that it better represents the chemistry.
              This will run in C++ to count the total number of OH and NH. (NH2 = 2 donors)
              
        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        hbd_series = df['mol'].apply(rdMolDescriptors.CalcNumLipinskiHBD)
        mask = Filters.filter_mask(hbd_series, filter_query, 'hbd')
        rejected_mask = ~mask
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['hbd'] = hbd_series.loc[rejected_mask]
            rejected_df['hbd'] = rejected_df['hbd'].apply(lambda x: f'hbd {x}')
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug: 
                logger.info(f"Removed {len(rejected_df)} molecules with number of H-bond donors requirements: {rejected_df['hbd'].values}")
                print(f"Removed {len(rejected_df)} molecules with number of H-bond donors requirements: {rejected_df['hbd'].values}")
        
        return df.loc[mask].copy()

    def filter_by_mw(df, filter_query, rejectedFile, debug = False) -> DataFrame:
        """Filter out molecules with required molecular weight using the RDKit GetMolWt().

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        mw_series = df['mol'].apply(rdMolDescriptors.CalcExactMolWt)
        mask = Filters.filter_mask(mw_series, filter_query, 'mw')
        rejected_mask = ~mask
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['mw'] = mw_series.loc[rejected_mask]
            rejected_df['mw'] = rejected_df['mw'].apply(lambda x: f'mw {x:.2f}')
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug: 
                logger.info(f"Removed {len(rejected_df)} molecules with molecular weight requirements: {rejected_df['mw'].values}")
                print(f"Removed {len(rejected_df)} molecules with molecular weight requirements: {rejected_df['mw'].values}")
        
        return df.loc[mask].copy()

    @staticmethod
    def filter_by_tpsa(df, filter_query, rejectedFile, debug = False) -> DataFrame:
        """Filter molecules by topological polar surface area (TPSA)."""
        tpsa_series = df['mol'].apply(rdMolDescriptors.CalcTPSA)
        mask = Filters.filter_mask(tpsa_series, filter_query, 'tpsa')
        rejected_mask = ~mask

        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['tpsa'] = tpsa_series.loc[rejected_mask]
            rejected_df['tpsa'] = rejected_df['tpsa'].apply(lambda x: f'tpsa {x:.2f}')
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug:
                logger.info(f"Removed {len(rejected_df)} molecules with TPSA requirements: {rejected_df['tpsa'].values}")
                print(f"Removed {len(rejected_df)} molecules with TPSA requirements: {rejected_df['tpsa'].values}")

        return df.loc[mask].copy()

    @staticmethod
    def filter_by_fsp3(df, filter_query, rejectedFile, debug = False) -> DataFrame:
        """Filter molecules by the fraction of sp3 carbon atoms (FSP3)."""
        fsp3_series = df['mol'].apply(rdMolDescriptors.CalcFractionCSP3)
        mask = Filters.filter_mask(fsp3_series, filter_query, 'fsp3')
        rejected_mask = ~mask

        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['fsp3'] = fsp3_series.loc[rejected_mask]
            rejected_df['fsp3'] = rejected_df['fsp3'].apply(lambda x: f'fsp3 {x:.3f}')
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug:
                logger.info(f"Removed {len(rejected_df)} molecules with FSP3 requirements: {rejected_df['fsp3'].values}")
                print(f"Removed {len(rejected_df)} molecules with FSP3 requirements: {rejected_df['fsp3'].values}")

        return df.loc[mask].copy()

    @staticmethod
    def count_unspecified_chiralcenters(mol):
        """Count the number of unspecified chiral centers in a molecule.

        Args:
            mol (rdkit mol object): The input molecule.

        Returns:
            int: The number of unspecified chiral centers.
        """
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        unassigned = [idx for idx, tag in centers if tag == '?']
        return len(unassigned)
    
    @staticmethod
    def filter_by_chiralcenters(df, filter_query, rejectedFile, debug = False) -> DataFrame:
        """Filter out molecules with required number of unspecified chiral centers.

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with molecules containing heavy atoms.
        """
        cc_series = df['mol'].apply(rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters)
        mask = Filters.filter_mask(cc_series, filter_query, 'chiralcenters')
        rejected_mask = ~mask
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['chiralcenters'] = cc_series.loc[rejected_mask]
            rejected_df['chiralcenters'] = rejected_df['chiralcenters'].apply(lambda x: f'chiralcenters {x}')
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug: 
                logger.info(f"Removed {len(rejected_df)} molecules with number of chiral centers requirements: {rejected_df['chiralcenters'].values}")
                print(f"Removed {len(rejected_df)} molecules with number of chiral centers requirements: {rejected_df['chiralcenters'].values}")
        return df.loc[mask].copy()
    
    @staticmethod
    def applyStandarizeFilters(mol):

        uncharged_parent_clean_mol = None

        try:

            clean_mol = rdMolStandardize.Cleanup(mol) 

            # if many fragments, get the "parent"
            parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)

            # try to neutralize molecule
            global uncharger
            uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
            

        except:

            logger.info(f'Molecule NOT processed: {Chem.MolToSmiles(mol)}')

        return uncharged_parent_clean_mol

    @staticmethod
    def standarizeFilters(df: DataFrame) -> DataFrame:


        # follows the steps in
        # https://github.com/greglandrum/RSC_OpenScience_Standardization_202104/blob/main/MolStandardize%20pieces.ipynb
        # https://www.youtube.com/watch?v=eWTApNX8dJQ
        # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule

        #Sc, Y, In, Sn, W, Ac are not included for now (need to check how they can hit other structures).
        organometallics={'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Ga', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',\
        'Cd', 'La', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Ce', 'Pr', 'Nd',\
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Ge', 'Sb'}


        # Filter out rows where 'smiles' contains any organometallics
        filtered_df = df[~df['smiles'].apply(lambda x: any(om in x for om in organometallics))].copy()
        if filtered_df.empty:
            return filtered_df

        filtered_df['mol'] = filtered_df['mol'].apply(lambda x:  Filters.applyStandarizeFilters(x))
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
    def painsFilter(df: DataFrame, rejectedFile: str, debug: bool = False,
                    numcores: int = 1) -> DataFrame:
        """
        Detect and filter out molecules with PAINS functional groups using the RDKit PAINS catalog.

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (str): Path to the file to save rejected molecules.
            debug (bool, optional): Debug mode. Defaults to False.
            numcores (int, optional): Number of RDKit catalog threads. Defaults to 1.

        Returns:
            DataFrame: A new DataFrame chunk with molecules that passed the PAINS filter.
        """
        # Set up the PAINS catalog
        global _pains_catalog
        if rdBase.rdkitVersion < '2025.09.3':
            print('\nThe warning is expected and can be ignored.\n')
        if _pains_catalog is None:
            params = RDFilterCatalog.FilterCatalogParams()
            params.AddCatalog(RDFilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
            _pains_catalog = RDFilterCatalog.FilterCatalog(params)

        if numcores == 1:
            # Avoid serialization and reparsing overhead when parallelism is
            # disabled, while retaining the exact original single-core path.
            reasons = df['mol'].apply(
                lambda mol: Filters.detect_and_label_pains(mol, _pains_catalog)
            )
        else:
            # Run the native catalog search in parallel. Generate the input from
            # the current molecules because transformations such as neutralization
            # may intentionally leave the original text in the 'smiles' column.
            current_smiles = [Chem.MolToSmiles(mol) for mol in df['mol'].array]
            matches = RDFilterCatalog.RunFilterCatalog(
                _pains_catalog, current_smiles, numThreads=numcores
            )
            reasons = Series(
                [
                    'OK' if not entries else
                    'PAINS violation: ' + str(entries[0].GetDescription().capitalize())
                    for entries in matches
                ],
                index=df.index,
            )

        # Separate rejected molecules
        mask_ok = (reasons == 'OK')
        rejected_mask = ~mask_ok
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['reason'] = reasons[rejected_mask]
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            if debug:
                print(f"Removed {len(rejected_df)} molecules due to PAINS violations.")

        # Filter passed molecules
        return df.loc[mask_ok].copy()
    
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
    def filterbysmarts(mol, smarts_df: DataFrame) -> str:
        """Return the first matching label from a catalog or SMARTS DataFrame."""
        if hasattr(smarts_df, 'GetFirstMatch'):
            entry = smarts_df.GetFirstMatch(mol)
            return entry.GetDescription() if entry is not None else 'OK'

        # Preserve the public helper API for callers passing a SMARTS DataFrame.
        for substructure in smarts_df.itertuples(index=False):
            if mol.HasSubstructMatch(substructure.mol):
                return substructure.label
        return 'OK'

    @staticmethod
    def buildSMARTScatalog(smarts_df: DataFrame):
        """Compile ordered SMARTS patterns into an RDKit FilterCatalog."""
        catalog = RDFilterCatalog.FilterCatalog()
        for substructure in smarts_df.itertuples(index=False):
            label = str(substructure.label)
            matcher = RDFilterCatalog.SmartsMatcher(label, substructure.mol)
            catalog.AddEntry(RDFilterCatalog.FilterCatalogEntry(label, matcher))
        return catalog
    
    @staticmethod
    def unwantedFilter(df: DataFrame, rejectedFile, unwanted_option = None,
                       unwanted_df: DataFrame = None, debug = False,
                       unwanted_catalog = None) -> DataFrame:
        """Filter out unwanted substructures using a default list of SMARTS patterns.

            Args:
                df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
                rejectedFile (str): Path to the file to save rejected molecules.
                unwanted_option (list): The mode input by the user.
                unwanted_df (DataFrame): DataFrame containing the SMARTS patterns and their corresponding RDKit molecule objects.
                debug (bool, optional): Debug mode. Defaults to False.
                unwanted_catalog: Precompiled RDKit catalog. Defaults to None.

            Returns:
                DataFrame: A new DataFrame chunk with molecules that passed the filter.
        """
        if unwanted_df is None:
            # Backward compatible, if the user already provide the DataFrame, then use it
            # Get the absolute path to the template SMARTS file using pathlib
            smartsFile = Path(__file__).parent.parent / 'Data' / 'filter_out.txt'

            # Load smarts to clean  from file
            unwanted_df = loadSMARTSdata(smartsFile.resolve(), unwanted_option)
            logger.info(f'Parsed {len(unwanted_df)} substructures from: {smartsFile}')

        if unwanted_catalog is None:
            unwanted_catalog = Filters.buildSMARTScatalog(unwanted_df)

        # Apply reactions to each SMILES in the DataFrame
        reasons = df['mol'].apply(lambda x: Filters.filterbysmarts(x, unwanted_catalog))
        
        mask_ok = (reasons == 'OK')
        rejected_mask = ~mask_ok
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['reason'] = reasons[rejected_mask]
            
            if debug:
                logger.info(f"Removed {len(rejected_df)} molecules with unwanted substructures")
                print(f"Removed {len(rejected_df)} molecules with unwanted substructures")
                
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            
        return df.loc[mask_ok].copy()


    @staticmethod
    def customFilter(df: DataFrame, rejectedFile, smartsFile, debug = False) -> DataFrame:
        """Filter out unwanted substructures using a customized list of SMARTS patterns.

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (_type_): Path to the file to save rejected molecules.
            smartsFile (_type_): Path to the file containing the customized list of SMARTS patterns.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with molecules that passed the filter.
        """
        # Load smarts to clean  from file
        unwanted_df = loadSMARTSdata(smartsFile.resolve())
        if debug: logger.info(f'Parsed {len(unwanted_df)} custom substructures from: {smartsFile}')
        unwanted_catalog = Filters.buildSMARTScatalog(unwanted_df)

        # Apply reactions to each SMILES in the DataFrame
        reasons = df['mol'].apply(lambda x: Filters.filterbysmarts(x, unwanted_catalog))
        mask_ok = (reasons == 'OK')
        rejected_mask = ~mask_ok
        
        if rejected_mask.any():
            rejected_df = df.loc[rejected_mask, ['smiles', 'ids']].copy()
            rejected_df['reason'] = reasons[rejected_mask]
            rejected_df.to_csv(rejectedFile, index=False, mode='a', sep=' ', header=False)
            
        return df.loc[mask_ok].copy()
    
    def filter_df(self, df: DataFrame, rejectedFile: str, debug: bool = False) -> DataFrame:
        """Apply the filters to the input DataFrame.

        Args:
            df (DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (str): Path to the file to save rejected molecules.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            DataFrame: A new DataFrame chunk with molecules that passed all the filters.
        """
        if self.removesalts:
            df = Filters.saltstripping(df, debug)
        if self.ha:
            df = Filters.filter_by_ha(df, self.ha, rejectedFile, debug)
        if self.hba:
            df = Filters.filter_by_hba(df, self.hba, rejectedFile, debug)
        if self.hbd:
            df = Filters.filter_by_hbd(df, self.hbd, rejectedFile, debug)
        if self.mw:
            df = Filters.filter_by_mw(df, self.mw, rejectedFile, debug)
        if self.tpsa:
            df = Filters.filter_by_tpsa(df, self.tpsa, rejectedFile, debug)
        if self.fsp3:
            df = Filters.filter_by_fsp3(df, self.fsp3, rejectedFile, debug)
        if self.chiral:
            df = Filters.filter_by_chiralcenters(df, self.chiral, rejectedFile, debug)
        if self.logp:
            df = Filters.filter_by_logp(df, self.logp, rejectedFile, debug)
        if self.unwanted_df is not None:
            df = Filters.unwantedFilter(
                df, 
                rejectedFile, 
                unwanted_option=None, 
                unwanted_df=self.unwanted_df, 
                debug=debug,
                unwanted_catalog=self.unwanted_catalog,
            )
        if self.pains:
            df = Filters.painsFilter(df, rejectedFile, debug, self.numcores)
        return df

def loadSMARTSdata(smartsFile: str, unwanted_option=None) -> DataFrame:
        """Load SMARTS patterns from a file and convert them to RDKit molecule objects.

        Args:
            smartsFile (str): Path to the file containing the SMARTS patterns.
            unwanted_option (list): The mode input by the user. Defaults to None.

        Returns:
            DataFrame: A DataFrame containing the SMARTS patterns and their corresponding RDKit molecule objects.
        """
        if unwanted_option is not None: 
            # Using default substructure file 
            smarts_df = read_csv(smartsFile, sep=r'\s+', header=0, names=['smarts','label', 'reason', 'mode', 'ref'])
            smarts_df = smarts_df[smarts_df["mode"].isin(unwanted_option)]
        else:
            # Using customized substructure file
            has_header = Filters.check_header(smartsFile)
            if has_header:
                logger.info(f'Found header in {smartsFile}')
                smarts_df = read_csv(smartsFile, sep=r'\s+', header=0, usecols=[0,1], names=['smarts','label'])
            else:
                smarts_df = read_csv(smartsFile, sep=r'\s+', header=None, usecols=[0,1], names=['smarts','label'])

        smarts_df['mol'] = smarts_df['smarts'].apply(lambda x: Chem.MolFromSmarts(x)) #do we need mergeHs here?
        
        return smarts_df

hba_smarts = Chem.MolFromSmarts("[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]")

def CalcNumHBA(mol):
    return len(mol.GetSubstructMatches(hba_smarts))
