import logging
import multiprocessing as mp
from functools import partial

from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from pandas import DataFrame, Series  # only what you use

from msani.filtering.filters import Filters, against_humanity, hold_up
from msani.moltransform.tautomerizer import Tautomerizer
from msani.moltransform.ionizer import Ionizer
from msani.moltransform.neutralizer import Neutralizer

logger = logging.getLogger('msani')
RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit

class Msani:
    """
    A class to store the filter options and conduct chemical modifications for MolSanitizer. 
    Initialize the class with the desired filter options and apply the filters to the input DataFrame.
    
    Example use:

        processor = Msani(
                    removesalts=True,
                    ha='>10',
                    logp='<5',
                    tautomers=True, 
                    protonation = True, 
                    debug=True)

        processed_df = processor.run(df)

    """
    tautomer_params = None  # Define as a class variable

    def __init__(self, 
                removesalts=False,
                custom = None,
                unwanted = None,
                pains = False,
                ha = None,
                logp = None,
                hba = None,
                hbd = None,
                mw = None,
                chiral = None,
                tautomers = False,
                neutralize = True,
                taurdkit = True,
                stereoisomers = False,
                max_stereoisomers = 8,
                protonation = False,
                pH = 7,
                pH_range = 0,
                numcores = 1,
                standardize = False,
                protonation_library = None,
                tautomer_library = None,
                debug = False):

        self.removesalts = removesalts
        self.ha = ha
        self.logp = logp
        self.hba = hba
        self.hbd = hbd
        self.mw = mw
        self.custom = custom
        self.unwanted = unwanted
        self.pains = pains
        self.chiral = chiral
        
        self.tautomers = tautomers
        self.neutralize = neutralize
        self.taurdkit = taurdkit
        self.stereoisomers = stereoisomers
        self.max_stereoisomers = max_stereoisomers
        self.protonation = protonation
        self.pH = pH
        self.pH_range = pH_range
        self.standardize = standardize    
        self.debug = debug
        self.numcores = numcores
        self.protonation_library = protonation_library
        self.tautomer_library = tautomer_library
    
    def __repr__(self):
        cls_name = self.__class__.__name__
        attrs = ', '.join(f'{k}={v!r}' for k, v in self.__dict__.items())
        return f'{cls_name}\n({attrs})'
    


    @staticmethod
    def _generate_stereoisomers(mol, max_isomers):
        """
        Generate stereoisomers for a given molecule.
        
        Args:
        mol (rdkit.Chem.Mol): RDKit molecule object.
        max_isomers (int): Maximum number of stereoisomers to generate.
        
        Returns:
        list: List of stereoisomer molecules.
        """
        opts = StereoEnumerationOptions(tryEmbedding=True, unique=True, maxIsomers=max_isomers)
        isomers = list(EnumerateStereoisomers(mol, options=opts))
        return isomers

    @staticmethod
    def _process_molecule_stereoisomer(row_data, max_isomers=8):
        """
        A wrapper to generate stereoisomers for a given molecule using multiprocessing.
        Args:
        row_data (Series): A row from the input DataFrame containing ['smiles', 'ids', 'mol'] and optionally 'highlights'.
        max_isomers (int): Maximum number of stereoisomers to generate.

        Returns:
        results (list): List of expanded stereoisomers dictionaries with 'smiles', 'ids', 'mol', and optionally 'highlights'.
        """
        mol = Chem.MolFromSmiles(row_data['smiles'])
        # If max_isomers is set to 1, return the original molecule and let the RDKit/CORINA guess it.
        if max_isomers == 1: return [row_data]
        try:
            isomers = Msani._generate_stereoisomers(mol, max_isomers=max_isomers)
        except Exception as e:
            logger.error(f"Error generating stereoisomers for compound {row_data['ids']}: {row_data['smiles']}")
            isomers = [mol]
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        unassigned = [idx for idx, tag in centers if tag == '?']
        num_possible_isomers = 2 ** len(unassigned)
        result = []
        # Get highlights if available
        highlights = row_data.get('highlights', None)
        if len(isomers) == 1:
            result.append({'smiles': Chem.MolToSmiles(isomers[0]),
                            'ids': row_data['ids'],
                            'mol': isomers[0],
                            'highlights': highlights})
        else:
            if max_isomers > 0 and num_possible_isomers > max_isomers:
                logger.info(f"{row_data['ids']}: Not all the stereoisomers are written out (capped at {max_isomers}/{num_possible_isomers}).")
                isomers = isomers[:max_isomers]
            two_digits = len(isomers) >= 10
            for i, isomer in enumerate(isomers):
                result.append({
                    'smiles': Chem.MolToSmiles(isomer, isomericSmiles=True),
                    'ids': row_data['ids'] + '.' + (f"{i+1:02}" if two_digits else f"{i+1}"),
                    'mol': isomer,
                    'highlights': highlights
                })
        
        return result

    @staticmethod
    def enum_stereoisomers(df: DataFrame, max_isomers=8, numcores=4, debug=False) -> DataFrame:
        """
        Generate stereoisomers for molecules in the 'smiles' column and expand the DataFrame using multiprocessing.
        
        Args:
        df (DataFrame): DataFrame with 'smiles' and 'ids' columns.
        max_isomers (int): Maximum number of stereoisomers to generate for each molecule.
        numcores (int): Number of processes to use. Default is 4.
        debug (bool): Enable debug messages.
        
        Returns:
        DataFrame: Expanded DataFrame with each stereoisomer as a separate row.
        """
        if df.empty:
            logger.warning("Empty DataFrame provided, skipping stereoisomer generation.")
            return df
        # Partial function to fix max_isomers as an argument
        process_func = partial(Msani._process_molecule_stereoisomer, max_isomers=max_isomers)
        results = []

        with mp.Pool(processes=numcores) as pool:
            # Submit all tasks to the pool in parallel, keeping track of the rows for error handling
            async_results = [(row, pool.apply_async(process_func, (row,))) for _, row in df.iterrows()]

            # Collect results as they complete
            for row_data, async_result in async_results:
                try:
                    chunk_result = async_result.get(timeout=60)
                    results.extend(chunk_result)
                except (mp.TimeoutError, Exception) as e:
                    # Differentiate the logging based on the type of exception if desired
                    if isinstance(e, mp.TimeoutError):
                        logger.warning(f"Timeout occurred for compound {row_data['ids']}. Using original molecule.")
                    else:
                        logger.error(f"Error processing compound {row_data['ids']}: {str(e)}. Using original molecule.")
                    
                    highlights = row_data.get('highlights', None)
                    results.append({
                        'smiles': row_data['smiles'],
                        'ids': row_data['ids'],
                        'mol': row_data['mol'],
                        'highlights': highlights
                    })

        return DataFrame(results)
    
    def expand_ids(self, group):
        if len(group) == 1:
            return group
        two_digits = len(group) >= 10
        new_id = [f"{group.iloc[0]}_{i+1:02}" if two_digits else f"{group.iloc[0]}_{i+1}" for i in range(len(group))]
        return Series(new_id, index=group.index)
    
    def run(self, df: DataFrame, rejected_file = None) -> DataFrame:
        """
        Perform preparation on the input DataFrame using the specified filters and rule-based chemical modifications.
        
        Args:
            df (DataFrame): Input DataFrame with ['smiles', 'ids'] columns. 
            
        Returns:
            DataFrame: A new DataFrame chunk with molecules that passed the filters.
        """
        df.loc[:, 'mol'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x))
        # df.loc[:, 'against_humanity'] = df['mol'].apply(lambda x: x.HasSubstructMatch(against_humanity) if x else False)
        # if len(df[df['against_humanity'] == True]) > 0: print(hold_up)
        df = Filters.remove_invalid_SMILES(df)
        if self.standardize: df = Filters.remove_exotic_chem_to_db2(df)
        if self.removesalts: df = Filters.saltstripping(df, debug=self.debug)
        if self.ha is not None: df = Filters.filter_by_ha(df, self.ha, rejectedFile=rejected_file, debug=self.debug)
        if self.logp is not None: df = Filters.filter_by_logp(df, self.logp, rejectedFile=rejected_file, debug=self.debug)
        if self.hba is not None: df = Filters.filter_by_hba(df, self.hba, rejectedFile=rejected_file, debug=self.debug)
        if self.hbd is not None: df = Filters.filter_by_hbd(df, self.hbd, rejectedFile=rejected_file, debug=self.debug)
        if self.mw is not None: df = Filters.filter_by_mw(df, self.mw, rejectedFile=rejected_file, debug=self.debug)
        if self.chiral is not None: df = Filters.filter_by_chiralcenters(df, self.chiral, rejectedFile=rejected_file, debug=self.debug)
        if self.standardize: 
            df = Filters.standarizeFilters(df)
            return df
        if self.tautomers or self.protonation:
            if self.neutralize: 
                df = Neutralizer.neutralize_df(df)
                df = Filters.remove_invalid_SMILES(df)
            df.loc[:, 'original_idx'] = df.index
            
        if self.tautomers: 
            tautomerizer = Tautomerizer(smartsFile = self.tautomer_library,
                                        taurdkit = self.taurdkit, 
                                        debug = self.debug, 
                                        neutralize = False,
                                        numcores = self.numcores) # Already neutralized
            df = tautomerizer.tautomerize_df(df)
        if self.pains: df = Filters.painsFilter(df,
                                                rejectedFile = rejected_file,
                                                debug = self.debug)
        if self.unwanted is not None: 
            df = Filters.unwantedFilter(df,
                                        rejectedFile = rejected_file,
                                        unwanted_option = self.unwanted,
                                        debug = self.debug)
            
        if self.custom is not None: 
            df = Filters.customFilter(df,
                                      rejectedFile = rejected_file,
                                      smartsFile = self.custom,
                                      debug = self.debug)
            
        if self.protonation: 
            ionizer = Ionizer(smartsFile = self.protonation_library,
                              pH = self.pH,
                              pH_range = self.pH_range,
                              numcores = self.numcores,
                              neutralize = False,
                              debug=self.debug)
            df = ionizer.ionize_df(df)
        if (self.tautomers or self.protonation) and not(df.empty):
            # Coalesce the 'smiles' column to ensure it is present
            # Drop duplicate rows based on 'smiles' and 'ids'
            df = df.drop_duplicates(subset=['smiles', 'ids'], keep='first')
            # Sort to match original order
            df = df.sort_values(by=['original_idx', 'ids']).reset_index(drop=True)
            # Find duplicate ids
            duplicated_ids = df['ids'][df['ids'].duplicated(keep=False)]

            if not duplicated_ids.empty:
                # Group duplicates and expand
                df.loc[duplicated_ids.index, 'ids'] = df.loc[duplicated_ids.index, 'ids'].groupby(df['ids']).transform(self.expand_ids)
            
        if self.stereoisomers: df = Msani.enum_stereoisomers(df, max_isomers=self.max_stereoisomers, debug=self.debug, numcores=self.numcores)
        if (not(self.protonation) and not(self.protonation) and not(self.stereoisomers)):
            # If no SMILES processing , just return a canonical SMILES of the input
            df.loc[:, 'smiles'] = df['mol'].apply(lambda x: Chem.MolToSmiles(x))
        return df
    

