import logging

from rdkit import Chem, RDLogger
from pandas import DataFrame, Series, concat  # only what you use
from pathlib import Path

from msani.filtering.filters import Filters, against_humanity, hold_up, loadSMARTSdata
from msani.moltransform.tautomerizer import Tautomerizer
from msani.moltransform.ionizer import Ionizer
from msani.moltransform.neutralizer import Neutralizer
from msani.moltransform.stereoisomers import Stereoisomerizer

logger = logging.getLogger('msani')
RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit

class Msani:
    """
    A class to store the filter options and conduct chemical modifications for MolSanitizer. 
    Initialize the class with the desired filter options and apply the filters to the input DataFrame.
    
    Examples
    ---------

    >>> processor = Msani(
                    removesalts=True,
                    ha='>10',
                    logp='<5',
                    tautomers=True, 
                    protonation = True, 
                    debug=True)
    >>> processed_df = processor.run(df)

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
                extended_tautomers = False,
                neutralize = False,
                taurdkit = True,
                stereoisomers = False,
                max_stereoisomers = 8,
                protonation = False,
                pH = 7,
                pH_range = 0,
                numcores = 1,
                randomSeed = 42,
                stereo_timeout = 60,
                standardize = False,
                protonation_library = None,
                tautomer_library = None,
                tpsa = None,
                fsp3 = None,
                debug = False
                ):

        self.removesalts = removesalts
        self.ha = str(ha) if ha is not None else None
        self.logp = str(logp) if logp is not None else None
        self.hba = str(hba) if hba is not None else None
        self.hbd = str(hbd) if hbd is not None else None
        self.mw = str(mw) if mw is not None else None
        self.tpsa = str(tpsa) if tpsa is not None else None
        self.fsp3 = str(fsp3) if fsp3 is not None else None

        self.custom = custom
        self.unwanted = [word.title() for word in unwanted if isinstance(word, str)] if unwanted is not None else None

        if self.custom is not None or self.unwanted is not None:
            if self.custom: temp_df_custom = loadSMARTSdata(self.custom)
            if self.unwanted:
                smartsFile = Path(__file__).parent / 'Data' / 'filter_out.txt'
                temp_df_unwanted = loadSMARTSdata(smartsFile.resolve(), self.unwanted)
            self.unwanted_df = concat([temp_df_custom, temp_df_unwanted]) if self.custom and self.unwanted \
                else temp_df_custom if self.custom else temp_df_unwanted
            logger.info(f"Loaded {len(self.unwanted_df)} SMARTS patterns for unwanted filtering.")
        else:
            self.unwanted_df = None

        self.pains = pains
        self.chiral = chiral
        
        self.tautomers = tautomers
        self.extended_tautomers = extended_tautomers
        self.neutralize = neutralize
        self.taurdkit = taurdkit
        self.stereoisomers = stereoisomers
        self.max_stereoisomers = max_stereoisomers
        self.protonation = protonation
        self.pH = pH
        self.pH_range = pH_range
        self.randomSeed = randomSeed
        self.stereo_timeout = stereo_timeout
        self.standardize = standardize    
        self.debug = debug
        self.numcores = numcores
        self.protonation_library = protonation_library
        self.tautomer_library = tautomer_library
    
    def __repr__(self):
        cls_name = self.__class__.__name__
        attrs = ', '.join(f'{k}={v!r}' for k, v in self.__dict__.items())
        return f'{cls_name}\n({attrs})'
    
    
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
        if self.tpsa is not None: df = Filters.filter_by_tpsa(df, self.tpsa, rejectedFile=rejected_file, debug=self.debug)
        if self.fsp3 is not None: df = Filters.filter_by_fsp3(df, self.fsp3, rejectedFile=rejected_file, debug=self.debug)
        if self.chiral is not None: df = Filters.filter_by_chiralcenters(df, self.chiral, rejectedFile=rejected_file, debug=self.debug)
        if self.standardize: 
            df = Filters.standarizeFilters(df)
            return df
        
        if self.neutralize: 
            df = Neutralizer.neutralize_df(df)
            df = Filters.remove_invalid_SMILES(df)

        if self.tautomers or self.protonation:
            df.loc[:, 'original_idx'] = df.index
        if self.tautomers: 
            tautomerizer = Tautomerizer(smartsFile = self.tautomer_library,
                                        taurdkit = self.taurdkit, 
                                        debug = self.debug, 
                                        neutralize = False,
                                        numcores = self.numcores,
                                        extended_tautomers = self.extended_tautomers) # Already neutralized
            df = tautomerizer.tautomerize_df(df)
        if self.pains: 
            if self.debug: print("Filtering PAINS")
            df = Filters.painsFilter(df,
                                                rejectedFile = rejected_file,
                                                debug = self.debug)
            
        # We already gathered the unwanted SMARTS patterns in the constructor
        if self.unwanted_df is not None: 
            if self.debug: print("Filtering Unwanted")
            df = Filters.unwantedFilter(df,
                                        rejectedFile = rejected_file,
                                        unwanted_df = self.unwanted_df,
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
            
        if self.stereoisomers: 
            stereoisomerizer = Stereoisomerizer(maxIsomers=self.max_stereoisomers,
                                  tryEmbedding=True,
                                  randomSeed=self.randomSeed,
                                  numcores=self.numcores,
                                  timeout=self.stereo_timeout,
                                  debug=self.debug)
            df = stereoisomerizer.enumerate_df(df)

        if (not(self.protonation) and not(self.tautomers) and not(self.stereoisomers)):
            # If no SMILES processing , just return a canonical SMILES of the input
            df.loc[:, 'smiles'] = df['mol'].apply(lambda x: Chem.MolToSmiles(x))
        return df
    
