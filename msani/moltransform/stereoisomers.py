import argparse
import multiprocessing as mp
import logging
import time
from functools import partial

from pandas import DataFrame, read_csv
from pathlib import Path
from rdkit import Chem, RDLogger

from msani.io.parsers import CustomHelpFormatter

RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit
logger = logging.getLogger('msani')

try:
    import msani_stereoisomers as ms
except ImportError:
    raise ImportError("Could not import the msani_stereoisomers module. Please ensure that the C++ extension has been built.")

class Stereoisomerizer:
    """
    A class for stereoisomerization of molecules.

    Parameters
    ----------
    maxIsomers : int, default = 0.
        The maximum number of isomers to yield. If the number of possible isomers exceeds this value,
        a random subset will be yielded. If 0, there is no maximum. 
    onlyUnassigned : bool, default = True.
        If true, stereocenters which have a specified stereochemistry will not be perturbed unless they are part of a relative stereo group.
    onlyStereoGroups : bool, default = False.
        If true, only find stereoisomers that differ at the StereoGroups associated with the molecule.
    unique : bool, default = True.
        If true, only stereoisomers that differ in canonical SMILES will be returned.
    tryEmbedding : bool, default = False.
        If true, the process attempts to generate a standard RDKit distance geometry conformation for
        the stereoisomer. If this fails, we assume that the stereoisomer is non-physical and don't return it.
    randomSeed: int, default = -1.
        Random seed for choosing a random subset of stereoisomers of a given compound.
    """
    def __init__(self, 
                maxIsomers: int = 0, 
                onlyUnassigned: bool = True, 
                onlyStereoGroups: bool = False, 
                unique: bool = True, 
                tryEmbedding: bool = False,
                randomSeed: int = -1,
                numcores: int = 1,
                debug: bool = False
                ):

        self.options = ms.StereoEnumerationOptions()
        self.options.maxIsomers = maxIsomers
        self.options.onlyUnassigned = onlyUnassigned
        self.options.onlyStereoGroups = onlyStereoGroups
        self.options.unique = unique
        self.options.tryEmbedding = tryEmbedding
        self.options.randomSeed = randomSeed
        self.numcores = numcores
        self.debug = debug
        


    def enumerate(self, smiles: str) -> list[str]:
        """
        Enumerate stereoisomers for a given SMILES string.
        Parameters
        ----------
        smiles : str
            The SMILES string of the molecule to enumerate stereoisomers for.
        Returns
        ------- 
        list[str]
            A list of SMILES strings representing the stereoisomers.
        """
        if self.options.maxIsomers == 1: return [smiles]

        stereoisomers = ms.enumerate_stereoisomers(smiles, self.options, self.debug)
        
        return stereoisomers
       
    def enumerate_df_mp(self, 
                          df: DataFrame,
                          smiles_column: str = 'smiles',
                          mol_column: str = 'mol',
                          name_column: str = 'ids') -> DataFrame:
        """
        """
        if mol_column not in df.columns:
            df.loc[:, mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        num_cores = min(self.numcores, len(df))  # Avoid using more cores than data chunks
        chunks = [df.iloc[i::num_cores] for i in range(num_cores)]

        # Extract parameters to pass instead of the C++ object
        stereo_params = {
            'maxIsomers': self.options.maxIsomers,
            'onlyUnassigned': self.options.onlyUnassigned,
            'onlyStereoGroups': self.options.onlyStereoGroups,
            'unique': self.options.unique,
            'tryEmbedding': self.options.tryEmbedding,
            'randomSeed': self.options.randomSeed,
            'debug': self.debug
        }

        process_func = partial(_process_stereoisomers_rows, 
                               stereo_params=stereo_params,
                               smiles_column=smiles_column,
                               mol_column=mol_column,
                               name_column=name_column)

        results = []
        with mp.Pool(processes = num_cores) as pool:
            async_results = [pool.apply_async(process_func, (chunk,)) for chunk in chunks]

            for async_result in async_results:
                try:
                    results.extend(async_result.get())  # Timeout for safety
                except Exception as e:
                    logger.error(f"Error processing a stereoisomer batch: {str(e)}")

        return DataFrame(results)

    def enumerate_df(self, 
                       df: DataFrame,
                       smiles_column: str = 'smiles',
                       mol_column: str = 'mol',
                       name_column: str = 'ids') -> DataFrame:
        """
        
        """
        if df.empty:
            return df
        if self.numcores > 1:
            return self.enumerate_df_mp(df, smiles_column, mol_column, name_column)
        else:
            # For single-core, pass the stereoisomerizer directly (no pickling needed)
            return DataFrame(_process_stereoisomers_rows(df,
                                                        stereoisomerizer=self,
                                                        smiles_column=smiles_column,
                                                        mol_column=mol_column,
                                                        name_column=name_column))

                       

def _process_stereoisomers_rows(df, stereoisomerizer=None, stereo_params=None, smiles_column='smiles', mol_column='mol', name_column='ids'):
    """
    Common function to process stereoisomerization for both single-core and multiprocessing.
    
    Parameters
    ----------
    df : DataFrame
        The dataframe chunk to process
    stereoisomerizer : Stereoisomerizer, optional
        The stereoisomerizer object (used for single-core processing)
    stereo_params : dict, optional
        Dictionary of parameters to reconstruct options (used for multiprocessing)
    """
    # If stereo_params is provided (multiprocessing), create a new stereoisomerizer
    if stereo_params is not None:
        options = ms.StereoEnumerationOptions()
        options.maxIsomers = stereo_params['maxIsomers']
        options.onlyUnassigned = stereo_params['onlyUnassigned']
        options.onlyStereoGroups = stereo_params['onlyStereoGroups']
        options.unique = stereo_params['unique']
        options.tryEmbedding = stereo_params['tryEmbedding']
        options.randomSeed = stereo_params['randomSeed']
        debug = stereo_params['debug']
    else:
        # Single-core processing: use the provided stereoisomerizer
        options = stereoisomerizer.options
        debug = stereoisomerizer.debug
    
    results = []
    if mol_column not in df.columns:
        df.loc[:, mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)
    for _, row in df.iterrows():
        smiles = row[smiles_column]
        mol = row[mol_column]
        highlights = row.get('highlights', None)
        original_idx = row.get('original_idx', None)  # For debugging purposes


        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        unassigned = [idx for idx, tag in centers if tag == '?']
        num_possible_isomers = 2 ** len(unassigned)
        max_isomers = options.maxIsomers


        stereoisomers_smiles = ms.enumerate_stereoisomers(smiles, options, debug)
        if len(stereoisomers_smiles) == 1:
            results.append({smiles_column: stereoisomers_smiles[0],
                            name_column: row[name_column],
                            mol_column: Chem.MolFromSmiles(stereoisomers_smiles[0]),
                            'highlights': highlights,
                            'original_idx': original_idx
                            })
        else:
            if max_isomers > 0 and num_possible_isomers > max_isomers:
                logger.info(f"{row[name_column]}: Not all the stereoisomers are written out (capped at {max_isomers}/{num_possible_isomers}).")
                stereoisomers_smiles = stereoisomers_smiles[:max_isomers]
            two_digits = len(stereoisomers_smiles) >= 10
            for _, stereoisomer in enumerate(stereoisomers_smiles):
                results.append({ smiles_column: stereoisomer,
                                mol_column: Chem.MolFromSmiles(stereoisomer),
                                name_column: row[name_column] + (f"_{_+1:02d}" if two_digits else f"_{_+1}"),                         
                                'highlights': highlights,
                                'original_idx': original_idx
                                })

    return results

def main():
    start = time.time()
    parser = argparse.ArgumentParser(
        description="Generate stereoisomers for a list of molecules.",
        formatter_class=CustomHelpFormatter
    )
    parser.add_argument(
        "-i", "--input", 
        type=str, 
        help="Input CSV file containing a column of SMILES strings."
    )
    parser.add_argument(
        '-s', '--smiles',
        type=str,
        help="SMILES string to process."
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        help="Output CSV file to write the stereoisomers to.",
       default="stereoisomers.csv"
    )
    parser.add_argument(
        '--maxIsomers',
        type=int,
        default=0,
        help="Maximum number of isomers to yield. If the number of possible isomers exceeds this value, a random subset will be yielded. If 0, there is no maximum."
    )
    parser.add_argument(
        '--onlyUnassigned',
        action='store_true',
        help="If true, stereocenters which have a specified stereochemistry will not be perturbed unless they are part of a relative stereo group. Default is True.",
        default=True
    )
    parser.add_argument(
        '--onlyStereoGroups',
        action='store_true',
        help="If true, only find stereoisomers that differ at the StereoGroups associated with the molecule. Default is False.",
        default=False
    )
    parser.add_argument(
        '--unique',
        action='store_true',
        help="If true, only stereoisomers that differ in canonical SMILES will be returned. Default is True.",
        default=True
    )
    parser.add_argument(
        '--tryEmbedding',
        action='store_true',
        help="If true, the process attempts to generate a standard RDKit distance geometry conformation for the stereoisomer. If this fails, we assume that the stereoisomer is non-physical and don't return it. Default is False.",
        default=False
    )
    parser.add_argument(
        '--randomSeed',
        type=int,
        default=-1,
        help="Random seed for reproducibility. Default is -1 (no seed)."
    )
    parser.add_argument(
        '--numcores',
        '-j',
        type=int,
        default=1,
        help="Number of CPU cores to use for parallel processing. Default is 1."
    )
    parser.add_argument(
        '--debug',
        '-d',
        action='store_true',
        help="If true, enable debug mode for more verbose output. Default is False.",
        default=False
    )
    args = parser.parse_args()
    if args.smiles is None and args.input is None:
        parser.error("Either --smiles or --input must be provided.")

    stereoisomerizer = Stereoisomerizer(
        maxIsomers=args.maxIsomers,
        onlyUnassigned=args.onlyUnassigned,
        onlyStereoGroups=args.onlyStereoGroups,
        unique=args.unique,
        tryEmbedding=args.tryEmbedding,
        randomSeed=args.randomSeed,
        numcores=args.numcores,
        debug=args.debug

    )

    if args.smiles is not None:
        # Process the single SMILES string
        stereoisomers = stereoisomerizer.enumerate(args.smiles)
        print(f"Generated {len(stereoisomers)} stereoisomers for {args.smiles}")
        for smi in stereoisomers:
            print(smi)
    else:
        # Process the input CSV file
        input_path = Path(args.input)
        if not input_path.is_file():
            raise FileNotFoundError(f"Input file {args.input} does not exist.")

        df = read_csv(input_path, sep = r'\s+', header = None, names = ['smiles', 'name'])
        if 'smiles' not in df.columns:
            raise ValueError("Input CSV must contain a 'smiles' column.")

        output_df = stereoisomerizer.enumerate_df(df, smiles_column='smiles', name_column='name', mol_column='mol')
        output_df.to_csv(args.output, index=False)
        print(f"Wrote {len(output_df)} stereoisomers to {args.output}")
    print(f"Time taken: {time.time() - start:.2f} seconds")
    
if __name__ == "__main__":
    main()