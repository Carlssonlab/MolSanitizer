from contextlib import nullcontext
import argparse
import multiprocessing as mp
import platform
import logging
import time

from pandas import DataFrame, read_csv
from pathlib import Path
from rdkit import Chem, RDLogger

from msani.io.parsers import CustomHelpFormatter
from msani.moltransform.utils import enumerate_stereoisomers_corina

# forkserver avoids inheriting parent locks (RDKit allocator, logging) on Linux.
# Windows only supports 'spawn'; macOS prefers 'spawn' too (fork is deprecated).
_MP_START_METHOD = 'forkserver' if platform.system() == 'Linux' else 'spawn'

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
    numcores: int, default = 1.
        Number of CPU cores to use for parallel processing.
    timeout: int, default = 60.
        Maximum number of seconds to spend enumerating stereoisomers for a single molecule. If the timeout is
        reached, a partial set of stereoisomers may be returned; if no stereoisomers are generated (e.g., due to
        timeout or an error), the original input SMILES is kept.
    useCorina: bool, default = False.
        Use CORINA rather than RDKit to enumerate stereoisomers.
    corinaPath: str, optional.
        Path to the CORINA executable. Required when ``useCorina`` is true.
    debug: bool, default = False.
        Enable verbose debug output.

     Examples
    ----------
    
    >>> from msani.moltransform.stereoisomers import Stereoisomerizer
    >>> stereoisomerizer = Stereoisomerizer(maxIsomers = 8,
    ...                                     tryEmbedding=True,
    ...                                     unique = True,
    ...                                     randomSeed = 42,
    ...                                     numcores= 4)

    Enumerate stereoisomers for a molecule from a SMILES

    >>> isomers = stereoisomerizer.enumerate(smiles='CCC(O)C(O)CCC')
    >>> len(isomers)
    4

    Enumerate stereoisomers for a dataframe of molecules

    >>> isomers_df = stereoisomerizer.enumerate_df(df, smiles_column='smiles', name_column='ids')

    """
    def __init__(self, 
                maxIsomers: int = 0, 
                onlyUnassigned: bool = True, 
                onlyStereoGroups: bool = False, 
                unique: bool = True, 
                tryEmbedding: bool = False,
                randomSeed: int = -1,
                numcores: int = 1,
                timeout: int = 60,
                useCorina: bool = False,
                corinaPath: str = None,
                debug: bool = False
                ):

        self.options = {
            'maxIsomers': maxIsomers,
            'onlyUnassigned': onlyUnassigned,
            'onlyStereoGroups': onlyStereoGroups,
            'unique': unique,
            'tryEmbedding': tryEmbedding,
            'randomSeed': randomSeed,
            'timeout': float(timeout),
            'useCorina': useCorina,
            'corinaPath': corinaPath
        }
        self.numcores = numcores
        self.timeout = timeout
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
        if self.options['maxIsomers'] == 1: return [smiles]

        if self.options['useCorina']:
            stereoisomers = enumerate_stereoisomers_corina(
                smiles=smiles,
                max_isomers=self.options['maxIsomers'],
                onlyUnassigned=self.options['onlyUnassigned'],
                corina_path=self.options['corinaPath'],
                timeout=self.timeout,
            )
        else:
            stereoisomers = ms.enumerate_stereoisomers(smiles, self.options, self.debug)
        
        return stereoisomers
       
    def enumerate_df_mp(self, 
                        df: DataFrame,
                        smiles_column: str = 'smiles',
                        mol_column: str = 'mol',
                        name_column: str = 'ids',
                        *, pool=None) -> DataFrame:
        """
        Enumerate stereoisomers for a dataframe of molecules using multiprocessing Pool.
        """
        if mol_column not in df.columns:
            df = df.copy()
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        num_cores = min(self.numcores, len(df))
        ctx = mp.get_context(_MP_START_METHOD)
        # Package the row data efficiently using python dicts to avoid memory bloat of large Pandas Series
        keys = df.columns.tolist()
        rows = ((dict(zip(keys, r)), smiles_column, mol_column, name_column) for r in df.itertuples(index=False, name=None))
        
        results = []
        pool_context = nullcontext(pool) if pool is not None else ctx.Pool(
            processes=num_cores,
            initializer=_init_stereoisomerizer_worker,
            initargs=(self,))
        with pool_context as pool:
            # Ordered results keep the output order reproducible.
            for res in pool.imap(_process_single_stereoisomer_row, rows, chunksize=10):
                if res:
                    results.extend(res)
                    
        return DataFrame(results)

    def enumerate_df(self, 
                       df: DataFrame,
                       smiles_column: str = 'smiles',
                       mol_column: str = 'mol',
                       name_column: str = 'ids',
                       *, pool=None, parallel=True) -> DataFrame:
        """
        Enumerate stereoisomers for a dataframe of molecules.

        Parameters
        ----------
        df : DataFrame
            The dataframe containing molecules to enumerate stereoisomers for.
        smiles_column : str, default = 'smiles'
            The name of the column containing SMILES strings.
        mol_column : str, default = 'mol'
            The name of the column containing RDKit Mol objects. If this column does not exist, it will be created from the SMILES column.
        name_column : str, default = 'ids'
            The name of the column containing molecule identifiers.
        Returns
        -------
        DataFrame
            A new dataframe with the enumerated stereoisomers.
        """
        if df.empty:
            return df

        if mol_column not in df.columns:
            df = df.copy()
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        if parallel and self.numcores > 1:
            return self.enumerate_df_mp(df, smiles_column, mol_column, name_column, pool=pool)
        else:
            keys = df.columns.tolist()
            rows = ((dict(zip(keys, r)), smiles_column, mol_column, name_column) for r in df.itertuples(index=False, name=None))
            results = []

            for row_args in rows:
                res = _process_single_stereoisomer_row(row_args, worker=self)
                if res:
                    results.extend(res)
            return DataFrame(results)


# ---------------------------------------------------------------------------
# Module-level worker helpers for multiprocessing
# ---------------------------------------------------------------------------

_stereoisomerizer_worker = None

def _init_stereoisomerizer_worker(stereoisomerizer):
    """Initializer run once per worker process to prevent memory spikes from pickling."""
    global _stereoisomerizer_worker
    _stereoisomerizer_worker = stereoisomerizer

def _process_single_stereoisomer_row(args, worker=None):
    """
    Worker function to process stereoisomerization for a single row.
    Caught exceptions simply log an error and return a fallback list so the job continues.
    """
    row, smiles_column, mol_column, name_column = args
    worker = _stereoisomerizer_worker if worker is None else worker
    results = []
    
    mol_name = row[name_column]
    smiles = row[smiles_column]
    
    try:
        mol = row.get(mol_column, None)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles)
        
        longname = row.get('longname', None)
        original_idx = row.get('original_idx', None)

        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        unassigned = [idx for idx, tag in centers if tag == '?']
        num_possible_isomers = 2 ** len(unassigned)
        max_isomers = worker.options['maxIsomers']

        stereoisomers_smiles = worker.enumerate(smiles)

        # Fall back to the input SMILES when the selected enumerator returns no
        # isomers, so the molecule is not silently dropped.
        if not stereoisomers_smiles:
            logger.warning(
                f"{mol_name}: Stereoisomerization timed out, the input SMILES is kept."
            )
            return [{
                smiles_column: smiles,
                name_column: mol_name,
                mol_column: mol,
                'longname': longname,
                'original_idx': original_idx,
            }]

        if len(stereoisomers_smiles) == 1:
            results.append({
                smiles_column: stereoisomers_smiles[0],
                name_column: mol_name,
                mol_column: Chem.MolFromSmiles(stereoisomers_smiles[0]),
                'longname': longname,
                'original_idx': original_idx
            })
        else:
            if max_isomers > 0 and num_possible_isomers > max_isomers:
                logger.info(f"{mol_name}: Not all the stereoisomers are written out (capped at {max_isomers}/{num_possible_isomers}).")
                stereoisomers_smiles = stereoisomers_smiles[:max_isomers]
            
            two_digits = len(stereoisomers_smiles) >= 10
            for idx, stereoisomer in enumerate(stereoisomers_smiles):
                results.append({
                    smiles_column: stereoisomer,
                    mol_column: Chem.MolFromSmiles(stereoisomer),
                    name_column: mol_name + (f".{idx+1:02d}" if two_digits else f".{idx+1}"),
                    'longname': longname,
                    'original_idx': original_idx
                })

    except Exception as e:
        logger.error(
            f"Error processing stereoisomer for {mol_name}: {str(e)}"
        )
        results.append({
            name_column: mol_name,
            smiles_column: smiles,
            mol_column: row.get(mol_column),
            'longname': row.get('longname', None),
            'original_idx': row.get('original_idx', None)
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
        default=4,
        help="Number of CPU cores to use for parallel processing. Default is 4."
    )
    parser.add_argument(
        '--timeout',
        type=int,
        default=60,
        help="Per-molecule timeout in seconds. Default is 60."
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
        timeout=args.timeout,
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
        output_df[['smiles', 'name']].to_csv(args.output, index=False, sep=' ', header=False)
        print(f"Wrote {len(output_df)} stereoisomers to {args.output}")
    print(f"Time taken: {time.time() - start:.2f} seconds")
    
if __name__ == "__main__":
    main()
