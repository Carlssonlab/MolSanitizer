"""
The MolSanitizer C++ Stereoisomer Enumerator

This module provides high-performance stereoisomer enumeration based on the rdEnumerateStereoisomers of RDKit version >= 2025.3.4.
"""

from typing import Dict, List

__author__: str

class StereoEnumerationOptions:
    """
    EnumerateStereoisomers options.

    NOTE: This class is kept for reference but is not directly used.
    Use Python dictionaries with enumerate_stereoisomers() instead.

    Attributes
    ----------
    tryEmbedding : bool
        If true, the process attempts to generate a standard RDKit distance geometry
        conformation for the stereoisomer. If this fails, we assume that the stereoisomer is
        non-physical and don't return it. NOTE that this is computationally expensive and is
        just a heuristic that could result in stereoisomers being lost. Default=False
    onlyUnassigned : bool
        If true, stereocenters which have a specified stereochemistry will not be
        perturbed unless they are part of a relative stereo group. Default=True.
    onlyStereoGroups : bool
        If true, only find stereoisomers that differ at the StereoGroups associated with
        the molecule. Default=False.
    unique : bool
        If true, only stereoisomers that differ in canonical SMILES will be
        returned. Default=True.
    maxIsomers : int
        The maximum number of isomers to yield. If the number of possible isomers
        is greater than maxIsomers, a random subset will be yielded. If 0, there
        is no maximum. Since every additional stereocenter doubles the number of
        results (and execution time) it's important to keep an eye on this.
    randomSeed : int
        Seed for random number generator. Default=-1 means no seed.
    timeout : float
        Wall-clock timeout in seconds. 0 = no limit. Default=0
    """

    tryEmbedding: bool
    onlyUnassigned: bool
    onlyStereoGroups: bool
    unique: bool
    maxIsomers: int
    randomSeed: int
    timeout: float

    def __init__(self) -> None: ...

def enumerate_stereoisomers(
    smiles: str,
    params: Dict[str, any],
    verbose: bool = False,
) -> List[str]:
    """
    Enumerate all stereoisomers of a molecule given its SMILES and a dictionary of options.

    Parameters
    ----------
    smiles : str
        The SMILES string of the molecule.
    params : dict
        Dictionary with the following optional keys:

        - maxIsomers (int): Maximum number of isomers to yield. Default=0 (no limit)
        - onlyUnassigned (bool): Only enumerate unassigned stereocenters. Default=True
        - onlyStereoGroups (bool): Only enumerate StereoGroups. Default=False
        - unique (bool): Return only unique stereoisomers. Default=True
        - tryEmbedding (bool): Validate stereoisomers via embedding. Default=False
        - randomSeed (int): Random seed for subset selection. Default=-1
        - timeout (float): Wall-clock timeout in seconds. 0 = no limit. Default=0

    verbose : bool, optional
        Enable verbose output. Default=False

    Returns
    -------
    list of str
        List of SMILES strings representing the stereoisomers.
    """
    ...
