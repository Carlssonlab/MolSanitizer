#!/usr/bin/env python3
"""
Connect dissociated mononuclear metal complexes -- a metal and its ligands
written as separate fragments -- into single structures with dative bonds.

Scope
    One metal centre per record, the metal being Sc-As, Y-Sb or La-Bi
    (METAL_ATOMIC_NUMS). Records without a metal and complexes already drawn
    as one fragment pass through unchanged. Polynuclear and multi-centre
    records, and records that cannot be built, are rejected with a reason
    (REJECT_*) rather than altered.

Method, per record (prepare_complex)
    1. Triage. A string prefilter (has_metal) skips organic SMILES without
       parsing them; per-fragment metal counts tell a cluster apart from
       several mononuclear centres in one row.
    2. Spectators. The listed spectator ions (SPECTATOR_SMILES: Li+, Na+,
       K+, Rb+, Cs+, NH4+, BF4-, PF6-) are removed. No other fragment is
       dropped, except by strict mode below.
    3. Ligand preparation (process_step2), per ligand: S/P ylides are
       redrawn with double bonds, charges are neutralized except deliberate
       ones (NEUTRALIZATION_EXCEPTIONS: hydroxide, oxide, charged carbon,
       charge-separated pairs), and COMPONENT_REACTIONS re-ionize the groups
       that bind as anions (acids, phenols, thiols, HX, HCN, acidic N-H).
       Every edit is recorded; metal charges are never inferred.
    4. Donors (find_donor_candidates). DONOR_SMARTS, in priority order;
       WEAK_DONOR_SMARTS (ether and carbonyl O) only when a ligand has no
       stronger donor.
    5. Conflicts (classify_donor_conflict). Two donors of one ligand
       conflict when they would close a chelate ring smaller than
       MIN_CHELATE_RING -- linkage alternatives such as S- versus N-bound
       thiocyanate -- or when more than one bond of the path between them
       lies in a small rigid ring (denticity). Bipyridine-type N-N paths,
       rings of MACROCYCLE_RING_CUTOFF atoms or more and eta-5
       cyclopentadienide rings do not count against binding together.
    6. Soft mode lists every chemically plausible structure. Per ligand it
       keeps each maximal conflict-free donor set: all linkage alternatives,
       but only the largest, best-priority set when denticity is ambiguous.
       Ligand choices are then combined (permutations of identical ligands
       skipped), capped by the coordination limit, max_combinations and
       max_variants; each combination is one output.
    7. Strict mode returns one structure. It fills the coordination sphere
       once for all ligands together, priority tier by priority tier, up to
       the largest reachable coordination number (from
       metal_coordination_numbers.txt, else coordination_limit). Ties are
       built and the lexicographically first canonical SMILES is kept;
       ligands left without a seat are dropped. Both are reported.
    8. Build (connect_ligands). Dative donor-to-metal bonds, sanitization
       and stereo re-perception; every SMILES must parse back as a single
       fragment. Candidates are then checked against the tabulated
       coordination numbers for the metal and charge, where there is an
       entry; a result that reaches none is kept but flagged.

Performance
    Conflict-free donor subsets are found on bitmasks (Bron-Kerbosch for
    maximal sets, a memoized exact search for the largest ones) instead of
    trying all 2**n subsets. Strict mode remembers each conflict
    classification per ordered atom pair and stops as soon as
    max_combinations is exceeded. Strict branches that a symmetry of the
    ligands maps onto each other -- identical ligand copies, or a ligand's
    own stereo-preserving symmetry -- are built once. SMILES are written and
    parsed only where the result is new.

Entry points
    prepare_complex for one record; partition_dataframe for Msani
    DataFrames, used by msani.api.Msani (--metal strict|soft|off);
    process_file and main for standalone use of this file.
"""

import itertools
from dataclasses import dataclass, field
import os
import re
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# ============================================================================
# CONFIGURATION
# ============================================================================

# Single source of truth: Sc-As, Y-Sb, La-Bi. METAL_SMARTS and the string
# prefilter below are both derived from it so they can never drift apart.
METAL_ATOMIC_NUMS = tuple(range(21, 34)) + tuple(range(39, 52)) + tuple(range(57, 84))

METAL_SMARTS = '[' + ','.join(f'#{n}' for n in METAL_ATOMIC_NUMS) + ']'
METAL_PATTERN = Chem.MolFromSmarts(METAL_SMARTS)

_PTABLE = Chem.GetPeriodicTable()
METAL_SYMBOLS = [_PTABLE.GetElementSymbol(n) for n in METAL_ATOMIC_NUMS]

# A metal can only ever appear inside square brackets: SMILES lets the organic
# subset (B, C, N, O, P, S, F, Cl, Br, I and the aromatic forms) be written
# bare, and nothing else. So the symbol always sits directly after '[' plus an
# optional isotope. The negative lookahead stops 'Y' matching the 'Y' of 'Yb'.
_METAL_RE = re.compile(
    # RDKit also accepts aromatic arsenic, e.g. [asH]1cccc1.
    r'\[\d*(?:' + '|'.join(METAL_SYMBOLS + ['as']) + r')(?![a-z])'
)


def has_metal(smiles):
    """Cheap, RDKit-free test for a metal in a SMILES string.

    Conservative by construction: it never returns False for a SMILES that
    really does contain one of METAL_ATOMIC_NUMS, so it is safe as a skip
    gate. It may return True for a string RDKit later rejects, which only
    costs one wasted parse.

    The bracket check avoids running a regex on ordinary organic SMILES.
    """
    return '[' in smiles and _METAL_RE.search(smiles) is not None


COMPONENT_REACTIONS = [
    '[AH1&+0;$([O,S]-A=[O,S])!$([O,S]~A~[-1]):1]>>[AH0&-:1]',
    '[#9H1&+0,#17H1&+0,#35H1&+0,#53H1&+0:1]>>[*H0&-:1]',
    '[AH1&+0;$([O,S]-[*;$(A#N),$(a)])!$([O,S]~*~[-1])!$([O,S]~[-1]):1]>>[AH0&-:1]',
    '[NH1&+0;$(N-,=[!#6]=,#[!#6]):1]>>[NH0&-:1]',
    '[CH1&+0;$(C#N):1]>>[CH0&-:1]',
    '[SH1X2&+0:1]>>[SH0&-:1]',
]

# Atoms a chemist would write charged on purpose, because the charged form
# cannot arise from ordinary aqueous (de)protonation -- so an explicit
# charge here is a deliberate statement, not a default-state shortcut, and
# neutralize_ligand must never touch them. Everything NOT covered by one of
# these (or restored by COMPONENT_REACTIONS below) gets neutralized: this is
# a real, growable exception list, not a closed one -- add to it as new
# genuinely-deliberate charged motifs turn up.
NEUTRALIZATION_EXCEPTIONS = [
    '[OH1&-1]',   # hydroxide: OH- cannot form by simple deprotonation of water
    '[O&-2]',     # oxide: same reasoning, one step further
    '[#6&!+0]',   # any charged carbon, cationic or anionic: too exotic to second-guess
    # Directly bonded opposite charges (nitro, N-oxide, azide, C#[O+]) are a
    # charge-separated representation of a neutral group, not a protonation
    # state: neutralizing only one side would leave a net-charged species.
    '[*&+{1-4}]~[*&-{1-4}]',
]

# S/P ylide forms are standardized to their double-bonded form BEFORE
# neutralization (sulfoxide, sulfone, sulfoximine, phosphine oxide...).
YLIDE_STANDARDIZATION_REACTIONS = [
    '[SX3&+,SX4&+,PX4&+:1]-[*-:2]>>[S&+0,P&+0:1]=[*&+0:2]',
    '[SX4&+2:1](-[*-:2])-[*-:3]>>[S&+0:1](=[*&+0:2])=[*&+0:3]',
]

DONOR_SMARTS = [
    '[O&-2]',
    '[#53&-]',
    '[#35&-]',
    '[S&-2]',
    '[#16&-]',
    '[#17&-]',
    '[N&-]',
    '[#9&-]',
    '[O&-]',
    '[n&+0;!$(n~*=[!#6])]',
    '[N&+0;!$(N~*=[!#6])!$(N#C)]',
    '[O&+0&H2]',
    '[OX2&+0&H1]',
    '[PX3&+0]',
    '[#6&-]',
    '[cH&+0;$(c1[c-][c+0][c+0][c+0]1),$(c1[c+0][c-][c+0][c+0]1)]',
    '[N&+0;$(N#C)]',
    '[S;$([SX3&+0]=O)]',
    '[O;$(O=[SX3&+0])]',
]

# Last-resort donors: used only when the fragment has no donor from the list
# above (see WEAK_DONOR_PRIORITY below) -- a neutral ether or ketone oxygen
# is a real ligand (THF, DMF, acetone) but far weaker than anything above,
# and without this split Auranofin's ring ether would compete with its
# thiolate. A SEPARATE list, not a priority number to remember: appending
# here, or to DONOR_SMARTS above, keeps the threshold correct automatically.
WEAK_DONOR_SMARTS = [
    '[O;$(O=[CX3])!$(O=C-O)]',
    '[O;$([OX2&+0&H0])!$(O-C=O)!$([OX2]-[SX4])!$([OX2]-[PX4])]',
]

# Limits bound both per-ligand subset search and whole-record enumeration.
MAX_DONORS = 16
MAX_COMBINATIONS = 128
MAX_VARIANTS = 16

# These are spectator ions, not a general salt/solvent removal catalogue.
SPECTATOR_SMILES = frozenset(Chem.CanonSmiles(s) for s in (
    '[Li+]', '[Na+]', '[K+]', '[Rb+]', '[Cs+]', '[NH4+]',
    'F[B-](F)(F)F', 'F[P-](F)(F)(F)(F)F',
))
# A SMILES spells out every atom of its fragment, so a fragment with more atoms
# than the largest spectator can never match; prepare_complex uses this to skip
# writing SMILES for every real ligand just to rule it out.
SPECTATOR_MAX_ATOMS = max(Chem.MolFromSmiles(s).GetNumAtoms() for s in SPECTATOR_SMILES)


class EnumerationLimit(ValueError):
    pass


class NeedsReview(str):
    """A note on a kept result that a person should check, such as an
    arbitrary tie-break. It is an ordinary string everywhere; Msani logs
    these notes as warnings and every other note as info."""


# ============================================================================
# STEP 1
# ============================================================================


def read_smiles_file(filename):
    molecules = []
    n_skipped = 0
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                print(
                    f"Warning: Line {line_num} has insufficient columns, skipping"
                )
                continue
            smiles, name = parts[0], parts[1]
            if not has_metal(smiles):
                n_skipped += 1
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Warning: Could not parse SMILES on line {line_num}")
                continue
            molecules.append((mol, name, smiles))

    print(f"Read {len(molecules)} molecules from '{filename}' "
          f"({n_skipped} skipped by string prefilter, never parsed)")
    return molecules


def metals_per_fragment(mol):
    """Metal count for each fragment, plus the flat list of metal matches.

    Counting metals across the whole record conflates a real cluster with two
    separate complexes that happen to share a row: [Mn][Mn] is polynuclear,
    but Cl[Pt]Cl.Cl[Pd]Cl is two mononuclear centres.
    """
    matches = mol.GetSubstructMatches(METAL_PATTERN)
    metal_idxs = {m[0] for m in matches}
    counts = [
        len(metal_idxs.intersection(frag))
        for frag in Chem.GetMolFrags(mol, asMols=False)
    ]
    return counts, matches


# ============================================================================
# STEP 2
# ============================================================================

COMPONENT_REACTION_OBJECTS = tuple(AllChem.ReactionFromSmarts(s) for s in COMPONENT_REACTIONS)
YLIDE_REACTION_OBJECTS = tuple(AllChem.ReactionFromSmarts(s) for s in YLIDE_STANDARDIZATION_REACTIONS)
NEUTRALIZATION_EXCEPTION_PATTERNS = tuple(Chem.MolFromSmarts(s) for s in NEUTRALIZATION_EXCEPTIONS)


def neutralize_ligand(comp, changes=None):
    """Protonate every charged atom back toward neutral, one unit of charge
    at a time, except atoms matching a NEUTRALIZATION_EXCEPTIONS pattern --
    those are left exactly as written.

    This replaces guessing at which specific notations (Pt ammine anions,
    or anything else) might mean "neutral really" on a case-by-case basis:
    every ligand gets the same treatment, and only the exception list -- a
    chemistry judgment about which charged forms cannot arise from ordinary
    (de)protonation -- decides what survives up front. Anything neutralized
    here that should genuinely be charged is expected to be restored by
    COMPONENT_REACTIONS immediately afterward; a charged motif covered by
    neither the exceptions nor a reaction rule will end up incorrectly
    neutral, which is a gap in that rule catalogue to fill, not a case for
    a new exception here.

    Charge-separated neutral groups (nitro, N-oxide, azide, C#[O+]) are
    protected as a pair by the directly-bonded-opposite-charges entry in
    NEUTRALIZATION_EXCEPTIONS; S/P ylide forms were already standardized
    upstream in process_step2.

    As a safety net, each remaining charged atom is still neutralized and
    validated independently (a trial edit, sanitized on its own, kept only
    if valid): a charge that cannot be removed by a valid neutral structure
    is not a protonation site and is left untouched.
    """
    # Nothing charged means no trial edit below could ever run.
    if not any(atom.GetFormalCharge() for atom in comp.GetAtoms()):
        return comp
    exception_atoms = set()
    for pattern in NEUTRALIZATION_EXCEPTION_PATTERNS:
        for match in comp.GetSubstructMatches(pattern):
            exception_atoms.update(match)

    working = Chem.Mol(comp)
    touched = False
    for atom_idx in range(working.GetNumAtoms()):
        if atom_idx in exception_atoms:
            continue
        if working.GetAtomWithIdx(atom_idx).GetFormalCharge() == 0:
            continue
        trial = Chem.RWMol(working)
        atom = trial.GetAtomWithIdx(atom_idx)
        while atom.GetFormalCharge() != 0:
            step = -1 if atom.GetFormalCharge() > 0 else 1
            atom.SetFormalCharge(atom.GetFormalCharge() + step)
            atom.SetNumExplicitHs(max(0, atom.GetTotalNumHs() + step))
        trial_mol = trial.GetMol()
        try:
            Chem.SanitizeMol(trial_mol)
        except Exception:
            continue
        working = trial_mol
        touched = True
    if not touched:
        return comp
    if changes is not None:
        # comp itself is never edited, so its SMILES can wait until needed.
        before, after = Chem.MolToSmiles(comp), Chem.MolToSmiles(working)
        if after != before:
            changes.append(f'neutralized (exceptions protected): {before} -> {after}')
    return working


def combine_components(components):
    result = Chem.Mol()
    for component in components:
        result = Chem.CombineMols(result, component)
    return result


def prepare_ligand(component):
    """Copy and canonically order a ligand without changing its chemistry."""
    component = Chem.Mol(component)
    # Canonical ordering makes reaction-site and donor tie choices reproducible.
    ranks = Chem.CanonicalRankAtoms(component, breakTies=True)
    return Chem.RenumberAtoms(component, sorted(range(len(ranks)), key=ranks.__getitem__))


def _apply_reactions(comp, reactions, label, changes=None):
    """Apply each reaction repeatedly until it no longer changes the ligand."""
    comp_smi = None  # SMILES of the current comp, written once per molecule
    for i, reaction in enumerate(reactions, 1):
        for _ in range(comp.GetNumAtoms() + 1):
            products = reaction.RunReactants((comp,), maxProducts=1)
            if not products:
                break
            new_comp = products[0][0]
            Chem.SanitizeMol(new_comp)
            if comp_smi is None:
                comp_smi = Chem.MolToSmiles(comp)
            new_smi = Chem.MolToSmiles(new_comp)
            if comp_smi == new_smi:
                break
            comp = new_comp
            if changes is not None:
                changes.append(f'{label} {i}: {comp_smi} -> {new_smi}')
            comp_smi = new_smi
    return comp


def process_step2(mol, changes=None):
    """Prepare free ligands and record every protonation edit.

    Per ligand: standardize S/P ylide forms, neutralize (exceptions
    protected), then re-ionize via COMPONENT_REACTIONS.

    This deliberately does not infer oxidation states from formal metal charges.
    Mixed covalent/dative charge conventions are flagged by prepare_complex.
    """
    final = _prepare_ligands(mol, changes)
    return final, Chem.MolToSmiles(final)


def _prepare_ligands(mol, changes=None):
    """process_step2 without the result's SMILES, which prepare_complex never reads."""
    processed = []
    for comp in Chem.GetMolFrags(mol, asMols=True):
        if comp.HasSubstructMatch(METAL_PATTERN):
            processed.append(comp)
            continue
        comp = prepare_ligand(comp)
        comp = _apply_reactions(comp, YLIDE_REACTION_OBJECTS, 'ylide rule', changes)
        comp = neutralize_ligand(comp, changes=changes)
        comp = _apply_reactions(comp, COMPONENT_REACTION_OBJECTS, 'ligand rule', changes)
        processed.append(comp)
    return combine_components(processed)


# ============================================================================
# STEP 3
# ============================================================================

# Priority is the position in DONOR_SMARTS + WEAK_DONOR_SMARTS (1-indexed);
# WEAK_DONOR_PRIORITY is wherever the weak section actually starts, so
# editing either list can never desync it from a hand-written number again.
WEAK_DONOR_PRIORITY = len(DONOR_SMARTS) + 1

DONOR_PATTERNS = []
for priority, smarts in enumerate(DONOR_SMARTS + WEAK_DONOR_SMARTS, start=1):
    q = Chem.MolFromSmarts(smarts)
    if q:
        DONOR_PATTERNS.append((priority, q))


def find_donor_candidates(component):
    """Find donors without changing hydrogen counts or ligand stereochemistry."""
    component = prepare_ligand(component)
    # Fragments from Chem.GetMolFrags(..., sanitizeFrags=False) do not have
    # RingInfo populated. bond.IsInRing() happens to trigger it lazily as a
    # side effect (which is the only reason Rule 2's ring-size logic has
    # ever worked), but anything that queries RingInfo directly -- like
    # _same_haptic_ring -- without going through IsInRing() first will
    # silently see no rings at all instead of failing loudly. Make this
    # explicit here, once, so every downstream conflict check can rely on
    # RingInfo being present regardless of call order.
    Chem.FastFindRings(component)
    candidate_priorities = {}

    for priority, query in DONOR_PATTERNS:
        for match in component.GetSubstructMatches(query):
            atom_idx = match[0]
            if atom_idx not in candidate_priorities:
                candidate_priorities[atom_idx] = priority

    # Last-resort donors are discarded outright if anything stronger is present.
    if any(p < WEAK_DONOR_PRIORITY for p in candidate_priorities.values()):
        candidate_priorities = {
            i: p for i, p in candidate_priorities.items()
            if p < WEAK_DONOR_PRIORITY
        }

    candidates = []
    for atom_idx, priority in sorted(
        candidate_priorities.items(), key=lambda x: (x[1], x[0])
    ):
        atom = component.GetAtomWithIdx(atom_idx)

        candidates.append({
            'atom_idx': atom_idx,
            'priority': priority,
            'atomic_num': atom.GetAtomicNum(),
            'formal_charge': atom.GetFormalCharge(),
            'num_h': atom.GetTotalNumHs(),
        })

    return component, candidates


# A conservative preference, not a claim that smaller chelate rings cannot
# exist. Small-ring binding (e.g. bidentate carboxylate) remains out of scope.
# No upper bound: distant donor pairs occur in valid polydentate ligands.
MIN_CHELATE_RING = 5

CONFLICT_LINKAGE = 'linkage'      # alternative candidate attachment modes
CONFLICT_DENTICITY = 'denticity'  # ring geometry: which subset binds is a guess


# Rings at or above this size are treated as flexible macrocyclic scaffolds
# (cyclen/cyclam/crown-ether/porphyrin-type ligands): donors sitting at
# different points around such a ring are exactly the ones meant to bind the
# same metal simultaneously, so bonds belonging to a ring this size or larger
# must not count as evidence that two donors are held apart. The cutoff is a
# judgment call, not something derivable from the graph -- 10 separates
# typical small chelate/aromatic rings (5-6-membered, occasionally fused into
# 7-9-membered bicyclics) from typical polyaza/polyoxa macrocycles.
MACROCYCLE_RING_CUTOFF = 10


def _is_constraining_ring_bond(mol, bond, cutoff=MACROCYCLE_RING_CUTOFF):
    """True if this bond's rigidity should count against simultaneous binding.

    A bond only disqualifies itself here if it belongs EXCLUSIVELY to rings
    smaller than the cutoff. If it is part of any macrocycle-sized ring --
    even a bond also shared with a smaller fused ring, such as a pyridine
    fused into a larger aza-macrocycle -- the macrocycle's own flexibility
    governs, so the bond is exempted. Using the largest of a bond's ring
    memberships (rather than its smallest, i.e. SMARTS "r") is what makes
    this work for fused rings.
    """
    if not bond.IsInRing():
        return False
    ring_info = mol.GetRingInfo()
    sizes = ring_info.BondRingSizes(bond.GetIdx())
    return max(sizes) < cutoff


def _same_haptic_ring(mol, idx_a, idx_b):
    """True if idx_a and idx_b both belong to the same aromatic, all-carbon,
    net-monoanionic ring -- a cyclopentadienide-type pi-donor system, where
    every ring member is meant to bind the metal simultaneously (eta-5
    hapticity) rather than compete for one coordination site the way
    ordinary chelate-ring geometry would otherwise suggest.

    Purely structural, not tied to which DONOR_SMARTS pattern matched
    either atom, so it stays correct if the donor priority list is
    reordered or extended later. Scoped tightly to the ferrocene/
    cyclopentadienide case as it exists today (5-membered, all-carbon,
    exactly one formal negative charge); a different haptic motif -- an
    eta-6 arene sandwich, say, which is net-neutral rather than anionic --
    would need its own, separately-reasoned check, not a widening of this
    one.
    """
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 5 or idx_a not in ring or idx_b not in ring:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetAtomicNum() == 6 and a.GetIsAromatic() for a in atoms):
            continue
        charges = [a.GetFormalCharge() for a in atoms]
        if sum(charges) == -1 and charges.count(-1) == 1 and all(c in (0, -1) for c in charges):
            return True
    return False


def classify_donor_conflict(comp, idx_a, idx_b):
    """Why two donors cannot bind the same metal, or None if they can.

    Returns CONFLICT_LINKAGE, CONFLICT_DENTICITY or None. The caller treats the
    two kinds differently: linkage conflicts are enumerated as separate output
    rows, denticity conflicts are resolved down to a single answer.
    """
    if _same_haptic_ring(comp, idx_a, idx_b):
        return None

    path = Chem.GetShortestPath(comp, idx_a, idx_b)
    if len(path) < 2:
        return None

    bonds = [
        comp.GetBondBetweenAtoms(path[i], path[i + 1])
        for i in range(len(path) - 1)
    ]

    # Rule 1: prefer alternative attachment for small prospective chelate rings.
    # The ring is the path atoms plus the metal, hence len(path) + 1.
    if len(path) + 1 < MIN_CHELATE_RING:
        return CONFLICT_LINKAGE

    # Recognize the short aromatic N-C-...-C-N paths of bipy/phenanthroline.
    # Ring membership alone must not exclude these common chelators.
    if (len(path) in (4, 5)
            and all(comp.GetAtomWithIdx(i).GetIsAromatic() for i in path)
            and comp.GetAtomWithIdx(idx_a).GetAtomicNum() == 7
            and comp.GetAtomWithIdx(idx_b).GetAtomicNum() == 7
            and all(comp.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in path[1:-1])):
        return None

    # Rule 2: more than one bond of the path lies in a ring small enough to
    # hold the donors rigidly apart. Bonds that only belong to a macrocycle
    # (see _is_constraining_ring_bond) are not counted: that ring is exactly
    # what lets several of its donors reach the same metal at once.
    if sum(1 for b in bonds if _is_constraining_ring_bond(comp, b)) > 1:
        return CONFLICT_DENTICITY

    return None


class _ConflictMemo(dict):
    """classify_donor_conflict for one ligand, remembered per atom pair.

    Strict mode asks about the same pairs again for every branch and tier,
    and once more for the unrestricted maximum. The key keeps the argument
    order: when several shortest paths exist, GetShortestPath can pick a
    different one for (b, a) than for (a, b), and the kind follows the path.
    """

    def __init__(self, comp):
        super().__init__()
        self.comp = comp

    def __missing__(self, pair):
        kind = self[pair] = classify_donor_conflict(self.comp, *pair)
        return kind


# Conflict-free donor subsets on bitmasks: candidates are numbered by list
# position, and bit j of adj[i] is set when candidates i and j conflict. Both
# searches return exactly what scanning every subset would, without visiting
# all 2**n of them -- sugars and cyclodextrins easily carry 20+ same-priority
# -OH donors, where that scan takes seconds to minutes.

def _bit_indices(mask):
    """Positions of the set bits of mask, lowest first."""
    return [i for i in range(mask.bit_length()) if mask >> i & 1]


def _maximal_conflict_free_masks(adj):
    """Every maximal conflict-free subset, as bitmasks, in no particular order.

    Bron-Kerbosch with pivoting on the compatibility graph, whose maximal
    cliques are exactly the maximal conflict-free subsets.
    """
    full = (1 << len(adj)) - 1
    compatible = [full & ~(conflicts | 1 << i) for i, conflicts in enumerate(adj)]
    found = []

    def expand(chosen, pool, excluded):
        if not pool:
            if not excluded:
                found.append(chosen)
            return
        pivot = max(_bit_indices(pool | excluded),
                    key=lambda i: (compatible[i] & pool).bit_count())
        for i in _bit_indices(pool & ~compatible[pivot]):
            expand(chosen | 1 << i, pool & compatible[i], excluded & compatible[i])
            pool &= ~(1 << i)
            excluded |= 1 << i

    expand(0, full, 0)
    return found


def _max_conflict_free_size(mask, adj, memo):
    """Size of the largest conflict-free subset of the candidates in mask."""
    if not mask:
        return 0
    if mask in memo:
        return memo[mask]
    # Grow the set of candidates linked to the lowest one through conflicts;
    # if that is not all of mask, the two parts are independent problems.
    part = frontier = mask & -mask
    while frontier:
        low = frontier & -frontier
        frontier ^= low
        linked = adj[low.bit_length() - 1] & mask & ~part
        part |= linked
        frontier |= linked
    if part != mask:
        size = (_max_conflict_free_size(part, adj, memo)
                + _max_conflict_free_size(mask & ~part, adj, memo))
    else:
        # Leave out the most-conflicted candidate, or take it and drop its conflicts.
        v = max(_bit_indices(mask), key=lambda i: (adj[i] & mask).bit_count())
        rest = mask & ~(1 << v)
        size = 1 if not rest else max(
            _max_conflict_free_size(rest, adj, memo),
            1 + _max_conflict_free_size(rest & ~adj[v], adj, memo))
    memo[mask] = size
    return size


def _conflict_free_subsets(adj, size, memo):
    """Yield every conflict-free subset of exactly `size` candidates, as index
    tuples in itertools.combinations order. Branches that cannot reach `size`
    are cut using _max_conflict_free_size, whose memo is shared.
    """
    chosen = []

    def extend(allowed):
        if len(chosen) == size:
            yield tuple(chosen)
            return
        while allowed and len(chosen) + _max_conflict_free_size(allowed, adj, memo) >= size:
            low = allowed & -allowed
            allowed ^= low
            i = low.bit_length() - 1
            chosen.append(i)
            yield from extend(allowed & ~adj[i])
            chosen.pop()

    yield from extend((1 << len(adj)) - 1)


def get_valid_donor_combinations(comp, candidates):
    """Split candidates into non-conflicting donor combinations for a ligand.

    Returns ``(combinations, conflict_kinds)``; the kinds drive the policy in
    select_donor_combinations.
    """
    if len(candidates) > MAX_DONORS:
        raise EnumerationLimit('too many donors in one ligand')
    if len(candidates) <= 1:
        return [candidates], set()

    adj = [0] * len(candidates)
    kinds = set()
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            kind = classify_donor_conflict(
                comp, candidates[i]['atom_idx'], candidates[j]['atom_idx']
            )
            if kind:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
                kinds.add(kind)

    if not any(adj):
        return [candidates], kinds

    # Every maximal conflict-free subset, largest first, then in
    # itertools.combinations order: the same list, in the same order, as
    # scanning all subsets largest-first and keeping each one that no
    # earlier pick contains.
    masks = sorted(_maximal_conflict_free_masks(adj),
                   key=lambda mask: (-mask.bit_count(), _bit_indices(mask)))
    return [[candidates[i] for i in _bit_indices(mask)] for mask in masks], kinds


def select_donor_combinations(combinations, kinds):
    """Apply the per-kind policy to a ligand's candidate combinations.

    Linkage conflicts only  -> keep every combination; they are genuine
                               alternative binding modes (cyanide C- vs N-bound).
    Any denticity conflict  -> keep one. Which subset of a polydentate ligand
                               binds is underdetermined without a 3D structure,
                               and enumerating the subsets mostly produces
                               compounds that do not exist. Prefer the highest
                               denticity (the chelate effect), then the
                               best-priority donor from DONOR_SMARTS.

    A ligand with both kinds collapses to one row; mixing the policies would
    multiply guesses by isomers for no gain.
    """
    if len(combinations) <= 1 or kinds == {CONFLICT_LINKAGE}:
        return combinations
    return [min(combinations, key=lambda c: (
        -len(c), tuple(sorted(x['priority'] for x in c)),
        tuple(sorted(x['atom_idx'] for x in c))))]


def map_components(mol, metal_idx, verbose=True):
    """Split into components and extract donor combinations for ligand components."""
    frag_atoms = Chem.GetMolFrags(mol, asMols=False)
    components = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

    metal_comp = None
    ligand_comps = []

    for atoms, comp in zip(frag_atoms, components):
        if metal_idx in atoms:
            metal_comp = comp
            if verbose:
                print(f"  Metal: {Chem.MolToSmiles(comp)}")
        else:
            mapped, cands = find_donor_candidates(comp)
            combos, kinds = get_valid_donor_combinations(mapped, cands)
            n_raw = len(combos)
            combos = select_donor_combinations(combos, kinds)
            ligand_comps.append((mapped, combos))
            if verbose:
                print(
                    f"  Ligand: {Chem.MolToSmiles(mapped)}"
                )
                print(f"    Raw Donors found: {len(cands)}")
                print(f"    Valid Combinations: {len(combos)}"
                      + (f" (from {n_raw}, {'/'.join(sorted(kinds))} conflict"
                         f"{'s' if len(kinds) > 1 else ''})" if kinds else ""))

    return metal_comp, ligand_comps


def _raw_ligand_fragments(mol, metal_idx):
    """Like map_components, but returns each ligand's own raw donor
    candidates, uncollapsed. Soft mode resolves conflicts per ligand via
    get_valid_donor_combinations/select_donor_combinations; strict mode
    instead resolves them itself, globally across every fragment at once,
    so it needs the raw candidate list rather than the pre-collapsed combos.
    """
    frag_atoms = Chem.GetMolFrags(mol, asMols=False)
    components = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    metal_comp = None
    fragments = []
    for atoms, comp in zip(frag_atoms, components):
        if metal_idx in atoms:
            metal_comp = comp
        else:
            mapped, cands = find_donor_candidates(comp)
            fragments.append((mapped, cands))
    return metal_comp, fragments


def _maximal_conflict_free(comp, candidates, max_size=None, conflicts=None):
    """All maximal (non-dominated), same-size conflict-free subsets of
    candidates from ONE fragment, capped at max_size donors.

    Called by strict_donor_selection one priority tier at a time, so
    "candidates" only ever holds same-priority donors here -- more than one
    subset coming back means a genuine tie that priority alone cannot break.

    conflicts is an optional _ConflictMemo for comp, shared between calls.
    """
    n = len(candidates)
    cap = n if max_size is None else max(0, min(n, max_size))
    if cap == 0:
        return [[]]
    if conflicts is None:
        conflicts = _ConflictMemo(comp)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if conflicts[candidates[i]['atom_idx'], candidates[j]['atom_idx']] is not None:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    # The largest size within the cap that any conflict-free subset reaches,
    # then every subset of that size, in itertools.combinations order.
    memo = {}
    size = min(cap, _max_conflict_free_size((1 << n) - 1, adj, memo))
    return [[candidates[i] for i in combo]
            for combo in _conflict_free_subsets(adj, size, memo)]


def _grow_strict_branches(ligand_fragments, branches, priority, conflicts):
    """Yield every extension of each branch by the donors of one priority
    tier, in order and with duplicates (strict_donor_selection drops those).
    """
    for accepted, remaining in branches:
        if remaining <= 0:
            yield accepted, remaining
            continue
        per_fragment_options = {}
        for fid, (mapped, cands) in enumerate(ligand_fragments):
            tier = [c for c in cands if c['priority'] == priority]
            if not tier:
                continue
            fixed = accepted.get(fid, [])
            tier = [c for c in tier if not any(
                conflicts[fid][c['atom_idx'], f['atom_idx']] is not None
                for f in fixed
            )]
            if tier:
                per_fragment_options[fid] = _maximal_conflict_free(
                    mapped, tier, conflicts=conflicts[fid])
        if not per_fragment_options:
            yield accepted, remaining
            continue
        # Fragments never conflict with each other, so their own tie
        # options combine freely; only the shared capacity constrains them.
        frag_ids = list(per_fragment_options)
        for combo in itertools.product(*(per_fragment_options[fid] for fid in frag_ids)):
            union = [(fid, c) for fid, opt in zip(frag_ids, combo) for c in opt]
            pieces = ([union] if len(union) <= remaining
                      else itertools.combinations(union, remaining))
            for piece in pieces:
                grown = {k: list(v) for k, v in accepted.items()}
                for fid, c in piece:
                    grown.setdefault(fid, []).append(c)
                yield grown, remaining - len(piece)


def strict_donor_selection(ligand_fragments, target, max_combinations=MAX_COMBINATIONS,
                           conflicts=None):
    """Global, donor-priority-ordered greedy fill of a metal's coordination
    sphere, to exactly 'target' donors where achievable.

    Processes DONOR_SMARTS priority tiers best-to-worst, across every ligand
    fragment at once, so a high-priority donor on one ligand can outrank a
    low-priority donor on a different ligand for one of the open seats.
    Within a tier, a candidate is dropped if it conflicts -- CONFLICT_LINKAGE
    or CONFLICT_DENTICITY alike, strict mode does not distinguish them the
    way soft mode does -- with something already accepted from the SAME
    fragment; conflict is only ever possible within one fragment, since
    classify_donor_conflict needs a shared bonded graph.

    Returns a list of branches, each a {fragment_index: [candidates]} dict.
    More than one branch means a genuine priority tie could not be broken:
    either several same-priority donors within one fragment are mutually
    compatible but there isn't room for all of them, or several same-
    priority donors from DIFFERENT fragments are competing for the last open
    seat(s). Both are treated the same way -- enumerate every distinct way
    to fill the tied seats -- capped by max_combinations like every other
    enumeration in this module. (prepare_complex then keeps one branch.)

    conflicts is an optional list of one _ConflictMemo per fragment, so that
    repeated calls on the same fragments classify each donor pair once.
    """
    if conflicts is None:
        conflicts = [_ConflictMemo(mapped) for mapped, _ in ligand_fragments]
    priorities = sorted({c['priority'] for _, cands in ligand_fragments for c in cands})
    branches = [({}, target)]
    for p in priorities:
        seen = set()
        deduped = []
        # Branches are grown lazily so that a runaway tier stops at the cap.
        # Only first occurrences are kept and the count never shrinks, so
        # once past the cap it would still be past it at the end of the tier.
        for accepted, remaining in _grow_strict_branches(ligand_fragments, branches, p,
                                                         conflicts):
            key = tuple(sorted(
                (fid, tuple(sorted(c['atom_idx'] for c in cs)))
                for fid, cs in accepted.items()
            ))
            if key not in seen:
                seen.add(key)
                deduped.append((accepted, remaining))
                if len(deduped) > max_combinations:
                    raise EnumerationLimit(
                        f'too many tied strict-mode candidate branches (> {max_combinations})')
        branches = deduped
    return [accepted for accepted, _ in branches]


def _strict_natural_max(ligand_fragments, max_combinations=MAX_COMBINATIONS, conflicts=None):
    """Size of the unrestricted maximum, via the real algorithm.

    This MUST use strict_donor_selection itself, not a cheaper simplified
    walk: a naive single-pass greedy pass (accept each candidate in
    priority order if it doesn't conflict with what's already accepted for
    that fragment) can permanently lock in a low-value choice made purely
    by atom-index tiebreak order -- e.g. picking one isolated carboxylate
    oxygen simply because it happens to have a lower atom index than a
    ring's own pair of mutually-compatible oxygens, permanently blocking
    that better pair from ever being reached. strict_donor_selection's own
    per-fragment, per-tier maximal-independent-set logic does not have this
    flaw, so reuse it here at an effectively unbounded target -- the total
    candidate count is a safe upper bound -- and take the largest resulting
    branch. (Branches are not guaranteed to be equally sized: different tie
    choices at one tier can leave different room at later tiers.)
    """
    unbounded = sum(len(cands) for _, cands in ligand_fragments)
    branches = strict_donor_selection(ligand_fragments, unbounded,
                                      max_combinations=max_combinations, conflicts=conflicts)
    return max((sum(len(v) for v in b.values()) for b in branches), default=0)


def connect_ligands(metal_comp, ligand_combo, metal_idx, verbose=True):
    """Connect a single specific combination of ligand donors to the metal.

    Ligands are merged with CombineMols so that chirality, isotopes, bond
    stereo and aromaticity survive; rebuilding atoms by atomic number loses
    all of it. The metal component is placed first, so ``metal_idx`` stays
    valid in the combined molecule.
    """
    combined = Chem.Mol(metal_comp)
    dative_bonds = []

    for lig_mol, cands in ligand_combo:
        if not cands:
            continue
        offset = combined.GetNumAtoms()
        combined = Chem.CombineMols(combined, lig_mol)
        dative_bonds.extend(c['atom_idx'] + offset for c in cands)

    if not dative_bonds:
        return combined

    rw = Chem.RWMol(combined)
    for donor_idx in dative_bonds:
        rw.AddBond(donor_idx, metal_idx, Chem.BondType.DATIVE)

    result = rw.GetMol()
    try:
        Chem.SanitizeMol(result)
    except Exception as e:
        if verbose:
            print(f"    Sanitization Warning: {e}")
        return None

    # Dative bonds were appended after the tags were set, so re-perceive
    # stereo against the final neighbour ordering.
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    return result


def _build_candidate(metal_comp, selection, metal_idx, known=()):
    """Build one donor selection: (canonical SMILES, molecule), or None if it fails.

    The SMILES is also parsed back, the serialization downstream tools will
    read, unless it is in `known`: those were checked already.
    """
    try:
        connected = connect_ligands(metal_comp, selection, metal_idx, verbose=False)
        if connected is None:
            return None
        connected = Chem.RemoveHs(connected)
        smi = Chem.MolToSmiles(connected)
        if smi not in known:
            parsed = Chem.MolFromSmiles(smi)
            if parsed is None or len(Chem.GetMolFrags(parsed)) != 1:
                return None
    except Exception:
        return None
    return smi, connected


# ----------------------------------------------------------------------------
# Strict mode builds every tied branch only to keep the lexicographically
# first SMILES, and most ties come from symmetry: identical ligand copies
# competing for the last seats, or equivalent donors on one ligand (the -OH
# groups of an inositol, the carbons of a Cp ring). Branches that a symmetry
# of the ligands maps onto each other build the same structure, so only the
# first of them needs building.
# ----------------------------------------------------------------------------

# Self-matches examined per ligand. The cap can only hide a symmetry, so
# fewer branches are recognized as equivalent, never a wrong pair.
MAX_SYMMETRIES = 256

_TETRAHEDRAL = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
_ORDINARY_BONDS = (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE,
                   Chem.BondType.AROMATIC)


def _ligand_symmetries(mol, donors):
    """How the symmetries of a ligand move its donors: one tuple per distinct
    action, giving the image of each atom in `donors`, identity included.

    A symmetry is a self-match of the ligand -- which already keeps elements,
    charges, isotopes, radicals and bond types -- that also keeps hydrogens,
    map numbers and the handedness of every tetrahedral centre. Mirror-related
    centres (a meso ligand) therefore do not count: they would build the
    enantiomer. Stereo and directional bonds, unusual bond types and other
    kinds of stereo centre must stay put together with their neighbours.
    """
    identity = tuple(donors)
    if len(donors) < 2:
        return [identity]
    probe = Chem.Mol(mol)  # ranking and matching must not touch the ligand itself
    ranks = list(Chem.CanonicalRankAtoms(probe, breakTies=False))
    if len({ranks[i] for i in donors}) == len(donors):
        return [identity]  # symmetries keep ranks, so no donor can move
    classes, atom_class, centres, fixed = {}, [], {}, set()
    for atom in probe.GetAtoms():
        idx, tag = atom.GetIdx(), atom.GetChiralTag()
        invariant = (atom.GetTotalNumHs(), atom.GetNumExplicitHs(), atom.GetNoImplicit(),
                     atom.GetIsAromatic(), atom.GetAtomMapNum(),
                     'tetrahedral' if tag in _TETRAHEDRAL else tag)
        atom_class.append(classes.setdefault(invariant, len(classes)))
        if tag in _TETRAHEDRAL:
            centres[idx] = (tag, [b.GetOtherAtomIdx(idx) for b in atom.GetBonds()])
        elif tag != Chem.ChiralType.CHI_UNSPECIFIED:
            fixed.add(idx)
            fixed.update(nb.GetIdx() for nb in atom.GetNeighbors())
    for bond in probe.GetBonds():
        bond_type = bond.GetBondType()
        if (bond.GetStereo() != Chem.BondStereo.STEREONONE
                or bond.GetBondDir() != Chem.BondDir.NONE
                or bond_type not in _ORDINARY_BONDS
                or bond.GetIsAromatic() != (bond_type == Chem.BondType.AROMATIC)):
            for end in (bond.GetBeginAtom(), bond.GetEndAtom()):
                fixed.add(end.GetIdx())
                fixed.update(nb.GetIdx() for nb in end.GetNeighbors())
    for group in probe.GetStereoGroups():
        fixed.update(a.GetIdx() for a in group.GetAtoms())

    def keeps_handedness(perm):
        # A centre keeps its handedness when its tag, read in the order of
        # its mapped neighbours, equals the tag of the atom it maps onto.
        for idx, (tag, neighbours) in centres.items():
            image_tag, image_neighbours = centres[perm[idx]]
            order = [image_neighbours.index(perm[n]) for n in neighbours]
            odd = sum(a > b for i, a in enumerate(order) for b in order[i + 1:]) % 2
            if (tag == image_tag) == bool(odd):
                return False
        return True

    actions = {identity}
    for perm in probe.GetSubstructMatches(probe, uniquify=False, useChirality=True,
                                          maxMatches=MAX_SYMMETRIES):
        if ([atom_class[j] for j in perm] == atom_class
                and all(perm[i] == i for i in fixed)
                and keeps_handedness(perm)):
            actions.add(tuple(perm[i] for i in donors))
    return sorted(actions)


def _branch_keys(metal_comp, metal_idx, ligand_fragments, branches):
    """One key per strict-mode branch, equal only when a symmetry of the
    ligands maps one branch's donors onto the other's: identical ligand
    copies may trade places, and each ligand may map onto itself (see
    _ligand_symmetries). Such branches build the same structure, or fail the
    same way. Without a usable symmetry, every key differs.
    """
    # Stereo at the metal would depend on the order its new bonds arrive in.
    metal = metal_comp.GetAtomWithIdx(metal_idx)
    if (len(branches) < 2 or metal.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            or any(b.GetStereo() != Chem.BondStereo.STEREONONE
                   or b.GetBondDir() != Chem.BondDir.NONE for b in metal.GetBonds())):
        return list(range(len(branches)))
    donor_sets = [[tuple(sorted(c['atom_idx'] for c in branch.get(fid, ())))
                   for fid in range(len(ligand_fragments))] for branch in branches]
    # Byte-identical pickles are the same molecule in the same atom order,
    # down to stereo, hydrogens and map numbers: copies of one kind.
    kinds, kind = {}, []
    for mapped, _ in ligand_fragments:
        kind.append(kinds.setdefault(mapped.ToBinary(), len(kinds)))
    members = {}
    for fid, k in enumerate(kind):
        members.setdefault(k, []).append(fid)
    # A lone ligand with the same donors in every branch cannot tell them apart.
    varying = [k for k, fids in members.items()
               if len(fids) > 1 or len({sets[fids[0]] for sets in donor_sets}) > 1]
    actions, canonical = {}, {}

    def canonical_donors(fid, chosen):
        # Symmetries map the donor candidates onto themselves, so none or
        # all of them has only one form.
        if len(chosen) in (0, len(ligand_fragments[fid][1])):
            return chosen
        k = kind[fid]
        if (k, chosen) not in canonical:
            if k not in actions:
                mapped, cands = ligand_fragments[fid]
                donors = [c['atom_idx'] for c in cands]
                actions[k] = (_ligand_symmetries(mapped, donors),
                              {atom: pos for pos, atom in enumerate(donors)})
            images, position = actions[k]
            canonical[k, chosen] = min(tuple(sorted(image[position[a]] for a in chosen))
                                       for image in images)
        return canonical[k, chosen]

    return [tuple((k, tuple(sorted(canonical_donors(fid, sets[fid]) for fid in members[k])))
                  for k in varying)
            for sets in donor_sets]


# ============================================================================
# SHARED PREPARATION / DATAFRAME / CLI
# ============================================================================

KEEP = 'keep'
PROCESS = 'process'
REJECT = 'reject'
REJECT_POLYNUCLEAR = 'organometallic_polynuclear'
REJECT_MULTIPLE_CENTRES = 'organometallic_multiple_centres'
REJECT_NO_DONOR = 'organometallic_no_donor'
REJECT_BUILD_FAILED = 'organometallic_build_failed'
REJECT_COVALENT = 'organometallic_connected'
REJECT_COORDINATION = 'organometallic_coordination_limit'
REJECT_ENUMERATION = 'organometallic_enumeration_limit'


@dataclass
class PreparationResult:
    action: str
    reason: str | None = None
    molecules: list[Chem.Mol] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    changes: list[str] = field(default_factory=list)
    removed_fragments: list[str] = field(default_factory=list)
    charge_delta: int = 0  # Protonation repairs only; excludes removed ions.


def coordination_limit(metal):
    """Generous guardrails, not predictions of geometry or oxidation state."""
    z = metal.GetAtomicNum()
    if 57 <= z <= 71:
        return 12
    if z in (21, 22, 39, 40, 72, 73):
        return 8
    if z in (47, 79):
        return 4
    return 6


# ----------------------------------------------------------------------------
# Precise, per-element/per-oxidation-state coordination numbers.
#
# This is chemistry judgment, not something derivable from the graph, so it
# lives in an external, human-editable file (default: metal_coordination_
# numbers.txt next to this script) in the project's standard tab-separated
# ion/valence-table format, rather than a dict in this module -- see that
# file's header for the format itself and for notes on what to enter.
#
# coordination_limit() above stays untouched: it is a cheap, generic guard
# against combinatorial blow-up, evaluated for every metal. This table is
# the opposite -- precise but incomplete, filled in only for metals/charges
# someone has actually vetted -- and is applied as a second pass, after
# candidates are already built, so it can prefer a good candidate, flag a bad
# one, or defer entirely when there is no entry to consult.
# ----------------------------------------------------------------------------

DEFAULT_COORDINATION_TABLE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'metal_coordination_numbers.txt'
)

# [Sym+n] with an explicit signed charge, e.g. [Bi+3], [Pd+0], [Fe-1].
# [Sym] with no charge at all, e.g. [Bi], is the element-only wildcard.
_ION_RE = re.compile(r'^\[([A-Za-z]{1,2})([+-]\d+)?\]$')


def load_coordination_table(path=DEFAULT_COORDINATION_TABLE_PATH):
    """Load {(symbol, charge_or_None): {allowed donor counts}} from the
    project's tab-separated ion/valence table (ION<TAB>COORDINATION_NUMBERS,
    '#' comments, blank lines skipped). charge_or_None is None for an
    element-only wildcard row, e.g. [Bi] with no charge sign.

    Missing file: warn once and return {}, so the rest of the pipeline
    behaves exactly as it did before this table existed (every metal takes
    the 'no data' branch of select_by_coordination_number). A malformed
    individual line is skipped with its own warning rather than failing the
    whole file, since one bad row (e.g. a stray element name with no
    brackets or value) should not silently disable every other entry.
    """
    if not os.path.exists(path):
        print(f"NOTE: no coordination-number table found at '{path}'; "
              f"proceeding without one (see coordination_limit() only).")
        return {}

    table = {}
    with open(path) as f:
        for line_num, raw_line in enumerate(f, 1):
            line = raw_line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) != 2:
                print(f"WARNING: '{path}' line {line_num}: expected "
                      f"'ION<TAB>COORDINATION_NUMBERS', got {line!r}; skipping")
                continue
            ion_field, values_field = fields[0].strip(), fields[1].strip()
            match = _ION_RE.match(ion_field)
            if not match:
                print(f"WARNING: '{path}' line {line_num}: '{ion_field}' is not "
                      f"a valid ion, e.g. [Bi+3] or [Bi]; skipping")
                continue
            symbol, charge_str = match.group(1), match.group(2)
            charge = int(charge_str) if charge_str is not None else None
            try:
                values = frozenset(int(v.strip()) for v in values_field.split(','))
            except ValueError:
                print(f"WARNING: '{path}' line {line_num}: could not parse "
                      f"coordination numbers {values_field!r}; skipping")
                continue
            table[(symbol, charge)] = values
    return table


COORDINATION_TABLE = load_coordination_table()


def coordination_numbers_for(metal_atom, table=COORDINATION_TABLE):
    """The allowed set for this atom's (element, formal charge).

    A specific-charge entry, e.g. [Bi+3], always takes precedence over the
    element-only wildcard, e.g. [Bi], for the same element; the wildcard is
    used only when no entry names this exact charge. Returns None if
    neither is present.
    """
    symbol = _PTABLE.GetElementSymbol(metal_atom.GetAtomicNum())
    charge = metal_atom.GetFormalCharge()
    if (symbol, charge) in table:
        return table[(symbol, charge)]
    return table.get((symbol, None))


def select_by_coordination_number(outputs, donor_counts, allowed, metal_atom, warnings):
    """Apply the table's verdict to a dict of {smiles: mol} candidates.

    donor_counts holds each candidate's TOTAL coordination number (bonds
    already drawn to the metal in the input plus new dative bonds).

    - No table entry for this metal/charge: return everything unchanged.
      (This is the common case until the table is filled in further.)
    - One or more candidates hit an allowed number: keep exactly those, so
      the person can pick between legitimate alternatives.
    - None do: keep only the candidate(s) closest to the allowed set (ties
      broken toward the fuller one) and add a warning -- the result is kept,
      not rejected, but flagged for manual review.
    """
    if not allowed:
        return outputs

    matched = {s: m for s, m in outputs.items() if donor_counts[s] in allowed}
    if matched:
        return matched

    def distance(smi):
        n = donor_counts[smi]
        return (min(abs(n - a) for a in allowed), -n)

    best = min(donor_counts, key=distance)
    best_n = donor_counts[best]
    symbol = _PTABLE.GetElementSymbol(metal_atom.GetAtomicNum())
    charge = metal_atom.GetFormalCharge()
    warnings.append(NeedsReview(
        f'no candidate reaches coordination number {sorted(allowed)} known for '
        f'{symbol}{charge:+d}; kept closest match (CN {best_n}) -- flagged, verify manually'))
    return {s: m for s, m in outputs.items() if donor_counts[s] == best_n}


def _donor_selections(ligands):
    """Lazy Cartesian enumeration, removing permutations of identical ligands.

    map_components canonicalizes each ligand before assigning donor indices.
    Consecutive copies with the same molecule and donor choices can therefore
    use nondecreasing mode indices without losing distinct connectivity graphs.
    """
    def signature(item):
        mol, combos = item
        return (Chem.MolToSmiles(mol),
                tuple(tuple(c['atom_idx'] for c in combo) for combo in combos))

    # Stable sort on the signature alone, computed once per ligand.
    keyed = sorted(((signature(item), item) for item in ligands), key=lambda pair: pair[0])
    ordered = [item for _, item in keyed]
    signatures = [sig for sig, _ in keyed]
    selection = []
    indices = []

    def visit(position):
        if position == len(ordered):
            yield tuple(selection)
            return
        mol, combos = ordered[position]
        start = (indices[-1] if position and signatures[position] == signatures[position - 1]
                 else 0)
        for index in range(start, len(combos)):
            selection.append((mol, combos[index]))
            indices.append(index)
            yield from visit(position + 1)
            indices.pop()
            selection.pop()

    yield from visit(0)


def prepare_complex(mol, smiles=None, *, max_variants=MAX_VARIANTS,
                    max_combinations=MAX_COMBINATIONS, reject_covalent=False,
                    mode='soft'):
    """Return connected candidates or an explicit unsupported/rejected result.

    Existing metal bonds are retained. More than one supported metal anywhere
    in a record is unsupported, including an already connected metal cluster.
    Only the explicitly listed spectator ions are removed. Unknown fragments
    reject the record; no ligand or solvent is silently discarded (strict
    mode may leave ligands uncoordinated when the sphere is full, and always
    says so in a warning).

    Formal metal charges are never guessed. Protonation repairs and mixed
    covalent/dative conventions remain visible in warnings and provenance.

    Coordination numbers from the external table are compared against the
    metal's TOTAL coordination: bonds already drawn to the metal in the input
    plus the new dative bonds.

    mode='soft' (default) enumerates every chemically plausible donor
    combination -- including genuine denticity ambiguity -- as separate
    output rows, for a person to pick between.

    mode='strict' instead fills the metal's coordination sphere
    deterministically by donor priority (DONOR_SMARTS order) and always
    returns exactly one structure: denticity and linkage conflicts collapse
    to the best answer, and a remaining exact priority tie is resolved by
    keeping the lexicographically first canonical SMILES, with a warning.
    See strict_donor_selection() for the algorithm.
    """
    if mode not in ('soft', 'strict'):
        raise ValueError("mode must be 'soft' or 'strict'")
    if max_variants < 1 or max_combinations < 1:
        raise ValueError('metal enumeration limits must be positive')
    if mol is None or mol.GetNumAtoms() == 0:
        return PreparationResult(REJECT, REJECT_BUILD_FAILED)
    counts, matches = metals_per_fragment(mol)
    if not matches:
        return PreparationResult(KEEP, molecules=[Chem.Mol(mol)])
    if len(matches) > 1:
        reason = REJECT_POLYNUCLEAR if any(c > 1 for c in counts) else REJECT_MULTIPLE_CENTRES
        return PreparationResult(REJECT, reason)
    metal = mol.GetAtomWithIdx(matches[0][0])
    # Early plausibility guard. The generic per-block cap is widened to the
    # table's largest value when a table entry exists, so the table -- not the
    # generic guardrail -- decides what is plausible for tabulated metals.
    known_early = coordination_numbers_for(metal)
    early_cap = (max(coordination_limit(metal), max(known_early))
                 if known_early else coordination_limit(metal))
    if metal.GetDegree() > early_cap:
        return PreparationResult(REJECT, REJECT_COORDINATION)
    # Connectivity alone does not distinguish covalent and dative complexes.
    if len(counts) == 1:
        if reject_covalent:
            return PreparationResult(REJECT, REJECT_COVALENT)
        return PreparationResult(KEEP, molecules=[Chem.Mol(mol)])

    result = PreparationResult(PROCESS)
    components = []
    for component in Chem.GetMolFrags(mol, asMols=True):
        if component.GetNumAtoms() <= SPECTATOR_MAX_ATOMS:
            smi = Chem.MolToSmiles(component)
            if smi in SPECTATOR_SMILES:
                result.removed_fragments.append(smi)
                continue
        components.append(component)
    if result.removed_fragments:
        result.warnings.append('removed recognized spectator ions: '
                               + '.'.join(result.removed_fragments))
    # Bound recursion and avoid any subset work on implausibly many fragments.
    if len(components) - 1 + metal.GetDegree() > early_cap:
        result.action, result.reason = REJECT, REJECT_COORDINATION
        return result

    try:
        working = combine_components(components)
        initial_charge = Chem.GetFormalCharge(working)
        final = _prepare_ligands(working, changes=result.changes)
        result.charge_delta = Chem.GetFormalCharge(final) - initial_charge
        if result.changes:
            result.warnings.append(
                f'protonation repairs applied; charge change {result.charge_delta:+d} '
                '(excluding removed spectators)')
        if metal.GetFormalCharge() < 0 and any(
                b.GetBondType() != Chem.BondType.DATIVE for b in metal.GetBonds()):
            result.warnings.append(NeedsReview(
                'mixed covalent/dative charge representation; metal oxidation state '
                'was not inferred, candidate charge needs review'))
        metal_idx = final.GetSubstructMatches(METAL_PATTERN)[0][0]

        if mode == 'soft':
            metal_comp, ligands = map_components(final, metal_idx, verbose=False)
            if any(not combos or not combos[0] for _, combos in ligands):
                result.action, result.reason = REJECT, REJECT_NO_DONOR
                return result
            metal_idx = metal_comp.GetSubstructMatches(METAL_PATTERN)[0][0]
            local_metal = metal_comp.GetAtomWithIdx(metal_idx)
            existing_bonds = local_metal.GetDegree()
            known_cns = coordination_numbers_for(local_metal)
            # coordination_limit() is a generic per-block guardrail against
            # combinatorial explosion, not a chemistry judgment (see its
            # docstring). When we have a precise, element/charge-specific table
            # entry instead, widen the cap to cover it, so a candidate the table
            # actually cares about is never silently discarded before the table
            # gets a chance to rule on it below.
            if known_cns:
                # A precise table entry exists: it is now the authority on "too
                # many donors" (via select_by_coordination_number below,
                # including its closest-match-but-flagged fallback), so do not
                # let the generic per-block guardrail prune candidates first --
                # that would silently discard exactly the over-coordinated
                # candidate the table needs to see and flag. MAX_DONORS /
                # MAX_COMBINATIONS / max_variants remain as the real backstops
                # against combinatorial blow-up.
                available = MAX_DONORS
            else:
                available = coordination_limit(local_metal) - existing_bonds
            outputs = {}
            donor_counts = {}
            failed = 0
            search_truncated = False
            for attempt, selection in enumerate(_donor_selections(ligands)):
                if attempt >= max_combinations:
                    result.warnings.append(NeedsReview(
                        f'candidate search truncated at {max_combinations} combinations'))
                    search_truncated = True
                    break
                n_donors = sum(len(cands) for _, cands in selection)
                if n_donors > available:
                    continue
                built = _build_candidate(metal_comp, selection, metal_idx, known=outputs)
                if built is None:
                    failed += 1
                    continue
                smi, connected = built
                if smi in outputs:
                    continue  # a duplicate: already checked and kept
                if len(outputs) >= max_variants:
                    result.warnings.append(NeedsReview(f'output truncated at {max_variants} variants'))
                    break
                outputs[smi] = connected
                # Total coordination: bonds already drawn to the metal
                # in the input plus the new dative bonds.
                donor_counts[smi] = n_donors + existing_bonds
            if failed:
                result.warnings.append(f'{failed} candidate builds failed sanitization')
            if not outputs:
                result.action = REJECT
                result.reason = (REJECT_ENUMERATION if search_truncated else
                                 REJECT_BUILD_FAILED if failed else REJECT_COORDINATION)
                return result
            outputs = select_by_coordination_number(
                outputs, donor_counts, known_cns, local_metal, result.warnings)
            result.molecules = [outputs[key] for key in sorted(outputs)]
            return result

        else:  # mode == 'strict'
            metal_comp, ligand_fragments = _raw_ligand_fragments(final, metal_idx)
            if any(not cands for _, cands in ligand_fragments):
                result.action, result.reason = REJECT, REJECT_NO_DONOR
                return result
            metal_idx = metal_comp.GetSubstructMatches(METAL_PATTERN)[0][0]
            local_metal = metal_comp.GetAtomWithIdx(metal_idx)
            existing_bonds = local_metal.GetDegree()
            known_cns = coordination_numbers_for(local_metal)
            # Largest number of NEW dative bonds reachable at all, computed
            # with the same algorithm as the final selection (see
            # _strict_natural_max).
            conflicts = [_ConflictMemo(mapped) for mapped, _ in ligand_fragments]
            n_max = _strict_natural_max(ligand_fragments, max_combinations=max_combinations,
                                        conflicts=conflicts)
            if known_cns:
                # Fill to the LARGEST allowed TOTAL coordination number
                # (bonds already drawn to the metal + new dative bonds) this
                # input can reach; if none is reachable, keep everything
                # reachable and flag it -- same policy and wording as soft
                # mode's coordination-number filter.
                reachable = [v for v in known_cns
                             if existing_bonds <= v <= existing_bonds + n_max]
                if reachable:
                    target = max(reachable) - existing_bonds
                else:
                    symbol = _PTABLE.GetElementSymbol(local_metal.GetAtomicNum())
                    charge = local_metal.GetFormalCharge()
                    result.warnings.append(NeedsReview(
                        f'no candidate reaches coordination number {sorted(known_cns)} '
                        f'known for {symbol}{charge:+d}; kept closest match '
                        f'(CN {existing_bonds + n_max}) -- flagged, verify manually'))
                    target = n_max
            else:
                target = coordination_limit(local_metal) - existing_bonds

            branches = strict_donor_selection(ligand_fragments, target,
                                              max_combinations=max_combinations,
                                              conflicts=conflicts)
            outputs = {}
            uncoordinated = {}
            failed = 0
            outcome = {}  # branch key -> SMILES it built, or None if the build failed
            keys = _branch_keys(metal_comp, metal_idx, ligand_fragments, branches)
            for branch, key in zip(branches, keys):
                if key in outcome:
                    # Symmetric to a branch built already: it would build the
                    # same SMILES, or fail in the same way.
                    if outcome[key] is None:
                        failed += 1
                    continue
                combo_input = [(ligand_fragments[fid][0], cands) for fid, cands in branch.items()]
                built = _build_candidate(metal_comp, combo_input, metal_idx, known=outputs)
                outcome[key] = built[0] if built else None
                if built is None:
                    failed += 1
                    continue
                smi, connected = built
                if smi in outputs:
                    continue  # a tie that built the same structure again
                if len(outputs) >= max_variants:
                    result.warnings.append(NeedsReview(f'output truncated at {max_variants} variants'))
                    break
                outputs[smi] = connected
                # Fragments that got no seat in this branch are not part of
                # the connected output; remembered so the kept result can
                # report them.
                uncoordinated[smi] = sorted(Chem.MolToSmiles(ligand_fragments[fid][0])
                                            for fid in range(len(ligand_fragments))
                                            if not branch.get(fid))
            if failed:
                result.warnings.append(f'{failed} candidate builds failed sanitization')
            if not outputs:
                result.action = REJECT
                result.reason = REJECT_BUILD_FAILED if failed else REJECT_COORDINATION
                return result
            ordered = sorted(outputs)
            if uncoordinated[ordered[0]]:
                dropped = uncoordinated[ordered[0]]
                result.warnings.append(NeedsReview(
                    f'{len(dropped)} ligand fragment(s) left uncoordinated (coordination '
                    f'sphere full) and not included in the output: ' + '.'.join(dropped)))
            if len(ordered) > 1:
                # A genuine priority tie (see strict_donor_selection) with no
                # way to prefer one option over another -- e.g. which of
                # several chemically equivalent pendant -OH groups fills a
                # leftover coordination slot. Keep one, deterministically
                # (lexicographically first by canonical SMILES), rather than
                # returning every indistinguishable variant.
                result.warnings.append(NeedsReview(
                    f'{len(ordered)} equivalent binding options tied on donor '
                    f'priority; arbitrarily kept the lexicographically first'))
                ordered = ordered[:1]
            result.molecules = [outputs[key] for key in ordered]
            return result
    except EnumerationLimit as exc:
        result.action, result.reason = REJECT, REJECT_ENUMERATION
        result.warnings.append(str(exc))
    except Exception as exc:
        result.action, result.reason = REJECT, REJECT_BUILD_FAILED
        result.warnings.append(f'{type(exc).__name__}: {exc}')
    return result


def triage(mol, smiles, reject_covalent=False, mode='soft'):
    """Compatibility helper; classification uses the same engine as the CLI."""
    result = prepare_complex(mol, smiles, reject_covalent=reject_covalent, mode=mode)
    return result.action, result.reason


def partition_dataframe(df, rejectedFile=None, debug=False, reject_covalent=False,
                        *, max_variants=MAX_VARIANTS, max_combinations=MAX_COMBINATIONS,
                        mode='soft', provenance=True):
    """Prepare and expand rows, returning (prepared, rejected).

    Only rows whose SMILES holds a metal (or whose mol is missing) are
    prepared; every other row, and every metal row that needs no change,
    passes through untouched and in order. All original columns survive.
    With provenance=True both frames also carry metal_* columns recording
    what happened to each row. With provenance=False rows keep only their
    own columns, the rejected frame gains metal_reason, and a frame without
    metal rows is returned as it is.

    The outputs of one record are numbered id_1, id_2, ... (two digits from
    ten outputs on, like Msani's own numbering). Notes are logged per
    record: as warnings when they ask for a manual look (NeedsReview), as
    info otherwise. A rejection file uses the existing Msani convention:
    original SMILES, ID, reason.
    """
    import logging
    import pandas as pd
    logger = logging.getLogger('msani')
    if max_variants < 1 or max_combinations < 1:
        raise ValueError('metal enumeration limits must be positive')
    # Organic rows bypass per-row RDKit work and reconstruction entirely. The
    # API requires smiles/mol to describe the same input, as the Msani parser does.
    candidate_mask = df['smiles'].map(has_metal) | df['mol'].isna()
    if not provenance and not candidate_mask.any():
        return df, df.iloc[:0].assign(metal_reason=None)
    base = df.reset_index(drop=True).copy()
    candidate_mask = candidate_mask.reset_index(drop=True)
    if provenance:
        defaults = dict(metal_original_smiles=base['smiles'], metal_parent_id=base['ids'],
                        metal_variant=0, metal_status=KEEP, metal_reason=None,
                        metal_warnings=[()] * len(base), metal_changes=[()] * len(base),
                        metal_removed_fragments=[()] * len(base), metal_charge_delta=0)
        for key, default in defaults.items():
            if key not in base:
                base[key] = default
    columns = list(base.columns)
    unchanged = ~candidate_mask
    kept, positions, rejected = [], [], []
    reserved_ids = set(base['ids'].astype(str))
    for position, row in base.loc[candidate_mask].iterrows():
        result = prepare_complex(row['mol'], row['smiles'], max_variants=max_variants,
                                 max_combinations=max_combinations,
                                 reject_covalent=reject_covalent, mode=mode)
        if result.action == KEEP:
            # Already connected, or no metal after all: nothing to correct.
            unchanged[position] = True
            continue
        entry = row.to_dict()
        if provenance:
            entry.update(metal_status=result.action, metal_warnings=tuple(result.warnings),
                         metal_changes=tuple(result.changes),
                         metal_removed_fragments=tuple(result.removed_fragments),
                         metal_charge_delta=result.charge_delta)
        if result.action == REJECT:
            entry['metal_reason'] = result.reason
            rejected.append(entry)
            logger.info('Metal preparation %s rejected (%s)%s', row['ids'], result.reason,
                        ''.join(f'; {note}' for note in result.warnings))
            continue
        for note in result.warnings:
            logger.log(logging.WARNING if isinstance(note, NeedsReview) else logging.INFO,
                       'Metal preparation %s: %s', row['ids'], note)
        count = len(result.molecules)
        width = 2 if count >= 10 else 1
        for i, mol in enumerate(result.molecules, 1):
            output = dict(entry)
            output.update(mol=mol, smiles=Chem.MolToSmiles(mol))
            if provenance:
                output['metal_variant'] = i
            if count > 1:
                candidate = f"{row['ids']}_{i:0{width}d}"
                while candidate in reserved_ids:
                    candidate += '_'
                reserved_ids.add(candidate)
                output['ids'] = candidate
            kept.append(output)
            positions.append(position)
    prepared_df = pd.DataFrame(kept, columns=columns, index=positions)
    passthrough = base.loc[unchanged]
    parts = [part for part in (passthrough, prepared_df) if not part.empty]
    kept_df = (pd.concat(parts).sort_index(kind='stable').reset_index(drop=True)
               if parts else base.iloc[:0].copy())
    rejected_df = pd.DataFrame(
        rejected, columns=columns if provenance else columns + ['metal_reason'])

    if not rejected_df.empty:
        if rejectedFile is not None:
            original = 'metal_original_smiles' if provenance else 'smiles'
            rejected_df[[original, 'ids', 'metal_reason']].to_csv(
                rejectedFile, index=False, mode='a', sep=' ', header=False)
        logger.info('Rejected %d metal records: %s', len(rejected_df),
                    rejected_df['metal_reason'].value_counts().to_dict())
    return kept_df, rejected_df


def process_file(input_file, verbose=True, *, rejected_file=None, output_file=None,
                 max_variants=MAX_VARIANTS, max_combinations=MAX_COMBINATIONS,
                 mode='soft'):
    """Prototype CLI using the same preparation and rejection policy as Msani.

    rejected_file/output_file, if given, are written to incrementally as each
    record is processed (and flushed after every write) rather than buffered
    in memory and written once at the end -- so a crash, timeout, or Ctrl-C
    partway through a large run does not lose the records already completed.
    """
    n_rejects = 0
    results = []
    reject_stream = open(rejected_file, 'w') if rejected_file is not None else None
    output_stream = open(output_file, 'w') if output_file is not None else None
    try:
        for mol, name, smi in read_smiles_file(input_file):
            prepared = prepare_complex(mol, smi, max_variants=max_variants,
                                       max_combinations=max_combinations, mode=mode)
            for warning in prepared.warnings:
                print(f'WARNING {name}: {warning}')
            if prepared.action == REJECT:
                print(f'REJECT {name}: {prepared.reason}')
                n_rejects += 1
                if reject_stream is not None:
                    reject_stream.write(f'{smi} {name} {prepared.reason}\n')
                    reject_stream.flush()
                continue
            for i, candidate in enumerate(prepared.molecules, 1):
                variant = f'{name}_{i}' if len(prepared.molecules) > 1 else name
                out = Chem.MolToSmiles(candidate)
                results.append((variant, smi, out))
                if verbose:
                    print(f'{variant}:\n  Input:  {smi}\n  Output: {out}')
                if output_stream is not None:
                    # Same whitespace-separated SMILES<TAB>NAME convention
                    # read_smiles_file reads, so this file round-trips
                    # straight back in as input elsewhere.
                    output_stream.write(f'{out}\t{variant}\n')
                    output_stream.flush()
    finally:
        if reject_stream is not None:
            reject_stream.close()
        if output_stream is not None:
            output_stream.close()
    print(f'Produced {len(results)} candidates; rejected {n_rejects} records')
    return results


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n\n')[0])
    parser.add_argument('input_file')
    parser.add_argument('--rejected-file',
                        help='write rejected records here as SMILES NAME REASON, '
                             'one per line (overwrites an existing file).')
    parser.add_argument('--output-file',
                        help='write accepted candidates here as SMILES<TAB>NAME, '
                             'one per line (overwrites an existing file) -- the same '
                             '.smi convention this script reads, so the file can be '
                             'fed back in as input.')
    parser.add_argument('--max-variants', type=int, default=MAX_VARIANTS)
    parser.add_argument('--max-combinations', type=int, default=MAX_COMBINATIONS)
    parser.add_argument('--mode', choices=('soft', 'strict'), default='soft',
                        help="soft: enumerate every plausible donor combination "
                             "(default). strict: fill the coordination sphere "
                             "deterministically by donor priority; always one answer "
                             "per input (exact priority ties are resolved "
                             "arbitrarily, with a warning).")
    args = parser.parse_args()
    if args.max_variants < 1 or args.max_combinations < 1:
        parser.error('enumeration limits must be positive')
    process_file(args.input_file, rejected_file=args.rejected_file,
                 output_file=args.output_file,
                 max_variants=args.max_variants, max_combinations=args.max_combinations,
                 mode=args.mode)


if __name__ == '__main__':
    main()
