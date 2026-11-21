"""Behaviour of the C++ torsion sampler: torsion setting, RMSD pruning, budgets and energy window."""
import contextlib
import io
import itertools
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import rdMolTransforms as T
from rdkit.Chem.Lipinski import RotatableBondSmarts

from msani.conformers import conformers
from msani.conformers.conformers import CPP_AVAILABLE, ConformerGenerator

ANGLES = [60.0, 180.0, 300.0]
RMSD = 0.5


def _rigid_rotor_input(smiles):
    """One MMFF-minimised conformer and a heavy-atom torsion per rotatable bond."""
    heavy = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(heavy)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    torsions = []
    for j, k in heavy.GetSubstructMatches(RotatableBondSmarts):
        i = next(a.GetIdx() for a in heavy.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx() != k)
        l = next(a.GetIdx() for a in heavy.GetAtomWithIdx(k).GetNeighbors() if a.GetIdx() != j)
        torsions.append([i, j, k, l])
    return mol, torsions


def _discrete_sample(mol, torsions, numConfs=200):
    angle_map = {n: ('test', t, ANGLES) for n, t in enumerate(torsions)}
    score_map = {n: [1.0] * len(ANGLES) for n in angle_map}
    return conformers.cpp_sampler.stochastic_sampling_discrete(
        mol=Chem.Mol(mol), angle_map=angle_map, score_map=score_map,
        possible_numConfs=len(ANGLES) ** len(torsions), importance_order=[1.0] * len(torsions),
        window=0.0, hetero_H_bonds=[], timeout_conf=60, rmsd=RMSD, numConfs=numConfs, randomSeed=42)


def _sample(smiles, numConfs, mode='fixed', energywindow=25):
    with contextlib.redirect_stdout(io.StringIO()):
        generator = ConformerGenerator(smiles, 'test', randomSeed=42, mode=mode)
        generator.conf_sampling(numConfs=numConfs, energywindow=energywindow, timeout_conf=1)
    return generator.ring_confs[0]


@unittest.skipUnless(CPP_AVAILABLE, 'requires the C++ sampling extension')
class TestDiscreteSampler(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # tert-butyl and para-phenylene give symmetric atoms, so the RMSD
        # pruning has non-trivial symmetry orbits to work with.
        cls.mol, cls.torsions = _rigid_rotor_input('CC(C)(C)c1ccc(OCCCN(C)C)cc1')
        cls.result = _discrete_sample(cls.mol, cls.torsions)

    def test_sampled_torsions_take_the_requested_angles(self):
        # The sampler only rotates torsions (no minimisation), so every
        # conformer must sit exactly on one of the requested angles.
        self.assertGreater(self.result.GetNumConformers(), 1)
        for conf in self.result.GetConformers():
            for torsion in self.torsions:
                angle = T.GetDihedralDeg(conf, *torsion) % 360.0
                deviation = min(abs(angle - target) for target in ANGLES + [420.0])
                self.assertLess(deviation, 1e-6, (conf.GetId(), torsion, angle))

    def test_torsion_rotation_keeps_the_geometry_rigid(self):
        reference = self.mol.GetConformer()
        for conf in self.result.GetConformers():
            for bond in self.mol.GetBonds():
                i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                self.assertAlmostEqual(T.GetBondLength(conf, i, j),
                                       T.GetBondLength(reference, i, j), places=5)

    def test_pruned_conformers_are_rmsd_distinct(self):
        # A lower bound that is too large would skip a real duplicate. Check
        # every kept pair with RDKit's symmetry-aware heavy-atom RMSD.
        heavy = Chem.RemoveHs(self.result)
        ids = [conf.GetId() for conf in heavy.GetConformers()]
        for i, j in itertools.combinations(ids, 2):
            rms = rdMolAlign.GetBestRMS(Chem.Mol(heavy), Chem.Mol(heavy), refId=i, prbId=j)
            self.assertGreaterEqual(rms, RMSD - 1e-6, (i, j))


@unittest.skipUnless(CPP_AVAILABLE, 'requires the C++ sampling extension')
class TestSamplingBudgets(unittest.TestCase):
    MALTOSE = 'OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O'

    def test_hydroxyl_budget_is_filled_from_the_cores_found(self):
        # With few distinct cores, the hydroxyl variants per core grow so a
        # small budget is still filled (75 of 100 before the change).
        self.assertGreaterEqual(_sample(self.MALTOSE, 100).GetNumConformers(), 90)

    def test_small_budget_is_kept(self):
        self.assertEqual(_sample(self.MALTOSE, 30).GetNumConformers(), 30)

    def test_random_mode_respects_the_energy_window(self):
        # Conformers accepted while the running minimum was still high must be
        # dropped once the final minimum is known.
        window = 3.0
        result = _sample('CCCCCCOc1ccc(CCN)cc1', 200, mode='random', energywindow=window)
        energies = [float(conf.GetProp('Energy')) for conf in result.GetConformers()]
        self.assertGreater(len(energies), 1)
        self.assertLessEqual(max(energies) - min(energies), window + 1e-6)


if __name__ == '__main__':
    unittest.main()
