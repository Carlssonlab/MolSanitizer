import pandas as pd
import argparse
from rdkit import Chem


class Neutralizer:
    """
    A class for neutralizing molecules.
    """

    @staticmethod
    def neutralize_mol(mol: Chem.Mol) -> Chem.Mol:
        """
        Neutralize the input molecule by balancing the charges on atoms.

        This function identifies charged atoms in the molecule and neutralizes them
        by adding or removing hydrogens as appropriate.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The input molecule to be neutralized

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The neutralized molecule, or the original molecule if no neutralization is needed

        Notes
        -----
        Adapted from RDKit Cookbook: https://rdkit.org/docs/Cookbook.html
        """
        # SMARTS pattern to find charged atoms that can be neutralized
        pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")

        # Find all matching atoms
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]

        # If there are charged atoms to neutralize
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                # Get the atom and its properties
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()

                # Neutralize the atom
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()

            # Convert to SMILES and back to ensure the molecule is valid
            smiles = Chem.MolToSmiles(mol)
            return Chem.MolFromSmiles(smiles)
        else:
            # No charged atoms to neutralize
            return mol

    @classmethod
    def neutralize_df(cls, df,
                      mol_column='mol',
                      smiles_column = 'smiles') -> pd.DataFrame:
        """
        Neutralize molecules in a DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing molecules to neutralize
        mol_column : str, optional
            Name of the column containing RDKit Mol objects, default is 'mol'
        smiles_column : str, optional
            Name of the column containing SMILES strings, default is 'smiles'
        Returns
        -------
        pandas.DataFrame
            DataFrame with neutralized molecules
        """
        df = df.copy()
        if mol_column not in df.columns:
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)
        df[mol_column] = df[mol_column].apply(cls.neutralize_mol)

        return df

    def neutralize_smiles(self, smiles: str) -> str:
        """
        Neutralize a SMILES string by converting it to a molecule and neutralizing it.

        Parameters
        ----------
        smiles : str
            The input SMILES string to neutralize

        Returns
        -------
        str
            The neutralized SMILES string
        """
        mol = Chem.MolFromSmiles(smiles)
        neutralized_mol = self.neutralize_mol(mol)
        return Chem.MolToSmiles(neutralized_mol)

def main():
    parser = argparse.ArgumentParser(description='Protonate a molecule or a file of SMILES')
    parser.add_argument('-i', '--input', default=None, type=str, help='Input file containing SMILES strings and names')
    parser.add_argument('-s', '--smiles', default=None, type=str, help='SMILES string of the molecule')
    parser.add_argument('-o', '--output', default='neutralized_molecules.smi', type=str, help='Output file to save protonated molecules')
    args = parser.parse_args()
    neutralizer = Neutralizer()
    if args.smiles and args.input:
        raise ValueError("Either SMILES or input file must be provided.")

    if args.smiles:
        results = neutralizer.neutralize_smiles(smiles = args.smiles)
        print(results)
    else:
        df = pd.read_csv(args.input, sep =r'\s+', header=None, names=['smiles', 'ids'])
        neutralized_df = neutralizer.neutralize_df(df)
        neutralized_df[['smiles', 'ids']].to_csv(args.output, index=False, header=False, sep =' ')

if __name__ == '__main__':
    main()