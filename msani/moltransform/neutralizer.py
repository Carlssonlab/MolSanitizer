import argparse
import logging

from pandas import DataFrame, read_csv
from rdkit import Chem
from rdkit.Chem import AllChem

from msani.conformers import utils

logger = logging.getLogger('msani')

neutralize_boronates = AllChem.ReactionFromSmarts('[BX4&-:1][OH]>>[B&+0:1]')

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
        try:
            smiles = Chem.MolToSmiles(mol)
            # First, neutralize boronates if present
            while mol.HasSubstructMatch(neutralize_boronates.GetReactantTemplate(0)):
                new_mol = neutralize_boronates.RunReactants((mol,))[0][0]
                error = Chem.SanitizeMol(new_mol, catchErrors=True)
                if error == 0:
                    mol = new_mol
                else:
                    logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)}")
                    return None
                
            # SMARTS pattern to find charged atoms that can be neutralized
            pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])!$([BX4&-;!$(B-[OH])])]")

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
                error = Chem.SanitizeMol(mol, catchErrors=True)
                if error: 
                    logger.info(f"Error neutralizing molecule: {error}")
                    return None
                smiles = Chem.MolToSmiles(mol)
                return Chem.MolFromSmiles(smiles)
            else:
                # No charged atoms to neutralize
                return mol
        except Exception as e:
            logger.info(f"Error in neutralizing molecule {smiles} ({e})")
            return None

    @classmethod
    def neutralize_df(cls, df,
                      mol_column: str = 'mol',
                      smiles_column: str = 'smiles',
                      name_column: str = 'ids') -> DataFrame:
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
        name_column : str, optional
            Name of the column containing molecule identifiers, default is 'ids'
        Returns
        -------
        pandas.DataFrame
            DataFrame with neutralized molecules
        """
        df = df.copy()
        if mol_column not in df.columns:
            df.loc[:, mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)
        if name_column not in df.columns:
            df.loc[:, name_column] = range(len(df)).astype(str)
        df.loc[:, mol_column] = df[mol_column].apply(cls.neutralize_mol)
        df = log_error_entries(df, mol_column, name_column, smiles_column)
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

def log_error_entries(df: DataFrame,
                      mol_column: str,
                      name_column: str,
                      smiles_column: str) -> DataFrame:
    for _, row in df[df[mol_column].isnull()].iterrows():
        utils.log_error(smiles = row[smiles_column], name = row[name_column])
    df = df[df[mol_column].notnull()]
    return df

def main():
    parser = argparse.ArgumentParser(description='Neutralize a molecule or a file of SMILES')
    parser.add_argument('-i', '--input', default=None, type=str, help='Input file containing SMILES strings and names')
    parser.add_argument('-s', '--smiles', default=None, type=str, help='SMILES string of the molecule')
    parser.add_argument('-o', '--output', default='neutralized_molecules.smi', type=str, help='Output file to save protonated molecules')
    args = parser.parse_args()
    neutralizer = Neutralizer()

    if (args.smiles and args.input) or (not args.smiles and not args.input):
        raise ValueError("Either SMILES or input file must be provided.")

    if args.smiles:
        results = neutralizer.neutralize_smiles(smiles = args.smiles)
        print(results)
    else:
        df = read_csv(args.input, sep =r'\s+', header=None, names=['smiles', 'ids'])
        neutralized_df = neutralizer.neutralize_df(df)
        neutralized_df[['smiles', 'ids']].to_csv(args.output, index=False, header=False, sep =' ')

if __name__ == '__main__':
    main()