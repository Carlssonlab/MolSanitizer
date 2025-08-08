from pandas import DataFrame

def log_error(smiles, name):
    with open('msani_error.err', 'a') as f:
        f.write(f"{smiles} \t {name}\n")

def log_error_mol2(current_comments_str, current_mol2_str):
    with open('strain_error.err', 'a') as f:
        f.write(current_comments_str)
        f.write('\n')
        f.write(current_mol2_str)

def log_error_entries(df: DataFrame,
                      mol_column: str,
                      name_column: str,
                      smiles_column: str) -> DataFrame:
    for _, row in df[df[mol_column].isnull()].iterrows():
        log_error(smiles = row[smiles_column], name = row[name_column])
    df = df[df[mol_column].notnull()]
    return df