from pandas import DataFrame

def log_error(smiles, name, reason=None):
    """Append a failed molecule and an optional explanation to the error log."""
    fields = [str(smiles), str(name)]
    if reason is not None:
        fields.append(" ".join(str(reason).splitlines()))
    with open('msani_error.err', 'a') as f:
        f.write("\t".join(fields) + "\n")

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
