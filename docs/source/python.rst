Python scripting
=================================
.. _python:

MolSanitizer is designed to be used as a Python package, allowing users to integrate its functionalities into their own scripts or applications. The package provides a set of functions and classes that can be imported and utilized for various molecular preparation tasks. There are several classes that can be used by the users, including:

Msani
-----

The **Msani** class is the main entry point for using MolSanitizer. It provides methods for preparing molecules, and applying filters. Users can create an instance of this class and call its methods to perform various tasks.

.. code-block:: python

    import pandas as pd
    from msani.api import Msani

    df = pd.read_csv('input.csv', names=['smiles', 'name']) 

    # Process with msani
    molsani = Msani(removesalts = True,
                    tautomers = True,
                    protonation = True,
                    pH = 7,
                    unwanted = 'regular',
                    pains = True, 
                    stereoisomers=stereoisomers)
    processed_df = molsani.run(df)


For chunked input, use a context manager to share one worker pool across
transformation stages and chunks. The CLI does this automatically. Rules are
cached on the processor and refreshed when its transformation settings or a
custom rule file changes. Workers are closed when the context exits, including
on errors; a standalone ``run()`` closes its workers before returning.

.. code-block:: python

    def main():
        with Msani(tautomers=True, protonation=True, numcores=2) as processor:
            for chunk in pd.read_csv('input.smi', sep=r'\s+',
                                     names=['smiles', 'ids'], chunksize=1000):
                result = processor.run(chunk)
                result[['smiles', 'ids']].to_csv(
                    'prepared.smi', sep=' ', mode='a', header=False, index=False)

    if __name__ == '__main__':
        main()

Without a warm pool, tautomer/protonation batches smaller than 256 rows after
initial filtering run serially to avoid process startup overhead. Set
``parallel_min_rows=0`` on ``Msani`` to always use the requested worker count.
Stereoisomer enumeration retains multiprocessing for small batches because its
per-molecule embedding cost can be high. Processor instances should not be
shared between concurrent threads.


Ionizer
--------

The **Ionizer** class is responsible for handling the ionization of molecules. It provides methods to ionize a molecule from a SMILES string, from RDKit Mol object, or a DataFrame. The ionization process can be customized by specifying the pH value and the ionization library to use.

.. code-block:: python
    
    from rdkit import Chem
    from msani.moltransform.ionizer import Ionizer

    ionizer = Ionizer(pH = 7, pH_range = 0, numcores = 4)
    #ionize a single molecule from SMILES string:
    results = ionizer.ionize(smiles = 'CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1')

    #ionize a single molecule from RDKit Mol object:
    mol = Chem.MolFromSmiles('CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1')
    results = ionizer.ionize(mol = mol)

    #ionize a DataFrame of molecules with SMILES strings:
    df = df = ionizer.ionize_df(mol_df)

Tautomerizer
----------------

The **Tautomerizer** class is responsible for generating tautomers of molecules. It provides methods to generate tautomers from a SMILES string, from RDKit Mol object, or a DataFrame. The tautomerization process can be customized by specifying the tautomer library to use.

.. code-block:: python

    from rdkit import Chem
    from msani.moltransform.tautomerizer import Tautomerizer

    tautomerizer = Tautomerizer(numcores= 4)
    #Tautomerize a molecule from a SMILES
    tautomers = tautomerizer.tautomerize(smiles='c1ccccc1O')

    #Tautomerize a molecule from an RDKit Mol object
    mol = Chem.MolFromSmiles('c1ccccc1O')
    tautomers = tautomerizer.tautomerize(mol = mol)

    #Tautomerize a DataFrame of molecules with SMILES strings:
    tautomers_df = tautomerize_df(df, smiles_column = 'smiles', name_column = 'ids')


Filters
--------

The **Filters** class provides methods for applying various filters to molecules. This includes filtering out unwanted substructures, customized substructures, and PAINS. Filtering on descriptors is also supported. Users can create an instance of this class and call its methods to apply the desired filters.

.. code-block:: python

    from msani.filtering.filters import Filters
    # Unused arguments could be removed, but kept for clarity
    filters = Filters(removesalts=True, ha='>=5',
                      logp='<=3.5', hba='1-3',
                      hbd='1-2', mw='200-500',
                      chiral='0-2', custom='path/to/custom.smarts',
                      unwanted=['regular'], pains=True,
                      rejectedFile='rejected.txt', debug=True)

    filtered_df = filters.filter_df(df)
    
the Filters class also has multiple static methods that can be used independently:

.. code-block:: python

    import pandas as pd
    from msani.filtering.filters import Filters

    df = pd.read_csv('input.csv', names=['smiles', 'name'])
    debug = True
    rejectedFile = 'rejected.txt'

    df = Filters.filter_by_ha(df, '<=30', rejectedFile, debug)
    df = Filters.filter_by_logp(df, '<350', rejectedFile, debug)
    df = Filters.filter_by_hba(df, '<=10', rejectedFile, debug)
    df = Filters.filter_by_hbd(df, '<=5', rejectedFile, debug)
    df = Filters.filter_by_mw(df, '300-500', rejectedFile, debug)
    df = Filters.filter_by_tpsa(df, '<=100', rejectedFile, debug)
    df = Filters.filter_by_fsp3(df, '>=0.25', rejectedFile, debug)
    df = Filters.filter_by_chiralcenters(df, '<2', rejectedFile, debug)
    df = Filters.unwantedFilter(df, rejectedFile, ['regular', 'optional'], debug)
    df = Filters.painsFilter(df, rejectedFile, debug)

ConformerGenerator
----------------------

The **ConformerGenerator** class is responsible for generating conformers of molecules. It provides methods to generate conformers from a SMILES string, from RDKit Mol object, or a DataFrame. The conformer generation process can be customized by specifying the number of conformers to generate and the conformer generation method to use.

Two main steps are involved, first, initial embedding of the molecule (generate 3D), then torsional sampling. Finally, the conformers can be saved in various formats such as PDBQT, SDF, MOL2, DB2, or OEB.

.. code-block:: python

    from msani.conformers.conformers import ConformerGenerator

    smiles = 'CC(=O)C1=CC=CC=C1C(=O)O'
    name = 'test'

    # Initial embedding of the molecule
    confgen = ConformerGenerator(smiles = smiles,
                                 name =  name, 
                                 method = 'rdkit', # Options: 'rdkit', 'openbabel', 'corina'
                                 )

    # Torsional driving
    confgen.conf_sampling(numConfs=2000, energywindow=25)

    # Save conformers to different formats
    confgen.to_pdbqt()
    confgen.to_sdf()
    confgen.to_mol2()
    confgen.to_db2()


With the optional OpenEye toolkits installed and licensed, save the sampled
ensemble as OEB. Create the output directory when using the Python API directly:

.. code-block:: python

    from pathlib import Path
    from contextlib import closing
    from msani.conformers import oebwriter

    Path("oeb").mkdir(exist_ok=True)
    confgen.to_oeb()  # oeb/test.oeb.gz

    # Open once, then append each sampled ConformerGenerator to the library.
    # Opening an existing file truncates it; this API does not resume a run.
    with closing(oebwriter.open_oeb("oeb/library.oeb.gz")) as stream:
        confgen.to_oeb(stream=stream)
