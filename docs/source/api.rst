MolSanitizer API Reference
==========================

Msani Module
-----------------

The Msani module is the main interface for users to access the functionalities of MolSanitizer. It provides methods for molecule preparation, filtering. Conformer generator is triggered separately via the ConformerGenerator class.

.. autoclass::  msani.api::Msani
   :members: 

Moltransform Module
-------------------

.. autoclass:: msani.moltransform.tautomerizer::Tautomerizer
   :members: 


.. autoclass:: msani.moltransform.ionizer::Ionizer
   :members:

.. autoclass:: msani.moltransform.neutralizer::Neutralizer
   :members:

Filtering Module
-----------------

.. autoclass:: msani.filtering.filters::Filters
   :members:

Conformers Module
-----------------

.. autoclass:: msani.conformers.conformers::ConformerGenerator
   :members:

.. autofunction:: msani.conformers.conformers::gen_conf_chunk

.. autoclass:: msani.conformers.mol2writer::Mol2Writer
   :members:

.. autoclass:: msani.conformers.torsions::TorsionLibrary
   :members:
