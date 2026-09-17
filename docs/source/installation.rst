Installation
=============================

Dependencies
------------

MolSanitizer is built upon the following packages:

- RDKit
- pandas
- numpy


Optional dependencies
---------------------

- `Meeko <https://pypi.org/project/meeko/>`_ is required only for PDBQT output.
  Install it with ``pip install "MolSanitizer[pdbqt]"``.

- `Open Babel <https://openbabel.org/docs/dev/Installation/install.html>`_ is
  required only when ``obabel`` is selected as the embedding method or used as
  the fallback after an RDKit embedding timeout. If it is unavailable after a
  timeout, MolSanitizer records the reason in ``msani_error.err`` and continues
  with the remaining molecules.

Installation
------------

Installation of MolSanitizer can be done either from PyPI or by building from source. The recommended way is to install it from PyPI.

From PyPI
^^^^^^^^^

.. code-block:: console

   $ pip install molsanitizer

In order to also generate PDBQT files, install the optional dependency Meeko:

.. code-block:: console

   $ pip install "molsanitizer[pdbqt]"

From Conda-forge
^^^^^^^^^

.. code-block:: console

   $ conda install -c conda-forge molsanitizer

Building from source
------------------------

.. warning::

    Building MolSanitizer from source is NOT meant to be done by regular users!


We will set up the environment using `Anaconda <https://docs.anaconda.com/anaconda/install/index.html>`_.


.. code-block:: console

   $ git clone https://github.com/carlssonlab/MolSanitizer.git
    

Mac OS and Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For UNIX-based systems (Mac OS and Linux), the installation and building process could be done quite straightforwardly. In the same folder as previous steps, use:

.. code-block:: console
   
   $ conda env create -f MolSanitizer/environment.yml
   $ conda activate msani
   $ pip install -e MolSanitizer

Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The installation in Windows requires the installation of Microsoft Visual Studio (VS) with C++ build tools, which can be downloaded from `here <https://visualstudio.microsoft.com/downloads/?cid=learn-onpage-download-install-visual-studio-page-cta>`_. Scroll down and only Download the option "Build Tools for Visual Studio 2022".

.. image:: _static/VisualStudio.png
   :width: 700px
   :align: center

Then, during the installation, make sure to select the "Desktop development with C++" workload, as shown below:

.. image:: _static/VisualStudio2.png
   :width: 700px
   :align: center

After installing Visual Studio, you can proceed with the installation of MolSanitizer. In the same folder as previous steps, use:

.. code-block:: console

   $ conda env create -f MolSanitizer/environment.yml
   $ conda activate msani
   $ pip install -e MolSanitizer

Testing
-------

MolSanitizer uses `unittest <https://docs.python.org/3/library/unittest.html>`_ for testing. To run the tests, use the following command:

In the same folder as previous steps, use:

.. code-block:: console

   $ python -m unittest MolSanitizer/test/test_msani.py

The test takes around 1-2 minutes to complete.
