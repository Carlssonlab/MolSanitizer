Early Evaluation
=================

MolSanitizer is currently in development and has not yet been publicly released.

We are excited to invite interested users to participate in an early evaluation of our software. Your feedback will play a crucial role in shaping the future of MolSanitizer and ensuring it meets the needs of our community.

Available options
-----------------------------
We are now currently offer working with MolSanitizer via two main sources:

- **Web application**: You can access MolSanitizer through our web application: `MolSanitizer Web App <https://carlssonlabtools.icm.uu.se/molsani>`_. This platform allows you to easily prepare molecules, generate conformers, and create multiple file formats but not capable of high-throughput processing.

- **Compiled software**: For those who prefer a more hands-on approach, we provide a compiled version of MolSanitizer. Please send an email to `phong.lam@icm.uu.se <mailto:phong.lam@icm.uu.se>`_ for access to the compiled software. This version is designed for high-throughput processing and can be integrated into your existing workflows.


We look forward to collaborating with you!

Working with the compiled version
-----------------------------

The compiled version contains of the full functionality of MolSanitizer, including the ability to prepare molecules, generate conformers, and create multiple file formats. We are mainly working forward to improve the predictivity of the protonation and tautomerization modules, therefore the rule libraries used could be generated and customized by the users. 

By default, the help message of MolSanitizer will only produce basic options, but could be extended by using the ``--advanced_help`` flag, where the user can find options to create and run customized rule libraries like:

.. code-block:: console

    $./msani --advanced_help

    Miscellaneous:
    --debug, -d           Enable debugging mode
    --lazy                Implement all the processing and preparation steps
    --numcores, -j        Number of cores to use for parallel processing (default: 4)
    --help, -h            Show this help message and exit
    --help_advanced, -xh  Show advanced help message with additional options
    --version, -v         Show the current version of MolSanitizer
    --create_protlib      Create a template for customized protonation scheme
    --create_taulib       Create a template for customized tautomerization scheme
    --protlib             Path to the protonation library file (default: msani/Data/ionizations_v3.txt).
    --taulib              Path to the tautomer library file (default:  msani/Data/tautomers_v3.txt).


The ``--create_protlib`` and ``--create_taulib`` options will extract the currently compiled protonation and tautomerization libraries, respectively, for the user to customize. The user can then modify the libraries and use them with the ``--protlib`` and ``--taulib`` options.

One other way is to provide us with the examples that the user think to be problematic, and would be best if the user can also provide the expected output and references. We will then use these examples to improve the libraries and the predictivity of the modules.