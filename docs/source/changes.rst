Changelog
=========

All notable changes to this project will be documented in this file.

[0.4.0] - 2025-08-18
--------------------

🚀 Features
~~~~~~~~~~~

-  Multiple input files can now be imported via a file containing paths
   to these input files (use –input_list or -il). -
   (`7de4546 <https://github.com/phonglam3103/MolSanitizer/commit/7de4546520e4213902310e333ad835a84a5680eb>`__)
-  SMILES strings with name are now supported (for example: msani -s ‘CN
   methylamine’ -prot). Multiple SMILES strings once is also supported
   (msani -s ‘CO methanol’ ‘CCCC’) -
   (`5ab2842 <https://github.com/phonglam3103/MolSanitizer/commit/5ab2842697d748cad58cc2698febee32e8132c0c>`__)
-  Supports the parsing of CLI arguments from a yaml file. This could be
   triggered using the ‘–config’ or ‘-c’ flag. The template could be
   obtained using the ‘–create_config’ flag. -
   (`ceff762 <https://github.com/phonglam3103/MolSanitizer/commit/ceff7623e96734124b2d88e49de261845ae2add6>`__)
-  The customized torsion definition will first try to replace the
   existing ones in the Torlib first before putting it in the top of the
   priority. This aims to ensure that the customized definitions would
   not overwrite the hierrarchy of the Torlib match. -
   (`a763293 <https://github.com/phonglam3103/MolSanitizer/commit/a76329330b127b66e9ab93040ddd737b11a5d426>`__)
-  Added new small ring entries and modified cycloheptadiene entries -
   (`a9823be <https://github.com/phonglam3103/MolSanitizer/commit/a9823be385d2751774eaf5fdd79a5b4fa0996e7a>`__)
-  Don’t enumerate sulfonamides with non-puckering N and merged the
   flippable Cs SMARTS. -
   (`56b518f <https://github.com/phonglam3103/MolSanitizer/commit/56b518f49a78416971e213764461ac6ceb818022>`__)
-  Allow the definition of custom torsions ``-tor`` or ``--torsion``
   flags. Use ``--create_torsion`` to create a template. -
   (`875d1e9 <https://github.com/phonglam3103/MolSanitizer/commit/875d1e928b6f59eccba8dc238132a37e30b58a03>`__)
-  New flags supported: –extended for extended SMILES supports, –torsion
   for customized torsion definitions. -
   (`7d6b4fd <https://github.com/phonglam3103/MolSanitizer/commit/7d6b4fd89cc0a950609e788aac95ba62ca65f351>`__)
-  Support the extraction of the protonation/tautomerization rules and
   allow the use of customized transformation rules. -
   (`85a10ba <https://github.com/phonglam3103/MolSanitizer/commit/85a10ba4fd682f6f88a3d858232da19d384d7710>`__)
-  Hide advanced options in the help text by default, enabled by the -xh
   or –help_advanced flag. -
   (`a11c1de <https://github.com/phonglam3103/MolSanitizer/commit/a11c1dec3b6fdafff227db3de23ea66c2c1340bd>`__)
-  Updated ionizer so that more probable (de)protonation sites get
   (de)protonated first (eg. guanidines > amines > imidazoles) -
   (`b66bed6 <https://github.com/phonglam3103/MolSanitizer/commit/b66bed6cec9001cfed513a09297b16a42f29a225>`__)
-  More general approach for reducing redundant symmetric rotations
   using graph traversing method. -
   (`263a452 <https://github.com/phonglam3103/MolSanitizer/commit/263a4521847b27dd32973c7f4e6ee15d9b97148f>`__)
-  Allow the modification of Dielectric constant from the CLI and the
   API -
   (`d978f58 <https://github.com/phonglam3103/MolSanitizer/commit/d978f58781fd03389f05d9111063340b8a082fb4>`__)
-  Updated new entries for small-ring library. Always enumerate the two
   puckering forms of cycloheptatriene and cyclohepta-1,4-diene-3-sp2 -
   (`d90b73d <https://github.com/phonglam3103/MolSanitizer/commit/d90b73d27b0eef74043ed607e97684e55467a0eb>`__)
-  New symmetric patterns -
   (`4aa9345 <https://github.com/phonglam3103/MolSanitizer/commit/4aa934508efd4786e9865dd52ea091471ab1f770>`__)
-  Improved implementation of OpenBabel-based embedding. CLI version is
   chosen as it is more accurate for embedding of very difficule
   compounds. -
   (`6e05081 <https://github.com/phonglam3103/MolSanitizer/commit/6e05081dd4f679ceff0401d56b70e805dc271a05>`__)

🐛 Bug Fixes
~~~~~~~~~~~~

-  Fix #20, “‘list’ object has no attribute ‘update’” in ionizer. This
   should not crash the whole program. -
   (`e8ddca9 <https://github.com/phonglam3103/MolSanitizer/commit/e8ddca91c79b7f7f9996f6e9562edcef88911d31>`__)
-  Trials to neutralize molecules breaking aromaticity should not fail
   now -
   (`5b07800 <https://github.com/phonglam3103/MolSanitizer/commit/5b0780016ff29bd0229ca654c6aa28b4f2f335d8>`__)
-  Unify the score of benzene to other aromatic rings. This should fix
   the prioritization of benzene over pyridine. -
   (`f5e60f1 <https://github.com/phonglam3103/MolSanitizer/commit/f5e60f166ac56eed9e56802e2e5214c7bd10669c>`__)
-  Fix a bug that the check of chiralities did not work -
   (`4a73775 <https://github.com/phonglam3103/MolSanitizer/commit/4a7377579114cdfc2df9e10b1e46dd0a5bb867d3>`__)
-  Updated p-EWG phenol rule, remove redundant duplicated rules. -
   (`dde92a8 <https://github.com/phonglam3103/MolSanitizer/commit/dde92a8b57dd3fda15e0093902f70e7d99a26ce7>`__)

🚜 Refactor
~~~~~~~~~~~

-  Refactored Torlib and SrLib into classes that can accept or add
   customized rules, supporting the customized rules. -
   (`c0fd9a4 <https://github.com/phonglam3103/MolSanitizer/commit/c0fd9a40ac4b53999c94e31f1aa5b95be2058482>`__)
-  Updated imports according to PEP8 -
   (`bc857fb <https://github.com/phonglam3103/MolSanitizer/commit/bc857fbf34d23f23e413e96bc606cd7083e5bf95>`__)
-  Change the name back to MolSanitizer -
   (`20b3af8 <https://github.com/phonglam3103/MolSanitizer/commit/20b3af8abaaac243d148abaebdfb04235ff68245>`__)

.. _section-1:

[0.3.0] - 2025-06-19
--------------------

.. _features-1:

🚀 Features
~~~~~~~~~~~

-  Special treatment for aliphatic Ns and Cs. Reinforcing the equatorial
   direction of the N- and C-substituted cyclohexane -
   (`c25fab4 <https://github.com/phonglam3103/MolSanitizer/commit/c25fab44926168c43025c0494eb1511dfbac4e92>`__)
-  Make CLI for confgen -
   (`741ec82 <https://github.com/phonglam3103/MolSanitizer/commit/741ec82b0753fd0985d690e74e78ed16337b52ef>`__)
-  Incorporate molvs tautomer enumerator into tautomerizer. This allows
   us more freedom to modify the scoring function. -
   (`0416ce3 <https://github.com/phonglam3103/MolSanitizer/commit/0416ce35caa865d0a4ce7023a68e0a23a06453f0>`__)
-  New tautomerization rules version 3 with improved performance. -
   (`cc9b224 <https://github.com/phonglam3103/MolSanitizer/commit/cc9b2244f4c307456fa0b1c0bafdf2cc6d67ee1c>`__)
-  New ‘extensive’ method, allowing for sampling of the same peak more
   than once. This is specially useful for molecules with small degree
   of freedom like fragments. Trigger via ‘-mode extensive’. -
   (`cd50ff9 <https://github.com/phonglam3103/MolSanitizer/commit/cd50ff984f32b5d26e51612f30d8fa778836f51b>`__)
-  New filter by the number of unspecified chiral centers. -
   (`f1654a8 <https://github.com/phonglam3103/MolSanitizer/commit/f1654a8a6c09984d2661c1bee82f586636572f5f>`__)
-  Expose numcores to smi2db2. This can help to speed up the initial
   embedding for a molecule if the user asks for more cores. -
   (`08ccdd7 <https://github.com/phonglam3103/MolSanitizer/commit/08ccdd79ff6746991e1d00088319c155ccf66523>`__)
-  Expose the options for the number of ring conformations (-nr) for
   Smi2DB2. -
   (`477c7a1 <https://github.com/phonglam3103/MolSanitizer/commit/477c7a152df9e65fc7a2be45f545887c5872632c>`__)
-  New –rigid flag, allowing the constrained docking concept. -
   (`f71b4aa <https://github.com/phonglam3103/MolSanitizer/commit/f71b4aa40f5788774290160631f8a0a8cfa7688f>`__)
-  Brings Smi2DB2 to Windows! Updated instructions on how to setup AMSOL
   for DB2 generation with Windows in the amsol folder. -
   (`c88df20 <https://github.com/phonglam3103/MolSanitizer/commit/c88df20e7a7c1875699fcec36a3fffddb8ef67a9>`__)
-  Updated the new Ionizer class. This allows the Ionizer module to be
   used as a Python API from other programs. -
   (`25e698a <https://github.com/phonglam3103/MolSanitizer/commit/25e698a89cb383309038e4cf6fa05011cac5dea1>`__)
-  Updated new Tautomerizer class and new tautomers_v2 library, allowing
   for the enumeration of tautomers. -
   (`60a4c1a <https://github.com/phonglam3103/MolSanitizer/commit/60a4c1abdbe035addbaa336c57726131fe2d182b>`__)

.. _bug-fixes-1:

🐛 Bug Fixes
~~~~~~~~~~~~

-  Fix a bug that cannot produce small aliphatic compounds. -
   (`4f67f20 <https://github.com/phonglam3103/MolSanitizer/commit/4f67f203dcbbed606cdbc160a532399cfa60a982>`__)
-  Corrections to barbiturates and hydantoin-like structures planarity -
   (`f204ccc <https://github.com/phonglam3103/MolSanitizer/commit/f204ccc4bd6923dee5cb0fa63d14020cd88531cd>`__)
-  Fallback to RDKit tautomerizer with the customized scoring function.
   Now it requires that the RDKit version over 2024.9.3 to run the
   tautomerizer. -
   (`1f9ae46 <https://github.com/phonglam3103/MolSanitizer/commit/1f9ae463241dc165f05fa62149f823f3b532504b>`__)
-  Better defined the planarity of non-aromatic 6-membered rings and
   cleaned up the SMARTS patterns -
   (`fa3d1c5 <https://github.com/phonglam3103/MolSanitizer/commit/fa3d1c50dc0e8870c9d4e7e1ca14ac1825ba92c9>`__)
-  Fix a bug in the batchmode happening when take input files from the
   other directory than the working one. -
   (`6043b3a <https://github.com/phonglam3103/MolSanitizer/commit/6043b3a5912e7f486c405a25618ce89cf2a83d9a>`__)
-  Canonicalize the given –rigid flag if possible to find the correct
   ring. -
   (`4063713 <https://github.com/phonglam3103/MolSanitizer/commit/40637135b3cac9ff5cdfd990abadfcb39df606d2>`__)
-  :bug: Fix a bug in tautomerizing accessing an infinite loop by limit
   to maximum 100 steps of standardization only. -
   (`a250cde <https://github.com/phonglam3103/MolSanitizer/commit/a250cded2d20908f8ca43a8f61a3dc85d1f5f034>`__)
-  Fix a bug so that in smi2pdbqt, it produces boat conformations. -
   (`30d30b9 <https://github.com/phonglam3103/MolSanitizer/commit/30d30b98429d735e99e70664db8de2d8831b9cec>`__)
-  Fix typos in constrained minimization so that the hydrogen attached
   to the aromatic nitrogen not being planarized correctly. -
   (`6f3c43f <https://github.com/phonglam3103/MolSanitizer/commit/6f3c43fb83a75d5dbb8399a8015cca0121414abc>`__)
-  Fix a bug in smi2db2 that in rare cases with AMSOL, the script runs
   into a ValueError and stop the whole script. This fix puts that part
   into a Try-except so that it would not crash the whole program. -
   (`4ee6675 <https://github.com/phonglam3103/MolSanitizer/commit/4ee66751b8c9b87b1136a82c99c86b9dd4ce328b>`__)
-  Fix a unexplainable bug in neutralizing molecules make the atoms
   cannot receive or delete substituents ([N+],…) -
   (`089ec52 <https://github.com/phonglam3103/MolSanitizer/commit/089ec52db8ff7540b5ade1f946dbdcf78708e229>`__)
-  Fix a bug in msani_batch that string arguments are not passed to the
   single script properly. -
   (`16df3eb <https://github.com/phonglam3103/MolSanitizer/commit/16df3eb8041000348e4a0de40da9471c4fcc26fb>`__)

.. _refactor-1:

🚜 Refactor
~~~~~~~~~~~

-  New Mol2Writer, better replicating the behaviour of CORINA. Supports
   aromaticity better. -
   (`432f8d9 <https://github.com/phonglam3103/MolSanitizer/commit/432f8d91e9b50c70fc26e6c548431bb8d979ece0>`__)
-  Smi2db2 and smi2pdbqt are now refactored in conformers.
   ConformerGenerator class allows the user to generate different file
   formats easily. -
   (`933eac2 <https://github.com/phonglam3103/MolSanitizer/commit/933eac287f6b4590cc5765a67eebb87a57f305cf>`__)
-  Reorganized the whole package. The organization now should be more
   intuitive. -
   (`6de0647 <https://github.com/phonglam3103/MolSanitizer/commit/6de0647233345c60791cda6ea07349a44d32921a>`__)

⚡ Performance
~~~~~~~~~~~~~~

-  Change the name of the package to EirVS. Once the branch is merged,
   the user is expected to use ``pip uninstall molsanitizer`` and
   ``pip install -e .`` again. -
   (`567a780 <https://github.com/phonglam3103/MolSanitizer/commit/567a78076dbb9858cc361c6aba7906e80ad03b6b>`__)

.. _section-2:

[0.2.3] - 2025-02-04
--------------------

.. _features-2:

🚀 Features
~~~~~~~~~~~

-  By default will filter out the DB2-incompatible atoms such as
   hypervalence S, P as they are not supported by Mol2, and MMFF94. This
   can reduce the power necessary to handle these entries. -
   (`a003301 <https://github.com/phonglam3103/MolSanitizer/commit/a0033010d266f676eeeef514c74c0ae3924b9c39>`__)
-  Disabled minor peaks and corrected the wrongly assigned peaks.
   Resolved TODO list: disable minor peaks of the Torsion Library #6 -
   (`d9bcfec <https://github.com/phonglam3103/MolSanitizer/commit/d9bcfecc2bc51eb52b3bfcc06b2aa087aaf75556>`__)
-  Updated three new descriptor-based filters, hba, hbd, and mw for
   better drug-like and lead-like characterization. -
   (`5015d6e <https://github.com/phonglam3103/MolSanitizer/commit/5015d6ebc35105a07eda3897db2bef7c9a1de63c>`__)
-  Bring PDBQT to MolSanitizer, thanks to the Meeko Library from the
   Forli’s Lab. The user now can trigger it by –pdbqt flag. -
   (`05904d2 <https://github.com/phonglam3103/MolSanitizer/commit/05904d2850df9ec04543af8b08671aa93cfff537>`__)
-  Added the filters by 2D descriptors: number of heavy atoms (HA; –ha)
   and cLogP (–logp) for MolSanitizer. -
   (`fd61e7e <https://github.com/phonglam3103/MolSanitizer/commit/fd61e7e9d52a779da54395143b0a19086540b753>`__)
-  Added the flag ‘–conformal’ (or ‘-cp’) for preparing the format for
   Conformal Predictor script. Only standardize molecules using RDKit
   original function and skip all other flags. -
   (`0eb2678 <https://github.com/phonglam3103/MolSanitizer/commit/0eb2678b8b6cc8cf1a23485a6707cde9d20c688b>`__)
-  Added the –version flag, as well as better versioning control
   strategy. -
   (`eae5fa0 <https://github.com/phonglam3103/MolSanitizer/commit/eae5fa0128a4fbcc64bdf66208e47d4f8d775a0a>`__)
-  Now the final job in the batch array will also check whether all the
   jobs are finished. A file called RESUBMIT_FAILED_JOBS.txt will be
   generated to guide the user resubmitting unfinished jobs. -
   (`e51f7ce <https://github.com/phonglam3103/MolSanitizer/commit/e51f7cefb88d6d6160bf6c9a62ff8fd0869eab00>`__)
-  Msani in batch mode will now check if the number of jobs to be
   submitted will exceed the quota and warn the user. This is to avoid
   potential loss of the DB2 generation. -
   (`4e9bf80 <https://github.com/phonglam3103/MolSanitizer/commit/4e9bf80aed657afaa11f7a3c09b68ce91f73d106>`__)
-  Remove salts option now also only remove small fragments in one entry
   of molecule and only retains the largest one. -
   (`857d30e <https://github.com/phonglam3103/MolSanitizer/commit/857d30e05e0ef60a973e66212682b69ef0c16e70>`__)
-  Initial commit of the new self-written Mol2Writer. The Mol2 file now
   should be intepretable by RDKit and should not confront problem with
   kekulization caused by the old OpenBabel-based approach. -
   (`c6cc898 <https://github.com/phonglam3103/MolSanitizer/commit/c6cc8988b5e8094d67905fd5e1836ee21790f8bd>`__)
-  New protonation rules and protonation method. Msani now could
   enumerate the protonation states at different pH values and within a
   range of pH. -
   (`9a779a2 <https://github.com/phonglam3103/MolSanitizer/commit/9a779a2214159a9d177491ca6b436356cfdb96cc>`__)

.. _bug-fixes-2:

🐛 Bug Fixes
~~~~~~~~~~~~

-  Fixed a bug that the customized filters and neutralize flags are not
   processed. This is an error from the newly refactored script
   SmilesSanitizer on Jan 21 and should not be a problem for the msani
   version before that commit. -
   (`d643b73 <https://github.com/phonglam3103/MolSanitizer/commit/d643b73fa7b648a597737aa950fb72cea0244b33>`__)
-  Fix a bug in applying reactions to modify molecules that once the
   sanitization fail, the molecule is omitted. Now the latest valid
   version of the molecule will be kept. -
   (`f32403b <https://github.com/phonglam3103/MolSanitizer/commit/f32403b918a93ee85dbb7c09348fe1804105ed8d>`__)
-  Only apply GetLargestFragment when there are more than 1 fragment in
   the RDKit Mol object. This could avoid the problem of detaching
   covalently bound metal atoms. -
   (`95862a6 <https://github.com/phonglam3103/MolSanitizer/commit/95862a6ce2c3d09bd9c1b3b58c424c1e0680a426>`__)
-  Fix a bug that CORINA generated conformer could not be processed
   properly by RDKit. Also make the amine in sulfonamide now pyramidal
   (more realistic) instead of planar. -
   (`640cc55 <https://github.com/phonglam3103/MolSanitizer/commit/640cc5561f3597ef0cb6dbb99b7448bc25a07076>`__)

.. _refactor-2:

🚜 Refactor
~~~~~~~~~~~

-  Improved argparser. Added help text for the new descriptor-based
   filters. Added more shortened keywords for the arguments. -
   (`aadcff5 <https://github.com/phonglam3103/MolSanitizer/commit/aadcff52a45184682176717d70ddf32bf8f8048a>`__)
-  Refactored the whole filters.py into the SmilesStandardizer class.
   Now MolSanitizer Smiles section could be used as a Python library. -
   (`17ee8dc <https://github.com/phonglam3103/MolSanitizer/commit/17ee8dcc3e223d98e2a8dff95c8b1008a58dee97>`__)
-  Refactored huge chunks of scripts. Now MolSanitizer can continue to
   handle unexpectedly failed jobs without restarting from the
   beginning. -
   (`0017709 <https://github.com/phonglam3103/MolSanitizer/commit/0017709b92af88a57aed7a13177bb9e6e5c118d0>`__)

📚 Documentation
~~~~~~~~~~~~~~~~

-  Improved logger. Removed information regarding the old ionization
   method. -
   (`9497d19 <https://github.com/phonglam3103/MolSanitizer/commit/9497d19224f416690974b99022d05d7caa31fbe7>`__)

.. _performance-1:

⚡ Performance
~~~~~~~~~~~~~~

-  Improved the rules for possible dihedrals of symmetric 6-membered
   aromatic rings. -
   (`c745454 <https://github.com/phonglam3103/MolSanitizer/commit/c745454c8c7093f5389ecde4b883f15bc2e22c3e>`__)

🎨 Styling
~~~~~~~~~~

-  Improved the format of the logger. -
   (`6ab03f7 <https://github.com/phonglam3103/MolSanitizer/commit/6ab03f7bc9dc2cf5b5e6f7cdcf99c7ce2b4d139a>`__)

🧪 Testing
~~~~~~~~~~

-  Fix a bug in unittest that triggered the wrong argument. -
   (`e67411e <https://github.com/phonglam3103/MolSanitizer/commit/e67411ee1775e72254b613f225d0b1773aeff642>`__)
-  Unittest for DB2 part added. -
   (`6648c36 <https://github.com/phonglam3103/MolSanitizer/commit/6648c3660bbb8ed536e2d2d2a94346f3e418565e>`__)
-  Update the new goldenData for the unittest.py -
   (`dfe57b8 <https://github.com/phonglam3103/MolSanitizer/commit/dfe57b879df9d245741f480df91298e4bc479e09>`__)

.. _section-3:

[0.2.2] - 2024-11-30
--------------------

.. _features-3:

🚀 Features
~~~~~~~~~~~

-  Neutralization filter applied before the tautomerization perception.
   This help to remove any possible charge-related error from the input.
   The user can turn off suchbehavior by -noneu or –noneutralize -
   (`4f7c53f <https://github.com/phonglam3103/MolSanitizer/commit/4f7c53fe299cba0d3fb522a8cb7597c5e41f8e1d>`__)
-  Added -igtor or –långben flags, that will apply stochastic but
   without adhere to any rules in TorLib. This will help to explore the
   conformational space without any constraints to the Torlib. -
   (`fbb57a9 <https://github.com/phonglam3103/MolSanitizer/commit/fbb57a9586866d4de486a9684c6427c49e4db576>`__)
-  Incorporated CORINA as an optional 3D initial embedding machine.
   Enable by -c or –corina, but the user is asked to set a path to
   CORINA in msani_configurations.yaml. -
   (`de95bf7 <https://github.com/phonglam3103/MolSanitizer/commit/de95bf7bfafde8c2306236c4a6990dd01cec0d97>`__)

.. _bug-fixes-3:

🐛 Bug Fixes
~~~~~~~~~~~~

-  Fix a bug that rdkit tautomer canonicalize lose the specified
   stereocenters. -
   (`907e799 <https://github.com/phonglam3103/MolSanitizer/commit/907e7994bb15de84401d6b06fae3f1b970d11d47>`__)

.. _performance-2:

⚡ Performance
~~~~~~~~~~~~~~

-  Removed some redundant rotations such as para-substituted phenyl or
   para-pyridine. This have shown an X2 faster and fewer conformers in
   some cases with symmetric ring-substituents. -
   (`435a0df <https://github.com/phonglam3103/MolSanitizer/commit/435a0df74dd47bec831bbd27aa74f92a59554652>`__)

.. _styling-1:

🎨 Styling
~~~~~~~~~~

-  Reorganize the argument parser of msani. The help text is now easier
   to navigate. -
   (`df23453 <https://github.com/phonglam3103/MolSanitizer/commit/df23453e9a118c47612e2e53240766fb6557823e>`__)
-  Configurations of both single and batch jobs are now in
   msani_configurations.yaml. The user can now change default values of
   Msani with this configuration file, and set the PATH to CORINA if the
   user want to use. -
   (`358e4dd <https://github.com/phonglam3103/MolSanitizer/commit/358e4dd0ce07ca6e3792eb8f4ea11945083555d5>`__)

.. _section-4:

[0.2.1] - 2024-11-18
--------------------

.. _features-4:

🚀 Features
~~~~~~~~~~~

-  Lower down the exhaustiveness as it takes unecessarily long. -
   (`dbb8aae <https://github.com/phonglam3103/MolSanitizer/commit/dbb8aaeb1d9ba8450f7221ecc9f69817d163990a>`__)
-  Now support enrichment mode, not putting the db2 files into db2.tgz.
   Maily for adapting with the internal enrichment pipeline. -
   (`004af69 <https://github.com/phonglam3103/MolSanitizer/commit/004af6938faef4507ee9c32f7229e78471e73b88>`__)
-  Trial of new stochastic sampling method that involved in increased
   sampling but keeping the failure count continously increase without
   being resetted after every good conformer. This shown a 2X
   performance on a small tricky subset. -
   (`b5f8e32 <https://github.com/phonglam3103/MolSanitizer/commit/b5f8e32d1608dc9de3e8ca7be67014f6e7691465>`__)
-  SMILES is now back to the DB2 format! -
   (`0da0468 <https://github.com/phonglam3103/MolSanitizer/commit/0da04682d7cea4588945ee4fddaf5e8f1fb4ff16>`__)
-  New implementation of stochastic sampling that can sample more
   exhaustively the conformational space. This involves sampling more
   and filter later at the end. Msani now doesn’t write out and read in
   intermediate files repeatedly but process in the memory to speed up
   the process. -
   (`6edbfad <https://github.com/phonglam3103/MolSanitizer/commit/6edbfadda576db3c4b819e88aa7881615fb84847>`__)
-  The default maximum stereoisomers to be expanded is now set to 8
   (previously 32). This could always be set specifically by flag
   –max_isomers. -
   (`536f5fe <https://github.com/phonglam3103/MolSanitizer/commit/536f5fe94af181e32a9d5b3ad7d54f11061b61df>`__)
-  Log files and input files now being organized to corresponding
   folders. New cleanup mechanism for msani_batch. -
   (`16c7f11 <https://github.com/phonglam3103/MolSanitizer/commit/16c7f111b43f67c7ec3b60844f89723a55180382>`__)
-  DB2 files now are put in tar.gz by default. The number of files in
   each tar.gz depends on the number of lines per job. -
   (`dbd74a4 <https://github.com/phonglam3103/MolSanitizer/commit/dbd74a42537fed8c4e123f1f255b3debbd57d958>`__)

.. _bug-fixes-4:

🐛 Bug Fixes
~~~~~~~~~~~~

-  The recent increased sampling seems to run too long on some case
   examples. This was the case of very flexbile but contain multiple
   repulsive parts. The current implementation should sampling faster
   for these cases. -
   (`aed5d98 <https://github.com/phonglam3103/MolSanitizer/commit/aed5d98369b116d8a084b01b8cd735802a45e2d7>`__)
-  Fix a bug inheriting from the recent mol2 implementation improvement.
   Mol2 objects are now deepcopy to avoid referencing issues. -
   (`58a0815 <https://github.com/phonglam3103/MolSanitizer/commit/58a081580eea581081b963e6b4512553a2a7eeac>`__)
-  Constraints to MMFF94s to make the N atoms in amide linkages planar.
   -
   (`f1f82b7 <https://github.com/phonglam3103/MolSanitizer/commit/f1f82b7b7705b1bb5e32a3624fa7890e49b5a773>`__)

.. _refactor-3:

🚜 Refactor
~~~~~~~~~~~

-  Refactor the organization of the Mol2 object. -
   (`77b6fed <https://github.com/phonglam3103/MolSanitizer/commit/77b6fed73721a91ff569e1808fba73e7ac03b6fe>`__)
-  Remove deprecated scripts. -
   (`4d82dfa <https://github.com/phonglam3103/MolSanitizer/commit/4d82dfa97a7bf0adb6a11f3c7d6656ad3cd12329>`__)
-  Remove deprecated scripts. -
   (`80f915c <https://github.com/phonglam3103/MolSanitizer/commit/80f915c7187d7d2b7b089f2885765b0f4d85c893>`__)
-  Remove deprecated script that is not used anymore. -
   (`3cfa3b8 <https://github.com/phonglam3103/MolSanitizer/commit/3cfa3b87c545e416eee007c0ca643b3a27e21246>`__)

.. _documentation-1:

📚 Documentation
~~~~~~~~~~~~~~~~

-  Redirect MolSanitizer README to the read-the-docs page. -
   (`6e6bc43 <https://github.com/phonglam3103/MolSanitizer/commit/6e6bc434bc69180c67b24950fb476b21898907ea>`__)

.. _performance-3:

⚡ Performance
~~~~~~~~~~~~~~

-  Modifications to mol2db2 allows the mol2 object being recorded
   on-the-fly rather than reading from the mol2 blocks. OpenBabel is now
   only being used for the initial conversion for AMSOL. -
   (`65eed12 <https://github.com/phonglam3103/MolSanitizer/commit/65eed12479d9d31fc11eeb31d0b40061f59fab5a>`__)

.. _styling-2:

🎨 Styling
~~~~~~~~~~

-  Solv files are now deleted even the nocleanup is used. The .solv file
   is still exist in the db2 folder if the user want to check for the
   partial charges and desolvation penalties. -
   (`b99efdf <https://github.com/phonglam3103/MolSanitizer/commit/b99efdf80ef94561b591f4b8bbd4bb107c33e8e8>`__)

.. _section-5:

[0.2.0] - 2024-11-06
--------------------

.. _features-5:

🚀 Features
~~~~~~~~~~~

-  Updated new rules for aromatic hydroxyls to make them coplanar with
   the aromatic rings. -
   (`b240a29 <https://github.com/phonglam3103/MolSanitizer/commit/b240a29fd03dde6ccd64da19dc1a7b79f86d7f0e>`__)
-  Initial implementation of OpenBabel 3D embedding for faster initial
   embedding process. Set the default timeout to 2 minutes as surveyed
   from the Tetralith clusters. -
   (`056270f <https://github.com/phonglam3103/MolSanitizer/commit/056270f5acb1205d84e10a81b87824e9fba80cf6>`__)
-  New default energywindow is 25 kcal/mol as we found that this window
   could compromise the accuracy in terms of both redocking and
   enrichment. -
   (`2241d1a <https://github.com/phonglam3103/MolSanitizer/commit/2241d1a0f34bdc7ec480f7b641c09adebdd14cb4>`__)
-  Implemented parallelization for tautomers and stereoisomers options -
   (`5996a32 <https://github.com/phonglam3103/MolSanitizer/commit/5996a3231cca650daa44fbe834fb8c9bceb80f5e>`__)
-  Improved the initial conformations of conjugated Ns in heterocyclics
   inherited from using MMFF94s forcefield. Now these heterocycles
   should be planar. -
   (`3660f8b <https://github.com/phonglam3103/MolSanitizer/commit/3660f8b30fdb1ca59bda1b24e2bf8f6f8f425b47>`__)
-  New mechanism of running AMSOL to avoid shell piping issues. -
   (`78f2176 <https://github.com/phonglam3103/MolSanitizer/commit/78f2176fd9c3c715ac9a6864a8a0ebbc0a55ce5c>`__)
-  New mechanism of calculating maximum possible stereoisomers based on
   unassigned chiral centers -
   (`142a3f6 <https://github.com/phonglam3103/MolSanitizer/commit/142a3f6ff7ab51e25455a069aaba6e7d8566d7ca>`__)
-  New cleanup method to support non-SLURM jobs -
   (`c89e127 <https://github.com/phonglam3103/MolSanitizer/commit/c89e127a5b301ce12c90311cf281b2aa82af86dd>`__)
-  Msani now supports the multithreading for stereoisomers enumeration
   and set the time out for this process for each entry to 1 minute
   only. -
   (`d5d4c9e <https://github.com/phonglam3103/MolSanitizer/commit/d5d4c9e7957ec31b386204894ef91d7b81285943>`__)
-  Msani now only allows up to 4 minutes in the initial embedding stage.
   This is to avoid compounds that take too long for embedding that are
   likely because of the error in the SMILES level. -
   (`7c66150 <https://github.com/phonglam3103/MolSanitizer/commit/7c6615084d948b6e2f2e362e8fc7d421ba7c2fdc>`__)
-  MolSanitizer now suggests the user to update rdkit to avoid known
   errors with stereoisomers and tautomers. -
   (`63750b3 <https://github.com/phonglam3103/MolSanitizer/commit/63750b3d52f3b12ac3a85f44ec7c1bfae015f2ae>`__)

.. _bug-fixes-5:

🐛 Bug Fixes
~~~~~~~~~~~~

-  Fix a bug in run_amsol that makes msani proceed although AMSOL
   failed. -
   (`ec4210c <https://github.com/phonglam3103/MolSanitizer/commit/ec4210cb76969f2cb021bd689893d954120f54d1>`__)
-  Fix a bug that the DB2 file loses the information about the input
   names -> make all the DB2 files have the same name as \****\* -
   (`755d696 <https://github.com/phonglam3103/MolSanitizer/commit/755d69641b1eb5df29a70b9d569e3b3a9c3f94d1>`__)

.. _refactor-4:

🚜 Refactor
~~~~~~~~~~~

-  Remove deprecated functions -
   (`9bc63b6 <https://github.com/phonglam3103/MolSanitizer/commit/9bc63b6fde4568f4e83a67823fe0177110cf4773>`__)

.. _section-6:

[0.1.3] - 2024-10-05
--------------------

.. _features-6:

🚀 Features
~~~~~~~~~~~

-  MolSanitizer now will skip generating DB2 file if the file already
   exist. -
   (`52d7a40 <https://github.com/phonglam3103/MolSanitizer/commit/52d7a4044d03276993b1e6061309f110d09606d4>`__)
-  Warn the user if not all the stereoisomers are written out. -
   (`1e56118 <https://github.com/phonglam3103/MolSanitizer/commit/1e561180b912a98af541163c07af701a011aea2e>`__)
-  New default values of energywindow=15 and max_isomers=32 (max
   stereoisomers to be enumerated) -
   (`d901665 <https://github.com/phonglam3103/MolSanitizer/commit/d901665b804bfb5e7fd0842b08731e7f6e483c38>`__)
-  :bug: New cleanup mechanism for sessions not running in a SLURM job.
   -
   (`2ae700a <https://github.com/phonglam3103/MolSanitizer/commit/2ae700a19d9141e15b9371f77a4fb8418ba5b6cf>`__)
-  Only commit CHANGELOG.md when CHANGELOG.md contains differences. -
   (`5f87498 <https://github.com/phonglam3103/MolSanitizer/commit/5f87498b2854b657766719a6a18162ad4ea97acd>`__)
-  New msani_batch interface, showing the user how many jobs prior to
   submission. -
   (`fcd9755 <https://github.com/phonglam3103/MolSanitizer/commit/fcd9755fc37a971785091defa73232fd3171a2d6>`__)
-  :bug: Update new stereoisomers and tautomers expansion name patterns.
   -
   (`239b92a <https://github.com/phonglam3103/MolSanitizer/commit/239b92aecf9f2146c151e0dab0d4ec0b9ec48133>`__)
-  New alignment rules for non-ring compounds -
   (`c2376ac <https://github.com/phonglam3103/MolSanitizer/commit/c2376acd3eb9c75e01787fa9d70c352c660e4907>`__)
-  Reduced sampling for non-ring-containing molecules to mimic the
   behavior of DB2Pipeline. -
   (`5c55c43 <https://github.com/phonglam3103/MolSanitizer/commit/5c55c433eb48cbbc77781758785105d727fef08a>`__)
-  New cleanup mechanism updated -
   (`727c5b6 <https://github.com/phonglam3103/MolSanitizer/commit/727c5b6c60c530da062b784a35e122f042417b82>`__)
-  New cleanup mechanism so one job should not interfere other parallel
   jobs (on SLURM system). -
   (`fbfe34a <https://github.com/phonglam3103/MolSanitizer/commit/fbfe34ab2c92a4d3d3b0f124c11a2498ccaca66f>`__)
-  Implementation of energy calculation for conformers and use
   energywindow to remove unfavorable conformers. -
   (`6fc4242 <https://github.com/phonglam3103/MolSanitizer/commit/6fc4242d83293dd18ba4456bc05a7526f4da6a7a>`__)
-  Added the new parameter: energywindow to avoid unreasonable
   conformations -
   (`658d08c <https://github.com/phonglam3103/MolSanitizer/commit/658d08ce81b9f8d25c530b6063bffb3d0f8388ad>`__)

.. _bug-fixes-6:

🐛 Bug Fixes
~~~~~~~~~~~~

-  New cleanup mechanism, which should now cleanup even with parallel
   jobs of different array_id being running simultaneously. -
   (`0bb2bc9 <https://github.com/phonglam3103/MolSanitizer/commit/0bb2bc9896907c3903425d11238429cdabd3fe68>`__)
-  Fix a bug in stereoisomers expansion -
   (`8f530c1 <https://github.com/phonglam3103/MolSanitizer/commit/8f530c1ee8bea97589514c48d1c077874805a863>`__)
-  Compounds that fail to tautomerize should not interrupt the whole
   msani for now. If error in generating stereoisomers or tautomers
   occurs, the smiles should be kept as input rather than skipping it in
   the earlier version. -
   (`e17a0a1 <https://github.com/phonglam3103/MolSanitizer/commit/e17a0a13189a3c17fcf0faf3000fd932e46dfc75>`__)

.. _refactor-5:

🚜 Refactor
~~~~~~~~~~~

-  Remove unused codes -
   (`8437f18 <https://github.com/phonglam3103/MolSanitizer/commit/8437f18d4afe59d018dc6b7d7a04f7e659898a1b>`__)

.. _section-7:

[0.1.2] - 2024-09-26
--------------------

.. _features-7:

🚀 Features
~~~~~~~~~~~

-  Msani not use the reset terminal hydrogen of mol2db2 anymore. -
   (`f4d2d6e <https://github.com/phonglam3103/MolSanitizer/commit/f4d2d6ec6b870f6a24fe4960c3622d983151de04>`__)

.. _bug-fixes-7:

🐛 Bug Fixes
~~~~~~~~~~~~

-  The enumerated stereoisomers in the db2 part should also be output to
   the \_clean.smi file. -
   (`1c12e74 <https://github.com/phonglam3103/MolSanitizer/commit/1c12e749b211869ca2b91267adde3906884e6251>`__)
-  Disable the default clash checking of mol2db2 program, which could
   make DOCK skips the potential conformations (msani already checked in
   the torsional sampling part). -
   (`09553b3 <https://github.com/phonglam3103/MolSanitizer/commit/09553b388f5567f22461360383aa1cbd96af55e3>`__)
-  Unspecified stereocenters now will be enumerated automatically before
   undergoing conformational embedding. -
   (`e04b6d6 <https://github.com/phonglam3103/MolSanitizer/commit/e04b6d6ff08692ad7c1f31d9fce1899531c81ac5>`__)
-  Fix a bug that generated compounds not containing the name -
   (`8618524 <https://github.com/phonglam3103/MolSanitizer/commit/86185246b4c3ba090ab5e6d08bdc0153a4a6b1de>`__)
-  Try to fix the weird behavior of SLURM where all the entries failed
   (worked with flag –debug) -
   (`069cf1f <https://github.com/phonglam3103/MolSanitizer/commit/069cf1f50736163512f3c4b2777d7595b8cab1a0>`__)
-  Failed initial embedding should not crash the whole session. -
   (`66c818b <https://github.com/phonglam3103/MolSanitizer/commit/66c818b88c7479d5e55d2ee20fada5cee9c03b02>`__)
-  Fix another bug so that the compounds with no Torlib-satisfied
   conformation should output at least one conformation (from rdkit). -
   (`d71ff37 <https://github.com/phonglam3103/MolSanitizer/commit/d71ff37cb3e94234edefbcdfc1f9d1786811b6a1>`__)
-  Fix a bug that make the molecules without any rotatable bonds failed
   to generate DB2 files. -
   (`4b0d04b <https://github.com/phonglam3103/MolSanitizer/commit/4b0d04b56ef7b87a7c799688dcc0201655c15d2f>`__)

.. _refactor-6:

🚜 Refactor
~~~~~~~~~~~

-  Make the script more pythonic, to avoid the speed inconsistent
   between subprocess and os/shutil of python. -
   (`db778dd <https://github.com/phonglam3103/MolSanitizer/commit/db778dd4ca7ab6fd75c488e14640eadc1c2cae6a>`__)
-  Rewrite the main script (molSanitizer.py) to increase readability and
   better timing logging. -
   (`225590d <https://github.com/phonglam3103/MolSanitizer/commit/225590da8d4a62f2b05366e077f935e60cc5f7ef>`__)
-  Refactor the script a little bit. Change rigid_part_rules so at least
   three atoms are matched. -
   (`e060c5a <https://github.com/phonglam3103/MolSanitizer/commit/e060c5aef3bae4e3bb2e259eba901d4232a25ebb>`__)

.. _section-8:

[0.1.1] - 2024-09-22
--------------------

.. _features-8:

🚀 Features
~~~~~~~~~~~

-  The msani_batch now allows setting up default settings using a yaml
   file (batch_configurations.yaml). -
   (`b2badad <https://github.com/phonglam3103/MolSanitizer/commit/b2badad1efad59673e41e9a9ee714824653a712d>`__)
-  Set initial embeddings to 100 to save time and computational cost -
   (`6e1a8b2 <https://github.com/phonglam3103/MolSanitizer/commit/6e1a8b234c7bb9ff689d9760d63817ce489c00be>`__)
-  Trial of using different alignment references and trial of 200
   initial conformations -
   (`ba4b8a1 <https://github.com/phonglam3103/MolSanitizer/commit/ba4b8a120fec799572e4fff6ec2c84aadc375fa2>`__)
-  Trial of using smaller initial embedding to speed up the process -
   (`85cf8e1 <https://github.com/phonglam3103/MolSanitizer/commit/85cf8e1e8a7c722e94f78d214fe022b93c5aa9c7>`__)
-  Trial of using smaller num_confs_ring (1 instead of 10) -
   (`725f2ff <https://github.com/phonglam3103/MolSanitizer/commit/725f2ffe659213e45c1488fa95b0f24a4db20f08>`__)

.. _bug-fixes-8:

🐛 Bug Fixes
~~~~~~~~~~~~

-  Fix an error that find_sulfonamide not function as expected -
   (`1818ea7 <https://github.com/phonglam3103/MolSanitizer/commit/1818ea71c6b8856d0603f125c5860639d09886ab>`__)

.. _refactor-7:

🚜 Refactor
~~~~~~~~~~~

-  Remove unused parameters (rmsd) -
   (`19bbd40 <https://github.com/phonglam3103/MolSanitizer/commit/19bbd4067fdd2ba918d7534c9eabacef23e9d00d>`__)
-  Remove unused files in the repository -
   (`744f694 <https://github.com/phonglam3103/MolSanitizer/commit/744f694c98720177145d3d3edeeefa29d729a7ae>`__)

.. _documentation-2:

📚 Documentation
~~~~~~~~~~~~~~~~

-  Update README to match the method implemented in smi2db2 -
   (`36270e6 <https://github.com/phonglam3103/MolSanitizer/commit/36270e61267e56bebb452c2231817d676cfead1a>`__)

◀️ Revert
~~~~~~~~~

-  Revert back to 300 initial conformations for better performance -
   (`31fabcb <https://github.com/phonglam3103/MolSanitizer/commit/31fabcb4e8f238f691c27a2cd518e653e37fb85f>`__)

.. _section-9:

[0.1.0] - 2024-09-17
--------------------

.. _features-9:

🚀 Features
~~~~~~~~~~~

-  Updated new rules and merged the SMARTS -
   (`217b61c <https://github.com/phonglam3103/MolSanitizer/commit/217b61cd2d65fbe1f3e8589c1d5f7c52208b7dc2>`__)
-  Try to implement rotating hydrogen within stochastic sampling to
   increase diversity and speed up the mol2db2 process -
   (`4c6d05a <https://github.com/phonglam3103/MolSanitizer/commit/4c6d05a3a5237f6cf85dbc7fcf66c1b4d454b42f>`__)
-  :zap: Boost the performance of stochastic sampling by switching
   between the two modes, based on the relationship between number of
   possible conformations and number of allowed conformations. -
   (`a4e7a57 <https://github.com/phonglam3103/MolSanitizer/commit/a4e7a57dcb828759d54c4178f044c15b1151f91b>`__)
-  Added timing feature for mol2db2 workflow -
   (`e38916e <https://github.com/phonglam3103/MolSanitizer/commit/e38916e5175263aa58123ff6703a4246baa73d3c>`__)
-  :sparkles: Small-ring Torlib updated! Msani should now produce up to
   10 (and favorable) rigid scaffolds based on the new SR-Torlib! -
   (`e33139e <https://github.com/phonglam3103/MolSanitizer/commit/e33139e1f5223c8a84c037b7cf252a621588b132>`__)
-  Small-ring Torlib updated! Msani should now produce up to 10 (and
   favorable) rigid scaffolds based on the new SR-Torlib! -
   (`fcad867 <https://github.com/phonglam3103/MolSanitizer/commit/fcad86777f0ef5bb3dc18c42d9723b88e96279e0>`__)
-  Now supports upto 8-membered ring as rigid part in smi2db2 part -
   (`de62a99 <https://github.com/phonglam3103/MolSanitizer/commit/de62a9940b30ba6d0e0770aee225ba3271933e7d>`__)
-  Added the debug mode for testing on large scale -
   (`7b304e9 <https://github.com/phonglam3103/MolSanitizer/commit/7b304e9bebf885c46f5f2158e75ae0df6947aaa3>`__)
-  Added an epsilon values so that angle scores at 0 can still have the
   possibility to sample -
   (`6afbc63 <https://github.com/phonglam3103/MolSanitizer/commit/6afbc638f73949e1cff8a9c2cff36a37c51eba4c>`__)
-  First effort to embed multiple ring conformations and cover multiple
   regioisomers of sulfonamide-like structures -
   (`afd59b1 <https://github.com/phonglam3103/MolSanitizer/commit/afd59b1294846c3346f77c0684d6a769a36075e1>`__)

.. _bug-fixes-9:

🐛 Bug Fixes
~~~~~~~~~~~~

-  Removed meaningless rules, updated timing and catch an exception
   where no good conformations could be found (fused-ring systems) -
   (`d73bc8e <https://github.com/phonglam3103/MolSanitizer/commit/d73bc8e3559175e3daa7130e53e54c6b80f7678e>`__)

.. _section-10:

[0.0.7] - 2024-09-01
--------------------

.. _features-10:

🚀 Features
~~~~~~~~~~~

-  *(install)* Added toml file and fixed null arguments -
   (`61c1380 <https://github.com/phonglam3103/MolSanitizer/commit/61c138077348b74af345a29aa34ef87613ce357f>`__)
-  :sparkles: Using srETKDGv3 (small-ring version) to hopefully reduce
   the failed cases with “boat” conformation of the rings with the
   previous ETKDGv3 (speciallized for macrocycles) -
   (`2970f10 <https://github.com/phonglam3103/MolSanitizer/commit/2970f10515dbf69565183e75660606d27683be44>`__)
-  Msani_batch will now ask the user to confirm to remove the folder
   before removing it + skip the jobs with more than 1000 subjobs -
   (`9a6b76c <https://github.com/phonglam3103/MolSanitizer/commit/9a6b76c9c52b4534a1dbfc8a168929b6915cbf86>`__)

.. _bug-fixes-10:

🐛 Bug Fixes
~~~~~~~~~~~~

-  Fix a bug so that MolSanitizer batch mode still runs although the
   user asked for not to. -
   (`b518b03 <https://github.com/phonglam3103/MolSanitizer/commit/b518b03479b7441ed41b1829e1c3a82849d57d11>`__)
-  :bug: Fix a typo in torsion scan that crash msani -
   (`4275824 <https://github.com/phonglam3103/MolSanitizer/commit/4275824384d8567703a5234da77e015561a69e17>`__)

.. _performance-4:

⚡ Performance
~~~~~~~~~~~~~~

-  :zap: Improved performance for the stochastic sampling, removed RMSD
   pruning dependent. -
   (`302e715 <https://github.com/phonglam3103/MolSanitizer/commit/302e7158a72527bd08ebb2f5c9b8240579c38bd6>`__)

.. _section-11:

[0.0.6] - 2024-08-22
--------------------

.. _features-11:

🚀 Features
~~~~~~~~~~~

-  Changing the default maxAttempts in stochastic sampling for more
   exhaustive sampling -
   (`aa88ccf <https://github.com/phonglam3103/MolSanitizer/commit/aa88ccfec57bb4dbc8a75d54f317b71168847069>`__)
-  Failed stereoisomers-enumerated compounds should now print to the
   screen to notify the user -
   (`36846e1 <https://github.com/phonglam3103/MolSanitizer/commit/36846e13334c7c290a6620aa16a0ec75f27602c0>`__)

.. _performance-5:

⚡ Performance
~~~~~~~~~~~~~~

-  :zap: Efforts to speed up the conformers generator of super-flexible
   and symmetrical compounds -
   (`b6a04ad <https://github.com/phonglam3103/MolSanitizer/commit/b6a04ad9adf4f988092b6c5af0eed96aede2deff>`__)

.. _styling-3:

🎨 Styling
~~~~~~~~~~

-  Fix typos -
   (`e51eefc <https://github.com/phonglam3103/MolSanitizer/commit/e51eefc47099fe49ccabe0598e260e4cc387de5d>`__)
-  :art: Improved logging of the time of running of each step of
   MolSanitizer (should now output hours:mins:secs) -
   (`a3ff715 <https://github.com/phonglam3103/MolSanitizer/commit/a3ff715dc9ed4b16f84a690d0751e954c74e24a3>`__)

.. _section-12:

[0.0.5] - 2024-08-21
--------------------

.. _features-12:

🚀 Features
~~~~~~~~~~~

-  Adopts the same technique of UCSF for rescaling the number of confs
   generated -
   (`01281aa <https://github.com/phonglam3103/MolSanitizer/commit/01281aa690dcca0b0e56ac19e83fbd8c3557ed09>`__)

.. _bug-fixes-11:

🐛 Bug Fixes
~~~~~~~~~~~~

-  :bug: Remove 5-membered ring as they are not working as expected.
   Added in CC bond as the last resort in case nothing else to align to.
   -
   (`1c9db8d <https://github.com/phonglam3103/MolSanitizer/commit/1c9db8d5fd254125b218aa0e97e783476c0c014f>`__)

.. _section-13:

[0.0.4] - 2024-08-21
--------------------

.. _features-13:

🚀 Features
~~~~~~~~~~~

-  *(smi2db2)* :sparkles: Rigid compounds without any rotatable bonds
   (or with only 1 conf during rotating rot bonds) will output all the
   3D conformations by Rdkit rather than only one like before. eg.
   steroids, morphine…🔥 -
   (`0ff023e <https://github.com/phonglam3103/MolSanitizer/commit/0ff023ed4ee262100fc8baa67865dd9346b457a4>`__)

.. _styling-4:

🎨 Styling
~~~~~~~~~~

-  :fire: Better logger for errorneous compounds -
   (`4627645 <https://github.com/phonglam3103/MolSanitizer/commit/4627645bd555a5b9ae51476762cde4c070003c61>`__)

.. _section-14:

[0.0.3] - 2024-08-20
--------------------

.. _features-14:

🚀 Features
~~~~~~~~~~~

-  *(Added the debug mode for strain_filter; The strained molecules now
   should be stored in another file.)* :zap: -
   (`921c6b9 <https://github.com/phonglam3103/MolSanitizer/commit/921c6b98ff2cbd4bbc3e93e008f8fa60c47f11fe>`__)

.. _bug-fixes-12:

🐛 Bug Fixes
~~~~~~~~~~~~

-  *(smi2db2)* :bug: Fix a bug so that rmsd only comparing between
   heavy_atoms –> boost the performance significantly -
   (`2ab67b2 <https://github.com/phonglam3103/MolSanitizer/commit/2ab67b2d4bc3269186fa2d70e55d860822439ff1>`__)

.. _section-15:

[0.0.2] - 2024-08-19
--------------------

.. _features-15:

🚀 Features
~~~~~~~~~~~

-  *(Strain_filter now has its own standalone script!)* :zap: The
   strain_filters now can be called by command ‘strain -i examples.mol2’
   -
   (`f05bf9b <https://github.com/phonglam3103/MolSanitizer/commit/f05bf9b754f0ce49d239e2f258f4284147dcdd73>`__)
-  *(Strain_filter now has its own standalone script!)* :zap: The
   strain_filters now can be called by command ‘strain -i examples.mol2’
   -
   (`60a7958 <https://github.com/phonglam3103/MolSanitizer/commit/60a795852eb6cea3283528b22d75dfb85f0e8b28>`__)

.. _bug-fixes-13:

🐛 Bug Fixes
~~~~~~~~~~~~

-  *(Fix an error in strain_filter doesnt have main attribute ‘main’)*
   :bug: Reorganizing the main script to the main() function and
   redefine the scope of the Torlib variable -
   (`d91868f <https://github.com/phonglam3103/MolSanitizer/commit/d91868f978de7fd777ff82fe008dec3506b871ba>`__)
-  *(Now MolSanitizer will try different conformations for desolvation
   with AMSOL.)* :sparkles: -
   (`e190e96 <https://github.com/phonglam3103/MolSanitizer/commit/e190e9675a87f9a13161586510ea5d43c0286529>`__)

.. _documentation-3:

📚 Documentation
~~~~~~~~~~~~~~~~

-  *(Better documentation for argparsers)* :memo: -
   (`844e4e3 <https://github.com/phonglam3103/MolSanitizer/commit/844e4e3b43a65af150b92fa95f4b8116a1e3f0b6>`__)
-  *(Better documentations for argsparser)* - Added more details to the
   documentation of the argsparser -
   (`7d81d74 <https://github.com/phonglam3103/MolSanitizer/commit/7d81d74df808404fd85a7a1862f57a4adfea4de2>`__)
-  *(Documentations for the new batch mode of MolSanitizer)* :fire: -
   (`abe3cfc <https://github.com/phonglam3103/MolSanitizer/commit/abe3cfc707dfb5d7e4e48f299080cf37f6d8c347>`__)

.. _styling-5:

🎨 Styling
~~~~~~~~~~

-  :construction: Fix Typos -
   (`e400636 <https://github.com/phonglam3103/MolSanitizer/commit/e400636ea89e660f98c2af31c17c779f0176ce75>`__)

.. _section-16:

[0.0.1] - 2024-08-16
--------------------

Updated
~~~~~~~

-  Stochastic sampling with probs; second tolerance sampling for clash
   compounds; RMSD clustering for stochastic sampling. -
   (`8e63d2c <https://github.com/phonglam3103/MolSanitizer/commit/8e63d2c3e98e268b6e3f3d4e32c0b7ae5cfa8b54>`__)

.. raw:: html

   <!-- generated by git-cliff -->
