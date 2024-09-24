# Changelog

All notable changes to this project will be documented in this file.

## [unreleased - most recent changes come first]

### 🐛 Bug Fixes

- Try to fix the weird behavior of SLURM where all the entries failed (worked with flag --debug) - ([069cf1f](https://github.com/Isra3l/MolSanitizer/commit/069cf1f50736163512f3c4b2777d7595b8cab1a0))
- Failed initial embedding should not crash the whole session. - ([66c818b](https://github.com/Isra3l/MolSanitizer/commit/66c818b88c7479d5e55d2ee20fada5cee9c03b02))
- Fix another bug so that the compounds with no Torlib-satisfied conformation should output at least one conformation (from rdkit). - ([d71ff37](https://github.com/Isra3l/MolSanitizer/commit/d71ff37cb3e94234edefbcdfc1f9d1786811b6a1))
- Fix a bug that make the molecules without any rotatable bonds failed to generate DB2 files. - ([4b0d04b](https://github.com/Isra3l/MolSanitizer/commit/4b0d04b56ef7b87a7c799688dcc0201655c15d2f))

### 🚜 Refactor

- Make the script more pythonic, to avoid the speed inconsistent between subprocess and os/shutil of python. - ([db778dd](https://github.com/Isra3l/MolSanitizer/commit/db778dd4ca7ab6fd75c488e14640eadc1c2cae6a))
- Rewrite the main script (molSanitizer.py) to increase readability and better timing logging. - ([225590d](https://github.com/Isra3l/MolSanitizer/commit/225590da8d4a62f2b05366e077f935e60cc5f7ef))
- Refactor the script a little bit. Change rigid_part_rules so at least three atoms are matched. - ([e060c5a](https://github.com/Isra3l/MolSanitizer/commit/e060c5aef3bae4e3bb2e259eba901d4232a25ebb))

## [0.1.1] - 2024-09-22

### 🚀 Features

- The msani_batch now allows setting up default settings using a yaml file (batch_configurations.yaml). - ([b2badad](https://github.com/Isra3l/MolSanitizer/commit/b2badad1efad59673e41e9a9ee714824653a712d))
- Set initial embeddings to 100 to save time and computational cost - ([6e1a8b2](https://github.com/Isra3l/MolSanitizer/commit/6e1a8b234c7bb9ff689d9760d63817ce489c00be))
- Trial of using different alignment references and trial of 200 initial conformations - ([ba4b8a1](https://github.com/Isra3l/MolSanitizer/commit/ba4b8a120fec799572e4fff6ec2c84aadc375fa2))
- Trial of using smaller initial embedding to speed up the process - ([85cf8e1](https://github.com/Isra3l/MolSanitizer/commit/85cf8e1e8a7c722e94f78d214fe022b93c5aa9c7))
- Trial of using smaller num_confs_ring (1 instead of 10) - ([725f2ff](https://github.com/Isra3l/MolSanitizer/commit/725f2ffe659213e45c1488fa95b0f24a4db20f08))

### 🐛 Bug Fixes

- Fix an error that find_sulfonamide not function as expected - ([1818ea7](https://github.com/Isra3l/MolSanitizer/commit/1818ea71c6b8856d0603f125c5860639d09886ab))

### 🚜 Refactor

- Remove unused parameters (rmsd) - ([19bbd40](https://github.com/Isra3l/MolSanitizer/commit/19bbd4067fdd2ba918d7534c9eabacef23e9d00d))
- Remove unused files in the repository - ([744f694](https://github.com/Isra3l/MolSanitizer/commit/744f694c98720177145d3d3edeeefa29d729a7ae))

### 📚 Documentation

- Update README to match the method implemented in smi2db2 - ([36270e6](https://github.com/Isra3l/MolSanitizer/commit/36270e61267e56bebb452c2231817d676cfead1a))

### ◀️ Revert

- Revert back to 300 initial conformations for better performance - ([31fabcb](https://github.com/Isra3l/MolSanitizer/commit/31fabcb4e8f238f691c27a2cd518e653e37fb85f))

## [0.1.0] - 2024-09-17

### 🚀 Features

- Updated new rules and merged the SMARTS - ([217b61c](https://github.com/Isra3l/MolSanitizer/commit/217b61cd2d65fbe1f3e8589c1d5f7c52208b7dc2))
- Try to implement rotating hydrogen within stochastic sampling to increase diversity and speed up the mol2db2 process - ([4c6d05a](https://github.com/Isra3l/MolSanitizer/commit/4c6d05a3a5237f6cf85dbc7fcf66c1b4d454b42f))
- :zap: Boost the performance of stochastic sampling by switching between the two modes, based on the relationship between number of possible conformations and number of allowed conformations. - ([a4e7a57](https://github.com/Isra3l/MolSanitizer/commit/a4e7a57dcb828759d54c4178f044c15b1151f91b))
- Added timing feature for mol2db2 workflow - ([e38916e](https://github.com/Isra3l/MolSanitizer/commit/e38916e5175263aa58123ff6703a4246baa73d3c))
- :sparkles: Small-ring Torlib updated! Msani should now produce up to 10 (and favorable) rigid scaffolds based on the new SR-Torlib! - ([e33139e](https://github.com/Isra3l/MolSanitizer/commit/e33139e1f5223c8a84c037b7cf252a621588b132))
- Small-ring Torlib updated! Msani should now produce up to 10 (and favorable) rigid scaffolds based on the new SR-Torlib! - ([fcad867](https://github.com/Isra3l/MolSanitizer/commit/fcad86777f0ef5bb3dc18c42d9723b88e96279e0))
- Now supports upto 8-membered ring as rigid part in smi2db2 part - ([de62a99](https://github.com/Isra3l/MolSanitizer/commit/de62a9940b30ba6d0e0770aee225ba3271933e7d))
- Added the debug mode for testing on large scale - ([7b304e9](https://github.com/Isra3l/MolSanitizer/commit/7b304e9bebf885c46f5f2158e75ae0df6947aaa3))
- Added an epsilon values so that angle scores at 0 can still have the possibility to sample - ([6afbc63](https://github.com/Isra3l/MolSanitizer/commit/6afbc638f73949e1cff8a9c2cff36a37c51eba4c))
- First effort to embed multiple ring conformations and cover multiple regioisomers of sulfonamide-like structures - ([afd59b1](https://github.com/Isra3l/MolSanitizer/commit/afd59b1294846c3346f77c0684d6a769a36075e1))

### 🐛 Bug Fixes

- Removed meaningless rules, updated timing and catch an exception where no good conformations could be found (fused-ring systems) - ([d73bc8e](https://github.com/Isra3l/MolSanitizer/commit/d73bc8e3559175e3daa7130e53e54c6b80f7678e))

## [0.0.7] - 2024-09-01

### 🚀 Features

- *(install)* Added toml file and fixed null arguments - ([61c1380](https://github.com/Isra3l/MolSanitizer/commit/61c138077348b74af345a29aa34ef87613ce357f))
- :sparkles: Using srETKDGv3 (small-ring version) to hopefully reduce the failed cases with "boat" conformation of the rings with the previous ETKDGv3 (speciallized for macrocycles) - ([2970f10](https://github.com/Isra3l/MolSanitizer/commit/2970f10515dbf69565183e75660606d27683be44))
- Msani_batch will now ask the user to confirm to remove the folder before removing it + skip the jobs with more than 1000 subjobs - ([9a6b76c](https://github.com/Isra3l/MolSanitizer/commit/9a6b76c9c52b4534a1dbfc8a168929b6915cbf86))

### 🐛 Bug Fixes

- Fix a bug so that MolSanitizer batch mode still runs although the user asked for not to. - ([b518b03](https://github.com/Isra3l/MolSanitizer/commit/b518b03479b7441ed41b1829e1c3a82849d57d11))
- :bug: Fix a typo in torsion scan that crash msani - ([4275824](https://github.com/Isra3l/MolSanitizer/commit/4275824384d8567703a5234da77e015561a69e17))

### ⚡ Performance

- :zap: Improved performance for the stochastic sampling, removed RMSD pruning dependent. - ([302e715](https://github.com/Isra3l/MolSanitizer/commit/302e7158a72527bd08ebb2f5c9b8240579c38bd6))

## [0.0.6] - 2024-08-22

### 🚀 Features

- Changing the default maxAttempts in stochastic sampling for more exhaustive sampling - ([aa88ccf](https://github.com/Isra3l/MolSanitizer/commit/aa88ccfec57bb4dbc8a75d54f317b71168847069))
- Failed stereoisomers-enumerated compounds should now print to the screen to notify the user - ([36846e1](https://github.com/Isra3l/MolSanitizer/commit/36846e13334c7c290a6620aa16a0ec75f27602c0))

### ⚡ Performance

- :zap: Efforts to speed up the conformers generator of super-flexible and symmetrical compounds - ([b6a04ad](https://github.com/Isra3l/MolSanitizer/commit/b6a04ad9adf4f988092b6c5af0eed96aede2deff))

### 🎨 Styling

- Fix typos - ([e51eefc](https://github.com/Isra3l/MolSanitizer/commit/e51eefc47099fe49ccabe0598e260e4cc387de5d))
- :art: Improved logging of the time of running of each step of MolSanitizer (should now output hours:mins:secs) - ([a3ff715](https://github.com/Isra3l/MolSanitizer/commit/a3ff715dc9ed4b16f84a690d0751e954c74e24a3))

## [0.0.5] - 2024-08-21

### 🚀 Features

- Adopts the same technique of UCSF for rescaling the number of confs generated - ([01281aa](https://github.com/Isra3l/MolSanitizer/commit/01281aa690dcca0b0e56ac19e83fbd8c3557ed09))

### 🐛 Bug Fixes

- :bug: Remove 5-membered ring as they are not working as expected. Added in CC bond as the last resort in case nothing else to align to. - ([1c9db8d](https://github.com/Isra3l/MolSanitizer/commit/1c9db8d5fd254125b218aa0e97e783476c0c014f))

## [0.0.4] - 2024-08-21

### 🚀 Features

- *(smi2db2)* :sparkles: Rigid compounds without any rotatable bonds (or with only 1 conf during rotating rot bonds) will output all the 3D conformations by Rdkit rather than only one like before.  eg. steroids, morphine...🔥 - ([0ff023e](https://github.com/Isra3l/MolSanitizer/commit/0ff023ed4ee262100fc8baa67865dd9346b457a4))

### 🎨 Styling

- :fire: Better logger for errorneous compounds - ([4627645](https://github.com/Isra3l/MolSanitizer/commit/4627645bd555a5b9ae51476762cde4c070003c61))

## [0.0.3] - 2024-08-20

### 🚀 Features

- *(Added the debug mode for strain_filter; The strained molecules now should be stored in another file.)* :zap: - ([921c6b9](https://github.com/Isra3l/MolSanitizer/commit/921c6b98ff2cbd4bbc3e93e008f8fa60c47f11fe))

### 🐛 Bug Fixes

- *(smi2db2)* :bug: Fix a bug so that rmsd only comparing between heavy_atoms --> boost the performance significantly - ([2ab67b2](https://github.com/Isra3l/MolSanitizer/commit/2ab67b2d4bc3269186fa2d70e55d860822439ff1))

## [0.0.2] - 2024-08-19

### 🚀 Features

- *(Strain_filter now has its own standalone script!)* :zap: The strain_filters now can be called by command 'strain -i examples.mol2' - ([f05bf9b](https://github.com/Isra3l/MolSanitizer/commit/f05bf9b754f0ce49d239e2f258f4284147dcdd73))
- *(Strain_filter now has its own standalone script!)* :zap: The strain_filters now can be called by command 'strain -i examples.mol2' - ([60a7958](https://github.com/Isra3l/MolSanitizer/commit/60a795852eb6cea3283528b22d75dfb85f0e8b28))

### 🐛 Bug Fixes

- *(Fix an error in strain_filter doesnt have main attribute 'main')* :bug: Reorganizing the main script to the main() function and redefine the scope of the Torlib variable - ([d91868f](https://github.com/Isra3l/MolSanitizer/commit/d91868f978de7fd777ff82fe008dec3506b871ba))
- *(Now MolSanitizer will try different conformations for desolvation with AMSOL.)* :sparkles: - ([e190e96](https://github.com/Isra3l/MolSanitizer/commit/e190e9675a87f9a13161586510ea5d43c0286529))

### 📚 Documentation

- *(Better documentation for argparsers)* :memo: - ([844e4e3](https://github.com/Isra3l/MolSanitizer/commit/844e4e3b43a65af150b92fa95f4b8116a1e3f0b6))
- *(Better documentations for argsparser)* - Added more details to the documentation of the argsparser - ([7d81d74](https://github.com/Isra3l/MolSanitizer/commit/7d81d74df808404fd85a7a1862f57a4adfea4de2))
- *(Documentations for the new batch mode of MolSanitizer)* :fire: - ([abe3cfc](https://github.com/Isra3l/MolSanitizer/commit/abe3cfc707dfb5d7e4e48f299080cf37f6d8c347))

### 🎨 Styling

- :construction: Fix Typos - ([e400636](https://github.com/Isra3l/MolSanitizer/commit/e400636ea89e660f98c2af31c17c779f0176ce75))

## [0.0.1] - 2024-08-16

### Updated

- Stochastic sampling with probs; second tolerance sampling for clash compounds; RMSD clustering for stochastic sampling. - ([8e63d2c](https://github.com/Isra3l/MolSanitizer/commit/8e63d2c3e98e268b6e3f3d4e32c0b7ae5cfa8b54))

<!-- generated by git-cliff -->
