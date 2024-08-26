# Changelog

All notable changes to this project will be documented in this file.

## [unreleased]

### 🚀 Features

- Msani_batch will now ask the user to confirm to remove the folder before removing it + skip the jobs with more than 1000 subjobs - ([9a6b76c](https://github.com/Isra3l/MolSanitizer/commit/9a6b76c9c52b4534a1dbfc8a168929b6915cbf86))
- :sparkles: Using srETKDGv3 (small-ring version) to hopefully reduce the failed cases with "boat" conformation of the rings with the previous ETKDGv3 (speciallized for macrocycles) - ([2970f10](https://github.com/Isra3l/MolSanitizer/commit/2970f10515dbf69565183e75660606d27683be44))

### 🐛 Bug Fixes

- :bug: Fix a typo in torsion scan that crash msani - ([4275824](https://github.com/Isra3l/MolSanitizer/commit/4275824384d8567703a5234da77e015561a69e17))

### ⚡ Performance

- :zap: Improved performance for the stochastic sampling, removed RMSD pruning dependent. - ([302e715](https://github.com/Isra3l/MolSanitizer/commit/302e7158a72527bd08ebb2f5c9b8240579c38bd6))

## [0.0.6] - 2024-08-22

### 🚀 Features

- Failed stereoisomers-enumerated compounds should now print to the screen to notify the user - ([36846e1](https://github.com/Isra3l/MolSanitizer/commit/36846e13334c7c290a6620aa16a0ec75f27602c0))
- Changing the default maxAttempts in stochastic sampling for more exhaustive sampling - ([aa88ccf](https://github.com/Isra3l/MolSanitizer/commit/aa88ccfec57bb4dbc8a75d54f317b71168847069))

### ⚡ Performance

- :zap: Efforts to speed up the conformers generator of super-flexible and symmetrical compounds - ([b6a04ad](https://github.com/Isra3l/MolSanitizer/commit/b6a04ad9adf4f988092b6c5af0eed96aede2deff))

### 🎨 Styling

- :art: Improved logging of the time of running of each step of MolSanitizer (should now output hours:mins:secs) - ([a3ff715](https://github.com/Isra3l/MolSanitizer/commit/a3ff715dc9ed4b16f84a690d0751e954c74e24a3))
- Fix typos - ([e51eefc](https://github.com/Isra3l/MolSanitizer/commit/e51eefc47099fe49ccabe0598e260e4cc387de5d))

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

- *(Strain_filter now has its own standalone script!)* :zap: The strain_filters now can be called by command 'strain -i examples.mol2' - ([60a7958](https://github.com/Isra3l/MolSanitizer/commit/60a795852eb6cea3283528b22d75dfb85f0e8b28))
- *(Strain_filter now has its own standalone script!)* :zap: The strain_filters now can be called by command 'strain -i examples.mol2' - ([f05bf9b](https://github.com/Isra3l/MolSanitizer/commit/f05bf9b754f0ce49d239e2f258f4284147dcdd73))

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
