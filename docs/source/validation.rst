Validation
############
.. _validation:

Tautomerizer Validation
========================

The tautomerization module of MolSanitizer is validated using a subset of 873 tautomer-pair entries in water of the TautoBase [1]_, always starting from the minor tautomeric form, to assess the predictive power of our workflow. For the detailed curation of the dataset, refer to the Supplementary Information. 

Four possible outcomes for each entry could be obtained, including:

- Correct: all the preferred tautomeric forms are captured.
- Incorrect: the predicted forms were not included in the expected results.
- Incomplete: the predicted forms were a part of the expected results.
- Excess: both, the preferred and incorrect tautomeric forms are presented.

.. image:: _static/tautomer_validation.png
   :width: 400px
   :align: center

The results indicate that MolSanitizer's tautomerization module achieves the highest accuracy with 69.4% of the entries being correctly predicted. The module also shows a lowest rate of incorrect predictions (20.4%), while the incomplete and excess predictions are at 8.0% and 2.2%, respectively.

Ionizer Validation
========================

Two datasets are used to validate the ionization module of MolSanitizer. The first dataset contains only molecules that have only one ionizable functional group (monoprotic set), from three sources (https://www.reaxys.com/), the IUPAC digitalized dataset [2]_ , and the monoprotic molecules database by Lee et al [3]_ . The second database compiled drug-like molecules from three pharmaceutical companies AstraZeneca (AZ), [4]_  Novartis (NV), [5]_ Vertex (VX). [6]_ The latter serves as a more complex yet more realistic test cases for all the benchmarked methods.

.. image:: _static/protonation_validation.png
   :width: 440px
   :align: center

For both the monoprotic and drug-like datasets, MolSanitizer reaches the accuracy of commercial software packages and outperforms free software packages. MolSanitizer reached 85.5% correct prediction rate in the monoprotic and 68.6-76.8% in the drug-like datasets.

Conformer Generator Validation
==============================

Validation of the conformational generation part of MolSanitizer has been conducted based on the two datasets. For bioactive pose reproduction, the Platinum Diverse Dataset [7]_ , and for enrichment capability of the active compounds, the DUDE-Z dataset [8]_ were used. 

Bioactive pose reproduction
---------------------------
Platinum Diverse Dataset contains 2859 high-quality ligand bioactive conformations from the Protein Data Bank (PDB). For the current stage of validation, the best aligned conformation from multiple conformer generators (MolSanitizer, RDKit, Conforge, Conformator) and the current DB2 pipeline employed in the ZINC-22 database was used. The number of conformers were capped at 600 conformers to reflect the current setting in the ZINC-22 database. The RMSD values were calculated based on the heavy atoms of the generated and the reference conformer using the RDKit's GetBestRMS function. 


.. image:: _static/Platinum.png
  :width: 800px
  :align: center

Among all the methods tested, RDKit's built-in conformer generator achieved the highest overall coverage, successfully retrieving a bioactive conformer within 2 Å RMSD for 99.7% of the molecules in the benchmark dataset. This strong performance is primarily due to its use of distance-geometry-based embedding, which allows flexible sampling of ring conformations and generates diverse conformers iteratively. However, this method also exhibited notable drawbacks: only 45.3% of the conformers achieved sub-0.5 Å RMSD accuracy, and the total runtime to process the dataset was approximately 34 hours, making it unsuitable for high-throughput or large-scale virtual screening campaigns.
Both Conforge and MolSanitizer achieved comparable accuracy and emerged as the second-best performers (success rate 99.1% and 98.7%, respectively). In the case of MolSanitizer, tuning the dielectric constant (ε) in MMFF94s increased the retrieval rate. A higher ε value reduced intramolecular electrostatic attractions, preventing the system from getting trapped in local minima due to salt bridges, thereby improving conformer diversity and retrieval of bioactive conformations. 

In contrast, the reference method used by ZINC-22 showed the lowest retrieval rate. While based on OMEGA, its performance is limited due to the undersampling of conformers for the rotation of polar hydrogens. Unlike this approach, MolSanitizer explicitly rotates polar hydrogens simultaneously with other rotatable bonds during stochastic sampling, providing greater structural diversity and a higher likelihood of recovering the bioactive state.
Notably, MolSanitizer often generated more conformers per molecule (subplot C). This is a result of its decision not to use RMSD-based clustering for deduplication. During development, we found that heavy-atom-only clustering failed to capture important dihedral variations for molecular docking (e.g., hydroxyl orientations), while hydrogen-inclusive clustering was computationally expensive. Instead, MolSanitizer reduces redundancy through a SMARTS-based pruning strategy and graph traversal algorithms for symmetry perception, effectively minimizing unnecessary conformers while maintaining accuracy.


Enrichment capability
---------------------

DUDE-Z is a comprehensive and challenging test set designed for evaluating molecular docking methods. It includes 2,312 ligands and 69,994 property-matched decoys, covering 43 diverse targets. For this benchmark, all methods were tested with a fixed number of 2,000 conformers. The evaluation metric used is the adjusted Log-AUC, which assesses early enrichment performance as recommended by Stein et al [8]_. 

.. figure:: _static/logauc-alltargets.png
   :width: 800px
   :align: center

In terms of logAUC, MolSanitizer consistently outperforms the reference method in 32 out of 43 targets. In those few cases where MolSanitizer performs worse, the logAUC values of the reference method are also low, except for the Fatty Acid Binding Protein Adipocyte (FABP4). Both the choice of initial embedder and the dielectric constant have a significant influence on enrichment performance. Overall, the RDKit embedder tends to yield better enrichment, albeit with a longer processing time. Nonetheless, for several targets such as Tyrosine-Protein Kinase ABL (ABL1), Coagulation Factor X (FA10), Factor VII (FA7), and Trypsin I (TRY1), CORINA-based conformers result in slightly better performance, although the improvements are not significant.

The use of a higher dielectric constant (ε = 20) during torsional sampling generally enhances enrichment, but not consistently across all targets or embedders. Notably, improvements of 1-3 logAUC units were observed for FABP4, Peptide Deformylase (DEF), Tryptase β1 (TRYB1), and Urokinase (UROK) when increasing ε from 1 to 20. Conversely, a performance drop was noted in targets such as AmpC β-Lactamase (AMPC), GAR Transformylase (PUR2), and Thrombin (THRB). Since increasing the dielectric constant typically results in the generation of more conformers, and therefore more molecules to dock, the computational cost for large-scale screening increases substantially. As such, MolSanitizer sets the dielectric constant to 1 by default, while still allowing users to adjust it via both the command-line interface and the Python API, depending on the desired purposes.



References
------------
.. [1]	Wahl, O.; Sander, T. Tautobase: An Open Tautomer Database. J. Chem. Inf. Model. 2020, 60 (3), 1085-1089. https://doi.org/10.1021/acs.jcim.0c00035.
.. [2] Zheng, J.; Chalk, S.; Lafontant-Joseph, O.; Li, Y. IUPAC Digitized pKa Dataset  v2.2, 2024. https://doi.org/10.5281/zenodo.13987352.
.. [3] Lee, A. C.; Yu, J.; Crippen, G. M. pKa Prediction of Monoprotic Small Molecules the SMARTS Way. J. Chem. Inf. Model. 2008, 48 (10), 2042-2053. https://doi.org/10.1021/ci8001815.
.. [4] Wenlock, M.; Tomkinson, N. Experimental in Vitro DMPK and Physicochemical Data on a Set of Publicly Disclosed Compounds (CHEMBL3301361). https://doi.org/10.6019/CHEMBL3301361.
.. [5] Liao, C.; Nicklaus, M. C. Comparison of Nine Programs Predicting pKa Values of Pharmaceutical Substances. J. Chem. Inf. Model. 2009, 49 (12), 2801-2812. https://doi.org/10.1021/ci900289x.
.. [6] Settimo, L.; Bellman, K.; Knegtel, R. M. A. Comparison of the Accuracy of Experimental and Predicted pKa Values of Basic and Acidic Compounds. Pharm Res 2014, 31 (4), 1082–1095. https://doi.org/10.1007/s11095-013-1232-z.

.. [7] Friedrich, N. O., de Bruyn Kops, C., Flachsenberg, F., Sommer, K., Rarey, M., & Kirchmair, J. (2017). Benchmarking commercial conformer ensemble generators. Journal of chemical information and modeling, 57(11), 2719-2728. Available at: https://pubs.acs.org/doi/10.1021/acs.jcim.7b00505

.. [8] Stein, R. M., Yang, Y., Balius, T. E., O’Meara, M. J., Lyu, J., Young, J., ... & Irwin, J. J. (2021). Property-unmatched decoys in docking benchmarks. Journal of chemical information and modeling, 61(2), 699-714. Available at: https://pubs.acs.org/doi/10.1021/acs.jcim.0c00598