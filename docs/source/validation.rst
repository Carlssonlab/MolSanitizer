Validation
############
.. _validation:

Tautomerizer Validation
========================

The Tautomerizer module of MolSanitizer was validated using a curated subset of 873 entries from TautoBase [1]_ describing tautomeric equilibria measured in aqueous solution. To reduce bias toward the experimentally preferred structure, the minor tautomeric form was used as the input whenever a dominant form could be assigned. Details of dataset curation and the definition of the expected tautomeric states are provided in the Supplementary Information.

Predictions were classified into four outcomes:

- Correct: all experimentally preferred tautomeric forms are captured.
- Incorrect: none of the predicted forms matches the expected result.
- Incomplete: only a subset of the expected tautomeric forms is predicted.
- Exceed: all expected forms are present, but additional non-reference tautomeric forms are also generated.

.. image:: _static/tautomer_validation.png
   :width: 400px
   :align: center

MolSanitizer achieved the highest proportion of correct predictions among the evaluated methods, with approximately 70% of the TautoBase entries classified as correct. ChemAxon cxcalc reached 48%, followed by Epik Classic (23%), Unicon (12%), Scrubber (10%), and Gypsum-DL (7%). MolSanitizer also showed a lower proportion of incorrect predictions (approximately 20%) than the other tested methods. In contrast, Gypsum-DL and Scrubber predominantly over-enumerated tautomeric states, with about 60% and 62% of their predictions, respectively, classified as exceed. These results indicate that MolSanitizer's two-layer tautomerization protocol provides a favorable balance between recovering experimentally preferred states and avoiding unnecessary tautomer enumeration.

Ionizer Validation
========================

The Ionizer module was evaluated using two sets of experimentally characterized compounds with increasing chemical complexity. The first was a monoprotic dataset containing molecules with a single ionizable functional group, assembled from Reaxys, the digitized IUPAC dataset [2]_, and the monoprotic molecule dataset reported by Lee et al. [3]_. Because each molecule contains only one relevant ionizable site, this benchmark allows the expected change in protonation state to be related directly to an experimental pKa value.

The second benchmark contains drug-like molecules from three public pharmaceutical-company datasets: AstraZeneca (AZ) [4]_, Novartis (NV) [5]_, and Vertex (VX) [6]_. These molecules contain more complex combinations of ionizable groups and therefore provide a more realistic test of protonation-state assignment in drug-like chemical space. Dataset preprocessing, filtering, and assignment of the expected protonation states are described in the Supplementary Information.

.. image:: _static/protonation_validation.png
   :width: 800px
   :align: center

Across both the monoprotic and drug-like datasets, MolSanitizer achieved prediction accuracies comparable to the commercial methods and higher than the evaluated freely available methods. MolSanitizer correctly predicted 86% of the monoprotic set and approximately 69--77% of the three drug-like datasets. A fraction of the remaining predictions was classified as incomplete, reflecting the design of MolSanitizer to prioritize a compact set of likely protomers rather than exhaustively enumerate all protonation states that may be populated within the investigated pH range.

Conformer Generator Validation
==============================

The conformer-generation module was evaluated in two complementary benchmarks. Recovery of experimentally observed bioactive conformations was assessed using the Platinum Diverse Dataset [7]_, while the impact of the generated conformers on retrospective virtual-screening performance was evaluated using the DUDE-Z dataset [8]_.

Bioactive pose reproduction
---------------------------

The Platinum Diverse Dataset contains 2,859 high-quality protein-bound ligand conformations derived from X-ray crystal structures in the Protein Data Bank. MolSanitizer was compared with RDKit, Conforge, Conformator, and the conformer-generation protocol used for the ZINC database (referred to here as the reference method) [9]_. For MolSanitizer, both RDKit and CORINA were evaluated as initial 3D embedders. To ensure comparability, the number of generated conformers was capped at 600 per molecule. Performance was measured using the symmetry-corrected heavy-atom RMSD between the generated conformers and the corresponding bioactive structure.

.. image:: _static/Platinum.png
  :width: 800px
  :align: center

RDKit's built-in conformer generator achieved the highest overall coverage, reproducing a bioactive conformer within 2 Å RMSD for 99.7% of the dataset. The differences among the highest-performing methods were small at this threshold: MolSanitizer-RDKit reached 99.5%, Conformator 99.3%, and Conforge 98.7%.

Larger differences were observed at the more stringent 0.5 Å threshold. Conforge achieved the highest recovery at 60%, followed by MolSanitizer-RDKit and Conformator at 57% each, MolSanitizer-CORINA at 50%, the reference method at 48%, and RDKit at 45%. Thus, although RDKit produced the highest overall coverage, MolSanitizer-RDKit more frequently generated conformers very close to the experimentally observed bioactive geometry.

MolSanitizer achieved this recovery with comparatively compact conformer ensembles. The median number of conformers generated per molecule was 115 for MolSanitizer-RDKit and 95 for MolSanitizer-CORINA, compared with 84 for RDKit, 72 for Conforge, 426 for Conformator, and 43 for the reference method.

Computational cost also differed substantially between methods. MolSanitizer-RDKit completed conformational sampling for the full Platinum dataset in 1.51 hours, compared with 2.13 hours for Conforge, 6.02 hours for Conformator, and 30.98 hours for RDKit. Using CORINA as the initial embedder reduced the MolSanitizer runtime further to 0.30 hours while retaining strong bioactive-conformation recovery. Overall, MolSanitizer provides a favorable balance between conformational accuracy, ensemble size, and computational efficiency.

Enrichment capability
---------------------

DUDE-Z is a retrospective molecular-docking benchmark containing 2,312 known active compounds and 69,904 property-matched decoys across 43 protein targets [8]_. The pre-generated DUDE-Z conformers were used as the reference ensembles. For MolSanitizer, conformational ensembles were generated from the corresponding DUDE-Z SMILES using either RDKit or CORINA as the initial embedder. To isolate the contribution of conformer generation, MolSanitizer tautomerization and protonation were not applied in this benchmark. All resulting conformer ensembles were converted to DB2 format and docked using the same receptor grids and DOCK3.8 parameters. Screening performance was quantified using the adjusted logAUC, which emphasizes early recovery of known active compounds.

.. figure:: _static/logauc-alltargets.png
   :width: 800px
   :align: center

MolSanitizer generated higher adjusted logAUC values than the reference conformers for 33 of 43 targets when RDKit was used as the initial embedder and for 31 of 43 targets when CORINA was used. Across all targets, MolSanitizer-RDKit achieved a median adjusted logAUC of 17 and a mean of 22, compared with a median of 15 and a mean of 18 for the reference method. MolSanitizer-CORINA showed similar performance, with a median of 18 and a mean of 21.

The magnitude of the improvement was target-dependent. MolSanitizer-RDKit showed particularly large increases for DEF, NRAM, and XIAP, whereas MolSanitizer-CORINA showed pronounced improvements for PUR2 and DRD4. Conversely, the reference conformers produced higher enrichment for several targets, including FABP4, FA10, and ROCK1. Differences between the two MolSanitizer initial embedders were comparatively small overall, and neither RDKit nor CORINA was consistently superior across all targets. These results show that the conformer-generation strategy can influence retrospective docking enrichment in a target-dependent manner and that MolSanitizer improves early enrichment relative to the reference conformers for the majority of the tested targets.


References
------------
.. [1] Wahl, O.; Sander, T. Tautobase: An Open Tautomer Database. J. Chem. Inf. Model. 2020, 60 (3), 1085-1089. https://doi.org/10.1021/acs.jcim.0c00035.
.. [2] Zheng, J.; Lafontant-Joseph, O.; Green, W. H. Digitized Dataset of Aqueous Acid Dissociation Constants. RSC Adv. 2026, 16 (23). https://doi.org/10.1039/d6ra02418a.
.. [3] Lee, A. C.; Yu, J.; Crippen, G. M. pKa Prediction of Monoprotic Small Molecules the SMARTS Way. J. Chem. Inf. Model. 2008, 48 (10), 2042-2053. https://doi.org/10.1021/ci8001815.
.. [4] Wenlock, M.; Tomkinson, N. Experimental in Vitro DMPK and Physicochemical Data on a Set of Publicly Disclosed Compounds (CHEMBL3301361). https://doi.org/10.6019/CHEMBL3301361.
.. [5] Liao, C.; Nicklaus, M. C. Comparison of Nine Programs Predicting pKa Values of Pharmaceutical Substances. J. Chem. Inf. Model. 2009, 49 (12), 2801-2812. https://doi.org/10.1021/ci900289x.
.. [6] Settimo, L.; Bellman, K.; Knegtel, R. M. A. Comparison of the Accuracy of Experimental and Predicted pKa Values of Basic and Acidic Compounds. Pharm. Res. 2014, 31 (4), 1082-1095. https://doi.org/10.1007/s11095-013-1232-z.
.. [7] Friedrich, N.-O.; de Bruyn Kops, C.; Flachsenberg, F.; Sommer, K.; Rarey, M.; Kirchmair, J. Benchmarking Commercial Conformer Ensemble Generators. J. Chem. Inf. Model. 2017, 57 (11), 2719-2728. https://doi.org/10.1021/acs.jcim.7b00505.
.. [8] Stein, R. M.; Yang, Y.; Balius, T. E.; O'Meara, M. J.; Lyu, J.; Young, J.; Tang, K.; Shoichet, B. K.; Irwin, J. J. Property-Unmatched Decoys in Docking Benchmarks. J. Chem. Inf. Model. 2021, 61 (2), 699-714. https://doi.org/10.1021/acs.jcim.0c00598.
.. [9] Xia, Q.; Fu, Q.; Shen, C.; Brenk, R.; Huang, N. Assessing Small Molecule Conformational Sampling Methods in Molecular Docking. J. Comput. Chem. 2025, 46 (1), e27516. https://doi.org/10.1002/jcc.27516.
