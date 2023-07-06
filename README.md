# MESOMICS data and phenotypic map

This repository contains data and processing scripts associated with the MESOMICS project from the Rare Cancers Genomics initiative (http://rarecancersgenomics.com).

The main analysis paper is available at: http://nature.com/articles/s41588-023-01321-1

The copy number changes, damaging small variants (SNVs and indels), and structural variants matrices are given in the `phenotypic_map/MESOMICS` folder, in the `TableS31-37_CNVs.xlsx`, `TableS44-46_SNVs.xlsx`, and `TableS41-42_SVs.xlsx` files respectively. The folder also contains a maf file with all somatic small variants (damaging and non damaging), var_annovar_maf_corr_allvariants.txt, the exact MOFA inputs (D_alt_MOFA.RData, D_cnv_MOFA.RData, D_exprB_MOFA.RData, D_loh_MOFA.RData, D_met.bodB_MOFA.RData, D_met.enhB_MOFA.RData, D_met.proB_MOFA.RData), and the expression matrices of the MESOMICS samples correspond to the `gene_count_matrix_1pass.csv` with the raw read counts and `vstexpr.zip` with the normalized read counts, in the `phenotypic_map/MESOMICS` folder. There is also a subfolder EGA_metadata with tables to match EGA IDs and file names with MESOMICS sample IDs.

The data note describing the data production, validation, and reuse is available at: https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac128/7007909

The `phenotypic_map` folder contains data and R scripts used to prepare matrices for each omic layer (`Preprocessing_*.r`), as well as scripts to run MOFA and the Pareto analysis (`PhenotypicMap_*.r`) for the three cohorts (MESOMICS, Bueno, TCGA). A R markdown document detailing the analyses for the MESOMICS cohort is available [here](phenotypic_map/MESOMICS/PhenotypicMap_MESOMICS.md), reproducing the discovery of the three MPM tumor phenotypes using MOFA and Pareto:

![MOFA-Pareto](/phenotypic_map/MESOMICS/PhenotypicMap_MESOMICS_files/figure-html/LFarcplot-1.png)

The interactive phenotypic map resulting from these analyses can be explored at https://tumormap.ucsc.edu/?bookmark=746c4bc0e8bc4eb5f280cdd81c7dcc783955faf2e2b493d0d205b7d1e92b98c4.

![tumormap](/tumormap-screen.png)

