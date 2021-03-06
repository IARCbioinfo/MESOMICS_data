# MESOMICS data and phenotypic map

This repository contains data and processing scripts associated with the MESOMICS project from the Rare Cancers Genomics initiative (http://rarecancersgenomics.com).

The main analysis paper preprint is available at: https://www.biorxiv.org/content/10.1101/2021.09.27.461908v1

The copy number changes, mutations, and structural variants matrices are given in the `phenotypic_map/MESOMICS` folder, in the `TableS9_CNVs.xlsx`, `TableS12_SNVs.xlsx`, and `TableS11_SVs.xlsx` files respectively.
The expression matrices of the MESOMICS samples correspond to the `gene_count_matrix_1pass.csv` with the raw read counts and `vstexpr.zip` with the normalized read counts, in the `phenotypic_map/MESOMICS` folder.

The data note describing the data production, validation, and reuse is available at: https://www.biorxiv.org/content/10.1101/2022.07.06.499003v1

The `phenotypic_map` folder contains data and R scripts used to prepare matrices for each omic layer (`Preprocessing_*.r`), as well as scripts to run MOFA and the Pareto analysis (`PhenotypicMap_*.r`) for the three cohorts (MESOMICS, Bueno, TCGA). A R markdown document detailing the analyses for the MESOMICS cohort is available [here](phenotypic_map/MESOMICS/PhenotypicMap_MESOMICS.md), reproducing the discovery of the three MPM tumor phenotypes using MOFA and Pareto:

![MOFA-Pareto](/phenotypic_map/MESOMICS/PhenotypicMap_MESOMICS_files/figure-html/LFarcplot-1.png)

The interactive phenotypic map resulting from these analyses can be explored at https://tumormap.ucsc.edu/?bookmark=746c4bc0e8bc4eb5f280cdd81c7dcc783955faf2e2b493d0d205b7d1e92b98c4.

![tumormap](/tumormap-screen.png)

