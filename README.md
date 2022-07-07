# MESOMICS data and phenotypic map

This repository contains data and processing scripts associated with the MESOMICS project from the Rare Cancers Genomics initiative (http://rarecancersgenomics.com).

The main analysis paper preprint is available at: https://www.biorxiv.org/content/10.1101/2021.09.27.461908v1

The data note describing the data production, validation, and reuse is available at:  

The `phenotypic_map` folder contains data and R scripts used to prepare matrices for each omic layer (`Preprocessing_*.r`), as well as scripts to run MOFA and the Pareto analysis (`PhenotypicMap_*.r`) for the three cohorts (MESOMICS, Bueno, TCGA). A mardkdown document detailing the analyses for the MESOMICS cohort is available [here](phenotypic_map/MESOMICS/PhenotypicMap_MESOMICS.md).

The interactive phenotypic map resulting from these analyses can be explored at https://tumormap.ucsc.edu/?bookmark=746c4bc0e8bc4eb5f280cdd81c7dcc783955faf2e2b493d0d205b7d1e92b98c4.
