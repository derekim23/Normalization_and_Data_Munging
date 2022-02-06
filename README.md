# Description of Study
The overarching goal of this study is in the realm of multiomics data and how they relate to certain phenotypic measures. All of the underlying data and code has been sanitized.

# Immunophenotypic Data Exploration

This is an exploratory analysis of immunophenotypic protein data across multiple subjects. Measurements were made before and after intense exercise (blood draws 1 vs. 2, 3, and 4). The results show no marked changes in the immunophenotypic profiles post-exercise (vs. pre-exercise).

# Metabolomics Data Merger
The core idea of these files is to merge open-source metabolomics dataset into a metabolomics dataframe of a DARPA study.
The athlete_norm_git.R file is an attempt at reverse engineering a batch normalization process conducted by Metabolon and replicating that on other raw datasets collected by the same company.
The corr_matrices_comp_git.R file is for comparing the empirical correlation matrices of metabolite abundances of two different datasets and highlighting correlations that are close.
