# Description of Study
The overarching goal of this study is in the realm of multiomics data and how they relate to certain phenotypic measures. All of the underlying data and code has been sanitized. Omics data were collected before and after intense exercise: blood draw 1 (pre-exercise) vs. blood draws 2, 3, and 4 (post-exercise). The experiment and data collection were repeated twice more after the very first, resulting in a total of three sessions of four blood draws.

# [Weighted Gene Co-expression Network Analysis](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/) (WGCNA)
This is a clustering technique that treats any given set of molecules as a fully connected network with each edge corresponding to a relationship score between two a pair of molecules. The score is calculated based on the absolute correlation of expression/abundance values. Likewise, each node corresponds to each molecule. A dissimilarity metric between every pair of nodes is calculated based on the edge scores and the surrounding first order neighbors. With this metric, hierarchical clustering is performed to bring together molecules that are highly correlated with each other. Finally, representatives of these clusters (eigenvectors) are correlated to target phenotypes to see whether each cluster has a meaningful relationship with the target phenotypes.

# Immunophenotypic Data Exploration
This is an exploratory analysis of immunophenotypic protein data across multiple subjects. The results show no marked changes in the immunophenotypic profiles post-exercise (vs. pre-exercise).

# Metabolomics Data Merger
The core idea of these files is to merge open-source metabolomics dataset into a metabolomics dataframe of a DARPA study.
The athlete_norm_git.R file is an attempt at reverse engineering a batch normalization process conducted by Metabolon and replicating that on other raw datasets collected by the same company.
The corr_matrices_comp_git.R file is for comparing the empirical correlation matrices of metabolite abundances of two different datasets and highlighting correlations that are close.
