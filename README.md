# Description of Study
The overarching goal of this study is in the realm of multiomics data and how they relate to certain phenotypic measures. All of the underlying data and code has been sanitized. Omics data were collected before and after intense exercise: blood draw 1 (pre-exercise) vs. blood draws 2, 3, and 4 (post-exercise). The experiment and data collection were repeated twice more after the very first, resulting in a total of three sessions of four blood draws.

# Immunophenotypic Data Exploration
This is an exploratory analysis of immunophenotypic protein data across multiple subjects. The results show no marked changes in the immunophenotypic profiles post-exercise (vs. pre-exercise).

# Metabolomics Data Merger
The core idea of these files is to merge open-source metabolomics dataset into a metabolomics dataframe of a DARPA study.
The athlete_norm_git.R file is an attempt at reverse engineering a batch normalization process conducted by Metabolon and replicating that on other raw datasets collected by the same company.
The corr_matrices_comp_git.R file is for comparing the empirical correlation matrices of metabolite abundances of two different datasets and highlighting correlations that are close.

# Metabolomics Analysis
The first step of the metabolomics analysis invovles temporal clustering of metabolite abundance levels across four blood draws within each session. The follow-up is an across-session clustering, which invovles fixing a certain blood draw and looking at it across sessions. At the core of these processes is the idea of vectorizing the slopes between two consecutive time points and then clustering these vectors with [fuzzy c-means clustering](https://www.bioconductor.org/packages/release/bioc/vignettes/Mfuzz/inst/doc/Mfuzz.pdf). The approach was inspired by the [tscR](http://www.bioconductor.org/packages/release/bioc/manuals/tscR/man/tscR.pdf) package in R. 

The second step invovles analyzing the correlations between certain metabolites and metabolic ratios pre- and post-exercise to see if there are any patterns that coincide with those recorded in biological literature. Also, several regressoion models are fit between VO2Max values and metabolomic abundance levels to see if any molecules may be related to VO2Max variations across individuals.

The last step invovles applying [sparse canonical correlation analysis](https://cran.r-project.org/web/packages/PMA/PMA.pdf) to the dataset to achieve dimensionality reduction, to see if there is any meaningful way to relate a linear combination of the molecular data to a linear combination of the target phenotypes. Also, a [LASSO stability selection](https://cran.r-project.org/web/packages/stabs/stabs.pdf) method is used to select features that may be fed into a machine learning pipeline for predictions of phenotypic targets.

# [Weighted Gene Co-expression Network Analysis](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/) (WGCNA)
This is a clustering technique that treats any given set of molecules as a fully connected network with each edge corresponding to a relationship score between two a pair of molecules. The score is calculated based on the absolute correlation of expression/abundance values. Likewise, each node corresponds to each molecule. A dissimilarity metric between every pair of nodes is calculated based on the edge scores and the surrounding first order neighbors. With this metric, hierarchical clustering is performed to bring together molecules that are highly correlated with each other. Finally, representatives of these clusters (eigenvectors) are correlated to target phenotypes to see whether each cluster has a meaningful relationship with the target phenotypes.
