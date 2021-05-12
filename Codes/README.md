# SING

This folder contains R codes that produce the simulation and real data analysis results presented in the paper *"Phylogenetically informed Bayesian truncated copula graphical models for microbial association networks"* by Hee Cheol Chung, Irina Gaynanova and Yang Ni.[arXiv link](https://github.com/irinagain/SING)

## Main codes

**functions/AlbertChib_exp.R** - Albert-Chip data augmentation for Dist model

**functions/AlbertChib_mrg.R** - Alber-Chip data augmentation for PhyloBCG

**functions/conditionalSample.R** - Sampling from the truncated normal conditional distribution

**functions/distGibbs.R** - for the gibbs sampling algorithm for Dist

**functions/genTree.R** - wrapper function random tree generation

**functions/library.R** - call libraries

**functions/orclGibbs.R** - for the gibbs sampling algorithm for Oracle

**functions/SSVS.R** - for the gibbs step for SSSL (Wang, 2015)

**functions/TreeDataGeneration.R** - random tree generation using the ape package

**functions/treeGibbs.R** - for the gibbs sampling algorithm for PhyloBCG

**fileGen_sim_ada.sh** - Shell script, which generates "Main codes" for replicated data sets.


## Supporting data

**simulation/QMPtree.RData**

**QMPdata/QMPtree.RData**

**QMPdata/tree_newick2.txt**




