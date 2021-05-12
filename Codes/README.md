# PhyloBCG

This folder contains R codes that produce the simulation and real data analysis results presented in the paper *"Phylogenetically informed Bayesian truncated copula graphical models for microbial association networks"* by Hee Cheol Chung, Irina Gaynanova and Yang Ni. [arXiv link](https://arxiv.org/pdf/2105.05082.pdf)

## Main functions

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


## Simulation codes

**simulation/m1.R** - Base simulation code for Oracle

**simulation/m2.R** - Base simulation code for Dist

**simulation/m3.R** - Base simulation code for PhyloBCG

**simulation/m4.R** - Base simulation code for Spiec-Easi

**simulation/m5.R** - Base simulation code for SPRING

**simulation/fileGen_sim_ada.sh** - Shell script, which generates simulation codes for each setting using a base simulation code

**simulation/subjob*.sh** - Base shell script that submits m*.R to the LSF cluster


## Supporting data

**simulation/QMPtree.RData**

**QMPdata/QMPtree.RData**

**QMPdata/tree_newick2.txt**




