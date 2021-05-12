--------------
README
--------------



This folder contains R codes that produce the simulation and real data analysis results presented in the paper "Phylogenetically informed Bayesian truncated copula graphical models for microbial association networks" by Hee Cheol Chung, Irina Gaynanova and Yang Ni.

1. functions/AlbertChib_exp.R: Albert-Chip data augmentation for Dist model

2. functions/AlbertChib_mrg.R: Alber-Chip data augmentation for PhyloBCG

3. functions/conditionalSample.R: Sampling from the truncated normal conditional distribution

4. functions/distGibbs.R: for the gibbs sampling algorithm for Dist

5. functions/genTree.R: wrapper function random tree generation

6. functions/library.R: call libraries

7. functions/orclGibbs.R: for the gibbs sampling algorithm for Oracle

8. functions/SSVS.R: for the gibbs step for SSSL (Wang, 2015)

9. functions/TreeDataGeneration.R: random tree generation using the ape package

10. functions/treeGibbs.R: for the gibbs sampling algorithm for PhyloBCG

10. fileGen_sim_ada: Shell script, which generates "Main codes" for replicated data sets.

The following directories need to exist and accessible by "Main codes":



To run :
1. Generate data by running 
2. Generate R codes and shell scripts using codeGen.sh
3. Submit **_sub.sh files to accessible LSF cluster.


