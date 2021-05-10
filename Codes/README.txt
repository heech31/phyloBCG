--------------
README
--------------



This folder contains R codes that produce the simulation and real data analysis results presented in the paper "Phylogenetically informed Bayesian truncated copula graphical models for microbial association networks" by Hee Cheol Chung, Irina Gaynanova and Yang Ni.

1. functions/AlbertChib_mrg.R: 

2. functions/conditionalSample.R: 

3. functions/distGibbs.R: 

4. functions/orclGibbs.R: 

5. functions/treeGibbs.R: 

6. functions/SSVS.R: 

7. functions/genTree.R:

8. functions/TreeDataGeneration.R:

9. codeGen: Shell script, which generates "Main codes," corresponding to replicated data sets.

The following directories need to exist and accessible by "Main codes":



To run :
1. Generate data by running 
2. Generate R codes and shell scripts using codeGen.sh
3. Submit **_sub.sh files to accessible LSF cluster.


