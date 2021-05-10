--------------
README
--------------



This folder contains R codes which reproduce the simulation and real data analysis results in the paper "Phylogenetically informed Bayesian truncated copula graphical models for microbial association networks" by Hee Cheol Chung, Irina Gaynanova and Yang Ni.

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

For the real data, one need to have folders named as ./RT_LDA_realData/"dataNames". For example, ./RT_LDA_realData/diab. Data names are listed in RealData/dataNames.csv.
The results (test error rates and basis of discriminant subspaces) will be saved to "resultpath". The "resultpath" is specified in the "Main codes." Please change the path accordingly.

For the simulated data sets, one need to have folders named as ./RT_LDA_simData/"n**s**". For example, ./RT_LDA_simData/n90s50, where n is the sample size and s is the number of variables that are relevant to classification. See Section 3.1 for more details.


To run :
1. Generate data using GenCVindex.m and GenData.m.
2. Generate Matlab codes and shell scripts using codeGen.sh
3. Submit **_sub.sh files to accessible HPC cluster.


