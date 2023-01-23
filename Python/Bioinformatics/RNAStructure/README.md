# RNA structure - Nussinov

RNA-RNA interactions can have a huge impact on the function of non-coding RNA. Therefore prediction of RNA-RNA interactions can be a powerful tool to analyse function of these non-coding sequences. 

By predicting the structure of the two sequences we can find free unpaired nucleotides in loops, that are able to interact with other RNA molecules. There are different ways to predict RNA structure. One of the simplest ways is to base the structure on the largest possible number of base pairs. The Nussinov algorithm will find the maximum number of base pairs an RNA molecule can make. The algorithm is a dynamic programming algorithm that will iteratively fill a scoring matrix based on previous entries. Another approach can be to predict the structure based on Gibbs free energy. Each base pair will contribute differently to the energy of the final model and therefore they should not all have the same score. Furthermore the positive energy of loops and bulges is also dependent on the size of these. Modelling the structure based on the minimum free energy is called energy folding. 

This analysis will explore the method of predicting RNA-RNA interactions using the Nussinov algorithm [^fn1]. 

### Analysis

In the notebook [RNA structure](Nussinov.ipynb) I have conducted a small analysis looking at structure of two RNA sequences and finding possible interactions between these two. 

### References

[^fn1]: "Nussinov, Ruth and Pieczenik, George and Griggs, Jerrold R. and Kleitman, Daniel J. (1978). Algorithms for Loop Matchings. *SIAM Journal on Applied Mathematics*, 35(1), 68-82. https://doi.org/10.1137/0135006"
