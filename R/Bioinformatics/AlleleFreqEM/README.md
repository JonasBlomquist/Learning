# Estimate allele frequency for the 7 pompe stains 

We want to estimate the allele frequecies for the four alleles (A, C, G, T) based on 7 different yeast strains. For each position in the genome we compute the allele frequencies of the four alleles. We need to take quality score for the bases into consideration. We will use a maximum likelihood approach and the EM-algorithm for estimation of the maximum likelihood.

**Write about EM and the likelihood function**


### Code

The algorithm is implemented in R. The implementation can be seen in the [markdown file](AlleleFreq_EM_algo.md).
