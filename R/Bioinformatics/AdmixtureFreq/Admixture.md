Ancestry based on genotype likelihoods
================

We start by loading the data and identifying structure of the files:

``` r
### Load data
pop <- read.table('pop.info', header=F)
colnames(pop) <- c('population', 'sample')
qopt <- read.table('assign3.qopt', header=F)

fopt <- read.table('assign3.fopt.gz', header=F)

input <- read.table('input.gz', header=T)

prob_geno <- input[,-c(1:3)]
#################### Identify columns ################################

# Identify populations
print(unique(pop$population))
```

    ## [1] "ASW" "MXL" "CHB" "CEU" "YRI"

``` r
# Identify columns
by(qopt, pop$population, summary)
```

    ## pop$population: ASW
    ##        V1               V2                V3          
    ##  Min.   :0.6104   Min.   :0.05399   Min.   :0.000000  
    ##  1st Qu.:0.7196   1st Qu.:0.12242   1st Qu.:0.000000  
    ##  Median :0.8094   Median :0.18091   Median :0.002583  
    ##  Mean   :0.7919   Mean   :0.19660   Mean   :0.011525  
    ##  3rd Qu.:0.8776   3rd Qu.:0.25185   3rd Qu.:0.005046  
    ##  Max.   :0.9460   Max.   :0.38508   Max.   :0.158170  
    ## ------------------------------------------------------------ 
    ## pop$population: CEU
    ##        V1              V2               V3           
    ##  Min.   :1e-09   Min.   :0.9760   Min.   :1.000e-09  
    ##  1st Qu.:1e-09   1st Qu.:1.0000   1st Qu.:1.000e-09  
    ##  Median :1e-09   Median :1.0000   Median :1.000e-09  
    ##  Mean   :1e-09   Mean   :0.9987   Mean   :1.298e-03  
    ##  3rd Qu.:1e-09   3rd Qu.:1.0000   3rd Qu.:6.000e-09  
    ##  Max.   :1e-09   Max.   :1.0000   Max.   :2.395e-02  
    ## ------------------------------------------------------------ 
    ## pop$population: CHB
    ##        V1              V2              V3   
    ##  Min.   :1e-09   Min.   :1e-09   Min.   :1  
    ##  1st Qu.:1e-09   1st Qu.:1e-09   1st Qu.:1  
    ##  Median :1e-09   Median :1e-09   Median :1  
    ##  Mean   :1e-09   Mean   :1e-09   Mean   :1  
    ##  3rd Qu.:1e-09   3rd Qu.:1e-09   3rd Qu.:1  
    ##  Max.   :1e-09   Max.   :1e-09   Max.   :1  
    ## ------------------------------------------------------------ 
    ## pop$population: MXL
    ##        V1                V2               V3        
    ##  Min.   :0.00000   Min.   :0.5524   Min.   :0.0000  
    ##  1st Qu.:0.00000   1st Qu.:0.6868   1st Qu.:0.1809  
    ##  Median :0.00000   Median :0.7613   Median :0.2387  
    ##  Mean   :0.01696   Mean   :0.7520   Mean   :0.2310  
    ##  3rd Qu.:0.03684   3rd Qu.:0.8017   3rd Qu.:0.3070  
    ##  Max.   :0.06817   Max.   :1.0000   Max.   :0.3951  
    ## ------------------------------------------------------------ 
    ## pop$population: YRI
    ##        V1          V2              V3       
    ##  Min.   :1   Min.   :1e-09   Min.   :1e-09  
    ##  1st Qu.:1   1st Qu.:1e-09   1st Qu.:1e-09  
    ##  Median :1   Median :1e-09   Median :1e-09  
    ##  Mean   :1   Mean   :1e-09   Mean   :1e-09  
    ##  3rd Qu.:1   3rd Qu.:1e-09   3rd Qu.:1e-09  
    ##  Max.   :1   Max.   :1e-09   Max.   :1e-09

We see that the first column corresponds to african, second to european
and third to chinese.

``` r
# Rename columns according to findings
colnames(qopt) <- c("African", "European", "Chinese")
```

Let’s take a look at the average chinese ancestry in all individuals
from the ASW population:

``` r
# Compute average chinese ancestry
mean(qopt[which(pop$population == "ASW"),"Chinese"])
```

    ## [1] 0.01152502

We are interested in testing the hypothesis that the chinese ancestry in
the african american individuals is due to statistical noise. For this
test we use the following hypothesis:

``` r
# Find admixture proportions for first ASW individual
qopt_null <- qopt
qopt_null[which(pop$population == "ASW"),"Chinese"] <- 0
qopt_null[which(pop$population == "ASW"),"European"] <- 
1 - qopt[which(pop$population == "ASW"),"African"]



print(c("Population admixture under alternative hypothesis for the first individual: ", 
        qopt[which(pop$population == "ASW")[1],]))
```

    ## [[1]]
    ## [1] "Population admixture under alternative hypothesis for the first individual: "
    ## 
    ## $African
    ## [1] 0.6104123
    ## 
    ## $European
    ## [1] 0.3850796
    ## 
    ## $Chinese
    ## [1] 0.004508112

``` r
print(c("Population admixture under null hypothesis for the first individual: ",
        qopt_null[which(pop$population == "ASW")[1],]))
```

    ## [[1]]
    ## [1] "Population admixture under null hypothesis for the first individual: "
    ## 
    ## $African
    ## [1] 0.6104123
    ## 
    ## $European
    ## [1] 0.3895877
    ## 
    ## $Chinese
    ## [1] 0

We set up the likelihood model to test our hypothesis:

``` r
# Function to compute h
compute_PGh <- function(fopt, qopt, index){
hj <- as.matrix(fopt) %*% t(qopt[index,])
PGh <- cbind((1-hj)**2, 2*hj*(1-hj), (hj)**2)
return(PGh)
}

# Function to compute log likelihood
loglike <- function(fopt, qopt, prob_geno, ind, pop_info){
index <- which(pop_info$sample == ind)
PGh <- compute_PGh(fopt, qopt, index)
prob_geno_ind <- prob_geno[,(index*3-2):(index*3)]
loglike <- sum(log(apply(prob_geno_ind*PGh, 1, sum)))
return(loglike)
}
```

Now we can apply the likelihood model for the first individual and
perform the likelihood ratio test to test our hypothesis:

``` r
# Compute logLikelihood for first ASW individual
ind <- "NA19818"
# For each model
logLikeAlt <- loglike(fopt, qopt, prob_geno, ind, pop)
logLikeNull <- loglike(fopt, qopt_null, prob_geno, ind, pop)
# Compute likelihood ratio test
LR <- -2 *(logLikeNull - logLikeAlt)
pvalue <- 1 - pchisq(LR,df=1)
# Print results
print(rbind(ind, LR, pvalue))
```

    ##        [,1]               
    ## ind    "NA19818"          
    ## LR     "0.5535428209987"  
    ## pvalue "0.456873701885743"

Now let us do this for all african american individuals to find the
number of individuals with significant amount of chinese ancestry:

``` r
# Now compute log likelihood for all ASW individuals
asw_df <- data.frame(row.names = c('LR', 'pvalue'))
i <- 1
for (index in which(pop$population == 'ASW')){
ind <- pop$sample[index]
logLikeAlt <- loglike(fopt, qopt, prob_geno, ind, pop)
logLikeNull <- loglike(fopt, qopt_null, prob_geno, ind, pop)

LR <- -2 *(logLikeNull - logLikeAlt)
pvalue <- 1 - pchisq(LR,df=1)

asw_df <- cbind(asw_df, rbind(LR, pvalue))
colnames(asw_df)[i] <- ind
i <- i + 1 
}
print(asw_df)
```

    ##          NA19818   NA19901   NA19700    NA19819  NA19625   NA19712
    ## LR     0.5535428 0.1166573 0.6317932 3.99317008 837.1511 0.3677067
    ## pvalue 0.4568737 0.7326886 0.4266983 0.04568504   0.0000 0.5442570
    ##              NA19916       NA19917     NA19701     NA19904   NA19835
    ## LR     -1.269713e-06 -7.197086e-07 7.778849583 7.736630839 0.4923060
    ## pvalue  1.000000e+00  1.000000e+00 0.005286145 0.005411167 0.4829004
    ##             NA19914       NA19711   NA19834   NA19900       NA19909
    ## LR     -0.000721465 -6.427828e-06 0.6018329 1.3357088 -8.061324e-07
    ## pvalue  1.000000000  1.000000e+00 0.4378795 0.2477922  1.000000e+00
    ##             NA19704       NA19703       NA19713       NA19908
    ## LR     -0.006520275 -4.565227e-07 -5.962793e-07 -6.422779e-07
    ## pvalue  1.000000000  1.000000e+00  1.000000e+00  1.000000e+00

``` r
sum(asw_df["pvalue",]<0.05)
```

    ## [1] 4
