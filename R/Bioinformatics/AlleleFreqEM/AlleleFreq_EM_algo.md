Allele Frequency using EM Algorithm
================

Load the data file:

``` r
# Open mpileup
mpileup <- read.delim('nice.mpileup', 
                      as.is=T, sep='\t', header=F, comment.char="")

n <- (length(mpileup) - 3)/3 # Find number of individuals
# Name columns in mpileup
colnames(mpileup) <- c("chr", "pos", "ref", paste(c("n", "base", "qual"), 
                                                  rep(1:n, each=3), sep = '_'))
```

Convert ascii Q-score encoding to probability:

``` r
# Define function to convert quality score to probability of error
QtoProb<-function(x){
  y <- utf8ToInt(x) ## ascii to integer
  Q <- y -33 #offset
  P <- 10**(-Q/10) # Convert to prob
  P
}
```

Now we can define our likelihood model and a function for our EM
algorithm:

``` r
# Define likelihood model
likelihood <- function(base, observed, quality){
  corr <- prod(1 - QtoProb(quality)[which(strsplit(observed, "")[[1]]==base)]) #If position equal to base
  other <- prod(1/3*QtoProb(quality)[which(strsplit(observed, "")[[1]]!=base)]) # If position not equal to base
  corr*other
}

# Define EM algorithm
EMalgo <- function(pk, like){
  ##P(Xi|G)theta_j
  Ab <- like[1,]*pk[1] 
  Cb <- like[2,]*pk[2]
  Gb <- like[3,]*pk[3]
  Tb <- like[4,]*pk[4]
  
  pX <- Ab + Cb + Gb + Tb
  
  #Q_i(G_j)
  pA <- Ab/pX
  pC <- Cb/pX
  pG <- Gb/pX
  pT <- Tb/pX
  
  all <- sum(pA + pC + pG + pT)
  
  # Compute theta new
  pk_new <- c(sum(pA)/all, sum(pC)/all, sum(pG)/all, sum(pT)/all)
  return(pk_new)
}
```

Finally we can run our EM algorithm to find likelihoods for each
possible base. We identify sites for which the most common base has a
probability $< 0.95$.

``` r
# Initiate count
n_sites <- 0

# loop over sites and then yeast strains
for (j in 1:nrow(mpileup)){
  like = c()
  for (i in 1:n){
    # Initial guess
    pk = c(.25,.25,.25,.25)
    # Extract the observed bases and quality scores
    observed = mpileup[j,paste("base", i, sep="_")]
    quality = mpileup[j,paste("qual", i, sep="_")]
    # Compute likelihoods for each possible base
    like_A = likelihood("A", observed, quality)
    like_C = likelihood("C", observed, quality)
    like_G = likelihood("G", observed, quality)
    like_T = likelihood("T", observed, quality)
    like <- cbind(like, c(like_A, like_C, like_G, like_T))
  }
  # Perform 10 iterations of the EM algorithm
  for (k in 1:10){
    pk <- EMalgo(pk, like)
  }
  # If the most common allele have a frequence < 0.95 print the site
  if (max(pk) < .95){
    cat("\n\nShowing probabilities for position: ", mpileup[j, "pos"], "in row", 
        j, "\n")
    print(data.frame("A"=pk[1], "C"=pk[2], "G"=pk[3], "T"=pk[4]))
    n_sites <- n_sites + 1 # Increase count by 1
  }
}
```

    ## 
    ## 
    ## Showing probabilities for position:  6378 in row 438 
    ##          A            C            G        T
    ## 1 0.133425 3.134878e-21 3.134878e-21 0.866575
    ## 
    ## 
    ## Showing probabilities for position:  6380 in row 439 
    ##             A         C         G           T
    ## 1 2.64081e-31 0.1419504 0.8580496 2.64081e-31
    ## 
    ## 
    ## Showing probabilities for position:  6381 in row 440 
    ##              A         C         G            T
    ## 1 2.963244e-36 0.8574283 0.1425717 2.963244e-36
    ## 
    ## 
    ## Showing probabilities for position:  6583 in row 456 
    ##           A             C         G             T
    ## 1 0.1974129 8.073915e-129 0.8025871 8.073915e-129
    ## 
    ## 
    ## Showing probabilities for position:  6727 in row 600 
    ##           A            C         G            T
    ## 1 0.4587739 1.360575e-49 0.5412261 1.360575e-49
    ## 
    ## 
    ## Showing probabilities for position:  6831 in row 605 
    ##          A            C            G        T
    ## 1 0.857251 7.140636e-39 7.140636e-39 0.142749
    ## 
    ## 
    ## Showing probabilities for position:  14931 in row 3627 
    ##           A           C           G         T
    ## 1 0.8571788 2.42199e-45 2.42199e-45 0.1428212
    ## 
    ## 
    ## Showing probabilities for position:  14932 in row 3628 
    ##              A         C         G            T
    ## 1 2.423006e-44 0.1428119 0.8571881 2.423006e-44
    ## 
    ## 
    ## Showing probabilities for position:  14933 in row 3629 
    ##              A            C         G         T
    ## 1 2.422688e-44 2.422688e-44 0.8571881 0.1428119
    ## 
    ## 
    ## Showing probabilities for position:  15780 in row 4308 
    ##              A         C            G         T
    ## 1 4.388745e-41 0.1427674 4.388745e-41 0.8572326
    ## 
    ## 
    ## Showing probabilities for position:  15781 in row 4309 
    ##           A           C         G           T
    ## 1 0.8572325 5.06747e-41 0.1427675 5.06747e-41
    ## 
    ## 
    ## Showing probabilities for position:  15782 in row 4310 
    ##           A            C            G         T
    ## 1 0.1427861 6.043355e-42 6.043355e-42 0.8572139
    ## 
    ## 
    ## Showing probabilities for position:  15783 in row 4311 
    ##           A            C            G         T
    ## 1 0.8572322 6.044065e-41 6.044065e-41 0.1427678
    ## 
    ## 
    ## Showing probabilities for position:  15785 in row 4313 
    ##              A            C        G        T
    ## 1 5.065042e-42 5.065042e-42 0.857214 0.142786

``` r
# Print number of sites with most common < 0.95
print(n_sites)
```

    ## [1] 14
