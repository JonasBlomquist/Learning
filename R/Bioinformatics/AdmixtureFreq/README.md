# Ancestry based on genotype likelihoods

We will reexamine the results from the exercises based on the 1000 Genomes human low depth sequencing data. We will try to answer the question of whether estimated small admixture proportions is evidence for admixture or if it is only statistic noise. 

We use NGSadmix which will estimate both admixture proportions and allele frequencies for each ancestral population.

Now we will focus only on the African Americans. Some of the African Americans show signs of having Chinese ancestry. However, many of them have only very low Chinese admixture proportions. We therefore want to test whether each individual have a significant amount of Chinese ancestry or if the low Chinese admixture proportions can be explained by statistical noise. To do this use the estimated ancestral allele frequencies and the estimated admixture proportions from the above NGSadmix run. In order to make the test we need to have two models. One where we don’t allow for Chinese ancestry (null model) and one where we allow for Chinese ancestry (alternative model). The admixture proportions for the alternative model is already estimate from the NGSadmix run. For the null model we will make the assumption that the estimated Chinese admixture proportions was actually European. For example if the estimated admixture proportions under the alternative model was European=15%, Chinese=5% and African=80% then the null model is European=20%, Chinese=0% and African=80%. This assumption can be somewhat reasonable because Chinese are more similar to Europeans than they are to Africans and thus NGSadmix is more likely to infer European ancestry as Chinese. 

The likelihood model of the data given the ancestral allele frequencies and the admixture proportions can be written as:

    \begin{equation}
        P(X | Q, F) = \prod_j^M P(X_j | Q, F_j) = \prod_j^M \sum_{G \in \{0,1,2\}} P(X_j | Q, F_j, G_j) P(G_j | Q, F_j)
    \end{equation}
    
    \begin{equation}
        \prod_j^M \sum_{G \in \{0,1,2\}} P(X_j | Q, F_j, G_j) P(G_j | Q, F_j) = \prod_j^M \sum_{G \in \{0,1,2\}} P(X_j | G_j) P(G_j | Q, F_j)
    \end{equation}
    with the notation:
    \begin{itemize}
        \item $X_j$: the sequencing data at site j
        
        \item $Q$: admixture frequencies
        
        \item $F_j$: ancestral allele frequencies for site j
        
        \item $M$: number of sites
        
        \item $G_j$: genotype for site j
    \end{itemize}
    
    $P(X_j|G_j)$ is the genotype likelihoods that we have in the input.gz file. 
    
    The probability of the genotype given the ancestral allele frequencies and the admixture proportions assuming Hardy-Weinberg equilibrium can be written as:
    
\begin{equation}
        P(G_j | Q, F_j) =  P(G_j | h^j) =
    \begin{cases} \mbox{$(1-h^j)^2$} & \text{if}\quad G_j = 0 \\ 
    \mbox{$2 h (1-h)$} & \text{if}\quad G_j = 1 \\
    \mbox{$(h^j)^2$} & \text{if}\quad G_j = 2\end{cases}
\end{equation}
With $h^j$ defined as:
\begin{equation}
    h^j = \sum_{k=1}^K f^{jk} q^k
\end{equation}
with notation:
\begin{itemize}
    \item $K$: number of ancestral populations

    \item $q^k$: admixture proportion for population $k$

    \item $f^{jk}$: allele frequency of allele A for population $k$ at position $j$
\end{itemize}


\item For the first african american individual the we have the following admixture proportions under null and alternative model:
\begin{table}[h!]
    \centering
    \begin{tabular}{|c|c|c|}
    \hline
         & Null model & Alternative model \\\hline
        African & 0.6104123 & 0.6104123 \\
        European & 0.3895877 & 0.3850796 \\
        Chinese & 0 & 0.004508112 \\\hline
    \end{tabular}
    \caption{Null model and alternative model for the first African American individual.}
    \label{tab:models}
\end{table}

The likelihood funtion is implemented using R and this is used to compute the likelihood of the two proposed models and they can then be tested using a likelihood ratio test. 

### Code

The analysis is performed in R and can be found in the [markdown file](Admixture.md).
