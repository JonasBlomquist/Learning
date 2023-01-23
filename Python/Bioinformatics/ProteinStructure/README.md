# Protein structure - Contact Number


The measure of solvent exposure is an important feature for protein folding. It provides a measure for how exposed vs buried the amino acid is and thus information about the placement in the folded protein. Solvent exposure can be measured in various ways such as accessible surface area, which describes the area of the protein, which is accessible for (usually) a water molecule. Other measures as contact number can also be used. This feature simply counts the number of $C^\alpha$ atoms within a radius (usually 10 Å) of the $C^\alpha$ atom of the amino acid. This is much easier to implement and faster to compute. 

Rosetta is a protein folding prediction program introduced in Simons et al., 1997 [^fn1], that tries to assemble proteins using fragments from similar local structures with no biological relevance to the target protein. The probability of observing the structure given the sequence is then used as scoring for the structure. This probability scoring is derived using Bayesian statistics and include various physical and geometric properties. The scoring function derived to score the fragments is: 

$$
 P(aa_1, aa_2, aa_3, \dots \vert structure) \simeq \prod_i P(aa_i \vert E_i) \times \prod_{i < j} \frac{P(aa_i, aa_j \vert r_{ij}, E_i, E_j)}{P(aa_i \vert r_{ij}, E_i, E_j)P(aa_j \vert r_{ij}, E_i, E_j)}
$$

Here $E_i$ can be different measures that describe the local structure environment. This is based on the number of neighbors to the amino acids or in other words the contact number. In the first term of the scoring function, the contact number is used directly to describe the environment, whereas in the second term it is discretized to a boolean value. Thus in the second term Rosetta only considers whether the amino acid has more than 16 $C^\beta$ atoms within a radius of 10 Å of its own $C^\beta$ atom. In the case, where it has, it will be described as buried and otherwise it is exposed Simons et al., 1997 [^fn1]. 

This scoring function thus successfully combine the residue-residue interactions with a measure of environment. This combination of features in the scoring function results in a prediction folding method, which was better than any other methods at the time. Rosetta were for many years the best folding prediction method until AlphaFold. The method presented in Simons et al., 1997 does however have problems with proteins containing β-strands. The paper mentions implementing more features to measure environment and moving away from the use of fixed number of neighbors to achieve the boolean value of buried, as this is not very efifcient for small proteins. 

The use of half sphere exposure could potentially help to solve some of these problems. Half sphere exposure (HSE) introduced in Hammelryck, 2005 [^fn2] is a simple but very effective extension of contact number. Here the sphere is divided into two half spheres dependent on the side chain. The contact number can then both be described in the direction of the side chain and in the opposite direction. This measure can thus better be used to describe side chain burial, thus helping with $\beta$-sheets and with small proteins as well. 

### References

[^fn1]: "Simmons, K., Kooperberg, C., Huang, E., & Baker, D. (1997). Assembly of protein tertiary structures from fragments with similar local sequences using simulated annealing and bayesian scoring functions. *Journal of molecular biology*, 268(1), 209–225. https://doi.org/10.1006/jmbi.1997.0959"
[^fn2]: "Hammelryck, T. (2005). An amino acid has two sides: A new 2d measure provides a different view of solvent exposure. PROTEINS: Structure, Function, and Bioinformatics, 59(1), 38–48. https://doi.org/10.1002/ prot.20379"



