# Shannon entropy-Biphasic system
Python script to calculate the Shannon entropy of a biphasic system using GROMACS structure files extracted over the trajectory.

**Shannon Entropy:** The entropy of mixing for a binary mixture is
calculated based on Shannon's entropy that 
in the mean field approximation is given by,

$$S = - \sum_x p(x) \log p(x)$$
where $p(x)$ is the probability of system being in a state $x$. As an example,
for magnetic lattice one might consider $x=\uparrow$ or $\downarrow$ as states
and find the probabilities of lattice sites (atoms) being in one of the two
states.

However, for composite systems, we use the improvement proposed by Camesasca et
al. that is to subdivide the space into $N$ regions, and then evaluate the Shannon entropy of each region $i$ as,

$$S_i = - \sum_x p_i(x) \log p_i(x)$$
where $p_i(x)$ is the number of molecules of type $x$ in region $i$ divided by
the total number of molecules in that region. The entropy of mixing of the
whole system is then estimated as the average of the entropies over all sub-regions.

$$S = -\frac{1}{N} \sum_i^N \sum_x p_i(x) \log p_i(x)$$

Ref: Saraf et al., Journal of Molecular Liquids, Volume 381, 121803 (2023) \
     https://doi.org/10.1016/j.molliq.2023.121803
