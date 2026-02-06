# VietorisRipsHypercube -- homology of VR complexes of hypercubes as representation

## Description

The goal is to compute the homology groups of the Vietoris-Rips complex
of the hypercube $Q_n$ at a given scale $r$ with respect to the
Hamming distance, and to describe them as representations of the
hyperoctahedral group $\mathfrak{H}\_n$, i.e., the automorphism group
of $Q_n$. We denote this simplicial complex as $X^{n,r}$. A simplex
$\Delta$ in $X^{n,r}$ has cubic dimension $p$ if it can be
isometrically embedded in $Q_p$ but not in $Q\_{p-1}$. The simplices
with cubic dimension at most $p$ form a subcomplex $X^{n,r}\_p$ of
$X^{n,r}$. In fact, these subcomplexes form a filtration
$\\X^{n,r}\_p\\\_{p\geqslant 0}$ of $X^{n,r}$ by simplicial
complexes, and the homology of $X^{n,r}$ can be computed via the
spectral sequence associated with this filtration. Let
$C^{n,r}\_{p+q,p}$ be the vector space spanned by all the simplices of
dimension $p+q$ and cubic dimension $p$. The terms in the zero page
of the spectral sequence are $E^0\_{p,q} = C^{n,r}\_{p+q,p}$. It is
enough to study the subcomplexes corresponding to the smaller simplicial
complexes $\\0\\^{n-p} \times X^{p,r}$. Then, the symmetry afforded by
the action of $\mathfrak{H}\_n$ allows us to recover the result for
the bigger complex. More formally, for every integer $i$, we have an
isomorphism $$C^{n,r}\_{i,p} \cong
\operatorname{Ind}^{\mathfrak{H}\_n}\_{\mathfrak{S}\_{n-p} \times
\mathfrak{H}\_p} \left( \\n-p\\ \otimes C^{p,r}\_{i,p} \right)$$ of
representations of $\mathfrak{H}\_n$, where $\\n-p\$ denotes the
1-dimensional trivial representation of $\mathfrak{S}\_{n-p}$. Since
induction between finite groups is an exact functor, we then get
$$E^1\_{p,q} \cong
\operatorname{Ind}^{\mathfrak{H}\_n}\_{\mathfrak{S}\_{n-p} \times
\mathfrak{H}\_p} \left( \\n-p\\ \otimes H\_{p+q}
(C^{p,r}\_{p+\bullet,p}) \right).$$ This reduces the problem to
computing the homology of the complex $C^{p,r}\_{p+\bullet,p}$, which
we refer to as a "small complex" for the purpose of this package. Note
that the small complexes do not depend on $n$ so they are finite,
which makes it possible to build them and compute their homology with
the use of software.
