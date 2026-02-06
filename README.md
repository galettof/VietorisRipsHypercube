<!DOCTYPE html>
<html lang="en">
  <body>
    <div>
      <h1>VietorisRipsHypercube -- homology of VR complexes of hypercubes as representation</h1>
      <div>
        <h2>Description</h2>
        <div>
          <p>This package implements computations described in the article &quot;TITLE&quot; by Federico Galetto, Jonathan Montaño, and Zoe Wellner. The most up-to-date version of this package can be found at <a href="https://github.com/galettof/VietorisRipsHypercube">https://github.com/galettof/VietorisRipsHypercube</a>.</p>
          <h2>Background and Scope</h2>
          <p>The goal is to compute the homology groups of the Vietoris-Rips complex of the hypercube $Q_n$ at a given scale $r$ with respect to the Hamming distance, and to describe them as representations of the hyperoctahedral group $\mathfrak{H}_n$, i.e., the automorphism group of $Q_n$. We denote this simplicial complex as $X^{n,r}$. A simplex $\Delta$ in $X^{n,r}$ has cubic dimension $p$ if it can be isometrically embedded in $Q_p$ but not in $Q_{p-1}$. The simplices with cubic dimension at most $p$ form a subcomplex $X^{n,r}_p$ of $X^{n,r}$. In fact, these subcomplexes form a filtration $\{X^{n,r}_p\}_{p\geqslant 0}$ of $X^{n,r}$ by simplicial complexes, and the homology of $X^{n,r}$ can be computed via the spectral sequence associated with this filtration. Let $C^{n,r}_{p+q,p}$ be the vector space spanned by all the simplices of dimension $p+q$ and cubic dimension $p$. The terms in the zero page of the spectral sequence are $E^0_{p,q} = C^{n,r}_{p+q,p}$. It is enough to study the subcomplexes corresponding to the smaller simplicial complexes $\{0\}^{n-p} \times X^{p,r}$. Then, the symmetry afforded by the action of $\mathfrak{H}_n$ allows us to recover the result for the bigger complex. More formally, for every integer $i$, we have an isomorphism

$$C^{n,r}_{i,p} \cong \text{Ind}^{\mathfrak{H}_n}_{\mathfrak{S}_{n-p} \times \mathfrak{H}_p} \left( \{n-p\} \otimes C^{p,r}_{i,p} \right)$$

of representations of $\mathfrak{H}_n$, where $\{n-p\}$ denotes the 1-dimensional trivial representation of $\mathfrak{S}_{n-p}$. Since induction between finite groups is an exact functor, we then get

$$E^1_{p,q} \cong \text{Ind}^{\mathfrak{H}_n}_{\mathfrak{S}_{n-p} \times \mathfrak{H}_p} \left( \{n-p\} \otimes H_{p+q} (C^{p,r}_{p+\bullet,p}) \right).$$

This reduces the problem to computing the homology of the complex $C^{p,r}_{p+\bullet,p}$, which we refer to as a &quot;small complex&quot; for the purpose of this package. Note that the small complexes do not depend on $n$ so they are finite, which makes it possible to build them and compute their homology with the use of software.</p>
          <h2>Basic Functionality</h2>
          <p>The methods in this package allow users to construct small complexes, to compute their homology groups, and to describe those homology groups as representations of a suitable hyperoctahedral group by computing their characters. We illustrate the workflow for $X^{n,2}$, the Vietoris-Rips complex at scale 2 of $Q_n$. First, we construct the simplicial complex for $n=4$.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i1 : X42 = VRHypercube(4,2)

o1 = Vietoris-Rips complex at scale 2 of Q4

o1 : VRData</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>Next, we form the small complex for cubic dimension 2 and compute its homology.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i2 : smallComplex(X42,2)

             2       4       1
o2 = 0 &lt;-- QQ  &lt;-- QQ  &lt;-- QQ  &lt;-- 0
                                    
     0     1       2       3       4

o2 : Complex</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i3 : homology(X42,2)

       1
o3 = QQ
      
     2

o3 : Complex</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>We see that the homology is a 1-dimensional vector space in homological degree 2. The following commands give a more detailed description of this homology group.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i4 : H22 = homology(X42,2,2)

o4 = subquotient (| -e_(0,0)e_(0,1)e_(1,1)+e_(0,0)e_(1,0)e_(1,1) |, | -e_(0,0)e_(0,1)e_(1,0)+e_(0,0)e_(0,1)e_(1,1)-e_(0,0)e_(1,0)e_(1,1)+e_(0,1)e_(1,0)e_(1,1) |)

                                                           1
o4 : QQ[e   ..e   ]-module, subquotient of (QQ[e   ..e   ])
         0,0   1,1                              0,0   1,1</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i5 : netList flatten entries generators H22

     +-----------------------------+
o5 = |- e   e   e    + e   e   e   |
     |   0,0 0,1 1,1    0,0 1,0 1,1|
     +-----------------------------+</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i6 : netList flatten entries relations H22

     +-----------------------------------------------------------+
o6 = |- e   e   e    + e   e   e    - e   e   e    + e   e   e   |
     |   0,0 0,1 1,0    0,0 0,1 1,1    0,0 1,0 1,1    0,1 1,0 1,1|
     +-----------------------------------------------------------+</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>We interpret these results. The simplicial complex $X^{2,2}$ has four vertices: $(0,0)$, $(0,1)$, $(1,0)$, and $(1,1)$. The vertex $u$ corresponds to a variable $e_u$ in a polynomial ring. The squarefree monomial $e_{u_1} \cdots e_{u_t}$ corresponds to the simplex $[u_1,\dots,u_t]$. Note that the monomial ordering on the polynomial ring determines the orientation of the simplex. Then, a linear combination of monomials corresponds to a formal linear combination of simplices. In our case, $-[(0,0),(0,1),(1,1)]+[(0,0),(1,0),(1,1)]$ is a linear combination of triangles that generates the homology of the small complex of cubic dimension 2. This linear combination corresponds to one of the two ways of dividing a square into two triangles by tracing a diagonal. The relation on the homology says that the other way of dividing the square is equivalent.</p>
          <p>Now, we compute the character of this homology group as a representation of $\mathfrak{H}_2$.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i7 : c22 = character(X42,2,2)

o7 = Character over QQ
      
     (2, {3})  |  -1  1  -1  1  1

o7 : Character</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>We proceed to decompose this into irreducible characters. This functionality requires at least version 2.6 of the <a title="finite group characters on free resolutions and graded modules" href="../../BettiCharacters/html/index.html">BettiCharacters</a> package, which can be downloaded at <a href="https://github.com/galettof/BettiCharacters">https://github.com/galettof/BettiCharacters</a>. First, we construct the character table of $\mathfrak{H}_2$. Then, we express the character above as a positive integer linear combination of irreducible characters.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i8 : T2 = hyperoctahedralGroupTable(2,QQ)

o8 = Character table over QQ
      
             |   2  1   2   2   1
     --------+-------------------
      (2;0)  |   1  1   1   1   1
       2     |                 
     (1 ;0)  |  -1  1   1  -1   1
      (1;1)  |   0  2   0   0  -2
      (0;2)  |   1  1  -1  -1   1
         2   |                 
     (0;1 )  |  -1  1  -1   1   1

o8 : CharacterTable</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i9 : c22 / T2

o9 = Decomposition table
      
               |      2
               |  (0;1 )
     ----------+----------
     (2, {3})  |       1

o9 : CharacterDecomposition</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>We see that our homology group is in fact isomorphic to the irreducible representation indexed by the pair of partitions $(\varnothing; 1^2)$.</p>
          <p>When we perform the same steps for cubic dimension 3, we obtain the following.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i10 : smallComplex(X42,3)

                    8       10
o10 = 0 &lt;-- 0 &lt;-- QQ  &lt;-- QQ   &lt;-- 0
                                    
      0     1     2       3        4

o10 : Complex</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i11 : homology(X42,3)

        2
o11 = QQ
       
      3

o11 : Complex</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i12 : H33 = homology(X42,3,3);</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i13 : netList flatten entries generators H33

      +--------------------------------------------------------------------------------------------------------------------------------------+
o13 = |- e     e     e     e      + e     e     e     e      - e     e     e     e      - e     e     e     e      + e     e     e     e     |
      |   0,0,0 0,0,1 0,1,0 1,0,0    0,0,1 0,1,0 0,1,1 1,1,1    0,0,1 0,1,0 1,0,0 1,1,1    0,0,1 1,0,0 1,0,1 1,1,1    0,1,0 1,0,0 1,1,0 1,1,1|
      +--------------------------------------------------------------------------------------------------------------------------------------+
      |- e     e     e     e      + e     e     e     e      + e     e     e     e      - e     e     e     e      + e     e     e     e     |
      |   0,0,0 0,0,1 0,1,1 1,0,1    0,0,0 0,1,0 0,1,1 1,1,0    0,0,0 0,1,1 1,0,1 1,1,0    0,0,0 1,0,0 1,0,1 1,1,0    0,1,1 1,0,1 1,1,0 1,1,1|
      +--------------------------------------------------------------------------------------------------------------------------------------+</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i14 : netList flatten entries relations H33

o14 = ++
      ++</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i15 : c33 = character(X42,3,3)

o15 = Character over QQ
       
      (3, {4})  |  2  -2  2  0  0  0  2  0  -2  0

o15 : Character</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i16 : T3 = hyperoctahedralGroupTable(3,QQ);</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i17 : c33 / T3

o17 = Decomposition table
       
                |    3         3
                |  (1 ;0)  (0;1 )
      ----------+----------------
      (3, {4})  |       1       1

o17 : CharacterDecomposition</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>We deduce that the homology of the small complex for cubic dimension 3 is a 2-dimensional vector space in homological degree 3. We get two explicit generators and there are no relations. As a represenatation of $\mathfrak{H}_3$ it is isomorphic to the direct sum of the irreducibles indexed by the pairs of partitions $(1^3;\varnothing)$ and $(\varnothing;1^3)$.</p>
          <p>Finally, we observe that the small complex for cubic dimension 4 has no homology.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i18 : smallComplex(X42,4)

                          16       16
o18 = 0 &lt;-- 0 &lt;-- 0 &lt;-- QQ   &lt;-- QQ
                                  
      0     1     2     3        4

o18 : Complex</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i19 : homology(X42,4)

o19 = 0
       
      0

o19 : Complex</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p></p>
          <h2>Additional Features</h2>
        </div>
        <div>
          <p>The <span class="tt">differential</span> method can be used to study certain maps in the first page of the spectral sequence. Specifically, consider the differentials $d^1_{i,0} \colon E^1_{i,0} \to E^1_{i-1,0}$ with $1\leqslant i\leqslant n$. The <span class="tt">differential</span> method takes the elements of $E^1_{i,0}$ that generate the homology of the $i$-th small complex (computed using <span class="tt">homology</span>) and returns their image under $d^1_{i,0}$. Once again, the result is expressed as a linear combination of squarefree monomials corresponding to simplices.</p>
          <p>For the example above, we obtain the following.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i20 : netList differential(X42,1)

      +-------------------+
o20 = |e        - e       |
      | 1 : (0)    1 : (1)|
      +-------------------+</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i21 : netList differential(X42,2)

      +-----------------------------------------+
o21 = |e   e    - e   e    + e   e    - e   e   |
      | 0,0 0,1    0,0 1,0    0,1 1,1    1,0 1,1|
      +-----------------------------------------+</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>We see that $d^1_{1,0}$ sends the edge $[0,1]$ to its vertices, and $d^1_{2,0}$ sends the square on the vertices $(00)$, $(01)$, $(10)$, and $(11)$ to its boundary.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i22 : L = differential(X42,3);</code></pre>
            </td>
          </tr>
          <tr>
            <td>
              <pre><code class="language-macaulay2">i23 : netList L

      +-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
o23 = |- e     e     e      + e     e     e      + e     e     e      - e     e     e      - e     e     e      + e     e     e      + e     e     e      - e     e     e      - e     e     e      + e     e     e      + e     e     e      - e     e     e     |
      |   0,0,0 0,0,1 0,1,1    0,0,0 0,1,0 0,1,1    0,0,0 0,0,1 1,0,1    0,0,0 1,0,0 1,0,1    0,0,0 0,1,0 1,1,0    0,0,0 1,0,0 1,1,0    0,0,1 0,1,1 1,1,1    0,1,0 0,1,1 1,1,1    0,0,1 1,0,1 1,1,1    1,0,0 1,0,1 1,1,1    0,1,0 1,1,0 1,1,1    1,0,0 1,1,0 1,1,1|
      +-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
      |- e     e     e      + e     e     e      + e     e     e      - e     e     e      - e     e     e      + e     e     e      + e     e     e      - e     e     e      - e     e     e      + e     e     e      + e     e     e      - e     e     e     |
      |   0,0,0 0,0,1 0,1,1    0,0,0 0,1,0 0,1,1    0,0,0 0,0,1 1,0,1    0,0,0 1,0,0 1,0,1    0,0,0 0,1,0 1,1,0    0,0,0 1,0,0 1,1,0    0,0,1 0,1,1 1,1,1    0,1,0 0,1,1 1,1,1    0,0,1 1,0,1 1,1,1    1,0,0 1,0,1 1,1,1    0,1,0 1,1,0 1,1,1    1,0,0 1,1,0 1,1,1|
      +-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>Now, the homology of the small complex of cubic dimension 3 has two generators, and $d^1_{3,0}$ sends them to the expressions above, which are the same. This means we can change basis in the homology by taking the sum and difference of the generators, so that $d^1_{3,0}$ sends one of these new basis elements to zero.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i24 : netList {L_0+L_1,L_0-L_1}

      +-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
o24 = |- 2e     e     e      + 2e     e     e      + 2e     e     e      - 2e     e     e      - 2e     e     e      + 2e     e     e      + 2e     e     e      - 2e     e     e      - 2e     e     e      + 2e     e     e      + 2e     e     e      - 2e     e     e     |
      |    0,0,0 0,0,1 0,1,1     0,0,0 0,1,0 0,1,1     0,0,0 0,0,1 1,0,1     0,0,0 1,0,0 1,0,1     0,0,0 0,1,0 1,1,0     0,0,0 1,0,0 1,1,0     0,0,1 0,1,1 1,1,1     0,1,0 0,1,1 1,1,1     0,0,1 1,0,1 1,1,1     1,0,0 1,0,1 1,1,1     0,1,0 1,1,0 1,1,1     1,0,0 1,1,0 1,1,1|
      +-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
      |0                                                                                                                                                                                                                                                                      |
      +-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p></p>
          <h2>Some Technical Details</h2>
          <p>All the information about the simplicial complex $X^{n,2}$ is stored in an object of type <span class="tt">VRData</span> which is a dedicated type of <a title="the class of all hash tables" href="/usr/share/doc/Macaulay2/Macaulay2Doc/html/___Hash__Table.html">HashTable</a>. Let us take a look at the object for the example above.</p>
        </div>
        <table class="examples">
          <tr>
            <td>
              <pre><code class="language-macaulay2">i25 : peek X42

o25 = VRData{character => MutableHashTable{...2...}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            }
             cubicDimension => MutableHashTable{...5...}
             homology => MutableHashTable{...4...}
             ideal => ideal (x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       , x       x       )
                              0,0,1,1 0,1,0,0   0,0,1,0 0,1,0,1   0,0,0,1 0,1,1,0   0,0,0,0 0,1,1,1   0,0,1,1 1,0,0,0   0,1,0,1 1,0,0,0   0,1,1,0 1,0,0,0   0,1,1,1 1,0,0,0   0,0,1,0 1,0,0,1   0,1,0,0 1,0,0,1   0,1,1,0 1,0,0,1   0,1,1,1 1,0,0,1   0,0,0,1 1,0,1,0   0,1,0,0 1,0,1,0   0,1,0,1 1,0,1,0   0,1,1,1 1,0,1,0   0,0,0,0 1,0,1,1   0,1,0,0 1,0,1,1   0,1,0,1 1,0,1,1   0,1,1,0 1,0,1,1   0,0,0,1 1,1,0,0   0,0,1,0 1,1,0,0   0,0,1,1 1,1,0,0   0,1,1,1 1,1,0,0   1,0,1,1 1,1,0,0   0,0,0,0 1,1,0,1   0,0,1,0 1,1,0,1   0,0,1,1 1,1,0,1   0,1,1,0 1,1,0,1   1,0,1,0 1,1,0,1   0,0,0,0 1,1,1,0   0,0,0,1 1,1,1,0   0,0,1,1 1,1,1,0   0,1,0,1 1,1,1,0   1,0,0,1 1,1,1,0   0,0,0,0 1,1,1,1   0,0,0,1 1,1,1,1   0,0,1,0 1,1,1,1   0,1,0,0 1,1,1,1   1,0,0,0 1,1,1,1
             presentation => MutableHashTable{...4...}
             ring => QQ[x       ..x       ]
                         0,0,0,0   1,1,1,1
             scale => 2
             simplicialComplex => simplicialComplex | x_(1,1,0,0)x_(1,1,0,1)x_(1,1,1,0)x_(1,1,1,1) x_(1,0,1,0)x_(1,0,1,1)x_(1,1,1,0)x_(1,1,1,1) x_(0,1,1,0)x_(0,1,1,1)x_(1,1,1,0)x_(1,1,1,1) x_(1,0,0,1)x_(1,0,1,1)x_(1,1,0,1)x_(1,1,1,1) x_(0,1,0,1)x_(0,1,1,1)x_(1,1,0,1)x_(1,1,1,1) x_(1,0,0,1)x_(1,0,1,0)x_(1,1,0,0)x_(1,1,1,1) x_(0,1,0,1)x_(0,1,1,0)x_(1,1,0,0)x_(1,1,1,1) x_(0,0,1,1)x_(0,1,1,1)x_(1,0,1,1)x_(1,1,1,1) x_(0,0,1,1)x_(0,1,1,0)x_(1,0,1,0)x_(1,1,1,1) x_(0,0,1,1)x_(0,1,0,1)x_(1,0,0,1)x_(1,1,1,1) x_(1,0,0,0)x_(1,0,1,1)x_(1,1,0,1)x_(1,1,1,0) x_(0,1,0,0)x_(0,1,1,1)x_(1,1,0,1)x_(1,1,1,0) x_(1,0,0,0)x_(1,0,1,0)x_(1,1,0,0)x_(1,1,1,0) x_(0,1,0,0)x_(0,1,1,0)x_(1,1,0,0)x_(1,1,1,0) x_(0,0,1,0)x_(0,1,1,1)x_(1,0,1,1)x_(1,1,1,0) x_(0,0,1,0)x_(0,1,1,0)x_(1,0,1,0)x_(1,1,1,0) x_(0,0,1,0)x_(0,1,0,0)x_(1,0,0,0)x_(1,1,1,0) x_(1,0,0,0)x_(1,0,0,1)x_(1,1,0,0)x_(1,1,0,1) x_(0,1,0,0)x_(0,1,0,1)x_(1,1,0,0)x_(1,1,0,1) x_(0,0,0,1)x_(0,1,1,1)x_(1,0,1,1)x_(1,1,0,1) x_(0,0,0,1)x_(0,1,0,1)x_(1,0,0,1)x_(1,1,0,1) x_(0,0,0,1)x_(0,1,0,0)x_(1,0,0,0)x_(1,1,0,1) x_(0,0,0,0)x_(0,1,1,0)x_(1,0,1,0)x_(1,1,0,0) x_(0,0,0,0)x_(0,1,0,1)x_(1,0,0,1)x_(1,1,0,0) x_(0,0,0,0)x_(0,1,0,0)x_(1,0,0,0)x_(1,1,0,0) x_(1,0,0,0)x_(1,0,0,1)x_(1,0,1,0)x_(1,0,1,1) x_(0,0,1,0)x_(0,0,1,1)x_(1,0,1,0)x_(1,0,1,1) x_(0,0,0,1)x_(0,0,1,1)x_(1,0,0,1)x_(1,0,1,1) x_(0,0,0,1)x_(0,0,1,0)x_(1,0,0,0)x_(1,0,1,1) x_(0,0,0,0)x_(0,0,1,1)x_(1,0,0,1)x_(1,0,1,0) x_(0,0,0,0)x_(0,0,1,0)x_(1,0,0,0)x_(1,0,1,0) x_(0,0,0,0)x_(0,0,0,1)x_(1,0,0,0)x_(1,0,0,1) x_(0,1,0,0)x_(0,1,0,1)x_(0,1,1,0)x_(0,1,1,1) x_(0,0,1,0)x_(0,0,1,1)x_(0,1,1,0)x_(0,1,1,1) x_(0,0,0,1)x_(0,0,1,1)x_(0,1,0,1)x_(0,1,1,1) x_(0,0,0,1)x_(0,0,1,0)x_(0,1,0,0)x_(0,1,1,1) x_(0,0,0,0)x_(0,0,1,1)x_(0,1,0,1)x_(0,1,1,0) x_(0,0,0,0)x_(0,0,1,0)x_(0,1,0,0)x_(0,1,1,0) x_(0,0,0,0)x_(0,0,0,1)x_(0,1,0,0)x_(0,1,0,1) x_(0,0,0,0)x_(0,0,0,1)x_(0,0,1,0)x_(0,0,1,1) x_(0,1,1,1)x_(1,0,1,1)x_(1,1,0,1)x_(1,1,1,0)x_(1,1,1,1) x_(0,1,1,0)x_(1,0,1,0)x_(1,1,0,0)x_(1,1,1,0)x_(1,1,1,1) x_(0,1,0,1)x_(1,0,0,1)x_(1,1,0,0)x_(1,1,0,1)x_(1,1,1,1) x_(0,0,1,1)x_(1,0,0,1)x_(1,0,1,0)x_(1,0,1,1)x_(1,1,1,1) x_(0,0,1,1)x_(0,1,0,1)x_(0,1,1,0)x_(0,1,1,1)x_(1,1,1,1) x_(0,1,0,0)x_(1,0,0,0)x_(1,1,0,0)x_(1,1,0,1)x_(1,1,1,0) x_(0,0,1,0)x_(1,0,0,0)x_(1,0,1,0)x_(1,0,1,1)x_(1,1,1,0) x_(0,0,1,0)x_(0,1,0,0)x_(0,1,1,0)x_(0,1,1,1)x_(1,1,1,0) x_(0,0,0,1)x_(1,0,0,0)x_(1,0,0,1)x_(1,0,1,1)x_(1,1,0,1) x_(0,0,0,1)x_(0,1,0,0)x_(0,1,0,1)x_(0,1,1,1)x_(1,1,0,1) x_(0,0,0,0)x_(1,0,0,0)x_(1,0,0,1)x_(1,0,1,0)x_(1,1,0,0) x_(0,0,0,0)x_(0,1,0,0)x_(0,1,0,1)x_(0,1,1,0)x_(1,1,0,0) x_(0,0,0,1)x_(0,0,1,0)x_(0,0,1,1)x_(0,1,1,1)x_(1,0,1,1) x_(0,0,0,0)x_(0,0,1,0)x_(0,0,1,1)x_(0,1,1,0)x_(1,0,1,0) x_(0,0,0,0)x_(0,0,0,1)x_(0,0,1,1)x_(0,1,0,1)x_(1,0,0,1) x_(0,0,0,0)x_(0,0,0,1)x_(0,0,1,0)x_(0,1,0,0)x_(1,0,0,0) |
             size => 4
             smallComplex => MutableHashTable{...5...}
                          1       16       80       160       120       16
             complex => QQ  &lt;-- QQ   &lt;-- QQ   &lt;-- QQ    &lt;-- QQ    &lt;-- QQ
                                                                       
                        -1      0        1        2         3         4</code></pre>
            </td>
          </tr>
        </table>
        <div>
          <p>Most of the data is self explanatory. The key <span class="tt">complex</span> stores the simplicial chain complex of the entire simplicial complex. Access to the simplicial complex and related data (e.g.: faces, chain complex, boundary map, homology, and Stanley-Reisner ideal) is provided by the Macaulay2 package <a title="exploring abstract simplicial complexes within commutative algebra" href="/usr/share/doc/Macaulay2/SimplicialComplexes/html/index.html">SimplicialComplexes</a>. The small complexes, their homology, and characters are stored in dedicated keys to avoid computing them again. The key <span class="tt">cubicDimension</span> contains a hash table indexed by integers. Each key $p$ corresponds to a table, indexed by dimension, indicating the position of the faces of $\{0\}^{n-p} \times X^{p,r}$ within the list of faces of $X^{n,r}$. This information is used to form the differentials in the small complexes by restricting the boundary map of $X^{n,r}$ to $\{0\}^{n-p} \times X^{p,r}$.</p>
        </div>
      </div>
      <div>
        <h2>Caveat</h2>
        <div>
          <p>The computers we have access to are able to perform the computations above for $X^{6,3}$, which take a few minutes, but crash Macaulay2 upon attempting $X^{6,4}$.</p>
        </div>
      </div>
      <div>
        <div>
          <div>
            <h2>Authors</h2>
            <ul>
              <li><a href="http://math.galetto.org">Federico Galetto</a><span> &lt;<a href="mailto:galetto.federico%40gmail.com">galetto.federico@gmail.com</a>></span></li>
              <li><a href="https://j-montano.github.io">Jonathan Montaño</a><span> &lt;<a href="mailto:montano%40asu.edu">montano@asu.edu</a>></span></li>
              <li><a href="https://zoewellner.com/">Zoe Wellner</a><span> &lt;<a href="mailto:zoe%40wellners.net">zoe@wellners.net</a>></span></li>
            </ul>
          </div>
          <div>
            <h2>Version</h2>
            <p>This documentation describes version <b>0.1</b> of VietorisRipsHypercube, released <b>December 17, 2025</b>.</p>
          </div>
          <div>
            <h2>Citation</h2>
            <p>If you have used this package in your research, please cite it as follows:</p>
            <table class="examples">
              <tr>
                <td>
                  <pre><code class="language-bib">@misc{VietorisRipsHypercubeSource,
  title = {{VietorisRipsHypercube: homology of VR complexes of hypercubes as group representations. Version~0.1}},
  author = {Federico Galetto and Jonathan Montaño and Zoe Wellner},
  howpublished = {A \emph{Macaulay2} package available at
    \url{https://github.com/galettof/VietorisRipsHypercube}}
}
</code></pre>
                </td>
              </tr>
            </table>
          </div>
          <div>
            <h2>Exports</h2>
            <div class="exports">
              <ul>
                <li>Functions and commands                  <ul>
                    <li><kbd>differential</kbd></li>
                    <li><kbd>smallComplex</kbd></li>
                    <li><kbd>VRHypercube</kbd></li>
                  </ul>
                </li>
                <li>Methods                  <ul>
                    <li><kbd>VRHypercube(ZZ,ZZ)</kbd></li>
                  </ul>
                </li>
                <li>Symbols                  <ul>
                    <li><kbd>cubicDimension</kbd></li>
                    <li><kbd>scale</kbd></li>
                  </ul>
                </li>
              </ul>
            </div>
          </div>
        </div>
        <div class="waystouse">
          <h2>For the programmer</h2>
          <p>The object <a title="homology of VR complexes of hypercubes as representation" href="index.html">VietorisRipsHypercube</a> is <span>a <a title="the class of all packages" href="/usr/share/doc/Macaulay2/Macaulay2Doc/html/___Package.html">package</a></span>, defined in <span class="tt">VietorisRipsHypercube.m2</span>.</p>
        </div>
        <hr>
        <div class="waystouse">
          <p>The source of this document is in <span class="tt">VietorisRipsHypercube.m2:526:0</span>.</p>
        </div>
      </div>
    </div>
  </body>

</html>
