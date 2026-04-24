newPackage(
     "VietorisRipsHypercube",
     Version => "1.0",
     Date => "April 23, 2026",
     AuxiliaryFiles => false,
     Authors => {
	 {Name => "Federico Galetto",
	     Email => "galetto.federico@gmail.com",
	     HomePage => "http://math.galetto.org"},
	 {Name => "Jonathan Montaño",
	     Email => "montano@asu.edu",
	     HomePage => "https://j-montano.github.io"},
	 {Name => "Zoe Wellner",
	     Email => "zoe@wellners.net",
	     HomePage => "https://zoewellner.com/"},
	   },
     Headline => "homology of VR complexes of hypercubes as group representations",
     HomePage => "https://github.com/galettof/VietorisRipsHypercube",
     DebuggingMode => false,
     PackageExports => {"SimplicialComplexes","BettiCharacters"},
     Keywords => {"Algebraic Topology","Representation Theory"},
     )

export {
    "VRHypercube",
    "scale",
    "differential"
    }

-- hamming distance between binary sequences
-- unexported
distance = memoize (
    (u,v) -> sum toList apply(u,v, (x,y) -> abs(x-y))
    )

-- type for collecting data on VR complex of hypercube
VRData = new Type of HashTable
-- and its printing method
net(VRData) := String => X -> (
    "Vietoris-Rips complex at scale " | toString(X#(symbol scale)) |
    " of Q" | toString(X#(symbol size))
    )

-- VR complex of hypercube
-- 1st parameter: n, for hypercube Q_n
-- 2nd parameter: r, VR scale
VRHypercube = method(Options => {MonomialSize=>16})
VRHypercube(ZZ,ZZ) := VRData => opts -> (n,r) -> (
    -- polynomial ring with variables indexed by binary sequences
    x := getSymbol "x";
    R := QQ[x_(n:0)..x_(n:1),MonomialSize=>opts.MonomialSize];
    if n<1 then error "first argument must be a positive integer";
    if r<0 or r>n then error "second argument out of range";
    -- Stanley-Reisner ideal of the VR simplicial complex
    I := ideal apply(
	-- pairs of variables with indices at distance > r
	select(subsets(gens R,2),
	    y -> distance(last baseName y_0,last baseName y_1) > r
	    ),
	-- multiply those variables
	product);
    Δ := simplicialComplex I;
    new VRData from (
	{
	    (symbol scale) => r,
	    (symbol size) => n,
	    (symbol ring) => R,
	    (symbol ideal) => I,
	    (symbol simplicialComplex) => Δ,
	    (symbol complex) => complex Δ,
	    (symbol MonomialSize) => opts.MonomialSize
	    } |
	for i to n list ( i => new MutableHashTable )
	)
    )

-- check all vertices of a simplex are indexed by binary
-- sequences with z initial zeros
-- the simplex is passed as a monomial m
-- unexported
hasInitialZeros = (m,z) -> (
    -- extract indexing sequences
    s := apply(support m, v -> last baseName v);
    -- check if all first z entries are zero
    heads := apply(s, u -> take(u,z));
    if any(heads, u -> isMember(1,u)) then false else true
    )

-- check that a simplex has cubic dimension c
-- the simplex is passed as a monomial m
-- unexported
hasCubicDimension = (m,c) -> (
    -- extract indexing sequences
    s := apply(support m, v -> last baseName v);
    -- find length of sequence
    n := length first s;
    rearranged := pack(#s, mingle s);
    -- count how many different entries are in each (1 or 2)
    counts := apply(rearranged, l -> length unique l);
    -- cubic dimension is the number of 2's
    if number(counts, i -> i==2) == c then true else false
    )

-- get the complex C^{c,r}_{c+\bullet,c}
-- throughout the comments, this is referred to as a "small" complex
-- this is the complex of the faces of cubic dimension c
-- in {0}^(n-c) x X^{c,r} inside X^{n,r}
-- 1st paramater: Δ, the VR complex of Q_R at scale r
-- 2nd parameter: c, cubic dimension
complex(VRData,ZZ) := {} >> opts -> (X,c) -> (
    -- find length of indexing sequences
    n := X#(symbol size);
    -- check parameter is viable
    if c<0 or c>n then error "integer index out of range";
    -- check if small complex is precomputed
    if not X#c#?(symbol complex) then (
	-- extract saved information
	Δ := X#(symbol simplicialComplex);
	C := X#(symbol complex);
	-- length of C
	l := length C;
	-- find the position of the faces of cubdim c inside
	-- the list of all faces of Δ hashed by dimension
	pos := hashTable parallelApply(toList(0..(l-1)), i-> (
		(i, positions(faces(i,Δ),
			m -> hasInitialZeros(m,n-c) and hasCubicDimension(m,c)))
		));
	-- extract the subcomplex of C corresponding to faces of
	-- cubdim c by restricting the boundary maps to those faces
	-- this gives a list of matrices
	mat := parallelApply(toList(1..l-1), i -> (C.dd_i)^(pos#(i-1))_(pos#i));
	-- store position of faces
	X#c#(symbol faces) = pos;
	-- assemble matrices into a complex, store and return
	X#c#(symbol complex) = complex mat;
	-- create mutable hash tables to store future computations
	-- such as homology modules, presentations, and characters
	X#c#(symbol homology) = new MutableHashTable;
	X#c#(symbol presentation) = new MutableHashTable;
	X#c#(symbol character) = new MutableHashTable;
	-- also store the map from the polynomial ring of the simplicial complex
	-- to the skew polynomial ring of the complex for cubic dimension c
	-- obtained by dropping initial n-c zeros of indexing binary sequences
	R := X#(symbol ring);
	e := getSymbol "e";
	S := QQ[e_(c:0)..e_(c:1),SkewCommutative=>true,MonomialSize=>X#(symbol MonomialSize)];
	phi := map(S,R,apply(gens R, u -> e_(take(last baseName u,-c))));
	X#c#(symbol map) = phi;
	);
    X#c#(symbol complex)
    )

-- homology of small complex of cubic dimension c
homology(VRData,ZZ) := opts -> (X,c) -> (
    prune HH complex(X,c)
    )

-- minimal presentation for the i-th homology of the small complex
-- of cubic dimension c as a module over the polynomial ring
-- with variables indexed by binary sequences
homology(VRData,ZZ,ZZ) := opts -> (X,c,i) -> (
    -- find length of indexing sequences
    n := X#(symbol size);
    -- check 1st parameter is viable
    if c<0 or c>n then error "1st integer index out of range";
    -- get small complex of cubic dimension c
    -- triggers computation if not previously done
    C := complex(X,c);
    if i<(min C) or i>(max C) then error "2nd integer index out of range";
    -- check if homology is precomputed
    if not X#c#(symbol homology)#?i then (
	-- compute homology of the small complex, which triggers
	-- creation of necessary data if missing
	X#c#(symbol presentation)#i = trim HH_i C;
	-- get simplicial complex, face indices, and ring map
	Δ := X#(symbol simplicialComplex);
	pos := X#c#(symbol faces);
	phi := X#c#(symbol map);
	-- get all i-faces in the small complex
	f := (faces(i,Δ))_(pos#i);
	-- row matrix of monomials corresponding to faces
	-- of dimension i and cubic dimension c
	m := phi( matrix{f} );
	-- convert generators and relations of homology to ideals
	I := ideal (m * (gens X#c#(symbol presentation)#i));
	J := ideal (m * (relations X#c#(symbol presentation)#i));
	-- in degree i+1 I/J is isomorphic to the homology of the
	-- complex of simplices of desired cubic dimension
	X#c#(symbol homology)#i = trim(I/J);
    );
    X#c#(symbol homology)#i
    )

-- set up Hn-action on the ring of the small complex
-- unexported
hActBinSeq = R -> (
    -- length of binary sequences
    n := length last baseName first gens R;
    -- create representatives of the conjugacy classes of Hn
    -- for each bipartition (p,q) of n we get a rep (s,sigma)
    -- p is the cycle type of the cycles without signs
    -- q is the cycle type of the cycles with a minus in 1st place
    -- reps are sorted with |p| from 0 to n, then q lexicographically,
    -- then p, lexicographically;
    -- this sorting makes the identity the last element
    reps := flatten for i to n list (
	flatten for p in partitions(n-i) list (
	    for q in partitions(i) list (
		L := toList(0..n-1);
		sigma := flatten for u in (toList p)|(toList q) list (
		    l := take(L,u);
		    L = drop(L,u);
		    rotate(1,l)
		    );
		s := toList(n-i:0);
		s = s | flatten for u in q list (
		    {1} | toList(u-1:0)
		    );
		(s,sigma)
		)
	    )
	);
    -- the pair (s,sigma) acts on a 01-sequence u by sending it to v
    -- where v_j = ( u_(sigma(j)) + s_j ) mod 2
    -- both s and sigma are lists, u is a sequence
    -- for each representative (s,sigma) send the variable e_u to e_v
    for g in reps list (
	s := g_0;
	sigma := g_1;
	matrix {
	    for f in gens R list (
		u := last baseName f;
		v := toSequence apply(u_sigma + s, i -> i % 2);
		value ( (first baseName f)_v )
		)
	    }
	)
    )

-- compute character of homology of small complex using BettiCharacters
character(VRData,ZZ,ZZ) := opts -> (X,c,i) -> (
    -- get i-th homology of small complex of cubic dimension c
    -- triggers computation if not previously done
    M := homology(X,c,i);
    -- check if character was precomputed
    if not X#c#(symbol character)#?i then (
	-- set up the action
	Hn := hActBinSeq ring M;
	A := action(M,Hn);
	-- compute the character in the right degree
	x := character(A,i+1);
	-- shift it to homological degree i and return it
	X#c#(symbol character)#i = x[-i];
	);
    X#c#(symbol character)#i
    )

-- compute the differential of the complex E^1_{*,0}, i.e.,
-- the horizontal complex at height 0 in the first page
-- of the spectral sequence
-- most of this complex should just be the cellular chain
-- complex of a cube, and computing the differential allows
-- us to confirm
differential = method();
differential(VRData,ZZ) := (X,t) -> (
    -- extract saved information
    Δ := X#(symbol simplicialComplex);
    C := X#(symbol complex);
    -- find length of indexing sequences
    n := X#(symbol size);
    if t<1 or t>n then error "integer index outside range";
    -- compute homology, in case it is missing
    homology(X,t,t);
    -- domain indices
    dom := X#t#(symbol faces)#t;
    -- codomain indices
    f := faces(t-1,Δ);
    cod := positions(f, m -> hasInitialZeros(m,n-t) and hasCubicDimension(m,t-1));
    -- matrix of differential
    mat := (C.dd_t)^cod_dom;
    -- get ring map
    phi := X#t#(symbol map);
    -- get generators of homology
    g := gens X#t#(symbol presentation)#t;
    -- get image of generators under differential
    im := mat*g;
    -- get image as matrix of polynomials
    M := (phi matrix{ f_cod }) * im;
    -- get relations on homology
    I := ideal relations homology(X,t-1,t-1);
    -- generate all relations with one longer binary sequences
    S := target phi;
    J := sum flatten for i to t-1 list (
	for j to 1 list (
	    L := apply(gens ring I,
		u -> (first baseName u)_(insert(i,j,last baseName u))
		);
	    psi := map(S,ring I,L);
	    psi I
	    )
	);
    -- map to the quotient
    psi := map(S/J,S);
    flatten entries psi M
    )

beginDocumentation()

doc ///

Node
    Key
	VietorisRipsHypercube
	VRHypercube
	(VRHypercube, ZZ, ZZ)
	scale
	(complex, VRData, ZZ)
	differential
	(differential, VRData, ZZ)
    Headline
	homology of VR complexes of hypercubes as representation
    Description
	Text
	    This Macaulay2 package implements computations described in the article
	    "TITLE"
	    by Federico Galetto, Jonathan Montaño, and Zoe Wellner.
	    The most up-to-date version of this package can be found at
	    @HREF "https://github.com/galettof/VietorisRipsHypercube"@.

	    @SUBSECTION "Background and Scope"@
	    
	    The goal is to compute the homology groups of the Vietoris-Rips complex of
	    the hypercube $Q_n$ at a given scale $r$ with respect to the Hamming
	    distance, and to describe them as representations of the hyperoctahedral
	    group $\mathfrak{H}_n$, i.e., the automorphism group of $Q_n$.
	    We denote this simplicial complex as $X^{n,r}$. A simplex $\Delta$ in $X^{n,r}$ has
	    cubic dimension $p$ if it can be isometrically embedded in $Q_p$ but not
	    in $Q_{p-1}$. The simplices with cubic dimension at most $p$ form a
	    subcomplex $X^{n,r}_p$ of $X^{n,r}$. In fact, these subcomplexes form a
	    filtration $\{X^{n,r}_p\}_{p\geqslant 0}$ of $X^{n,r}$ by simplicial
	    complexes, and the homology of $X^{n,r}$ can be computed via the
	    spectral sequence associated with this filtration. Let $C^{n,r}_{p+q,p}$
	    be the vector space spanned by all the simplices of dimension $p+q$ and
	    cubic dimension $p$. The terms in the zero page of the spectral sequence
	    are $E^0_{p,q} = C^{n,r}_{p+q,p}$. It is enough to study the subcomplexes
	    corresponding to the smaller simplicial complexes
	    $\{0\}^{n-p} \times X^{p,r}$. Then, the symmetry afforded by the action of
	    $\mathfrak{H}_n$ allows us to recover the result for the bigger complex.
	    More formally, for every integer $i$, we have an isomorphism
	    $$C^{n,r}_{i,p} \cong
	    \operatorname{Ind}^{\mathfrak{H}_n}_{\mathfrak{S}_{n-p} \times \mathfrak{H}_p}
	    \left( \{n-p\} \otimes C^{p,r}_{i,p} \right)$$
	    of representations of $\mathfrak{H}_n$, where $\{n-p\}$ denotes the
	    1-dimensional trivial representation of $\mathfrak{S}_{n-p}$. Since
	    induction between finite groups is an exact functor, we then get
	    $$E^1_{p,q} \cong
	    \operatorname{Ind}^{\mathfrak{H}_n}_{\mathfrak{S}_{n-p} \times \mathfrak{H}_p}
	    \left( \{n-p\} \otimes H_{p+q} (C^{p,r}_{p+\bullet,p}) \right).$$
	    This reduces the problem to computing the homology of the complex
	    $C^{p,r}_{p+\bullet,p}$, which we refer to as the complex for cubic
	    dimension $p$. Note that $C^{p,r}_{p+\bullet,p}$ does not depend on
	    $n$, so it is finite, which makes it possible to build it and compute
	    its homology using software.

	    @SUBSECTION "Basic Functionality"@

	    The methods in this package allow users to construct the complexes
	    $C^{p,r}_{p+\bullet,p}$, to
	    compute their homology groups, and to describe those homology groups as
	    representations of a suitable hyperoctahedral group by computing their
	    characters. We illustrate the workflow for $X^{n,2}$, the Vietoris-Rips
	    complex at scale $r=2$ of $Q_n$. First, we construct the simplicial
	    complex for $n=4$.
	Example
	    X42 = VRHypercube(4,2)
	Text
	    Next, we form the complex $C^{2,2}_{2+\bullet,2}$ for cubic dimension 2
	    and compute its homology.
	Example
	    complex(X42,2)
	    homology(X42,2)
	Text
	    We see that the homology is a 1-dimensional vector space in homological
	    degree 2. The following commands give a more detailed description of this
	    homology group.
	Example
	    H22 = homology(X42,2,2)
	    netList flatten entries generators H22
	    netList flatten entries relations H22
	Text
	    We interpret these results. The simplicial complex $X^{2,2}$ has four
	    vertices: $(0,0)$, $(0,1)$, $(1,0)$, and $(1,1)$. The vertex $u$ corresponds
	    to a variable $e_u$ in a polynomial ring. The squarefree monomial
	    $e_{u_1} \cdots e_{u_t}$ corresponds to the simplex $[u_1,\dots,u_t]$.
	    Note that the monomial ordering on the polynomial ring determines the
	    orientation of the simplex. Then, a linear combination of monomials
	    corresponds to a formal linear combination of simplices. In our case,
	    $-[(0,0),(0,1),(1,1)]+[(0,0),(1,0),(1,1)]$ is a linear combination of
	    triangles that generates $H_2 (C^{2,2}_{2+\bullet,2})$.
	    This linear combination corresponds to one of the two ways of
	    dividing a square into two triangles by tracing a diagonal. The relation
	    on the homology says that the other way of dividing the square is equivalent.

	    Now, we compute the character of this homology group as a representation
	    of $\mathfrak{H}_2$.
	Example
	    c22 = character(X42,2,2)
	Text
	    We proceed to decompose this into irreducible characters.
	    This functionality requires at least version 2.6
	    of the @TT "BettiCharacters"@ package, which can be downloaded at
	    @HREF "https://github.com/galettof/BettiCharacters"@.
	    First, we construct the character table of $\mathfrak{H}_2$. Then, we
	    express the character above as a positive integer linear combination
	    of irreducible characters.
	Example
	    T2 = hyperoctahedralGroupTable(2,QQ)
	    c22 / T2
	Text
	    We see that our homology group is in fact isomorphic to the irreducible
	    representation indexed by the pair of partitions $(\varnothing; 1^2)$.

	    When we perform the same steps for cubic dimension 3, we obtain the
	    following.
	Example
	    complex(X42,3)
	    homology(X42,3)
	    H33 = homology(X42,3,3);
	    netList flatten entries generators H33
	    netList flatten entries relations H33
	    c33 = character(X42,3,3)
	    T3 = hyperoctahedralGroupTable(3,QQ);
	    c33 / T3
	Text
	    We deduce that the homology of $C^{3,2}_{3+\bullet,3}$
	    is a 2-dimensional vector space in homological degree 3. We get two
	    explicit generators and there are no relations. As a represenatation of
	    $\mathfrak{H}_3$ it is isomorphic to the direct sum of the irreducibles
	    indexed by the pairs of partitions $(1^3;\varnothing)$ and
	    $(\varnothing;1^3)$.

	    Finally, we observe that the complex $C^{4,2}_{4+\bullet,4}$ for
	    cubic dimension 4 has no homology.
	Example
	    complex(X42,4)
	    homology(X42,4)

	Text
	    
	    @SUBSECTION "Additional Features"@
	Text
	    The @TT "differential"@ method can be used to study certain maps in the
	    first page of the spectral sequence. Specifically, consider the differentials
	    $d^1_{i,0} \colon E^1_{i,0} \to E^1_{i-1,0}$ with $1\leqslant i\leqslant n$.
	    The @TT "differential"@ method takes the elements of $E^1_{i,0}$ that
	    generate the homology group $H_i (C^{i,r}_{i+\bullet,i})$ (as computed
		using @TT "homology"@) and returns their image under $d^1_{i,0}$.
	    Once again, the result is expressed as a linear combination of squarefree
	    monomials corresponding to simplices.

	    For the example above, we obtain the following.
	Example
	    netList differential(X42,1)
	    netList differential(X42,2)
	Text
	    We see that $d^1_{1,0}$ sends the edge $[0,1]$ to its vertices, and
	    $d^1_{2,0}$ sends the square on the vertices $(00)$, $(01)$, $(10)$, and
	    $(11)$ to its boundary.
	Example
	    L = differential(X42,3);
	    netList L
	Text
	    Now, the homology group $E^1_{3,0} = H_3 (C^{3,2}_{3+\bullet,3})$ has two
	    generators, and $d^1_{3,0}$ sends them to the expressions above,
	    which are equal. This means we can change basis in the homology
	    by taking the sum and difference of the generators, so that $d^1_{3,0}$
	    sends one of these new basis elements to zero. The other basis element
	    is sent to twice the boundary of the cube.
	Example
	    netList {L_0+L_1,L_0-L_1}
	Text
	    
	    @SUBSECTION "Technical Details"@

	    All the information about the simplicial complex $X^{n,2}$ is stored in
	    an object of type @TT "VRData"@ which is a dedicated type of
	    @TT "HashTable"@. Let us take a look at the object for the example above.
	Example
	    peek X42
	Text
	    The keys @TT "size"@ and @TT "scale"@ store the parameters $n$ and $r$.
	    We store the simplicial complex, its Stanley-Reisner ideal and ambient
	    polynomial ring, and its simplicial chain complex in dedicated keys.
	    This data, along with other functionality for working with simplicial
	    complexes such as computing faces, boundary maps, and homology, is
	    provided by the Macaulay2 package @TT "SimplicialComplexes"@.
	    
	    The integer key $p$ stores the information about the complex for cubic
	    dimension $p$, after it is computed. For example, let us look at $p=2$.
	Example
	    peek X42#2
	Text
	    The complex $C^{p,r}_{p+\bullet,p}$ is stored under @TT "complex"@.
	    The @TT "faces"@ key holds a hash table indexed by integers.
	    Each key $i$ corresponds to a list holding the position of the $i$-faces
	    of cubic dimension $p$
	    of $\{0\}^{n-p} \times X^{p,r}$ within the list of faces of $X^{n,r}$.
	    This data is used to form the differentials in $C^{p,r}_{p+\bullet,p}$
	    by restricting the boundary map of $X^{n,r}$ to the appropriate faces.
	    The @TT "map"@ key holds a ring homomorphism between the polynomial rings
	    of $X^{n,r}$ and $X^{p,r}$. The ring of $X^{p,r}$ is skew-commutative so
	    that a squarefree monomial represents an oriented simplex; this is crucial
	    for the correct computation of characters.
	    The @TT "homology"@ and @TT "character"@ keys contain mutable hash tables
	    with integer keys corresponding to homological degrees, and values
	    used to store homology (as modules over polynomial rings) and characters
	    after they are computed.
	    The @TT "presentation"@ key holds equivalent information about
	    homology but as a vector space.

	    Finally, the @TT "MonomialSize"@ option specifies the minimum number of bits
	    to be used for storing each exponent in a monomial. Its default value for
	    this package is 16, which speeds up some computations. Users can change this
	    when constructing the VR complex (e.g. @TT "VRHypercube(4,2,MonomialSize=>32)"@).
///

end
