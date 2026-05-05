from itertools import combinations, product
import numpy as np
from scipy.sparse import coo_matrix, save_npz
import time

# This code is designed to create matrices for certain differentials
# in the zeroth-page of the spectral sequence for scale 3

# TO USE: set a value of p from 5 to 9 for the cubic dimension
c=8
# then run "python matrix_generators.py" from a terminal
# opened to the directory containing this file
# the matrices will be saved in the same directory
# By default, the matrices are saved in SciPy NPZ format
# set the next variable to true to out in Magma format
magma=False

# method that converts a set of tuples of binary sequences into a tuple of decimal numbers.
def convert(bag):
    return {
        tuple(
            sum(bit << i for i, bit in enumerate(reversed(inner)))
            for inner in outer)
        for outer in bag}


# method to generate the special permutations used in method act
# these permutations represent distinc cosets of Hn modulo
# a subgroup of the stabilizer of the simplices as in Prop. 6.1
def block_perms(n, k):
    A = tuple(range(2, n - k + 1))
    B = (0,) + tuple(range(n - k + 1, n))

    for pos_A in combinations(range(n), len(A)):
        rem1 = [i for i in range(n) if i not in pos_A]
        for pos_B in combinations(rem1, len(B)):
            perm = [None] * n

            for i, p in enumerate(pos_A):
                perm[p] = A[i]

            for i, p in enumerate(pos_B):
                perm[p] = B[i]

            for p in range(n):
                if perm[p] is None:
                    perm[p] = 1

            yield tuple(perm)


# This method returns a set with the action of the hyperoctahedral group on an individual simplex
def act(simplex, k):
    n = len(simplex[0])
    return {
        tuple(sorted(tuple(inner[i] ^ mask[i] for i in perm) for inner in simplex))
        for mask in product((0, 1), repeat=n)
        for perm in block_perms(n, k)}


# The methods below produce an accurate answer only for p <= 10
# where p is the cubic dimension

#  method to generate the initial bag of simplices of dimension p-2
def initial_bag(p=5):
    # the vertex with all zero entries
    def zeros():
        return (0,) * p

    # the vertices v_i of Prop. 6.1
    def middle_pattern(i):
        return (1, 1) + (0,) * i + (1,) + (0,) * (p - 3 - i)

    # the vertices w_i of Prop 6.1
    def edge_pattern(k):
        return (1,) + (0,) * (p - k - 1) + (1,) + (0,) * (k - 1)

    bags = [
        (zeros(),) + tuple(middle_pattern(i) for i in range(p - 2)),
        (zeros(), edge_pattern(1)) + tuple(middle_pattern(i) for i in range(p - 3)),
        (zeros(), edge_pattern(2), edge_pattern(1)) + tuple(middle_pattern(i) for i in range(p - 4)),
    ]

    if p >= 9:
        bags.append(
            (zeros(), edge_pattern(3), edge_pattern(2), edge_pattern(1))
            + tuple(middle_pattern(i) for i in range(p - 5))
        )

    return tuple(bags)


# method to generate the candidate vertices to add to the dim p-2 simplices
# the candidate vertices are obtained by inspection
# note that these vertices change the stabilizer but we skirt the issue by
# introducing all their allowable permutations as additional candidates
def temp_bags(p=5):
    temp_bag1 = {
        ((1,) + (0,) * (p - 1),),
        ((0, 1) + (0,) * (p - 2),),
        ((1, 1) + (0,) * (p - 2),)}
    for i in range(p - 2):
        temp_bag1.update({
            ((1, 0) + (0,) * i + (1,) + (0,) * (p - 3 - i),),
            ((0, 1) + (0,) * i + (1,) + (0,) * (p - 3 - i),)})

    temp_bag2 = {
        ((1,) + (0,) * (p - 1),),
        ((0, 1) + (0,) * (p - 2),),
        ((1, 1) + (0,) * (p - 2),),
        ((0, 1) + (0,) * (p - 3) + (1,),),
        ((1, 1) + (0,) * (p - 3) + (1,),)}
    for i in range(p - 3):
        temp_bag2.add(((1, 0) + (0,) * i + (1,) + (0,) * (p - 3 - i),))

    temp_bag3 = {
        ((1,) + (0,) * (p - 1),),
        ((0, 1) + (0,) * (p - 2),),
        ((1, 1) + (0,) * (p - 2),),
        ((1, 1) + (0,) * (p - 4) + (1, 0),),
        ((1, 1) + (0,) * (p - 3) + (1,),)}
    for i in range(p - 4):
        temp_bag3.add(((1, 0) + (0,) * i + (1,) + (0,) * (p - 3 - i),))

    temp_bag4 = {
        ((1,) + (0,) * (p - 1),),
        ((0, 1) + (0,) * (p - 2),),
        ((1, 1) + (0,) * (p - 2),),
        ((1, 1) + (0,) * (p - 5) + (1, 0, 0),),
        ((1, 1) + (0,) * (p - 4) + (1, 0),),
        ((1, 1) + (0,) * (p - 3) + (1,),)}
    for i in range(p - 5):
        temp_bag4.add(((1, 0) + (0,) * i + (1,) + (0,) * (p - 3 - i),))

    if p >= 9:
        return temp_bag1, temp_bag2, temp_bag3, temp_bag4
    else:
        return temp_bag1, temp_bag2, temp_bag3


#  method to generate the simplices of dimension p-2
def pminus2(p=5):
    return convert({t for n in range(4 if p >= 9 else 3) for t in act(initial_bag(p)[n], n + 1)})


#  method to generate the simplices of dimension p-1
def pminus1(p=5):
    oo = initial_bag(p)
    temps = temp_bags(p)
    return convert({t for n in range(4 if p >= 9 else 3) for a in temps[n] for t in act(oo[n] + a, n + 1)})


#  method to generate the simplices of dimension p
def dimp(p=5):
    oo = initial_bag(p)
    temps = temp_bags(p)
    return convert({t
                    for n in range(4 if p >= 9 else 3)
                    for a, b in combinations(temps[n], 2)
                    # check we are not adding two vertices at distance >3
                    if sum(abs(x - y) for x, y in zip(a[0], b[0])) < 4
                    for t in act(oo[n] + a + b, n + 1)})

#  method to generate the simplices of dimension p+1
def pplus1(p=5):
    oo = initial_bag(p)
    temps = temp_bags(p)
    return convert({t
                    for n in range(4 if p >= 9 else 3)
                    for a, b, c in combinations(temps[n], 3)
                    # check we are not adding two vertices at distance >3
                    if all(sum(abs(x - y) for x, y in zip(u[0], v[0])) < 4 for u, v in ((a, b), (a, c), (b, c)))
                    for t in act(oo[n] + a + b + c, n + 1)})


# method to built the differentials of the zeroth-page over ZZ/2
# by restricting the usual boundary map of simplicial complexes
# working over ZZ/2 means we can ignore signs
# matrices are saved in SciPy NPZ file
def make_boundary(p,cols_faces, rows_faces, k, name):
    row_index = {face: i for i, face in enumerate(rows_faces)}

    rows, cols, data = [], [], []

    for c, face in enumerate(cols_faces):
        for i in range(k):
            t = face[:i] + face[i+1:]
            if t in row_index:
                rows.append(row_index[t])
                cols.append(c)
                data.append(1)

    rows = np.array(rows, dtype=np.int64)
    cols = np.array(cols, dtype=np.int64)
    data = np.array(data, dtype=np.int8)

    shape = (len(rows_faces), len(cols_faces))
    matrix = coo_matrix((data, (rows, cols)), shape=shape).tocsr()

    # the matrices are saved in a sparse format
    # see the docs of scipy.sparse for more info
    save_npz(f"matrix{p}{name}.npz", matrix)


# as above but matrices are saved in Magma sparse matrix format
# Warning: this function gives the transpose of the actual matrices
def make_boundary_magma(p,cols_faces, rows_faces, k, name):
    row_index = {face: i for i, face in enumerate(rows_faces)}

    # create a list to be converted to a flat sequence of integers
    # following the second Magma specification for sparse matrices
    # we construct the transpose of the boundary map since it's easier
    seq = []
    middle = ""
    for c, face in enumerate(cols_faces):
        temp_list = []
        for i in range(k):
            t = face[:i] + face[i+1:]
            if t in row_index:
                temp_list.append(row_index[t]+1)
        temp_list.sort()
        temp_str = ',1,'.join(map(str,temp_list))
        middle += str(len(temp_list)) + "," + temp_str + ",1,"

    # NOT DONE!
    # Magma string
    start = f"M{name}:=SparseMatrix(GF(2),{len(cols_faces)},{len(rows_faces)},["
    middle = middle[:-1]
    finish = "]);"
    magma_string = start + middle + finish
    with open(f"matrix{p}{name}.txt", "w") as text_file:
        print(magma_string, file=text_file)


# method to generate the matrices of differentials. It generates three matrices for p=8 and two otherwise
def build_matrices(magma,p=5):

    # find faces
    print(f"Looking for faces of diameter 3 and cubic dimension {c}.\n")
    start_time = time.time()
    facesm2 = list(pminus2(p))
    print(f"Found {len(facesm2)} faces of dimension {c-2} in {round(time.time() - start_time,5)} seconds.")
    start_time = time.time()
    facesm1 = list(pminus1(p))
    print(f"Found {len(facesm1)} faces of dimension {c-1} in {round(time.time() - start_time,5)} seconds.")
    start_time = time.time()
    facesp  = list(dimp(p))
    print(f"Found {len(facesp)} faces of dimension {c} in {round(time.time() - start_time,5)} seconds.")
    if p < 9:
        start_time = time.time()
        facesplus1 = list(pplus1(p))
        print(f"Found {len(facesplus1)} faces of dimension {c+1} in {round(time.time() - start_time,5)} seconds.")

    # construct differentials
    print(f"\nComputing differentials for cubic dimension {c} over ZZ/2.\n")
    start_time = time.time()
    if magma:
        make_boundary_magma(p, facesm1, facesm2, p,   "1")
    else:
        make_boundary(p, facesm1, facesm2, p,   "1")
    print(f"Computed first differential in {round(time.time() - start_time,5)} seconds.")
    start_time = time.time()
    if magma:
        make_boundary_magma(p, facesp,  facesm1, p+1, "2")
    else:
        make_boundary(p, facesp,  facesm1, p+1, "2")
    print(f"Computed second differential in {round(time.time() - start_time,5)} seconds.")
    if p < 9:
        start_time = time.time()
        if magma:
            make_boundary_magma(p, facesplus1, facesp, p+2, "3")
        else:
            make_boundary(p, facesplus1, facesp, p+2, "3")
        print(f"Computed third differential in {round(time.time() - start_time,5)} seconds.")

    # final warnings
    print("\nWarning! Faces of diameter 2 are NOT included in the above.")
    if magma:
        print("Warning! The Magma matrices computed here are transpose of the actual differentials.")

# this runs the computations
build_matrices(magma,c)
