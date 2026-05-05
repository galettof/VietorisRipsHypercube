from scipy.sparse import load_npz
#from matrix_generator import build_matrices
import time

# this code computes the rank/homology of the differentials
# in the zeroth-page of the spectral sequence for scale 3
# there are no parameters and ranks/homology are only computed
# for cubic dimension p=8, dimension 6, 7, 8; and
# for cubic dimension p=9, dimension 7, 8,
# which are the only terms needed for the arguments of the paper

# IMPORTANT: before running this, you need to run matrix_generator.py
# once with c=8, magma=False
# and once with c=9, magma=False
# to create the necessary matrix files

start_time = time.time()

# this function implements Gaussian elimination over ZZ/2
# notice that addition over ZZ/2 operates as XOR (^=)
# currently, this does not take advantage of sparseness
def rank_sparse_Z2(M):
    pivots = {}
    rank = 0

    for i in range(M.shape[0]):

        row = {c for c in M.indices[M.indptr[i]:M.indptr[i + 1]]}

        while row:
            j = min(row)

            if j not in pivots:
                pivots[j] = row
                rank += 1
                break

            row ^= pivots[j]

    return rank


# Fred: by commeting this out, the assumption is that files exist on drive
#matrices are saved as matripi.npz, were i=1,2,3 if p=8, and i=1,2 otherwise
#build_matrices(magma=False,p=8)
#build_matrices(magma=False,p=9)


# rank and homology computations
M = {k: load_npz(f"matrix{k}.npz").tocsr() for k in (81, 82, 83, 91, 92)}
shape = {k: M[k].shape for k in (81, 82, 83, 91, 92)}
rank = {k: rank_sparse_Z2(M[k]) for k in (81, 82, 83, 91, 92)}

info = [
    (8, 81, 7, 6, None),
    (8, 82, 8, 7, 81),
    (8, 83, 9, 8, 82),
    (9, 91, 8, 7, None),
    (9, 92, 9, 8, 91),
]

for p, k, d1, d0, prev in info:
    print(f"p = {p}, matrix from dimension {d1} to {d0}")
    print("matrix size:", shape[k])
    print("matrix rank:", rank[k])

    if prev is not None:
        print(f"homology at dim {d0}", shape[k][0] - rank[prev] - rank[k])

    print()


print("--- %s seconds ---" % (time.time() - start_time))