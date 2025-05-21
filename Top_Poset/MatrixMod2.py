# some functions for doing matrix operations over F_2 using True/False
import numpy as np

def rref2(matrix):
    ''' Returns (A, p), where A is the reduced row echelon form over F_2
    and p is the indices of the pivot columns'''
    A = np.array(matrix, dtype=bool)
    n, m = A.shape
    A_rref = A.copy()

    row = 0
    pivot_cols = []

    # Row reduction (mod 2, using XOR logic)
    for col in range(m):
        # Find pivot
        pivot_row = None
        for r in range(row, n):
            if A_rref[r, col]:
                pivot_row = r
                break
        if pivot_row is None:
            continue

        # Swap pivot row into position
        if pivot_row != row:
            A_rref[[row, pivot_row]] = A_rref[[pivot_row, row]]

        # Eliminate other rows
        for r in range(n):
            if r != row and A_rref[r, col]:
                A_rref[r] ^= A_rref[row]  # XOR for row elimination

        pivot_cols.append(col)
        row += 1
        if row == n:
            break
        
    return A_rref, np.array(pivot_cols)

def im2(matrix):
    ''' Returns a basis for the image of a matrix over F2'''
    A, pivots = rref2(matrix)
    return np.transpose(np.array(matrix)[:, pivots])

def ker2(matrix):
    ''' Returns a basis for the kernel of a matrix over F2'''
    A, pivots = rref2(matrix)
    m = A.shape[1] # number of columns
    
    # Kernel basis: construct for each free variable
    free_vars = [j for j in range(m) if j not in pivots]
    kernel_basis = []

    for free in free_vars:
        vec = np.zeros(m, dtype=bool)
        vec[free] = True
        for i in reversed(range(len(pivots))):
            col = pivots[i]
            row_vals = A[i]
            sum_ = False
            for j in free_vars:
                sum_ ^= row_vals[j] and vec[j]
            vec[col] = sum_
        kernel_basis.append(vec)
    return np.array(kernel_basis).astype(bool)

def extendBasis2(subspace, space):
    '''Given a basis b_1 of subspace and a basis b_2 of the whole space,
    returns vectors b_3=[v_1, ..., v_n] such that b_1 U b_3 is a basis
    of the whole space. Equivalently, b_3 is representatives for a basis
    of the quotient space, which is isomorphic to the orthogonal complement
    subspace^perp

    The bases should be given as [v_j] where each v_j is a row vector
      '''
    # first make an augmented matrix whose columns are A = [ b_1 | b_2 ]
    U = np.transpose(np.array(subspace))
    V = np.transpose(np.array(space))
    dim_subspace = U.shape[1]
    A = np.hstack((U, V))
    pivots = rref2(A)[1]
    orth_pivots = pivots[pivots >= dim_subspace]
    return np.transpose(A[:, orth_pivots])

def homology2(complex : dict[np.ndarray], basis : dict = None,
              as_set : bool = True) -> dict[dict[np.ndarray]]:
    ''' complex is a dict {n: d_n}, where d_n is the boundary map from
    C_n -> C_{n-1}. complex.items should be a contingous set of integers,
    starting at 0 for regular homology or -1 for reduced homology.
    
    Returns a dict {n: {'generators': basis, 'boundaries': basis}}.
    The approach is:
    0. Start at the highest dimension N and compute a basis of ker d_N
    1. At spot n < N, extend a basis of im(d_{n+1}) to ker(d_n)
    '''
    m, M = min(complex), max(complex)
    homology = {
                M: {
                    'generators': ker2(complex[M]), 
                    'boundaries': np.array([[False] * complex[M].shape[1]])
                    }
                }
    for n in reversed(range(m, M)):
        assert n in complex, 'You need to specify all boundary maps in complex'
        d_n , d_n1 = complex[n], complex[n+1]
        ker , im = ker2(d_n), im2(d_n1)
        homology[n] = {'generators': extendBasis2(im, ker), 'boundaries': im}
    
    if as_set and basis:
        def boolToSet(arr, basis):
            return {b for i, b in enumerate(basis) if arr[i]}
        for i in homology:
            homology[i] = {j: [boolToSet(k, basis[i]) 
                               for k in homology[i][j] ] 
                               for j in homology[i]}
    return homology

