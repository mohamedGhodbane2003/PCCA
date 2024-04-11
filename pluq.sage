########################
#  PLUQ factorization  #
#  (compact storage)   #
########################

def PLUQ(A):
    ## input: matrix A over a field, dimensions m x n.
    #
    ## output: (LU,P,Q,rank) such that 
    #  - LU is L and U compactly stored ## <-- give more precise definition
    #  - A = P L U Q is the PLUQ decomposition of A
    #  - with P a permutation (list in [0...m])
    #  - and Q a permutation (list in [0...n])
    #  - rank is the rank of A
    #
    ## Rotations are performed to preserve the row rank profile, so that in the
    # end we have P = rrp+nrrp, where rrp is the row rank profile of A and nrrp
    # is the list of rows not in rrp (both rrp and nrrp are increasing)
    #
    ## Using lexico order + row rotations + column rotations
    # ==> preserves rank profile matrix

    # matrix dimensions
    m,n = A.dimensions()

    # initialize LU with copy of A
    LU = copy(A)

    # initialize permutations with identity permutation
    Q = list(range(n))
    P = list(range(m))

    # currently discovered rank and dimension of left nullspace
    rank = 0
    nullity = 0

    # at each iteration, either rank or nullity is incremented by 1
    while (rank + nullity < m):
        #find column with pivot element on row `rank`, if there is some
        pivot = rank # pivot is at column index >= rank; take first one
        while (pivot < n and LU[rank,pivot] == 0):
            pivot += 1
        if pivot == n:
            P = row_rotation(LU,rank,P)
            #P = row_transposition(LU,rank,m-1-nullity,P)
            nullity += 1
        else:
            #Q = column_rotation(LU,rank,pivot+1,Q)
            Q = column_transposition(LU,rank,pivot,Q)
            for k in range(rank+1,m):
                LU[k,rank] = LU[k,rank]/LU[rank,rank]
                for j in range(rank+1,n):
                    LU[k,j] = LU[k,j] - LU[k,rank]*LU[rank,j]
            rank += 1
    return LU,P,Q,rank


def PLUQ_Crout(A):
    ## same as above, with Crout scheduling of updates
    m,n = A.dimensions()
    LU = copy(A)
    Q = [i for i in range(n)]
    P = [i for i in range(m)]
    nullity = 0
    rank = 0
    while (rank + nullity < m):
        # update row
        LU[rank, rank:] = LU[rank, rank:] - LU[rank, :rank] * LU[:rank, rank:]
        #find column with pivot element on row rank, if there is some
        pivot = rank
        while (pivot < n and LU[rank,pivot] == 0):
            pivot += 1
        if pivot == n:  # no pivot: put row in last position
            #print(f"---rank{rank}--nullity---\n{LU}\n")
            P = row_rotation(LU,rank,P)
            #P = row_transposition(LU,rank,m-1-nullity,P)
            #print(f"---rank{rank}--nullity---\n{LU}\n")
            nullity += 1
        else:  # found pivot
            invpiv = LU[rank, pivot].inverse()
            for i in range(rank+1,m-nullity):
                LU[i,pivot] = LU[i,pivot] - (LU[i, :rank] * LU[:rank, pivot])[0,0]
                LU[i,pivot] = LU[i,pivot] * invpiv
            #print(f"---rank{rank}--rank---\n{LU}\n")
            Q = column_rotation(LU,rank,pivot+1,Q)
            #Q = column_transposition(LU,rank,pivot,Q)
            #print(f"---rank{rank}--rank---\n{LU}\n")
            rank += 1
    return LU,P,Q,rank



###########
#  Utils  #
###########

def column_rotation(A,cstart,cend,C):
    # columns 0, ..., cstart, .., cend-1, ..., n-1
    # are permuted into
    # columns 0, ..., cstart-1, cstart+1, .., cend-1, cstart, cend, ..., n-1
    # If a list C (of length n) is given, it is also
    # permuted accordingly
    n = A.ncols()
    perm_list = list(range(cstart)) + list(range(cstart+1,cend)) + [cstart] + list(range(cend,n))
    perm = Permutation([i+1 for i in perm_list]).inverse()
    A.permute_columns(perm)
    return perm.action(C)

def column_transposition(A,src,tgt,C):
    # columns 0, ..., src, ..., tgt, ..., n-1
    # are permuted into
    # columns 0, ..., tgt, src+1 ..., tgt-1, src, tgt+1, ... n-1
    # If a list R (of length m) is given, it is also
    # permuted accordingly
    A.swap_columns(src,tgt)
    tmp = C[src]
    C[src] = C[tgt]
    C[tgt] = tmp
    return C

def row_rotation(A,row,R):
    # rows 0, ..., row, ..., m-1
    # are permuted into
    # rows 0, ..., row-1, row+1, ..., m-1, row
    # If a list R (of length m) is given, it is also
    # permuted accordingly
    m = A.nrows()
    perm = Permutation([i+1 for i in range(m) if i != row] + [row+1])
    A.permute_rows(perm)
    return perm.action(R)

def row_transposition(A,src,tgt,R):
    # rows 0, ..., src, ..., tgt, ..., m-1
    # are permuted into
    # rows 0, ..., tgt, src+1 ..., tgt-1, src, tgt+1, ... m-1
    # If a list R (of length m) is given, it is also
    # permuted accordingly
    A.swap_rows(src,tgt)
    tmp = R[src]
    R[src] = R[tgt]
    R[tgt] = tmp
    return R

def expand_PLUQ(LU,P,Q,rank):
    m,n = LU.dimensions()
    # retrieve L
    L = Matrix(LU.base_ring(), m, m)
    for j in range(rank):
        for i in range(j+1,m):
            L[i,j] = LU[i,j]
    for i in range(m):
        L[i,i] = 1
    # retrieve U
    U = Matrix(LU.base_ring(), m, n)
    U[:rank,:] = LU[:rank,:]
    for i in range(rank):
        for j in range(i):
            U[i,j] = 0
    return L,U


#############
#  Testing  #
#############

def check_triL(L):
    m,n = L.dimensions()
    if m != n:
        return False
    unit = all([L[i,i] == 1 for i in range(m)]) 
    if not unit:
        return False
    tri = all([L[i,j] == 0 for i in range(m) for j in range(i+1,m)])
    if not tri:
        return False
    return True

def check_triU(U,rank):
    m,n = U.dimensions()
    if rank > m or rank > n:
        return False
    inv = all([U[i,i] != 0 for i in range(rank)]) 
    if not inv:
        return False
    tri = all([U[i,j] == 0 for i in range(m) for j in range(min(i,n))])
    if not tri:
        return False
    return True

def check_many_PLUQ(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        LU,P,Q,rank = PLUQ(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("tall rectangular: ok")
    else:
        print("tall rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,3,6)
        LU,P,Q,rank = PLUQ(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        LU,P,Q,rank = PLUQ(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

    return True

def check_many_PLUQ_Crout(field_prime=2,max_iter=1000):
    field = GF(field_prime)
    i = 0
    correct = True
    while correct and i < max_iter:
        A = Matrix.random(field,6,3)
        LU,P,Q,rank = PLUQ_Crout(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("tall rectangular: ok")
    else:
        print("tall rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,3,6)
        LU,P,Q,rank = PLUQ_Crout(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("wide rectangular: ok")
    else:
        print("wide rectangular: wrong")
        return A

    i = 0
    while correct and i < max_iter:
        A = Matrix.random(field,4,4)
        LU,P,Q,rank = PLUQ_Crout(A)
        L,U = expand_PLUQ(LU,P,Q,rank)
        correct = rank==A.rank()             \
                      and check_triL(L)      \
                      and check_triU(U,rank)
        L.permute_rows(Permutation([i+1 for i in P]).inverse())
        #L.permute_columns(Permutation([i+1 for i in P]).inverse())
        #U.permute_rows(Permutation([i+1 for i in P]).inverse())
        U.permute_columns(Permutation([i+1 for i in Q]).inverse())
        correct = correct and L*U == A
        i += 1
    if correct:
        print("square: ok")
    else:
        print("square: wrong")
        return A

    return True
