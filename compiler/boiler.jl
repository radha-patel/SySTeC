# intitialize symbols
A = :A
B = :B
C = :C
i = index(:i)
j = index(:j)
k = index(:k)
l = index(:l)

# expressions to test with
ex = @finch_program C[i] += A[i, j] * B[j]
symmetrize(ex, [A])
ex = @finch_program C[] += B[i] * A[i, j] * B[j]
symmetrize(ex, [A])
ex = @finch_program C[i, j] += A[i, k] * B[k, j]
symmetrize(ex, [A])
ex = @finch_program C[i, j] += A[i, k] * A[k, j]
symmetrize(ex, [A])
ex = @finch_program C[i, j, l] += A[k, j, l] * B[k, i]
symmetrize(ex, [A])
ex = @finch_program C[i, j] += A[i, k, l] * B[l, j] * B[k, j]
symmetrize(ex, [A])
ex = @finch_program C[i, j, l] += A[i, k, j, l] * B[k, i]
symmetrize(ex, [A])