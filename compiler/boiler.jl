# intitialize symbols
A = :A
B = :B
C = :C
D = :D
i = index(:i)
j = index(:j)
k = index(:k)
l = index(:l)
m = index(:m)
n = index(:n)

y = :y
x = :x

# ssymv
ex = @finch_program C[i] += A[i, j] * B[j]
symmetrize(ex, [A], [i, j])

ex = @finch_program y[i] += A[i, j] * x[j]
symmetrize(ex, [A], [i, j])

# syprd
ex = @finch_program C[] += B[i] * A[i, j] * B[j]
symmetrize(ex, [A], [i, j])

# ssymm
ex = @finch_program C[i, j] += A[i, k] * B[k, j]
symmetrize(ex, [A], [i, k, j])

# ssyrk (A symmetric)
ex = @finch_program C[i, j] += A[i, k] * A[k, j]
symmetrize(ex, [A], [i, k, j])

# ssyrk
ex = @finch_program C[i, j] += A[i, k] * A[j, k]
symmetrize(ex, [], [])


# ttm
ex = @finch_program C[i, j, l] += A[k, j, l] * B[k, i]
symmetrize(ex, [A], [i, j, k, l])

# 3d mttkrp
ex = @finch_program C[i, j] += A[i, k, l] * B[l, j] * B[k, j]
symmetrize(ex, [A], [j, i, k, l])

# 4d mttkrp
ex = @finch_program C[i, j] += A[i, k, l, m] * B[l, j] * B[k, j] * B[m, j]
symmetrize(ex, [A], [j, i, k, l, m])

# 5d mttkrp
ex = @finch_program C[i, j] += A[i, k, l, m, n] * B[l, j] * B[k, j] * B[m, j] * B[n, j]
symmetrize(ex, [A], [j, i, k, l, m, n], true, true)

ex = @finch_program C[i, j] += A[i, k] * B[k, j]
symmetrize(ex, [A, B], [i, k, j])
