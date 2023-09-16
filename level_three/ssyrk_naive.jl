using Finch
using BenchmarkTools

A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1))
C = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))

# C = AA^T (A symmetric)
# TIME: ~80s
eval(@finch_kernel function ssyrk_naive(C, A) 
    for k = _, j = _, i = _
        C[i, j] += A[i, k] * A[j, k]
    end
end)

@btime ssyrk_naive($C, $A)