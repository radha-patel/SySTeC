using Finch
using BenchmarkTools


n = 100
A = Fiber!(Dense(Dense(SparseList(Element(0.0)))), fsprand((n, n, n), 0.1))
C = Fiber!(Dense(Dense(Dense(Element(0.0)))), zeros((n, n, n)))
X = Fiber!(Dense(Dense(Element(0.0))), rand(n, n))
# TIME: ~35ms
eval(@finch_kernel mode=fastfinch function n_mode_product_naive(C, A, X)
    for k=_, l=_, j=_, i=_ 
        C[i, j, l] += A[i, k, l] * X[k, j] 
    end
end)

@btime n_mode_product_naive($C, $A, $X)