using Finch
using BenchmarkTools

A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1))
B = Fiber!(Dense(Dense(Element(0.0))), rand(10000, 10000))
C = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))

# C = AB (A symmetric)
# TIME: ~90ms
eval(@finch_kernel function elementwise_naive(C, A, B)
           for j=_, i=_
               C[i, j] = A[i, j] * B[i, j]
           end
       end)

@btime elementwise_naive($C, $A, $B)