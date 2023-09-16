using Finch
using BenchmarkTools

A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1))
B = Fiber!(Dense(Dense(Element(0.0))), rand(10000, 10000))
C = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))

# C = AB (A symmetric)
# Optimization: 1x flops, 2x memory bandwidth
# TIME: 205.812 ms
eval(@finch_kernel function elementwise(C, A, B)
           for j=_, i=_
               if uptrimask[i, j]
                   temp = A[i, j]
                   C[i, j] = temp * B[i, j]
                   C[j, i] = temp * B[j, i]
               end
           end
       end)

@btime elementwise($C, $A, $B)