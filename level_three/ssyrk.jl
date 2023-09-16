using Finch
using BenchmarkTools

A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1))
C = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))

# C = AA^T (A symmetric)
# Optimization: 1x flops, 2x memory bandwidth
# TIME: ~90s
eval(@finch_kernel function ssyrk(C, A) 
    for l = _, j = _
        temp = A[j, l]
        for i = _
            C[i, j] += temp * A[i, l]
        end
    end
end)

# C = AA^T (A symmetric)
# Optimization: 2x flops, 2x memory bandwidth
# TIME: ~100s
# eval(@finch_kernel function ssyrk(C, A) 
#     for l = _, j = _
#         temp1 = A[j, l]
#         for i = _
#             if i > j
#                 temp2 = temp1 * A[i, l]
#                 C[i, j] += temp2
#                 C[j, i] += temp2
#             end 
#             if i == j
#                 C[i, j] += temp1 * A[i, l]
#             end
#         end
#     end
# end)

@btime ssyrk($C, $A)