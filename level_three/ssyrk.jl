using Finch
using BenchmarkTools

# A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1))
# C = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))

# C = AA^T (A symmetric)
# Optimization: 1x flops, 2x memory bandwidth
# TIME: ~90s
# eval(@finch_kernel function ssyrk(C, A) 
#     for l = _, j = _
#         temp = A[j, l]
#         for i = _
#             C[i, j] += temp * A[i, l]
#         end
#     end
# end)

# C = AA^T (A symmetric)
# Optimization: 2x flops, 2x memory bandwidth
# TIME: ~55s
A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1))
C = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))
w = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))
# w = Fiber!(SparseHash{2}(Element(0.0)))
eval(@finch_kernel function ssyrk(C, A, w) 
    for l = _, j = _
        temp1 = A[j, l]
        for i = _
            if uptrimask[i, j]
                C[i, j] += temp1 * A[i, l]
            end 
        end
    end
    # copyto!(C^T, swizzle(C, 2, 1))
    for j = _, i = _ 
        w[j, i] = C[i, j] 
    end
    for i = _, j = _
        if uptrimask[i+1, j]
            C[j, i] = w[j, i]
        end
    end
end)

@btime ssyrk($C, $A, $w)