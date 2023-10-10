using Finch
using BenchmarkTools

# C = AA^T (A symmetric)
# Optimization: 2x flops, 2x memory bandwidth
# TIME: ~40s
A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1))
C = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))
w = Fiber!(Dense(Dense(Element(0.0))), zeros(10000, 10000))
eval(@finch_kernel mode=fastfinch function ssyrk_dense(C, A, w) 
    for l = _, j = _
        let temp1 = A[j, l]
            for i = _
                if uptrimask[i, j]
                    C[i, j] += temp1 * A[i, l]
                end 
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

@btime ssyrk_dense($C, $A, $w)

# C = AA^T (A symmetric)
# Optimization: 2x flops, 2x memory bandwidth
# TIME: -- (too long)
A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10, 10), 0.1))
C = Fiber!(Dense(SparseList(Element(0.0))), zeros(10, 10))
c = Fiber!(SparseHash{2}(Element(0.0)), zeros(10, 10))
w = Fiber!(SparseHash{2}(Element(0.0)), zeros(10, 10))
eval(@finch_kernel mode=fastfinch function ssyrk_sparse(C, A, w, c) 
    c .= 0
    for l = _, j = _
        let temp1 = A[j, l]
            for i = _
                if uptrimask[i, j]
                    c[i, j] += temp1 * A[i, l]
                end 
            end
        end
    end
    w .= 0
    for j = _, i = _ 
        w[j, i] = c[i, j] 
    end
    C .= 0
    for j = _, i = _
        C[i, j] = ifelse(i <= j, c[i, j], w[i, j])
    end
end)

@btime ssyrk_sparse($C, $A, $w, $c)