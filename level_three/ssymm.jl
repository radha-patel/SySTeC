using Finch
using BenchmarkTools

A = Tensor(Dense(SparseList(Element(0.0))), fsprand(10000, 10000, 0.1))
B = Tensor(Dense(Dense(Element(0.0))), rand(10000, 10000))
C = Tensor(Dense(Dense(Element(0.0))), zeros(10000, 10000))
temp2 = Scalar(0.0)

# C = AB (A symmetric)
# Optimization: 1x flops, 2x memory bandwidth
# TIME: ~165s
eval(@finch_kernel function ssymm(C, A, B, temp2)
    for j=_, i =_
        temp2 .= 0
        let temp1 = B[i, j]
            for k =_
                let temp3 = A[k, i]
                    if uptrimask[k+1, i]
                        C[k, j] += temp1 * temp3
                    end
                    if uptrimask[k, i]
                        temp2[] += B[k, j] * temp3
                    end
                end
            end
        end
        C[i, j] += temp2[]
    end
end)

@btime ssymm($C, $A, $B, $temp2)