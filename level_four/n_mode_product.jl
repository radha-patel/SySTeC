using Finch
using BenchmarkTools
# using SymmetricTensors
# using Random
# M = rand(SymmetricTensor{Float64, 3}, n)


n = 100
A = Fiber!(Dense(Dense(SparseList(Element(0.0)))), fsprand((n, n, n), 0.1))
C = Fiber!(Dense(Dense(Dense(Element(0.0)))), zeros((n, n, n)))
X = Fiber!(Dense(Dense(Element(0.0))), rand(n, n))
temp2 = Scalar(0.0)
# TIME: ~12ms
eval(@finch_kernel mode=fastfinch function n_mode_product(C, A, X, temp2) 
    for l=_, j=_, i=_
        temp2 .= 0
        let temp1 = X[i, j]
            for k=_ 
                let temp3 = A[k, i, l]
                    if uptrimask[k+1, i] 
                        C[k, j, l] += temp1 * temp3
                    end
                    if uptrimask[k, i]
                        temp2[] += X[k, j] * temp3
                    end
                end
            end
        end
        C[i, j, l] += temp2[]
    end
end)

@btime n_mode_product($C, $A, $X, $temp2)