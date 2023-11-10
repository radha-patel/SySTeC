using Finch
using BenchmarkTools
using SymmetricTensors
using Random

n = 3
M = rand(SymmetricTensor{Float64, 3}, n)
# A = Fiber!(Dense(Dense(SparseList(Element(0.0)))), fsprand((n, n, n), 0.1))
A = Fiber!(Dense(Dense(Dense(Element(0.0)))), Array(M))
C = Fiber!(Dense(Dense(Dense(Element(0.0)))), zeros((n, 1, n)))
X = Fiber!(Dense(Dense(Element(0.0))), rand(n, 1))
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
println("compiled n_mode_product")

_C = Fiber!(Dense(Dense(Dense(Element(0.0)))), zeros((n, 1, n)))
# eval(@finch_kernel mode=fastfinch function n_mode_product2(C, A, X, temp2) 
#     for l=_, j=_, i=_
#         temp2 .= 0
#         let temp1 = X[i, j]
#             for k=_ 
#                 if i >= k && k >= l
#                     let temp3 = A[k, i, l]
#                         if uptrimask[k+1, i] 
#                             C[k, j, l] += temp1 * temp3
#                         end
#                         if uptrimask[k, i]
#                             temp2[] += X[k, j] * temp3
#                         end
#                         if k != l
#                             C[l, j, k] += X[i, j] * temp3
#                         end
#                         if i != l
#                             C[k, j, i] += X[l, j] * temp3 
#                         end
#                         if k != i && k != l && i != l
#                             C[i, j, k] += X[l, j] * temp3
#                             C[l, j, i] += X[k, j] * temp3 
#                         end
#                     end
#                 end
#             end
#         end
#         C[i, j, l] += temp2[]
#     end
# end)
# println("compiled n_mode_product 2")

# C[i, j, l] = A[i, k, l] * X[k, j]
# C[l, j, i] = A[l, k, i] * X[k, j]
# C[i, j, l] = C[l, j, i]
temp4 = Fiber!(Dense(Element(0.0)), zeros(n))
eval(@finch_kernel mode=fastfinch function n_mode_product3(C, A, X, temp2, temp4) 
    for j=_, i=_
        temp4 .= 0
        for l=_
            if i <= l
                temp2 .= 0
                let temp1 = X[i, j]
                    for k=_ 
                        let temp3 = A[k, i, l]
                            if k <= i
                                C[k, j, l] += temp1 * temp3
                            end
                            if k < i
                                temp2[] += X[k, j] * temp3
                            end
                        end
                    end
                end
                C[i, j, l] += temp2[]
            end
            if i < l
                for k=_
                    if k <= i
                        temp4[k] += A[k, i, l] * X[l, j]
                    end
                end
            end
        end
        for k=_
            if k <= i
                C[k, j, i] += temp4[k]
            end
        end
    end
end)

n_mode_product(C, A, X, temp2)
println("run n_mode_product")
n_mode_product3(_C, A, X, temp2, temp4)
println("run n_mode_product 2")
C == _C