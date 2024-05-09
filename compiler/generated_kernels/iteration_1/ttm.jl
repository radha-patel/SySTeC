using Finch
using BenchmarkTools

# n = 100
# triA = fsprand(Int, n, n, n, 0.1)
# symA = [triA[sort([i, j, k])...] for i = 1:n, j = 1:n, k = 1:n]
# b = rand(Int, n, n)

# A = Tensor(Dense(SparseList(SparseList(Element(0)))), symA)
# B = Tensor(Dense(Dense(Element(0))), b)    
# B_T = Tensor(Dense(Dense(Element(0))), transpose(b)) 
# C = Tensor(Dense(Dense(Dense(Element(0)))), zeros(Int, n, n, n))

n = 100
triA = fsprand(Int, n, n, n, 0.1)
symA = [triA[sort([i, j, k])...] for i = 1:n, j = 1:n, k = 1:n]
nondiagA = zeros(Int, n, n, n)
diagA = zeros(Int, n, n, n)
b = rand(Int, n, n)

for k=1:n, j=1:n, i=1:n
    if i != j && j != k && i != k
        nondiagA[i, j, k] = symA[i, j, k]
    end
    if i == j || j == k || i == k
        diagA[i, j, k] = symA[i, j, k]
    end
end

A = Tensor(Dense(Dense(Dense(Element(0)))), symA)
A_nondiag = Tensor(Dense(SparseList(SparseList(Element(0)))), nondiagA)
A_diag = Tensor(Dense(SparseList(SparseList(Element(0)))), diagA)
B = Tensor(Dense(Dense(Element(0))), b)    
B_T = Tensor(Dense(Dense(Element(0))), transpose(b)) 
C = Tensor(Dense(Dense(Dense(Element(0)))), zeros(Int, n, n, n))
C_nondiag = Tensor(Dense(Dense(Dense(Element(0)))), zeros(Int, n, n, n))
C_diag = Tensor(Dense(Dense(Dense(Element(0)))), zeros(Int, n, n, n))

# ~7.9ms
eval(@finch_kernel mode=:fast function ttm_ref(C, A, B_T)
    C .= 0
    for l=_, j=_, k=_, i=_
        C[i, j, l] += A[k, j, l] * B_T[i, k]
    end
end)

# ~3.6ms
# eval(@finch_kernel mode=:fast function ttm_gen(C, A, B)
#     C .= 0
#     for l=_, k=_, j=_, i=_
#         let jk_leq = (j <= k), kl_leq = (k <= l)
#             let A_jkl = A[j, k, l]
#                 if jk_leq && kl_leq
#                     C[i, j, k] += A_jkl * B[l, i]
#                 end
#                 if j < k && kl_leq
#                     C[i, k, l] += A_jkl * B[j, i]
#                 end
#                 if jk_leq && k < l
#                     C[i, j, l] += A_jkl * B[k, i]
#                 end
#             end
#         end
#     end
# end)

# ~3.0ms
# eval(@finch_kernel mode=:fast function ttm_opt(C, A, B_T)
#     C .= 0
#     for l=_, k=_, j=_, i=_
#         let jk_leq = (j <= k), kl_leq = (k <= l)
#             let A_jkl = A[j, k, l]
#                 if jk_leq && kl_leq
#                     C[i, j, k] += A_jkl * B_T[i, l]
#                 end
#                 if j < k && kl_leq
#                     C[i, k, l] += A_jkl * B_T[i, j]
#                 end
#                 if jk_leq && k < l
#                     C[i, j, l] += A_jkl * B_T[i, k]
#                 end
#             end
#         end
#     end
# end)

eval(@finch_kernel mode=:fast function ttm_opt_1(C, A_nondiag, B_T)
    C .= 0
    for l=_, k=_, j=_, i=_
        let A_jkl = A_nondiag[j, k, l], B_ik = B_T[i, k], B_il = B_T[i, l], B_ij = B_T[i, j]
            if j < k && k < l
                C[i, j, l] += A_jkl * B_ik
                C[i, j, k] += B_il * A_jkl
                C[i, k, l] += B_ij * A_jkl
            end
        end
    end
end)

eval(@finch_kernel mode=:fast function ttm_opt_2(C, A_diag, B_T)
    C .= 0
    for l=_, k=_, j=_, i=_
        if j <= k && k <= l
            let jk_eq = (identity(j) == identity(k)), kl_eq = (identity(k) == identity(l))
                let A_jkl = A_diag[j, k, l], B_ik = B_T[i, k], B_il = B_T[i, l], B_ij = B_T[i, j]
                    if (jk_eq && !kl_eq ) || (!jk_eq  && kl_eq)
                        C[i, j, k] += B_il * A_jkl
                        C[i, l, j] += A_jkl * B_ik
                        C[i, k, l] += B_ij * A_jkl
                    end
                    if jk_eq && kl_eq
                        C[i, l, j] += A_jkl * B_ik
                    end
                end
            end
        end
    end
end)

function main()
    ref = Tensor(Dense(Dense(Dense(Element(0)))), zeros(Int, n, n, n))
    @btime(ttm_ref($ref, $A, $B_T))

    # @btime(ttm_gen($C, $A, $B))
    # check = Scalar(true)
    # @finch for l=_, j=_, i=_
    #     if j <= l
    #         check[] &= C[i, j, l] == ref[i, j, l]
    #     end
    # end
    # @info "check generated code" check[]


    # @btime(ttm_opt($C, $A, $B_T))
    # check = Scalar(true)
    # @finch for l=_, j=_, i=_
    #     if j <= l
    #         check[] &= C[i, j, l] == ref[i, j, l]
    #     end
    # end
    # @info "check generated code" check[]

    @btime(ttm_opt_1($C_nondiag, $A_nondiag, $B_T))
    @btime(ttm_opt_2($C_diag, $A_diag, $B_T))
    check = Scalar(true)
    @finch for l=_, j=_, i=_
        if j <= l
            check[] &= (C_nondiag[i, j, l] + C_diag[i, j, l]) == ref[i, j, l]
        end
    end
    @info "check generated code" check[]
end

main()