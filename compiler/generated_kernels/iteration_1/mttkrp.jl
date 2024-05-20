using Finch
using BenchmarkTools

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

A = Tensor(Dense(SparseList(SparseList(Element(0)))), symA)
A_nondiag = Tensor(Dense(SparseList(SparseList(Element(0)))), nondiagA)
A_diag = Tensor(Dense(SparseList(SparseList(Element(0)))), diagA)
B = Tensor(Dense(Dense(Element(0))), b)   
B_T = Tensor(Dense(Dense(Element(0))), transpose(b)) 
C = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
C_nondiag = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
C_diag = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))

# ~11ms
eval(@finch_kernel mode=:fast function mttkrp_ref(C, A, B)
    C .= 0
    for l=_, j=_, k=_, i=_
        C[i, j] += A[i, k, l] * B[l, j] * B[k, j]
    end
end)

# DID NOT COMPILE
# println("before eval mttkrp_gen")
# eval(@finch_kernel mode=:fast function mttkrp_gen(C, A, B)
#     C .= 0
#     for l=_, j=_, k=_, i=_
#         if i <= k && k <= l 
#             let ik_eq = (identity(i) == identity(k)), kl_eq = (identity(k) == identity(l))
#                 let B_ij = B[i, j], A_ikl = A[i, k, l], B_kj = B[k, j], B_lj = B[l, j]
#                     if i != k && k != l
#                         C[l, j] += 2 * B_kj * B_ij * A_ikl
#                         C[k, j] += 2 * B_lj * B_ij * A_ikl
#                         C[i, j] += 2 * B_kj * B_lj * A_ikl
#                     end
#                     if (ik_eq && identity(k) != identity(l)) || (identity(i) != identity(k) && kl_eq)
#                         C[i, j] += B_kj * B_lj * A_ikl
#                         C[k, j] += B_lj * B_ij * A_ikl
#                         C[l, j] += B_kj * B_ij * A_ikl
#                     end
#                     if ik_eq && kl_eq 
#                         C[i, j] += B_kj * B_lj * A_ikl
#                     end
#                 end
#             end
#         end
#     end
# end)
# println("after eval mttkrp_gen")

# println("before eval mttkrp_opt1")
# # ~4.5ms
# eval(@finch_kernel mode=:fast function mttkrp_opt1(C, A, B_T)
#     C .= 0
#     for l=_, k=_, i=_, j=_
#         if i <= k && k <= l 
#             let ik_eq = (identity(i) == identity(k)), kl_eq = (identity(k) == identity(l))
#                 let B_ij = B_T[j, i], A_ikl = A[i, k, l], B_kj = B_T[j, k], B_lj = B_T[j, l]
#                     if !ik_eq && !kl_eq
#                         C[l, j] += 2 * B_kj * B_ij * A_ikl
#                         C[k, j] += 2 * B_lj * B_ij * A_ikl
#                         C[i, j] += 2 * B_kj * B_lj * A_ikl
#                     end
#                     if (ik_eq && !kl_eq) || (!ik_eq && kl_eq)
#                         C[i, j] += B_kj * B_lj * A_ikl
#                         C[k, j] += B_lj * B_ij * A_ikl
#                         C[l, j] += B_kj * B_ij * A_ikl
#                     end
#                     if ik_eq && kl_eq 
#                         C[i, j] += B_kj * B_lj * A_ikl
#                     end
#                 end
#             end
#         end
#     end
# end)
# println("after eval mttkrp_opt1")

println("before eval mttkrp_opt2_1")
# ~3.5ms
eval(@finch_kernel mode=:fast function mttkrp_opt2_1(C, A_nondiag, B_T)
    C .= 0
    for l=_, k=_, i=_, j=_
        if i < k && k < l 
            let B_ij = B_T[j, i], A_ikl = A_nondiag[i, k, l], B_kj = B_T[j, k], B_lj = B_T[j, l]
                C[l, j] += 2 * B_kj * B_ij * A_ikl
                C[k, j] += 2 * B_lj * B_ij * A_ikl
                C[i, j] += 2 * B_kj * B_lj * A_ikl
            end
        end
    end
end)
println("after eval mttkrp_opt2_1")

# println("before eval mttkrp_opt2_2")
# # ~240μs
# eval(@finch_kernel mode=:fast function mttkrp_opt2_2(C, A_diag, B_T)
#     C .= 0
#     for l=_, k=_, i=_, j=_
#         if identity(i) <= identity(k) && identity(k) <= identity(l) 
#             let ik_eq = (identity(i) == identity(k)), kl_eq = (identity(k) == identity(l))
#                 let B_ij = B_T[j, i], A_ikl = A_diag[i, k, l], B_kj = B_T[j, k], B_lj = B_T[j, l]
#                     if (ik_eq && !kl_eq) || (!ik_eq && kl_eq)
#                         C[i, j] += B_kj * B_lj * A_ikl
#                         C[k, j] += B_lj * B_ij * A_ikl
#                         C[l, j] += B_kj * B_ij * A_ikl
#                     end
#                     if ik_eq && kl_eq 
#                         C[i, j] += B_kj * B_lj * A_ikl
#                     end
#                 end
#             end
#         end
#     end
# end)
# println("after eval mttkrp_opt2_2")

println("before eval mttkrp_opt2_3")
# ~220μs
eval(@finch_kernel mode=:fast function mttkrp_opt2_3(C, A_diag, B_T)
    C .= 0
    for l=_, k=_, i=_, j=_
        if identity(i) <= identity(k) && identity(k) <= identity(l) 
            let ik_eq = (i == k), kl_eq = (k == l)
                let B_ij = B_T[j, i], A_ikl = A_diag[i, k, l], B_kj = B_T[j, k], B_lj = B_T[j, l]
                    if (ik_eq && !kl_eq) || (!ik_eq && kl_eq)
                        C[i, j] += B_kj * B_lj * A_ikl
                        C[k, j] += B_lj * B_ij * A_ikl
                        C[l, j] += B_kj * B_ij * A_ikl
                    end
                    if ik_eq && kl_eq 
                        C[i, j] += B_kj * B_lj * A_ikl
                    end
                end
            end
        end
    end
end)
println("after eval mttkrp_opt2_3")

# println("before eval mttkrp_opt2_4")
# # ~220μs
# eval(@finch_kernel mode=:fast function mttkrp_opt2_4(C, A_diag, B_T)
#     C .= 0
#     for l=_, k=_, i=_, j=_
#         if i <= k && k <= l
#             let ik_eq = (identity(i) == identity(k)), kl_eq = (identity(k) == identity(l))
#                 let B_ij = B_T[j, i], A_ikl = A_diag[i, k, l], B_kj = B_T[j, k], B_lj = B_T[j, l]
#                     if (ik_eq && !kl_eq) || (!ik_eq && kl_eq)
#                         C[i, j] += B_kj * B_lj * A_ikl
#                         C[k, j] += B_lj * B_ij * A_ikl
#                         C[l, j] += B_kj * B_ij * A_ikl
#                     end
#                     if ik_eq && kl_eq 
#                         C[i, j] += B_kj * B_lj * A_ikl
#                     end
#                 end
#             end
#         end
#     end
# end)
# println("after eval mttkrp_opt2_4")

function main()
    ref = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
    @btime(mttkrp_ref($ref, $A, $B))

    # @btime(mttkrp_gen($C, $A, $B))
    # @info C == ref

    # @btime(mttkrp_opt1($C, $A, $B_T))
    # @info C == ref

    # @btime(mttkrp_opt2_1($C_nondiag, $A_nondiag, $B_T))
    # @btime(mttkrp_opt2_2($C_diag, $A_diag, $B_T))
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= C_diag[i, j] + C_nondiag[i, j] == ref[i, j]
    # end
    # @info "check" check[]

    @btime(mttkrp_opt2_1($C_nondiag, $A_nondiag, $B_T))
    @btime(mttkrp_opt2_3($C_diag, $A_diag, $B_T))
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= C_diag[i, j] + C_nondiag[i, j] == ref[i, j]
    end
    @info "check" check[]

    # @btime(mttkrp_opt2_1($C_nondiag, $A_nondiag, $B_T))
    # @btime(mttkrp_opt2_4($C_diag, $A_diag, $B_T))
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= C_diag[i, j] + C_nondiag[i, j] == ref[i, j]
    # end
    # @info "check" check[]
end

main()