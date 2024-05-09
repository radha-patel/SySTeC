using Finch
using BenchmarkTools

n = 1000
triA = fsprand(Int, n, n, 0.1)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
nondiagA = zeros(Int, n, n)
diagA = zeros(Int, n, n)

for j=1:n, i=1:n
    if i != j
        nondiagA[i, j] = symA[i, j]
    end
    if i == j
        diagA[i, j] = symA[i, j]
    end
end

A = Tensor(Dense(SparseList(Element(0))), symA)
A_T = Tensor(Dense(SparseList(Element(0))), transpose(symA))
A_dense = Tensor(Dense(Dense(Element(0))), symA)
A_diag = Tensor(Dense(SparseList(Element(0))), diagA)
A_nondiag = Tensor(Dense(SparseList(Element(0))), nondiagA)
C = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
ref = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
temp = Scalar(0)

# ~12ms
eval(@finch_kernel mode=:fast function syrk_ref(C, A)
    C .= 0
    for k=_, j=_, i=_
        C[i, j] += A[i, k] * A[j, k]
    end
end)

# ~290ms
# println("before eval gen")
# eval(@finch_kernel mode=:fast function syrk_gen(C, A)
#     C .= 0
#     for j=_, k=_, i=_
#         let jk_leq = (j <= k), ij_leq = (i <= j)
#             let A_ij = A[i, j], A_ik = A[i, k], A_jk = A[j, k]
#                 if ij_leq && jk_leq
#                     C[i, j] += A_jk * A_ik
#                 end
#                 if ij_leq && j < k
#                     C[i, k] += A_jk * A_ij
#                 end
#                 if i < j && jk_leq 
#                     C[j, k] += A_ik * A_ij
#                 end
#             end
#         end
#     end
# end)
# println("after eval gen")

# println("before eval opt1")
# # ~1.6s
# eval(@finch_kernel mode=:fast function syrk_opt1(C, A)
#     C .= 0
#     for k=_, j=_, i=_
#         let jk_leq = (identity(j) <= identity(k)), ij_leq = (identity(i) <= identity(j))
#             let A_ij = A[i, j], A_ik = A[i, k], A_jk = A[j, k]
#                 if ij_leq && jk_leq
#                     C[i, j] += A_jk * A_ik
#                 end
#                 if ij_leq && identity(j) < identity(k)
#                     C[i, k] += A_jk * A_ij
#                 end
#                 if identity(i) < identity(j) && jk_leq 
#                     C[j, k] += A_ik * A_ij
#                 end
#             end
#         end
#     end
# end)
# println("after eval opt1")

# println("before eval opt2")
# # ~190ms
# eval(@finch_kernel mode=:fast function syrk_opt2(C, A_nondiag)
#     C .= 0
#     for k=_, j=_, i=_
#         let A_jk = A_nondiag[j, k], A_ik = A_nondiag[i, k], A_ij = A_nondiag[i, j]
#             if i < j && j < k
#                 C[i, j] += A_jk * A_ik
#                 C[i, k] += A_jk * A_ij
#                 C[j, k] += A_ik * A_ij
#             end
#         end
#     end
# end)
# println("after eval opt2")

# println("before eval opt2_1")
# # ~2.4s
# eval(@finch_kernel mode=:fast function syrk_opt2_1(C, A_diag)
#     for k=_, j=_, i=_
#         let A_jk = A_diag[j, k], A_ik = A_diag[i, k], A_ij = A_diag[i, j]
#             let ij_neq = (identity(i) != identity(j)), jk_neq = (identity(j) != identity(k)), ij_eq = (identity(i) == identity(j)), jk_eq = (identity(j) == identity(k))
#                 if identity(i) <= identity(j) && identity(j) <= identity(k)
#                     if ij_eq && jk_neq
#                         C[i, k] += A_jk * A_ij
#                         C[i, j] += A_jk * A_ik
#                     end
#                     if ij_neq && jk_eq
#                         C[i, j] += A_jk * A_ik
#                         C[j, k] += A_ik * A_ij
#                     end
#                     if ij_eq && jk_eq
#                         C[i, j] += A_jk * A_ik
#                     end
#                 end
#             end
#         end
#     end
# end)
# println("after eval opt2_1")

# # ~7.0ms
println("before eval opt3")
eval(@finch_kernel mode=:fast function syrk_opt3(C, A)
    C .= 0
    for k=_, j=_, i=_
        if i <= j
            C[i, j] += A[i, k] * A[j, k]
        end
    end
end)
println("after eval opt3")

# ~750ms
# println("before eval opt3.5")
# eval(@finch_kernel mode=:fast function syrk_opt3_5(C, A, A_T)
#     C .= 0
#     for j=_, i=_, k=_
#         if i <= j
#             C[i, j] += A_T[k, i] * A[k, j]
#         end
#     end
# end)
# println("after eval opt3.5")

# # ~280ms
# println("after eval opt4")
# eval(@finch_kernel mode=:fast function syrk_opt4(C, A_dense)
#     C .= 0
#     for j=_, k=_, i=_
#         let jk_leq = (j <= k), ij_leq = (i <= j)
#             let A_ij = A_dense[i, j], A_ik = A_dense[i, k], A_jk = A_dense[j, k]
#                 if ij_leq && jk_leq
#                     C[i, j] += A_jk * A_ik
#                 end
#                 if ij_leq && j < k
#                     C[i, k] += A_jk * A_ij
#                 end
#                 if i < j && jk_leq 
#                     C[j, k] += A_ik * A_ij
#                 end
#             end
#         end
#     end
# end)
# println("after eval opt4")

# # ~200ms
# println("after eval opt4_1")
# eval(@finch_kernel mode=:fast function syrk_opt4_1(C, A, A_dense)
#     C .= 0
#     for j=_, k=_, i=_
#         let jk_leq = (j <= k), ij_leq = (i <= j)
#             let A_ij = A[i, j], A_ik = A_dense[i, k], A_jk = A_dense[j, k]
#                 if ij_leq && jk_leq
#                     C[i, j] += A_jk * A_ik
#                 end
#                 if ij_leq && j < k
#                     C[i, k] += A_jk * A_ij
#                 end
#                 if i < j && jk_leq 
#                     C[j, k] += A_ik * A_ij
#                 end
#             end
#         end
#     end
# end)
# println("after eval opt4_1")

function main()
    # ref = Tensor(Dense(Dense(Element(0))), zeros(n, n))
    @btime(syrk_ref($ref, $A))

    # @btime(syrk_gen($C, $A))
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # @info "check generated code" check[]

    # @btime(syrk_opt1($C, $A))
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # Sometimes doesn't match ref because of overflow errors?
    # @info "check optimized code (1)" check[]

    # @btime(syrk_opt2($C, $A_nondiag))
    # @btime(syrk_opt2_1($C, $A_diag))
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # Sometimes doesn't match ref because of overflow errors?
    # @info "check optimized code (2)" check[]

    @btime(syrk_opt3($C, $A))
    check = Scalar(true)
    @finch for j=_, i=_
        if i <= j
            check[] &= C[i, j] == ref[i, j]
        end
    end
    @info "check optimized code (3)" check[]

    # @btime(syrk_opt3_5($C, $A, $A_T))
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # @info "check optimized code (3.5)" check[]

    # @btime(syrk_opt4($C, $A_dense))
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # @info "check optimized code (4)" check[]

    # @btime(syrk_opt4_1($C, $A, $A_dense))
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # @info "check optimized code (4.1)" check[]
end

main()