using Finch
using BenchmarkTools

n = 1000
triA = fsprand(Int, n, n, 0.1)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
nondiagA = zeros(Int, n, n)
diagA = zeros(Int, n, n)
triB = fsprand(Int, n, n, 0.1)
symB = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
nondiagB = zeros(Int, n, n)
diagB = zeros(Int, n, n)

for j=1:n, i=1:n
    if i != j
        nondiagA[i, j] = symA[i, j]
        nondiagB[i, j] = symB[i, j]
    end
    if i == j
        diagA[i, j] = symA[i, j]
        diagB[i, j] = symB[i, j]
    end
end

A = Tensor(Dense(SparseList(Element(0))), symA)
A_diag = Tensor(Dense(SparseList(Element(0))), diagA)
A_nondiag = Tensor(Dense(SparseList(Element(0))), nondiagA)
B = Tensor(Dense(SparseList(Element(0))), triB)
B_diag = Tensor(Dense(SparseList(Element(0))), diagB)
B_nondiag = Tensor(Dense(SparseList(Element(0))), nondiagB)
C = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))

eval(@finch_kernel mode=fastfinch function symsym_ref(C, A, B)
    # C .= 0
    for j=_, k=_, i=_
        C[i, j] += A[i, k] * B[k, j]
    end
end)

# DID NOT COMPILE
# println("before eval symsym_gen")
# eval(@finch_kernel mode=fastfinch function symsym_gen(C, A, B)
#     C .= 0
#     for k=_, j=_, i=_
#         let jk_lt = (j < k), ij_lt = (i < j), jk_leq = (j <= k), ij_leq = (i <= j)
#             let A_ik = A[i, k], A_ij = A[i, j], B_ij = B[i, j], B_ik = B[i, k], A_jk = A[j, k], B_jk = B[j, k]
#                 if ij_leq && jk_leq
#                     C[i, j] += B_jk * A_ik 
#                 end
#                 if ij_leq && jk_lt
#                     C[k, i] += A_jk * B_ij
#                     C[i, k] += B_jk * A_ij
#                 end
#                 if ij_lt && jk_lt
#                     C[k, j] += A_ik * B_ij
#                 end
#                 if ij_lt && jk_leq
#                     C[j, k] += A_ij * B_ik
#                     C[j, i] += A_jk * B_ik
#                 end
#             end
#         end
#     end
# end)
# println("after eval symsym_gen")

println("before eval symsym_opt1_1")
eval(@finch_kernel mode=fastfinch function symsym_opt1_1(C, A, B)
    C .= 0
    for k=_, j=_, i=_
        if i < k && k < j
            let A_ik = A[i, k], A_ij = A[i, j], B_ij = B[i, j], B_ik = B[i, k], A_jk = A[j, k], B_jk = B[j, k]
                C[i, j] += A_ik * B_jk
                C[k, j] += A_ik * B_ij
                C[i, k] += A_ij * B_jk
                C[k, i] += A_jk * B_ij
                C[j, k] += A_ij * B_ik
                C[j, i] += A_jk * B_ik
            end
        end
    end
end)
println("after eval symsym_opt1_1")

println("before eval symsym_opt1_2")
eval(@finch_kernel mode=fastfinch function symsym_opt1_2(C, A, B)
    for k=_, j=_, i=_
        if i <= k && k <= j
            let ik_eq = (identity(i) == identity(k)) && kj_eq = (identity(k) == identity(j))
                let A_ik = A[i, k], A_ij = A[i, j], B_ij = B[i, j], B_ik = B[i, k], A_jk = A[j, k], B_jk = B[j, k]
                    if ik_eq && !kj_eq
                        C[i, j] += A_ik * B_jk
                        C[i, k] += A_ij * B_jk
                        C[j, k] += A_ij * B_ik
                    end
                    if !ik_eq && kj_eq
                        C[i, j] += A_ik * B_jk
                        C[k, j] += A_ik * B_ij
                        C[k, i] += A_jk * B_ij
                    end
                    if ik_eq && kj_eq
                        C[i, j] += A_ik * B_jk
                    end
                end
            end
        end
    end
end)
println("after eval symsym_opt1_2")

function main()
    ref = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
    @btime(symsym_ref($ref, $A, $B))

    # @btime(symsym_gen($C, $A, $B))
    # @info C == ref

    @btime(symsym_opt1_1($C, $A_nondiag, $B_nondiag))
    @btime(symsym_opt1_2($C, $A_diag, $B_nondiag))
    @btime(symsym_opt1_2($C, $A_nondiag, $B_diag))
    @btime(symsym_opt1_2($C, $A_diag, $B_diag))
    @info C == ref
end

main()