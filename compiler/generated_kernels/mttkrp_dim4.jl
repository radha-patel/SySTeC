using Finch
using BenchmarkTools

n = 10
triA = fsprand(Int, n, n, n, n, 1.0)
symA = [triA[sort([i, j, k, l])...] for i = 1:n, j = 1:n, k = 1:n, l = 1:n]
nondiagA = zeros(Int, n, n, n, n)
diagA = zeros(Int, n, n, n, n)
b = rand(Int, n, n)

for l=1:n, k=1:n, j=1:n, i=1:n
    if i != j && j != k && k != l && i != k && i != l && j != l
        nondiagA[i, j, k, l] = symA[i, j, k, l]
    end
    if i == j || j == k || k == l || i == k || i == l || j == l
        diagA[i, j, k, l] = symA[i, j, k, l]
    end
end

A = Tensor(Dense(SparseList(SparseList(SparseList(Element(0))))), symA)
A_nondiag = Tensor(Dense(SparseList(SparseList(SparseList(Element(0))))), nondiagA)
A_diag = Tensor(Dense(SparseList(SparseList(SparseList(Element(0))))), diagA)
B = Tensor(Dense(Dense(Element(0))), b)   
B_T = Tensor(Dense(Dense(Element(0))), transpose(b)) 
C_T = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
C_nondiag = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
C_diag = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))

eval(@finch_kernel mode=:fast function mttkrp_ref(C_T, A, B_T)
    C_T .= 0
    for m=_, l=_, k=_, i=_, j=_
        C_T[j, i] += A[i, k, l, m] * B_T[j, l] * B_T[j, k] * B_T[j, m]
    end
    return C_T
end)

eval(@finch_kernel mode=:fast function mttkrp_opt_1(C_T, A_nondiag, B_T)
    C_T .= 0
    for m=_, l=_, k=_, i=_, j=_
        if i < k && k < l && l < m
            let A_iklm = A_nondiag[i, k, l, m], B_T_jl = B_T[j, l], B_T_jk = B_T[j, k], B_T_ji = B_T[j, i], B_T_jm = B_T[j, m]
                C_T[j, m] += 6 * B_T_jl * B_T_jk * B_T_ji * A_iklm
                C_T[j, l] += 6 * B_T_jk * B_T_ji * A_iklm * B_T_jm
                C_T[j, k] += 6 * B_T_jl * B_T_ji * A_iklm * B_T_jm
                C_T[j, i] += 6 * B_T_jl * B_T_jk * A_iklm * B_T_jm 
            end
        end
    end
    return C_T
end)

eval(@finch_kernel mode=:fast function mttkrp_opt_2(C_T, A_diag, B_T)
    C_T .= 0
    for m=_, l=_, k=_, i=_, j=_
        if identity(i) <= identity(k) && identity(k) <= identity(l) && identity(l) <= identity(m)
            let ik_eq = (i == k), kl_eq = (k == l), lm_eq = (l == m)
                let A_iklm = A_diag[i, k, l, m], B_T_jl = B_T[j, l], B_T_jk = B_T[j, k], B_T_ji = B_T[j, i], B_T_jm = B_T[j, m]
                    if (ik_eq && !kl_eq && !lm_eq) || (!ik_eq && kl_eq && !lm_eq) || (!ik_eq && !kl_eq && lm_eq)
                        C_T[j, m] += 3 * B_T_jl * B_T_jk * B_T_ji * A_iklm
                        C_T[j, l] += 3 * B_T_jk * B_T_ji * A_iklm * B_T_jm
                        C_T[j, k] += 3 * B_T_jl * B_T_ji * A_iklm * B_T_jm
                        C_T[j, i] += 3 * B_T_jl * B_T_jk * A_iklm * B_T_jm 
                    end
                    if (ik_eq && !kl_eq && lm_eq)
                        C_T[j, m] += 3 * B_T_jl * B_T_jk * B_T_ji * A_iklm
                        C_T[j, k] += 3 * B_T_jl * B_T_ji * A_iklm * B_T_jm
                    end
                    if (ik_eq && kl_eq && !lm_eq) || (!ik_eq && kl_eq && lm_eq)
                        C_T[j, m] += B_T_jl * B_T_jk * B_T_ji * A_iklm
                        C_T[j, l] += B_T_jk * B_T_ji * A_iklm * B_T_jm
                        C_T[j, k] += B_T_jl * B_T_ji * A_iklm * B_T_jm
                        C_T[j, i] += B_T_jl * B_T_jk * A_iklm * B_T_jm 
                    end
                    if ik_eq && kl_eq && lm_eq
                        C_T[j, m] += B_T_jl * B_T_jk * B_T_ji * A_iklm
                    end
                end
            end
        end
    end
    return C_T
end)

function main()
    ref = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
    @btime(mttkrp_ref($ref, $A, $B_T))

    @btime(mttkrp_opt_1($C_nondiag, $A_nondiag, $B_T))
    @btime(mttkrp_opt_2($C_diag, $A_diag, $B_T))
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= C_diag[i, j] + C_nondiag[i, j] == ref[i, j]
    end
    @info "check" check[]
end

main()

