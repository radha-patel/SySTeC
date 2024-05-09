using Finch
using BenchmarkTools

n = 10
triA = fsprand(Int, n, n, n, n, n, 1.0)
symA = [triA[sort([i, j, k, l, m])...] for i = 1:n, j = 1:n, k = 1:n, l = 1:n, m = 1:n]
nondiagA = zeros(Int, n, n, n, n, n)
diagA = zeros(Int, n, n, n, n, n)
b = rand(Int, n, n)

for m=1:n, l=1:n, k=1:n, j=1:n, i=1:n
    if i != j && j != k && k != l && l != m && i != k && i != l && i != m && j != l && j != m && k != m
        nondiagA[i, j, k, l, m] = symA[i, j, k, l, m]
    end
    if i == j || j == k || k == l || l == m || i == k || i == l || i == m || j == l || j == m || k == m
        diagA[i, j, k, l, m] = symA[i, j, k, l, m]
    end
end

A = Tensor(Dense(SparseList(SparseList(SparseList(SparseList(Element(0)))))), symA)
A_nondiag = Tensor(Dense(SparseList(SparseList(SparseList(SparseList(Element(0)))))), nondiagA)
A_diag = Tensor(Dense(SparseList(SparseList(SparseList(SparseList(Element(0)))))), diagA)
B = Tensor(Dense(Dense(Element(0))), b)   
B_T = Tensor(Dense(Dense(Element(0))), transpose(b)) 
C_T = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
C_nondiag = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
C_diag = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))

eval(@finch_kernel mode=:fast function mttkrp_ref(C_T, A, B_T)
    C_T .= 0
    for n=_, m=_, l=_, k=_, i=_, j=_
        C_T[j, i] += A[i, k, l, m, n] * B_T[j, l] * B_T[j, k] * B_T[j, m] * B_T[j, n]
    end
    return C_T
end)

eval(@finch_kernel mode=:fast function mttkrp_opt_1(C_T, A_nondiag, B_T)
    C_T .= 0
    for n=_, m=_, l=_, k=_, i=_, j=_
        if i < k && k < l && l < m && m < n
            let A_iklmn = A_nondiag[i, k, l, m, n], B_T_jl = B_T[j, l], B_T_jk = B_T[j, k], B_T_ji = B_T[j, i], B_T_jm = B_T[j, m], B_T_jn = B_T[j, n]
                C_T[j, i] += 24 * B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                C_T[j, l] += 24 * A_iklmn * B_T_jk * B_T_ji * B_T_jm * B_T_jn
                C_T[j, k] += 24 * B_T_jl * A_iklmn * B_T_ji * B_T_jm * B_T_jn
                C_T[j, n] += 24 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jm
                C_T[j, m] += 24 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jn
            end
        end
    end
    return C_T
end)

eval(@finch_kernel mode=:fast function mttkrp_opt_2(C_T, A_nondiag, B_T)
    C_T .= 0
    for m=_, l=_, k=_, i=_, j=_
        if identity(i) <= identity(k) && identity(k) <= identity(l) && identity(l) <= identity(m) && identity(m) <= identity(n)
            let ik_eq = (i == k), mn_eq = (m == n), kl_eq = (k == l), lm_eq = (l == m)
                let A_iklmn = A_nondiag[i, k, l, m, n], B_T_jl = B_T[j, l], B_T_jk = B_T[j, k], B_T_ji = B_T[j, i], B_T_jm = B_T[j, m], B_T_jn = B_T[j, n]
                    if (ik_eq && !kl_eq && !lm_eq && !mn_eq) || (!ik_eq && kl_eq && !lm_eq && !mn_eq) || (!ik_eq && !kl_eq && lm_eq && !mn_eq) || (!ik_eq && !kl_eq && !lm_eq && mn_eq)
                        C_T[j, i] += 12 * B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                        C_T[j, l] += 12 * A_iklmn * B_T_jk * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, k] += 12 * B_T_jl * A_iklmn * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, n] += 12 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jm
                        C_T[j, m] += 12 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jn
                    end
                    if (ik_eq && !kl_eq && lm_eq && !mn_eq) || (!ik_eq && kl_eq && !lm_eq && mn_eq) || (ik_eq && !kl_eq && !lm_eq && mn_eq)
                        C_T[j, i] += 6 * B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                        C_T[j, l] += 6 * A_iklmn * B_T_jk * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, k] += 6 * B_T_jl * A_iklmn * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, n] += 6 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jm
                        C_T[j, m] += 6 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jn
                    end
                    if (ik_eq && kl_eq && !lm_eq && !mn_eq) || (!ik_eq && kl_eq && lm_eq && !mn_eq) || (!ik_eq && !kl_eq && lm_eq && mn_eq)
                        C_T[j, i] += 4 * B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                        C_T[j, l] += 4 * A_iklmn * B_T_jk * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, k] += 4 * B_T_jl * A_iklmn * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, n] += 4 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jm
                        C_T[j, m] += 4 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jn
                    end
                    if (ik_eq && !kl_eq && lm_eq && mn_eq) || (ik_eq && kl_eq && !lm_eq && mn_eq)
                        C_T[j, i] += 2 * B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                        C_T[j, l] += 2 * A_iklmn * B_T_jk * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, k] += 2 * B_T_jl * A_iklmn * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, n] += 2 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jm
                        C_T[j, m] += 2 * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jn
                    end
                    if (ik_eq && kl_eq && lm_eq && !mn_eq) || (!ik_eq && kl_eq && lm_eq && mn_eq)
                        C_T[j, i] += B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                        C_T[j, l] += A_iklmn * B_T_jk * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, k] += B_T_jl * A_iklmn * B_T_ji * B_T_jm * B_T_jn
                        C_T[j, n] += B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jm
                        C_T[j, m] += B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jn
                    end
                    if ik_eq && kl_eq && lm_eq && mn_eq
                        C_T[j, i] += B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                    end
                end
            end
        end
    end
    return C_T
end)

function main()
    ref = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
    @btime(mttkrp_ref($ref, $A_nondiag, $B_T))

    @btime(mttkrp_opt_1($C_nondiag, $A_nondiag, $B_T))
    @btime(mttkrp_opt_2($C_diag, $A_diag, $B_T))
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= C_diag[i, j] + C_nondiag[i, j] == ref[i, j]
    end
    @info "check" check[]
    @info "check" ref == C_nondiag
end

main()

