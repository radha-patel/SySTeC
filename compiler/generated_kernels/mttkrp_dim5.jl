using Finch
using BenchmarkTools

n = 10
triA = fsprand(Int, n, n, n, n, n, 1.0)
symA = [triA[sort([i, j, k, l, m])...] for i = 1:n, j = 1:n, k = 1:n, l = 1:n, m = 1:n]
nondiagA = zeros(Int, n, n, n, n, n)
diagA = zeros(Int, n, n, n, n, n)
b = rand(Int, n, n)

lookup = Tensor(Dense(Element(0)), zeros(Int, 211))
lookup[2 + 1] = 12
lookup[3 + 1] = 12
lookup[5 + 1] = 12
lookup[7 + 1] = 12
lookup[10 + 1] = 6
lookup[21 + 1] = 6
lookup[14 + 1] = 6
lookup[6 + 1] = 4
lookup[15 + 1] = 4
lookup[35 + 1] = 4
lookup[42 + 1] = 2
lookup[70 + 1] = 2
lookup[30 + 1] = 1
lookup[105 + 1] = 1

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
A_diag = Tensor(Dense(Dense(SparseList(SparseList(SparseList(Element(0)))))), diagA)
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

println("begin eval mttkrp_opt_1")
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
println("after eval mttkrp_opt_1")

println("begin eval mttkrp_opt_2")
eval(@finch_kernel mode=:fast function mttkrp_opt_2(C_T, A_diag, B_T, lookup)
    C_T .= 0
    for n=_, m=_, l=_, k=_, i=_, j=_
        if i <= k && k <= l && l <= m && m <= n
        # if identity(i) <= identity(k) && identity(k) <= identity(l) && identity(l) <= identity(m) && identity(m) <= identity(n)
            let ik_eq = (i == k), mn_eq = (m == n), kl_eq = (identity(k) == identity(l)), lm_eq = (identity(l) == identity(m))
                let A_iklmn = A_diag[i, k, l, m, n], B_T_jl = B_T[j, l], B_T_jk = B_T[j, k], B_T_ji = B_T[j, i], B_T_jm = B_T[j, m], B_T_jn = B_T[j, n]
                    let idx = (ik_eq) * 2 + (kl_eq) * 3 + (lm_eq) * 5 + (mn_eq) * 7 + 1
                        let factor = lookup[idx]
                            C_T[j, i] += factor * B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                            C_T[j, l] += factor * A_iklmn * B_T_jk * B_T_ji * B_T_jm * B_T_jn
                            C_T[j, k] += factor * B_T_jl * A_iklmn * B_T_ji * B_T_jm * B_T_jn
                            C_T[j, n] += factor * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jm
                            C_T[j, m] += factor * B_T_jl * A_iklmn * B_T_jk * B_T_ji * B_T_jn
                        end
                    end
                    if ik_eq && mn_eq && kl_eq && lm_eq
                        C_T[j, i] += B_T_jl * A_iklmn * B_T_jk * B_T_jm * B_T_jn
                    end
                end
            end
        end
    end
    return C_T
end)
println("after eval mttkrp_opt_2")


function main()
    ref = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
    @btime(mttkrp_ref($ref, $A, $B_T))

    @btime(mttkrp_opt_1($C_nondiag, $A_nondiag, $B_T))
    @btime(mttkrp_opt_2($C_diag, $A_diag, $B_T, $lookup))
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= C_diag[i, j] + C_nondiag[i, j] == ref[i, j]
    end
    @info "check" check[]
end

main()

# ik_eq, mn_eq, kl_eq, lm_eq
# 2 * 3 * 5 * 7

# 1 equal -> 2, 3, 5, 7 -> 12
# 2 equal -> 10, 21, 14 -> 6
# 2 equal tg -> 6, 15, 35 -> 4
# 3 equal -> 42, 70 -> 2
# 3 equal tg -> 30, 105 -> 1
# (ik_eq + kl_eq + lm_eq + mn_eq) -> 2 * 3 * 5 * 7 = 210 = 





# lookup = []
