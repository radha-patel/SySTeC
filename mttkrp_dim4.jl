eval(@finch_kernel mode=:fast function mttkrp_dim4_finch_opt_helper_base(A_nondiag, B_T, C_T)
    C_T .= 0
    for m = _, l = _, k = _, i = _, j = _
        if and((i < k), (k < l), (l < m))
            let B_T_ji = B_T[j, i], B_T_jk = B_T[j, k], A_iklm = A_nondiag[i, k, l, m], B_T_jm = B_T[j, m], B_T_jl = B_T[j, l]
                C_T[j, m] += *(6, *(B_T_jk, A_iklm, B_T_jl, B_T_ji))
                C_T[j, k] += *(6, *(A_iklm, B_T_jl, B_T_ji, B_T_jm))
                C_T[j, i] += *(6, *(B_T_jk, A_iklm, B_T_jl, B_T_jm))
                C_T[j, l] += *(6, *(B_T_jk, A_iklm, B_T_ji, B_T_jm))
            end
        end
    end
    return C_T
end)

eval(@finch_kernel mode=:fast function mttkrp_dim4_finch_opt_helper_edge(A_diag, B_T, C_T)
    C_T .= 0
    for m = _, l = _, k = _, i = _, j = _
        if and((identity(i) <= identity(k)), (identity(k) <= identity(l)), (identity(l) <= identity(m)))
            let ik_eq = (i == k), lm_eq = (l == m), kl_eq = (k == l), B_T_jm = B_T[j, m], B_T_jk = B_T[j, k], B_T_jl = B_T[j, l], B_T_ji = B_T[j, i], A_iklm = A_diag[i, k, l, m]
                if or(and(ik_eq, !kl_eq, !lm_eq), and(!ik_eq, kl_eq, !lm_eq), and(!ik_eq, !kl_eq, lm_eq))
                    C_T[j, m] += *(3, *(B_T_jk, A_iklm, B_T_jl, B_T_ji))
                    C_T[j, k] += *(3, *(A_iklm, B_T_jl, B_T_ji, B_T_jm))
                    C_T[j, i] += *(3, *(B_T_jk, A_iklm, B_T_jl, B_T_jm))
                    C_T[j, l] += *(3, *(B_T_jk, A_iklm, B_T_ji, B_T_jm))
                end
                if or(and(ik_eq, kl_eq, !lm_eq), and(!ik_eq, kl_eq, lm_eq))
                    C_T[j, l] += *(B_T_jk, A_iklm, B_T_ji, B_T_jm)
                    C_T[j, m] += *(B_T_jk, A_iklm, B_T_jl, B_T_ji)
                    C_T[j, i] += *(B_T_jk, A_iklm, B_T_jl, B_T_jm)
                    C_T[j, k] += *(A_iklm, B_T_jl, B_T_ji, B_T_jm)
                end
                if and(ik_eq, !kl_eq, lm_eq)
                    C_T[j, i] += *(3, *(B_T_jk, A_iklm, B_T_jl, B_T_jm))
                    C_T[j, l] += *(3, *(B_T_jk, A_iklm, B_T_ji, B_T_jm))
                end
                if and(ik_eq, kl_eq, lm_eq)
                    C_T[j, i] += *(B_T_jk, A_iklm, B_T_jl, B_T_jm)
                end
            end
        end
    end
    return C_T
end)

