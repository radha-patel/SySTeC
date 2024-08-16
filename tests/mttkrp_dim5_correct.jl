eval(@finch_kernel mode=:fast function mttkrp_dim5_finch_opt_helper_base(A, B_T, C_T)
    C_T .= 0
    for n = _, m = _, l = _, k = _, i = _, j = _
        if and((i < k), (k < l), (l < m), (m < n))
            let A_iklmn = A[i, k, l, m, n], B_T_ji = B_T[j, i], B_T_jk = B_T[j, k], B_T_jm = B_T[j, m], B_T_jn = B_T[j, n], B_T_jl = B_T[j, l]
                C_T[j, i] += *(24, *(B_T_jk, A_iklmn, B_T_jn, B_T_jl, B_T_jm))
                C_T[j, n] += *(24, *(B_T_jk, A_iklmn, B_T_jl, B_T_ji, B_T_jm))
                C_T[j, l] += *(24, *(B_T_jk, A_iklmn, B_T_jn, B_T_ji, B_T_jm))
                C_T[j, k] += *(24, *(A_iklmn, B_T_jn, B_T_jl, B_T_ji, B_T_jm))
                C_T[j, m] += *(24, *(B_T_jk, A_iklmn, B_T_jn, B_T_jl, B_T_ji))
            end
        end
    end
    return C_T
end)

eval(@finch_kernel mode=:fast function mttkrp_dim5_finch_opt_helper_edge(A, B_T, C_T, lookup)
    C_T .= 0
    for n = _, m = _, l = _, k = _, i = _, j = _
        if and((identity(i) <= identity(k)), (identity(k) <= identity(l)), (identity(l) <= identity(m)), (identity(m) <= identity(n)))
            let mn_eq = (m == n), ik_eq = (i == k), lm_eq = (l == m), kl_eq = (k == l), B_T_jm = B_T[j, m], B_T_jk = B_T[j, k], B_T_jl = B_T[j, l], A_iklmn = A[i, k, l, m, n], B_T_ji = B_T[j, i], B_T_jn = B_T[j, n]
                if and(ik_eq, kl_eq, lm_eq, mn_eq)
                    C_T[j, m] += *(B_T_jk, A_iklmn, B_T_jn, B_T_jl, B_T_ji)
                end
                let idx = +(+(+(*(2, ik_eq), *(3, kl_eq)), *(5, lm_eq)), *(7, mn_eq)), factor = lookup[idx]
                    C_T[j, i] += *(factor, *(B_T_jk, A_iklmn, B_T_jn, B_T_jl, B_T_jm))
                    C_T[j, n] += *(factor, *(B_T_jk, A_iklmn, B_T_jl, B_T_ji, B_T_jm))
                    C_T[j, m] += *(factor, *(B_T_jk, A_iklmn, B_T_jn, B_T_jl, B_T_ji))
                    C_T[j, l] += *(factor, *(B_T_jk, A_iklmn, B_T_jn, B_T_ji, B_T_jm))
                    C_T[j, k] += *(factor, *(A_iklmn, B_T_jn, B_T_jl, B_T_ji, B_T_jm))
                end
            end
        end
    end
    return C_T
end)

lookup = [0.0, 12.0, 12.0, 0.0, 12.0, 4.0, 12.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 6.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]