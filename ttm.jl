eval(@finch_kernel mode=:fast function ttm_finch_opt_helper(A, B_T, _C1)
    _C1 .= 0
    for l = _, k = _, j = _, i = _
        if and((j <= k), (k <= l))
            let jk_eq = (j == k), kl_eq = (k == l), A_jkl = A[j, k, l], B_T_ij = B_T[i, j], B_T_ik = B_T[i, k], B_T_il = B_T[i, l]
                if and(!jk_eq, !kl_eq)
                    _C1[i, k, l] += *(A_jkl, B_T_ij)
                    _C1[i, j, l] += *(A_jkl, B_T_ik)
                    _C1[i, j, k] += *(A_jkl, B_T_il)
                end
                if and(jk_eq, !kl_eq)
                    _C1[i, l, k] += *(A_jkl, B_T_ij)
                    _C1[i, j, l] += *(A_jkl, B_T_ik)
                    _C1[i, j, k] += *(A_jkl, B_T_il)
                end
                if and(!jk_eq, kl_eq)
                    _C1[i, k, j] += *(A_jkl, B_T_il)
                    _C1[i, k, l] += *(A_jkl, B_T_ij)
                    _C1[i, j, l] += *(A_jkl, B_T_ik)
                end
                if and(jk_eq, kl_eq)
                    _C1[i, j, l] += *(A_jkl, B_T_ik)
                end
            end
        end
    end
    return _C1
end)

