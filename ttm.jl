eval(@finch_kernel mode=:fast function ttm_finch_opt_helper(A, B_T, _C1)
    _C1 .= 0
    for l = _, k = _, j = _, i = _
        let jk_eq = (j == k), kl_eq = (k == l), A_jkl = A[j, k, l]
            if and(jk_eq, (k < l))
                _C1[i, j, l] += *(A_jkl, B_T[i, k])
            end
            if and((j < k), kl_eq)
                _C1[i, j, k] += *(A_jkl, B_T[i, l])
            end
            if and(jk_eq, kl_eq)
                _C1[i, k, l] += *(A_jkl, B_T[i, j])
            end
        end
    end
    return _C1
end)

