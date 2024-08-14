eval(@finch_kernel mode=:fast function ttm_finch_opt_helper(A, B_T, _C1)
    _C1 .= 0
    for l = _, k = _, j = _, i = _
        let jk_leq = (j <= k), kl_leq = (k <= l), A_jkl = A[j, k, l]
            if (j < k) && kl_leq
                _C1[i, k, l] += A_jkl * B_T[i, j]
            end
            if jk_leq && (k < l)
                _C1[i, j, l] += B_T[i * k], A_jkl
            end
            if jk_leq && kl_leq
                _C1[i, j, k] += A_jkl * B_T[i, l]
            end
        end
    end
    return _C1
end)