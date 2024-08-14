eval(@finch_kernel mode=:fast function ssymm_finch_opt_helper(A, B_T, C_T)
    C .= 0
    for k = _, i = _, j = _
        let A_ik = A[i, k]
            if (i <= k)
                C_T[j, i] += A_ik * B_T[j, k]
            end
            if (i < k)
                C_T[j, k] += A_ik * B_T[j, i]
            end
        end
    end
    return C
end)