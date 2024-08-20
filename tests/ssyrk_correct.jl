eval(@finch_kernel mode=:fast function ssyrk_finch_opt_helper(A, C)
    C .= 0
    for k = _, j = _, i = _
        if (i <= j)
            C[i, j] += *(A[i, k], A[j, k])
        end
    end
    return C
end)

