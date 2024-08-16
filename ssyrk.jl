eval(@finch_kernel mode=:fast function ssyrk_finch_opt_helper(A, C)
    C .= 0
    if (i <= j)
        C[i, j] += *(A[i, k], A[j, k])
    end
    return C
end)

