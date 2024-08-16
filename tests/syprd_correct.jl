eval(@finch_kernel mode=:fast function syprd_finch_opt_helper(A, diag, x, y)
    y .= 0
    for j = _
        for i = _
            y[] += *(2, *(x[j], A[i, j], x[i]))
        end
        y[] += *(x[j], diag[j], x[j])
    end
    return y
end)

