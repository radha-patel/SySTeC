eval(@finch_kernel mode=:fast function ssymv_finch_opt_helper(A, diag, temp, x, y)
    y .= 0
    for j = _
        temp .= 0
        for i = _
            let A_ij = A[i, j]
                y[i] += *(x[j], A_ij)
                temp[] += *(A_ij, x[i])
            end
        end
        y[j] += *(diag[j], x[j])
        y[j] += temp[]
    end
    return y
end)

