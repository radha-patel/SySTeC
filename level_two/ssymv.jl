using Finch
using BenchmarkTools

n = 10000
triA = fsprand(Int, n, n, 0.1)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]

A = Tensor(Dense(SparseList(Element(0))), symA)
x = Tensor(Dense(Element(0)), rand(Int, n))
y = Tensor(Dense(Element(0)), zeros(n))
diag = Tensor(Dense(Element(0)), zeros(n))
temp = Scalar(0)

eval(@finch_kernel mode=:fast function ssymv_ref(y, A, x) 
    y .= 0
    for j = _, i = _
        y[i] += A[i, j] * x[j]
    end
    y
end)

eval(@finch_kernel mode=:fast function ssymv_gen(y, A, x)
    y .= 0
    for j=_, i=_
        let A_ij = A[i, j]
            if i <= j
                y[i] += A_ij * x[j]
            end
            if i < j
                y[j] += A_ij * x[i]
            end
        end
    end
end)

eval(@finch_kernel mode=:fast function ssymv_opt1(y, A, x, temp)
    y .= 0
    for j = _
        let x_j = x[j]
            temp .= 0
            for i = _
                let A_ij = A[i, j]
                    if i <= j
                        y[i] += x_j * A_ij
                    end
                    if i < j
                        temp[] += A_ij * x[i]
                    end
                end
            end
            y[j] += temp[]
        end
    end
    y
end)

eval(@finch_kernel mode=:fast function ssymv_opt2(y, A, x, diag, temp)
    y .= 0
    for j = _
        let x_j = x[j]
            temp .= 0
            for i = _
                let A_ij = A[i, j]
                    y[i] += x_j * A_ij
                    temp[] += A_ij * x[i]
                end
            end
            y[j] += temp[] + diag[j] * x_j
        end
    end
    y
end)

function main()
    ref = Tensor(Dense(Element(0)), zeros(n))
    @btime ssymv_ref($ref, $A, $x)

    @btime ssymv_gen($y, $A, $x)
    @info "compiler generated code" ref == y

    @btime ssymv_opt1($y, $A, $x, $temp)
    @info "hand-optimized code #1" ref == y

    _A = Tensor(Dense(SparseList(Element(0))))
    _d = Tensor(Dense(Element(0)))
    @finch mode=:fast begin
        _A .= 0
        _d .= 0
        for j = _, i = _
            if i < j
                _A[i, j] = A[i, j]
            end
            if i == j
                _d[i] = A[i, j]
            end
        end
    end
    @btime ssymv_opt2($y, $_A, $x, $_d, $temp)
    @info "hand-optimized code #2" ref == y
end

main()