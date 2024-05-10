using Finch
using BenchmarkTools
using SparseArrays

n = 10000
triA = fsprand(Int, n, n, 0.1)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]

A = Tensor(Dense(SparseList(Element(0))), symA)
x = Tensor(Dense(Element(0)), rand(Int, n))
y = Scalar(0)
diag = Scalar(0)

eval(@finch_kernel mode=:fast function syprd_ref(y, A, x)
    y .= 0
    for j=_, i=_ 
        y[] += x[i] * A[i, j] * x[j] 
    end
end)

eval(@finch_kernel mode=:fast function syprd_gen(y, A, x)
    y .= 0
    for j=_, i=_
        let x_i = x[i], A_ij = A[i, j], x_j = x[j]
            if i < j
                y[] += 2 * A_ij * x_i * x_j
            end
            if i == j
                y[] += A_ij * x_i * x_j
            end
        end
    end
end)

eval(@finch_kernel mode=:fast function syprd_opt1(y, diag, A, x)
    y .= 0
    diag .= 0
    for j=_ 
        let x_j = x[j]
            for i=_ 
                let A_ij = A[i, j]
                    if i < j 
                        y[] += x[i] * A_ij * x_j
                    end
                    if i == j
                        diag[] += x[i] * A_ij * x_j
                    end
                end
            end
        end
    end
end)



function main()
    ref = Scalar(0)
    @btime syprd_ref($ref, $A, $x)

    @btime syprd_gen($y, $A, $x)
    @info "compiler generated code" ref[] == y[]

    d = Scalar(0)
    @btime syprd_opt1($y, $d, $A, $x)
    @info "hand-optimized code #1" ref[] == 2 * y[] + d[]
end

main()