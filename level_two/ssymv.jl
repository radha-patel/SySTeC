using Finch
using BenchmarkTools

n = 10000
A = Tensor(Dense(SparseList(Element(0))), fsprand(Int, n, n, 0.1))
x = Tensor(Dense(Element(0)), rand(Int, n))
y = Tensor(Dense(Element(0)), zeros(n))
temp = Scalar(0)

# A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((n, n), 0.1));
# x = Fiber!(Dense(Element(0.0)), rand(n));
# y = Fiber!(Dense(Element(0.0)), zeros(n));
# temp = Scalar(0.0)

eval(@finch_kernel mode=:fast function ssymv_ref(y, A, x) 
    y .= 0
    for j = _, i = _
        y[i] += A[i, j] * x[j]
    end
    y
end)

function ssymv_ref_helper(y, A, x) 
    time = @belapsed ssymv_ref($y, $A, $x)
    return (;time = time)
end

# y = Ax (A symmetric)
# Optimization: 1x flops, 2x memory bandwidth
# TIME: ~15ms
eval(@finch_kernel mode=:fast function ssymv_opt1(y, A, x, temp)
    y .= 0
    for j = _
        let temp1 = x[j]
            temp .= 0
            for i = _
                let temp3 = A[i, j]
                    if i <= j
                        y[i] += temp1 * temp3
                    end
                    if i < j
                        temp[] += temp3 * x[i]
                    end
                end
            end
            y[j] += temp[]
        end
    end
    y
end)

# eval(@finch_kernel mode=:fast function ssymv_opt2(y, A, x)
#     y .= 0
#     for j=_, i=_
#         let A_ij = A[i, j], x_j = x[j]
#             if i != j
#                 y[i] += A_ij * x_j 
#                 y[j] += A_ij * x[i]
#             end
#             if i == j
#                 y[i] += A_ij * x_j
#             end
#         end
#     end
# end)

eval(@finch_kernel mode=:fast function ssymv_opt3(y, A, x, temp)
    y .= 0
    for j=_
        let x_j = x[j]
            temp .= 0
            for i=_
                if i <= j
                    let A_ij = A[i, j], y_i = A_ij * x_j 
                        if i != j
                            y[i] += y_i
                            temp[] += A_ij * x[i]
                        end
                        # TODO: benchmark this but with following segment w runtime compilation
                        if i == j
                            y[i] += y_i
                        end
                    end
                end
            end
            y[j] += temp[]
        end
    end
end)

function ssymv_opt_helper(y, A, x) 
    temp = Scalar(0)
    time = @belapsed ssymv_opt($y, $A, $x, $temp)
    return (;time = time)
end

function main()
    ref = Tensor(Dense(Element(0)), zeros(n))
    @btime ssymv_ref($ref, $A, $x)
    @btime ssymv_opt1($y, $A, $x, $temp)
    # @btime ssymv_opt2($y, $A, $x)
    @btime ssymv_opt3($y, $A, $x, $temp)
end

main()