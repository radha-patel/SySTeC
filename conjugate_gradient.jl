using Finch
using BenchmarkTools

n = 10
println("n = ", n)
A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((n, n), 0.1))
t = Fiber!(Dense(Element(0.0)), zeros(n))
x = Fiber!(Dense(Element(0.0)), zeros(n))
_x = Fiber!(Dense(Element(0.0)), zeros(n))
t = Fiber!(Dense(Element(0.0)), zeros(n))
r = Fiber!(Dense(Element(0.0)), rand(n))
_r = Fiber!(Dense(Element(0.0)), zeros(n))
p = r
_p = Fiber!(Dense(Element(0.0)), zeros(n))
alpha = Scalar(0.0)
beta = Scalar(0.0)
num1 = Scalar(0.0)
denom1 = Scalar(0.0)
num2 = Scalar(0.0)
denom2 = Scalar(0.0)

# @finch begin 
#     for k = 1:100
#         t .= 0
#         for j = _, i = _
#             t[i] += A[i, j] * p[j]
#         end

#         num1 .= 0
#         denom1 .= 0
#         for i = _
#             num1[] += r[i]^2
#             denom1[] += (p[i] * t[i])
#         end

#         alpha .= 0
#         alpha[] = num[] / denom[]

#         _x .= 0
#         _r .= 0
#         for i = _
#             _x[i] = x[i] + alpha[] * p[i]
#             _r[i] = r[i] - alpha[] * t[i]
#         end

#         num2 .= 0
#         for i = _
#             num2[] += _r[i]^2
#         end

#         beta .= 0
#         beta[] = num2[] / num1[]

#         _p .= 0
#         for i = _
#             _p[i] = _r[i] + beta[] * p[i]
#         end

#         p .= 0
#         x .= 0
#         r .= 0
#         for i = _
#             p[i] = _p[i]
#             x[i] = _x[i]
#             r[i] = _r[i]
#         end
#     end
# end

# eval(@finch_kernel function cg_naive(A, t, p, _p, r, _r, x, _x, alpha, beta, num1, denom1, num2) 
#     t .= 0
#     # SpMV
#     for j = _, i = _
#         t[i] += A[i, j] * p[j]
#     end

#     num1 .= 0
#     denom1 .= 0
#     for i = _
#         num1[] += r[i]^2
#         denom1[] += (p[i] * t[i])
#     end

#     alpha .= 0
#     alpha[] = num1[] / denom1[]

#     _x .= 0
#     _r .= 0
#     for i = _
#         _x[i] = x[i] + alpha[] * p[i]
#         _r[i] = r[i] - alpha[] * t[i]
#     end

#     num2 .= 0
#     for i = _
#         num2[] += _r[i]^2
#     end

#     beta .= 0
#     beta[] = num2[] / num1[]

#     _p .= 0
#     for i = _
#         _p[i] = _r[i] + beta[] * p[i]
#     end
# end)

# @btime cg_naive($A, $t, $p, $_p, $r, $_r, $x, $_x, $alpha, $beta, $num1, $denom1, $num2)

temp2 = Scalar(0.0)
eval(@finch_kernel function cg(A, t, p, _p, r, _r, x, _x, alpha, beta, num1, denom1, num2, temp2) 
    for l = 1:10
        t .= 0
        # SpMV - A * p_k
        for j = _
            let temp1 = x[j]
                temp2 .= 0
                for i = _
                    let temp3 = A[i, j]
                        if uptrimask[i, j]
                            t[i] += temp1 * temp3
                        end
                        if uptrimask[i, j - 1]
                            temp2[] += temp3 * x[i]
                        end
                    end
                end
            end
            t[j] += temp2[]
        end

        # alpha_k = r_k^T * r_k / (p_k^T * A * p_k)
        num1 .= 0
        denom1 .= 0
        for i = _
            num1[] += r[i]^2
            denom1[] += (p[i] * t[i])
        end
        alpha .= 0
        alpha[] = num1[] / denom1[]

        # x_k+1 = x_k + alpha_k * p_k
        # r_k+1 = r_k - alpha_k * A * p_k
        _x .= 0
        _r .= 0
        for i = _
            _x[i] = x[i] + alpha[] * p[i]
            _r[i] = r[i] - alpha[] * t[i]
        end

        # beta_k = r_k+1^T * r_k+1 / (r_k^T * r_k)
        num2 .= 0
        for i = _
            num2[] += _r[i]^2
        end
        beta .= 0
        beta[] = num2[] / num1[]

        # p_k+1 = r_k+1 * beta_k * p_k
        _p .= 0
        for i = _
            _p[i] = _r[i] + beta[] * p[i]
        end
    end
end)

@btime cg($A, $t, $p, $_p, $r, $_r, $x, $_x, $alpha, $beta, $num1, $denom1, $num2, $temp2)
