using Finch
using BenchmarkTools

n = 100
j = n

    triA = rand(Int, n, n, n)
    symA = [triA[sort([i, j, k])...] for i = 1:n, j = 1:n, k = 1:n]
    x = rand(Int, n, j)

    A = Fiber!(Dense(SparseList(SparseList(Element(0)))), symA)
    X = Fiber!(Dense(Dense(Element(0))), x)
    X_T = Fiber!(Dense(Dense(Element(0))), transpose(x))
    C = Fiber!(Dense(Dense(Dense(Element(0)))), zeros(n, j, n))
    _C = Fiber!(Dense(Dense(Dense(Element(0)))), zeros(j, n, n))
    temp2 = Scalar(0)
    temp4 = Fiber!(Dense(Element(0)), rand(Int, n))

    eval(@finch_kernel mode=fastfinch function mode2_product_ref(C, A, X)
        C .= 0
        for l=_, j=_, k=_, i=_ 
            C[i, j, l] += A[i, k, l] * X[k, j] 
        end
    end)

    eval(@finch_kernel mode=fastfinch function mode1_product_ref(C, A, X)
        C .= 0
        for l=_, j=_, k=_, i=_ 
            C[i, j, l] += A[k, j, l] * X[k, i] 
        end
    end)

    # eval(@finch_kernel mode=fastfinch function mode3_product_ref(C, A, X)
    #     C .= 0
    #     for l=_, k=_, j=_, i=_ 
    #         C[i, j, l] += A[i, j, k] * X[k, l] 
    #     end
    # end)

    # eval(@finch_kernel mode=fastfinch function mode2_product_opt1(C, A, X, temp2) 
    #     C .= 0
    #     for l=_, j=_, i=_
    #         temp2 .= 0
    #         let temp1 = X[i, j]
    #             for k=_ 
    #                 let temp3 = A[k, i, l]
    #                     if k <= i 
    #                         C[k, j, l] += temp1 * temp3
    #                     end
    #                     if k < i
    #                         temp2[] += X[k, j] * temp3
    #                     end
    #                 end
    #             end
    #         end
    #         C[i, j, l] += temp2[]
    #     end
    # end)

    # eval(@finch_kernel mode=fastfinch function mode2_product_opt2(C, A, X, temp2, temp4) 
    #     C .= 0
    #     for j=_, i=_
    #         temp4 .= 0
    #         for l=_
    #             if i <= l
    #                 temp2 .= 0
    #                 let temp1 = X[i, j]
    #                     for k=_ 
    #                         let temp3 = A[k, i, l]
    #                             if k <= i
    #                                 C[k, j, l] += temp1 * temp3
    #                             end
    #                             if k < i
    #                                 temp2[] += X[k, j] * temp3
    #                             end
    #                         end
    #                     end
    #                 end
    #                 C[i, j, l] += temp2[]
    #             end
    #             if i < l
    #                 for k=_
    #                     if k <= i
    #                         temp4[k] += A[k, i, l] * X[l, j]
    #                     end
    #                 end
    #             end
    #         end
    #         for k=_
    #             if k <= i
    #                 C[k, j, i] += temp4[k]
    #             end
    #         end
    #     end
    # end)

    # eval(@finch_kernel mode=fastfinch function mode2_product_opt3(C, A, X, temp2) 
    #     C .= 0
    #     for l=_, j=_, i=_
    #         temp2 .= 0
    #         let temp1 = X[i, j]
    #             for k=_ 
    #                 if k <= l
    #                     let temp3 = A[k, i, l]
    #                         if k <= i 
    #                             C[k, j, l] += temp1 * temp3
    #                         end
    #                         if k < i
    #                             temp2[] += X[k, j] * temp3
    #                         end
    #                     end
    #                 end
    #             end
    #         end
    #         C[i, j, l] += temp2[]
    #     end
    # end)

    # eval(@finch_kernel mode=fastfinch function mode2_product_opt4(C, A, X, temp2) 
    #     C .= 0
    #     for j=_, l=_, i=_
    #         if i <= l
    #             temp2 .= 0
    #             for k=_ 
    #                 let temp3 = A[k, i, l]
    #                     if k <= i
    #                         C[k, j, l] += X[i, j] * temp3
    #                         if i < l
    #                             C[k, j, i] += X[l, j] * temp3
    #                         end
    #                     end
    #                     if k < i
    #                         temp2[] += X[k, j] * temp3
    #                     end
    #                 end
    #             end
    #             C[i, j, l] += temp2[]
    #         end
    #     end
    # end)

    # eval(@finch_kernel mode=fastfinch function mode1_product_opt4(C, A, X, temp2) 
    #     C .= 0
    #     for j=_, l=_, i=_
    #         if i <= l
    #             temp2 .= 0
    #             for k=_ 
    #                 let temp3 = A[k, i, l]
    #                     if k <= i
    #                         C[k, j, l] += X[i, j] * temp3
    #                         if i < l
    #                             C[k, j, i] += X[l, j] * temp3
    #                         end
    #                     end
    #                     if k < i
    #                         temp2[] += X[k, j] * temp3
    #                     end
    #                 end
    #             end
    #             C[i, j, l] += temp2[]
    #         end
    #     end
    # end)

    # eval(@finch_kernel mode=fastfinch function mode3_product_opt4(C, A, X, temp2) 
    #     C .= 0
    #     for j=_, l=_, i=_
    #         if i <= l
    #             temp2 .= 0
    #             for k=_ 
    #                 let temp3 = A[k, i, l]
    #                     if k <= i
    #                         C[k, j, l] += X[i, j] * temp3
    #                         if i < l
    #                             C[k, j, i] += X[l, j] * temp3
    #                         end
    #                     end
    #                     if k < i
    #                         temp2[] += X[k, j] * temp3
    #                     end
    #                 end
    #             end
    #             C[i, j, l] += temp2[]
    #         end
    #     end
    # end)

    # will need to do a permute afterwards to switch i and j in C[j, i, l]
    # eval(@finch_kernel mode=fastfinch function mode2_product_opt5(C, _C, A, X, X_T, temp2) 
    #     C .= 0
    #     _C .= 0
    #     for l=_, i=_
    #         if i <= l
    #             for j=_
    #                 temp2 .= 0
    #                 for k=_ 
    #                     let temp3 = A[k, i, l]
    #                         if k <= i
    #                             C[k, j, l] += X_T[j, i] * temp3
    #                             if i < l
    #                                 C[k, j, i] += X_T[j, l] * temp3
    #                             end
    #                         end
    #                         if k < i
    #                             temp2[] += X[k, j] * temp3
    #                         end
    #                     end
    #                 end
    #                 _C[j, i, l] += temp2[]
    #             end
    #         end
    #     end
    # end)

function main()
    # MODE 2
    # ref = Fiber!(Dense(Dense(Dense(Element(0), n), n), n))
    # @btime mode2_product_ref($ref, $A, $X)

    # @btime mode2_product_opt1($C, $A, $X, $temp2)

    # @info "check" C == ref

    # @btime mode2_product_opt2($C, $A, $X, $temp2, $temp4)

    # check = Scalar(true)
    # @finch for l=_, j=_, i=_
    #     if i <= l
    #         check[] &= C[i, j, l] == ref[i, j, l]
    #     end
    # end
    # @info "check" check[]

    # @btime mode2_product_opt3($C, $A, $X, $temp2)

    # check = Scalar(true)
    # @finch for l=_, j=_, i=_
    #     if i <= l
    #         check[] &= C[i, j, l] == ref[i, j, l]
    #     end
    # end
    # @info "check" check[]

    # @btime mode2_product_opt4($C, $A, $X, $temp2)

    # check = Scalar(true)
    # @finch for l=_, j=_, i=_
    #     if i <= l
    #         check[] &= C[i, j, l] == ref[i, j, l]
    #     end
    # end
    # @info "check" check[]

    # @btime mode2_product_opt5($C, $_C, $A, $X, $X_T, $temp2)

    # check = Scalar(true)
    # @finch for l=_, j=_, i=_
    #     if i <= l
    #         check[] &= (C[i, j, l] + _C[j, i, l]) == ref[i, j, l]
    #     end
    # end
    # @info "check" check[]

    # MODE 1
    ref = Fiber!(Dense(Dense(Dense(Element(0), n), n), n))
    @btime mode1_product_ref($ref, $A, $X)

    @btime mode1_product_opt1($C, $A, $X, $temp2)

    @info "check" C == ref
end

main()