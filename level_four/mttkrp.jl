using Finch
using BenchmarkTools

n = 10

    triA = rand(Int, n, n, n)
    # triA = fill(1, (n, n, n))
    symA = [triA[sort([i, j, k])...] for i = 1:n, j = 1:n, k = 1:n]
    # clear diagonals
    for k=1:n, j=1:n, i=1:n 
        # if i == j || j == k || i == k
        # if j == k
        # if i == j == k
        #     symA[i, j, k] = 0
        # end
    end

    C = Tensor(Dense(Dense(Element(0))), zeros(n, n))
    T = Tensor(Dense(Dense(Element(0))), zeros(n, n))
    # _C = Tensor(Dense(Dense(Dense(Dense(Element(0))))), zeros(n, n, n, n))
    # _T = Tensor(Dense(Dense(Dense(Dense(Element(0))))), zeros(n, n, n, n))
    A = Tensor(Dense(Dense(Dense(Element(0)))), symA)
    B = Tensor(Dense(Dense(Element(0))), rand(Int, n, n))
    # B = Tensor(Dense(Dense(Element(0))), fill(1, (n, n)))

    eval(@finch_kernel mode=fastfinch function mttkrp_ref(C, A, B)
        C .= 0
        for l=_, j=_, k=_, i=_
            C[i, j] += A[i, k, l] * B[l, j] * B[k, j]
        end
    end)

    # eval(@finch_kernel mode=fastfinch function _mttkrp_ref(_C, A, B)
    #     _C .= 0
    #     for l=_, j=_, k=_, i=_
    #         _C[i, k, l, j] += A[i, k, l] * B[l, j] * B[k, j]
    #     end
    # end)

    # eval(@finch_kernel mode=fastfinch function mttkrp_opt1(C, A, B)
    #     C .= 0
    #     for j=_, l=_, k=_, i=_
    #         if k <= l # can replace <= with < if A[i, k, l] = 0 when i == j || j == k || i == k
    #             let a = A[i, k, l]
    #                 if i < k
    #                     C[i, j] += a * B[l, j] * B[k, j]
    #                 end
    #                 if i <= k # here too
    #                     C[k, j] += a * B[l, j] * B[i, j]
    #                     if k < l
    #                         C[l, j] += a * B[k, j] * B[i, j]
    #                     end
    #                 end
    #             end
    #         end
    #     end
    # end)

    # eval(@finch_kernel mode=fastfinch function mttkrp_opt2(_C, _T, A, B)
    #     _C .= 0
    #     # _T .= 0
    #     for j=_, l=_, k=_, i=_
    #         if k <= l
    #             let a = A[i, k, l]
    #                 if i < k
    #                     _C[i, k, l, j] += a * B[l, j] * B[k, j]
    #                 end
    #                 if i <= k
    #                     _C[k, i, l, j] += a * B[l, j] * B[i, j]
    #                     if k < l
    #                         _C[l, i, k, j] += a * B[k, j] * B[i, j]
    #                     end
    #                 end
    #             end
    #         end
    #         # if k <= l 
    #         #     let a = A[i, k, l]
    #         #         if i < k
    #         #             _T[i, k, l, j] += a * B[l, j] * B[k, j]
    #         #         end
    #         #         if i <= k 
    #         #             _T[k, i, l, j] += a * B[l, j] * B[i, j]
    #         #             if k < l
    #         #                 _T[l, i, k, j] += a * B[k, j] * B[i, j]
    #         #             end
    #         #         end
    #         #     end
    #         # end
    #     end
    # end)

    println("evaluating opt3")
    eval(@finch_kernel mode=fastfinch function mttkrp_opt3(C, A, B)
        C .= 0
        for j=_, l=_, k=_, i=_
            let A_ikl = A[i, k, l], B_kj = B[k, j], B_lj = B[l, j], B_ij = B[i, j]
                let ik_equal = (i == k), kl_equal = (k == l), il_equal = (i == l)
                let any_equal = ik_equal || kl_equal || il_equal, all_equal = ik_equal && kl_equal
                    if i <= k && k <= l
                        if !any_equal
                            C[i, j] += 2 * B_kj * A_ikl * B_lj
                            C[k, j] += 2 * B_ij * A_ikl * B_lj
                            C[l, j] += 2 * B_kj * B_ij * A_ikl
                        end
                        if any_equal && !all_equal
                            C[i, j] += B_kj * A_ikl * B_lj
                            C[k, j] += B_ij * A_ikl * B_lj
                            C[l, j] += B_kj * B_ij * A_ikl
                        end
                        if all_equal
                            C[i, j] += B_kj * A_ikl * B_lj
                        end
                    end
                end
                end
            end
        end
    end)
    println("done evaluating opt3")

    # eval(@finch_kernel mode=fastfinch function mttkrp_opt4(C, A, B)
    #     C .= 0
    #     for j=_, l=_, k=_, i=_
    #         if i <= k && k <= l 
    #             C[i, j] += B[k, j] * A[i, k, l] * B[l, j]
    #             if i < k && k < l 
    #                 C[i, j] += B[k, j] * A[i, k, l] * B[l, j]
    #             end
    #         end
    #         if i < k && k <= l 
    #             C[k, j] += B[i, j] * A[i, k, l] * B[l, j]
    #             if k < l 
    #                 C[k, j] += B[i, j] * A[i, k, l] * B[l, j]
    #             end
    #         end
    #         if i <= k && k < l 
    #             C[l, j] += B[k, j] * B[i, j] * A[i, k, l] 
    #             if i < k
    #                 C[l, j] += B[k, j] * B[i, j] * A[i, k, l] 
    #             end
    #         end
    #     end
    # end)

# begin
#     if and(<(i, k), <(k, l))
#         begin
#         C[i, j] <<+>>= *(2, *(B[k, j], A[i, k, l], B[l, j]))
#         C[k, j] <<+>>= *(2, *(B[i, j], A[i, k, l], B[l, j]))
#         C[l, j] <<+>>= *(2, *(B[k, j], B[i, j], A[i, k, l]))
#         end
#     end
#     if and(==(i, k), ==(k, l))
#         begin
#         C[i, j] <<+>>= *(B[k, j], A[i, k, l], B[l, j])
#         end
#     end
#     if and(!=(i, k), ==(k, l))
#         begin
#         C[k, j] <<+>>= *(B[i, j], A[i, k, l], B[l, j])
#         end
#     end
#     if and(==(i, k), !=(k, l))
#         begin
#         C[l, j] <<+>>= *(B[k, j], B[i, j], A[i, k, l])
#         end
#     end
# end

    ref = Tensor(Dense(Dense((Element(0)))), zeros(n, n))
    # _ref = Tensor(Dense(Dense(Dense(Dense((Element(0)))))), zeros(n, n, n, n))
    @btime mttkrp_ref($ref, $A, $B)
    # @btime _mttkrp_ref($_ref, $A, $B)

    # @btime mttkrp_opt1($C, $A, $B)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= (2 * C[i, j]) == ref[i, j]
    # end
    # @info "check" check[]

    # @btime mttkrp_opt2($_C, $_T, $A, $B)
    # check = Scalar(true)
    # # _res = Tensor(Dense(Dense(Dense(Dense(Element(0))))), zeros(n, n, n, n))
    # @finch for j=_, l=_, k=_, i=_
    #     if k <= l 
    #         check[] &= 2 * _C[i, k, l, j] == _ref[i, k, l, j]
    #     end
    #     # _res[i, k, l, j] = ifelse(k > l, _T[i, l, k, j], _T[i, k, l, j])
    # end
    # @finch for j=_, l=_, k=_, i=_
    #     check[] &= (_C[i, k, l, j] + _T[i, k, l, j]) == ref[i, k, l, j]
    # end
    # @info "check" check[]

    @btime mttkrp_opt3($C, $A, $B)
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= C[i, j] == ref[i, j]
    end
    @info "check" check[]

    # @btime mttkrp_opt4($C, $A, $B)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= C[i, j] == ref[i, j]
    # end
    # @info "check" check[]
