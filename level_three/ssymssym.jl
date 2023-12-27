using Finch
using BenchmarkTools
using SparseArrays

n = 10

# triA = rand(Int, n, n)
triA = fill(1, (n, n))
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
# triB = rand(Int, n, n)
triB = fill(1, (n, n))
symB = [triB[sort([i, j])...] for i = 1:n, j = 1:n]

# A = Fiber!(Dense(SparseList(Element(0))), symA)    
# B = Fiber!(Dense(SparseList(Element(0))), symB)    
A = Fiber!(Dense(Dense(Element(0))), symA)    
B = Fiber!(Dense(Dense(Element(0))), symB) 
C = Fiber!(Dense(Dense(Element(0))), zeros(n, n))
_C = Fiber!(Dense(Dense(Element(0))), zeros(n, n))

temp1 = Scalar(0)
temp2 = Scalar(0)
temp3 = Scalar(0)

eval(@finch_kernel mode=fastfinch function ssymssym_ref(C, A, B)
    C .= 0
    for j = _, k = _, i = _
        C[i, j] += A[i, k] * B[k, j]
    end
end)

eval(@finch_kernel mode=fastfinch function ssymssym_opt1(C, A, B, temp2)
    C .= 0
    for j=_, i =_
        temp2 .= 0
        let temp1 = B[i, j]
            for k =_
                let temp3 = A[k, i]
                    if uptrimask[k+1, i]
                        C[k, j] += temp1 * temp3
                    end
                    if uptrimask[k, i]
                        temp2[] += B[k, j] * temp3
                    end
                end
            end
        end
        C[i, j] += temp2[]
    end
end)

eval(@finch_kernel mode=fastfinch function ssymssym_opt2(C, A, B, temp2)
    C .= 0
    for j=_, k=_, i=_
        if i <= j
            C[i, j] += A[i, k] * B[k, j]
        end
    end
end)

eval(@finch_kernel mode=fastfinch function ssymssym_opt3(C, A, B, temp2)
    C .= 0
    for j=_, i=_, k=_
        if i <= j
            if k < i 
                C[i, j] += A[k, i] * B[k, j]
            end 
            if k <= i 
                C[k, j] += A[k, i] * B[i, j]
                if i < j
                    C[k, i] += A[k, j] * B[i, j]
                end
            end
        end
    end
end)

eval(@finch_kernel mode=fastfinch function ssymssym_opt4(C, A, B, temp2)
    C .= 0
    for j=_, i=_, k=_
        if i <= j
            let a = A[k, i]
                let b = B[i, j]
                    if k < i 
                        C[i, j] += a * B[k, j]
                    end 
                    if k <= i 
                        C[k, j] += a * b
                        if i < j
                            C[k, i] += A[k, j] * b
                        end
                    end
                end
            end
        end
    end
end)

# eval(@finch_kernel mode=fastfinch function ssymssym_opt5(C, A, B, temp2)
#     C .= 0
#     for j=_, i=_, k=_
#         if k <= i 
#             if j <= k
#                 C[i, j] += A[k, i] * B[k, j]
#                 if k < i 
#                     C[k, j] += A[k, i] * B[i, j]
#                 end
#             end
#             if j < k
#                 C[i, k] += A[i, j] * B[k, j]
#             end
#         end
#     end
# end)

eval(@finch_kernel mode=fastfinch function ssymssym_opt5(C, _C, A, B, temp2)
    C .= 0
    _C .= 0
    for j=_, i=_, k=_
        if k <= i 
            if j <= k
                C[i, j] += A[k, i] * B[k, j]
                if k < i 
                    C[k, j] += A[k, i] * B[i, j]
                    C[k, i] += A[k, j] * B[i, j] # i<->j from prev line to account for implicit i <= j  
                end
            end
            if j < k
                C[i, k] += A[i, j] * B[k, j]
                _C[k, j] += A[i, j] * B[k, i] # i<->j from prev line to account for implicit i <= j
                # we place this condition here to prevent double counting
                if k < i 
                    _C[i, j] += A[k, j] * B[k, i] # i<->k for the missing update + i<->j to account for implicit i <= j 
                end
            end
        end
    end
end)

# eval(@finch_kernel mode=fastfinch function ssymssym_opt6(C, _C, A, B, temp1, temp2, temp3)
#     C .= 0
#     _C .= 0
#     for j=_, i=_
#         temp1 .= 0
#         temp2 .= 0
#         for k=_
#             let A_ki = A[k, i]
#                 let A_kj = A[k, j]
#                     let A_ij = A[i, j]
#                         let B_kj = B[k, j]
#                             let B_ij = B[i, j]
#                                 let B_ki = B[k, i]
#                                     if k <= i 
#                                         if j <= k
#                                             temp1[] += A_ki * B_kj # C[i, j] -> temp1
#                                             if k < i 
#                                                 C[k, j] += A_ki * B_ij
#                                                 C[k, i] += A_kj * B_ij
#                                             end
#                                         end
#                                         if j < k
#                                             C[i, k] += A_ij * B_kj
#                                             _C[k, j] += A_ij * B_ki
#                                             if k < i 
#                                                 temp2[] += A_kj * B_ki # _C[i, j] -> temp2
#                                             end
#                                         end
#                                     end
#                                 end
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#         C[i, j] += temp1[]
#         C[j, i] += temp2[]
#     end
# end)

function main()
    ref = Fiber!(Dense(Dense(Element(0))), zeros(n, n))
    @btime ssymssym_ref($ref, $A, $B)

    # @btime ssymssym_opt1($C, $A, $B, $temp2)
    # @info "check opt1" C == ref

    # @btime ssymssym_opt2($C, $A, $B, $temp2)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # @info "check opt2" check[]

    # @btime ssymssym_opt3($C, $A, $B, $temp2)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # @info "check opt3" check[]

    # @btime ssymssym_opt4($C, $A, $B, $temp2)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     if i <= j
    #         check[] &= C[i, j] == ref[i, j]
    #     end
    # end
    # @info "check opt4" check[]

    @btime ssymssym_opt5($C, $_C, $A, $B, $temp2)
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= (C[i, j] + _C[j, i]) == ref[i, j]
    end
    @info "check opt5" check[]

    # @btime ssymssym_opt6($C, $_C, $A, $B, $temp1, $temp2, $temp3)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= (C[i, j] + _C[j, i]) == ref[i, j]
    # end
    # @info "check opt6" check[]
end

main()