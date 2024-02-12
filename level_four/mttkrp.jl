using Finch
using BenchmarkTools

n = 100

    triA = rand(Int, n, n, n)
    # triA = fill(1, (n, n, n))
    symA = [triA[sort([i, j, k])...] for i = 1:n, j = 1:n, k = 1:n]
    nondiagA = zeros(Int, n, n, n)
    diagA = zeros(Int, n, n, n)
    # clear diagonals
    for k=1:n, j=1:n, i=1:n 
        if i != j && j != k && i != k
            nondiagA[i, j, k] = symA[i, j, k]
        end
        if i == j || j == k || i == k
            diagA[i, j, k] = symA[i, j, k]
        end
    end

    C = Tensor(Dense(Dense(Element(0))), zeros(n, n))
    T = Tensor(Dense(Dense(Element(0))), zeros(n, n))
    # _C = Tensor(Dense(Dense(Dense(Dense(Element(0))))), zeros(n, n, n, n))
    # _T = Tensor(Dense(Dense(Dense(Dense(Element(0))))), zeros(n, n, n, n))
    A_nondiag = Tensor(Dense(SparseList(SparseList(Element(0)))), nondiagA)
    A_diag = Tensor(Dense(SparseList(SparseList(Element(0)))), diagA)
    A = Tensor(Dense(Dense(Dense(Element(0)))), symA)
    B = Tensor(Dense(Dense(Element(0))), rand(Int, n, n))
    # B = Tensor(Dense(Dense(Element(0))), fill(1, (n, n)))

    # temporaries
    C_jk = Scalar(0)
    C_jl = Scalar(0)

    eval(@finch_kernel mode=fastfinch function mttkrp_ref(C, A, B)
        C .= 0
        for l=_, j=_, k=_, i=_
            C[i, j] += A[i, k, l] * B[l, j] * B[k, j]
        end
    end)

    eval(@finch_kernel mode=fastfinch function mttkrp_ref(C, A_diag, B)
        C .= 0
        for l=_, j=_, k=_, i=_
            C[i, j] += A_diag[i, k, l] * B[l, j] * B[k, j]
        end
    end)

    # eval(@finch_kernel mode=fastfinch function _mttkrp_ref(_C, A, B)
    #     _C .= 0
    #     for l=_, j=_, k=_, i=_
    #         _C[i, k, l, j] += A[i, k, l] * B[l, j] * B[k, j]
    #     end
    # end)

    # assumes diagonals have been zeroed (~3.8x faster than ref) - A is Dense(Sparse(Sparse))
    eval(@finch_kernel mode=fastfinch function mttkrp_opt1(C, A_nondiag, B)
        C .= 0
        for l=_, k=_, i=_, j=_
            if i < k && k < l
                let a = A_nondiag[i, k, l], B_lj = B[j, l], B_kj = B[j, k], B_ij = B[j, i]
                    C[j, i] += a * B_lj * B_kj
                    C[j, k] += a * B_lj * B_ij
                    C[j, l] += a * B_kj * B_ij
                end
            end
        end
    end)

    # assumes diagonals have been zeroed (~2.6x faster than ref)
    # eval(@finch_kernel mode=fastfinch function mttkrp_opt1_1(C, C_jk, C_jl, A, B)
    #     C .= 0
    #     for l=_, k=_, j=_
    #         C_jk .= 0
    #         C_jl .= 0
    #         for i=_
    #             if i < k && k < l
    #                 let a = A[i, k, l], B_lj = B[j, l], B_kj = B[j, k], B_ij = B[j, i]
    #                     C[j, i] += a * B_lj * B_kj
    #                     C_jk[] += a * B_lj * B_ij
    #                     C_jl[] += a * B_kj * B_ij
    #                 end
    #             end
    #         end
    #         C[j, k] += C_jk[]
    #         C[j, l] += C_jl[]
    #     end
    # end)

    # assumes diagonals have been zeroed (~2.2x faster than ref)
    # eval(@finch_kernel mode=fastfinch function mttkrp_opt1_2(C, C_jk, C_jl, A, B)
    #     C .= 0
    #     for j=_, l=_, k=_
    #         C_jk .= 0
    #         C_jl .= 0
    #         for i=_
    #             if i < k && k < l
    #                 let a = A[i, k, l], B_lj = B[l, j], B_kj = B[k, j], B_ij = B[i, j]
    #                     C[i, j] += a * B_lj * B_kj
    #                     C_jk[] += a * B_lj * B_ij
    #                     C_jl[] += a * B_kj * B_ij
    #                 end
    #             end
    #         end
    #         C[k, j] += C_jk[]
    #         C[l, j] += C_jl[]
    #     end
    # end)

    # assumes diagonals have been zeroed (~2.8x faster than ref)
    # eval(@finch_kernel mode=fastfinch function mttkrp_opt1_3(C, C_jk, C_jl, A, B)
    #     C .= 0
    #     for l=_, j=_, k=_
    #         C_jk .= 0
    #         C_jl .= 0
    #         for i=_
    #             if i < k && k < l
    #                 let a = A[i, k, l], B_lj = B[l, j], B_kj = B[k, j], B_ij = B[i, j]
    #                     C[i, j] += a * B_lj * B_kj
    #                     C_jk[] += a * B_lj * B_ij
    #                     C_jl[] += a * B_kj * B_ij
    #                 end
    #             end
    #         end
    #         C[k, j] += C_jk[]
    #         C[l, j] += C_jl[]
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

    # println("evaluating opt3")
    # eval(@finch_kernel mode=fastfinch function mttkrp_opt3(C, A, B)
    #     C .= 0
    #     for l=_, k=_, i=_, j=_
    #         let A_ikl = A[i, k, l], B_kj = B[j, k], B_lj = B[j, l], B_ij = B[j, i]
    #             let ik_equal = (i == k), kl_equal = (k == l)
    #             let any_equal = ik_equal || kl_equal, all_equal = ik_equal && kl_equal
    #                 if i <= k && k <= l
    #                     let C_ij = B_kj * A_ikl * B_lj, C_kj = B_ij * A_ikl * B_lj, C_lj = B_kj * B_ij * A_ikl
    #                         if !any_equal
    #                             C[j, i] += 2 * C_ij
    #                             C[j, k] += 2 * C_kj
    #                             C[j, l] += 2 * C_lj
    #                         end
    #                         if any_equal && !all_equal
    #                             C[j, i] += C_ij
    #                             C[j, k] += C_kj
    #                             C[j, l] += C_lj
    #                         end
    #                         if all_equal
    #                             C[j, i] += C_ij
    #                         end
    #                     end
    #                 end
    #             end
    #             end
    #         end
    #     end
    # end)
    # println("done evaluating opt3")

    # Only evaluates diagonals
    println("evaluating opt3_1")
    eval(@finch_kernel mode=fastfinch function mttkrp_opt3_1(C, A_diag, B)
        C .= 0
        for l=_, k=_, i=_, j=_
            let A_ikl = A_diag[i, k, l], B_kj = B[j, k], B_lj = B[j, l], B_ij = B[j, i]
                let ik_equal = (identity(i) == identity(k)), kl_equal = (identity(k) == identity(l))
                let any_equal = ik_equal || kl_equal, all_equal = ik_equal && kl_equal
                    if identity(i) <= identity(k) && identity(k) <= identity(l)
                        let C_ij = B_kj * A_ikl * B_lj, C_kj = B_ij * A_ikl * B_lj, C_lj = B_kj * B_ij * A_ikl
                            if any_equal && !all_equal
                                C[j, i] += C_ij
                                C[j, k] += C_kj
                                C[j, l] += C_lj
                            end
                            if all_equal
                                C[j, i] += C_ij
                            end
                        end
                    end
                end
                end
            end
        end
    end)
    println("done evaluating opt3")

    # println("evaluating opt4")
    # eval(@finch_kernel mode=fastfinch function mttkrp_opt4(C, A, B)
    #     C .= 0
    #     for l=_, k=_, i=_
    #         let A_ikl = A[i, k, l]
    #             let ik_equal = (i == k), kl_equal = (k == l)
    #             let any_equal = ik_equal || kl_equal, all_equal = ik_equal && kl_equal
    #                 if i <= k && k <= l
    #                     for j=_
    #                         let B_kj = B[j, k], B_lj = B[j, l], B_ij = B[j, i]
    #                             let C_ij = B_kj * A_ikl * B_lj, C_kj = B_ij * A_ikl * B_lj, C_lj = B_kj * B_ij * A_ikl
    #                                 if !any_equal
    #                                     C[j, i] += 2 * C_ij
    #                                     C[j, k] += 2 * C_kj
    #                                     C[j, l] += 2 * C_lj
    #                                 end
    #                                 if any_equal && !all_equal
    #                                     C[j, i] += C_ij
    #                                     C[j, k] += C_kj
    #                                     C[j, l] += C_lj
    #                                 end
    #                                 if all_equal
    #                                     C[j, i] += C_ij
    #                                 end
    #                             end
    #                         end
    #                     end
    #                 end
    #             end
    #             end
    #         end
    #     end
    # end)
    # println("done evaluating opt4")


    ref = Tensor(Dense(Dense((Element(0)))), zeros(n, n))
    ref_nondiag = Tensor(Dense(Dense((Element(0)))), zeros(n, n))
    ref_diag = Tensor(Dense(Dense((Element(0)))), zeros(n, n))
    # _ref = Tensor(Dense(Dense(Dense(Dense((Element(0)))))), zeros(n, n, n, n))
    @btime mttkrp_ref($ref, $A, $B)
    @btime mttkrp_ref($ref_nondiag, $A_nondiag, $B)
    @btime mttkrp_ref($ref_diag, $A_diag, $B)

    B_T = Tensor(Dense(Dense(Element(0))), zeros(Int, n, n))
    @finch for j=_, i=_
        B_T[i, j] = B[j, i]
    end
    @btime mttkrp_opt1($C, $A_nondiag, $B_T)
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= (2 * C[j, i]) == ref_nondiag[i, j]
    end
    @info "check" check[]

    # C_jk = Scalar(0)
    # C_jl = Scalar(0)
    # @btime mttkrp_opt1_1($C, $C_jk, $C_jl, $A, $B_T)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= (2 * C[j, i]) == ref[i, j]
    # end
    # @info "check" check[]

    # @btime mttkrp_opt1_2($C, $C_jk, $C_jl, $A, $B)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= (2 * C[i, j]) == ref[i, j]
    # end
    # @info "check" check[]

    # @btime mttkrp_opt1_3($C, $C_jk, $C_jl, $A, $B)
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

    # @btime mttkrp_opt3($C, $A, $B_T)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= C[j, i] == ref[i, j]
    # end
    # @info "check" check[]

    C_diagonals = Tensor(Dense(Dense(Element(0))), zeros(n, n))
    @btime mttkrp_opt3_1($C_diagonals, $A_diag, $B_T)
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= C_diagonals[j, i] == ref_diag[i, j]
    end
    @info "check" check[]

    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= ref_diag[i, j] + ref_nondiag[i, j] == ref[i, j]
    end
    @info "check" check[]

    # @btime mttkrp_opt4($C, $A, $B_T)
    # check = Scalar(true)
    # @finch for j=_, i=_
    #     check[] &= C[j, i] == ref[i, j]
    # end
    # @info "check" check[]
