using Finch
using BenchmarkTools

n = 1000
triA = fsprand(Int, n, n, 0.1)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
b = rand(Int, n, n)

A = Tensor(Dense(SparseList(Element(0))), symA)
B = Tensor(Dense(Dense(Element(0))), b)
B_T = Tensor(Dense(Dense(Element(0))), transpose(b))
C = Tensor(Dense(Dense(Element(0))), zeros(n, n))
C_T = Tensor(Dense(Dense(Element(0))), zeros(n, n))
temp = Scalar(0)

# ~100ms
eval(@finch_kernel mode=:fast function symm_ref(C, A, B)
    C .= 0
    for j=_, k=_, i=_
        C[i, j] += A[i, k] * B[k, j]
    end
end)

# ~90ms
eval(@finch_kernel mode=:fast function symm_gen(C, A, B)
    C .= 0
    for j=_, k=_, i=_
        let A_ik = A[i, k]
            if i <= k
                C[i, j] += A_ik * B[k, j]
            end
            if i < k
                C[k, j] += A_ik * B[i, j]
            end
        end
    end
end)

# ~70ms
eval(@finch_kernel mode=:fast function symm_opt1(C, A, B, temp)
    C .= 0
    for j=_, k=_
        temp .= 0
        for i=_
            let A_ik = A[i, k]
                if i <= k
                    C[i, j] += A_ik * B[k, j]
                end
                if i < k
                    temp[] += A_ik * B[i, j]
                end
            end
        end
        C[k, j] += temp[]
    end
end)

# ~40ms
eval(@finch_kernel mode=:fast function symm_opt2(C_T, A, B_T)
    C_T .= 0
    for k=_, i=_, j=_
        let A_ik = A[i, k]
            if i <= k
                C_T[j, i] += A_ik * B_T[j, k]
            end
            if i < k
                C_T[j, k] += A_ik * B_T[j, i]
            end
        end
    end
end)

function main()
    ref = Tensor(Dense(Dense(Element(0))), zeros(n, n))
    @btime(symm_ref($ref, $A, $B))

    @btime(symm_gen($C, $A, $B))
    @info "check generated code" C == ref

    @btime(symm_opt1($C, $A, $B, $temp))
    @info "check optimized code (1)" C == ref

    @btime(symm_opt2($C_T, $A, $B_T))
    check = Scalar(true)
    @finch for j=_, i=_
        check[] &= C_T[j, i] == ref[i, j]
    end
    @info "check optimized code (2)" check[]
end

main()