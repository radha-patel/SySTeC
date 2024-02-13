using Finch
using BenchmarkTools

n = 1000
triA = fsprand(Int, n, n, 0.1)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]

A = Tensor(Dense(SparseList(Element(0))), symA)
B = Tensor(Dense(Dense(Element(0))), rand(Int, n, n))
C = Tensor(Dense(Dense(Element(0))), zeros(n, n))
temp = Scalar(0)

# ~125ms
eval(@finch_kernel mode=fastfinch function symm_ref(C, A, B)
    C .= 0
    for j=_, k=_, i=_
        C[i, j] += A[i, k] * B[k, j]
    end
end)

# ~100ms
eval(@finch_kernel mode=fastfinch function symm_gen(C, A, B)
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
eval(@finch_kernel mode=fastfinch function symm_opt(C, A, B, temp)
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

function main()
    ref = Tensor(Dense(Dense(Element(0))), zeros(n, n))
    @btime(symm_ref($ref, $A, $B))

    @btime(symm_gen($C, $A, $B))
    @info "check generated code" C == ref

    @btime(symm_opt($C, $A, $B, $temp))
    @info "check optimized code" C == ref
end

main()