using Finch
using BenchmarkTools

n = 100000
triA = fsprand(Int, n, n, 0.1)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]

A = Tensor(Dense(SparseList(Element(0))), symA)
B = Tensor(Dense(Element(0)), rand(Int, n))
C = Tensor(Dense(Element(0)), zeros(n))
temp = Scalar(0)

# ~16ms
eval(@finch_kernel mode=fastfinch function symv_ref(C, A, B)
    C .= 0
    for j=_, i=_
        C[i] += A[i, j] * B[j]
    end
end)

# ~11ms
eval(@finch_kernel mode=fastfinch function symv_gen(C, A, B)
    C .= 0
    for j=_, i=_
        let A_ij = A[i, j]
            if i <= j
                C[i] += A_ij * B[j]
            end
            if i < j
                C[j] += A_ij * B[i]
            end
        end
    end
end)

# ~10ms
eval(@finch_kernel mode=fastfinch function symv_opt(C, A, B, temp)
    C .= 0
    for j=_
        temp .= 0
        for i=_
            let A_ij = A[i, j]
                if i <= j
                    C[i] += A_ij * B[j]
                end
                if i < j
                    temp[] += A_ij * B[i]
                end
            end
        end
        C[j] += temp[]
    end
end)

function main()
    ref = Tensor(Dense(Element(0)), zeros(n))
    @btime(symv_ref($ref, $A, $B))

    @btime(symv_gen($C, $A, $B))
    @info "check generated code" C == ref

    @btime(symv_opt($C, $A, $B, $temp))
    @info "check optimized code" C == ref
end

main()
