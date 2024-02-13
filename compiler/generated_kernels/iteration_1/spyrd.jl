using Finch
using BenchmarkTools

n = 10000
triA = fsprand(Int, n, n, 0.1)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]

A = Tensor(Dense(SparseList(Element(0))), symA)
B = Tensor(Dense(Element(0)), rand(Int, n))
C = Scalar(0)

# ~18ms
eval(@finch_kernel mode=fastfinch function syprd_ref(C, A, B)
    C .= 0
    for j=_, i=_
        C[] += B[i] * A[i, j] * B[j]
    end
end)

# ~8.5ms
eval(@finch_kernel mode=fastfinch function syprd_gen(C, A, B)
    C .= 0
    for j=_, i=_
        let B_i = B[i], A_ij = A[i, j], B_j = B[j]
            if i < j
                C[] += 2 * A_ij * B_i * B_j
            end
            if i == j
                C[] += A_ij * B_i * B_j
            end
        end
    end
end)

# ~8.5ms
eval(@finch_kernel mode=fastfinch function syprd_opt(C, A, B)
    C .= 0
    for j=_, i=_
        let B_i = B[i], A_ij = A[i, j], B_j = B[j]
            let prod = A_ij * B_i * B_j
                if i < j
                    C[] += 2 * prod
                end
                if i == j
                    C[] += prod
                end
            end
        end
    end
end)

function main()
    ref = Scalar(0)
    @btime(syprd_ref($ref, $A, $B))

    @btime(syprd_gen($C, $A, $B))
    @info "check generated code" C[] == ref[]

    @btime(syprd_opt($C, $A, $B))
    @info "check optimized code" C[] == ref[]
end

main()
