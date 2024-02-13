using Finch
using BenchmarkTools

n = 100
triA = fsprand(Int, n, n, n, 0.1)
symA = [triA[sort([i, j, k])...] for i = 1:n, j = 1:n, k = 1:n]
b = rand(Int, n, n)

A = Tensor(Dense(SparseList(SparseList(Element(0)))), symA)
B = Tensor(Dense(Dense(Element(0))), b)    
B_T = Tensor(Dense(Dense(Element(0))), transpose(b)) 
C = Tensor(Dense(Dense(Dense(Element(0)))), zeros(Int, n, n, n))

# ~7.9ms
eval(@finch_kernel mode=fastfinch function ttm_ref(C, A, B)
    C .= 0
    for l=_, j=_, k=_, i=_
        C[i, j, l] += A[k, j, l] * B[k, i]
    end
end)

# ~3.6ms
eval(@finch_kernel mode=fastfinch function ttm_gen(C, A, B)
    C .= 0
    for l=_, k=_, j=_, i=_
        let jk_leq = (j <= k), kl_leq = (k <= l)
            let A_jkl = A[j, k, l]
                if jk_leq && kl_leq
                    C[i, j, k] += A_jkl * B[l, i]
                end
                if j < k && kl_leq
                    C[i, k, l] += A_jkl * B[j, i]
                end
                if jk_leq && k < l
                    C[i, j, l] += A_jkl * B[k, i]
                end
            end
        end
    end
end)

# ~3.0ms
eval(@finch_kernel mode=fastfinch function ttm_opt(C, A, B_T)
    C .= 0
    for l=_, k=_, j=_, i=_
        let jk_leq = (j <= k), kl_leq = (k <= l)
            let A_jkl = A[j, k, l]
                if jk_leq && kl_leq
                    C[i, j, k] += A_jkl * B_T[i, l]
                end
                if j < k && kl_leq
                    C[i, k, l] += A_jkl * B_T[i, j]
                end
                if jk_leq && k < l
                    C[i, j, l] += A_jkl * B_T[i, k]
                end
            end
        end
    end
end)

function main()
    ref = Tensor(Dense(Dense(Dense(Element(0)))), zeros(Int, n, n, n))
    @btime(ttm_ref($ref, $A, $B))

    @btime(ttm_gen($C, $A, $B))
    check = Scalar(true)
    @finch for l=_, j=_, i=_
        if j <= l
            check[] &= C[i, j, l] == ref[i, j, l]
        end
    end
    @info "check generated code" check[]


    @btime(ttm_opt($C, $A, $B_T))
    check = Scalar(true)
    @finch for l=_, j=_, i=_
        if j <= l
            check[] &= C[i, j, l] == ref[i, j, l]
        end
    end
    @info "check generated code" check[]
end

main()