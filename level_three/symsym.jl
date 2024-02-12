using Finch
using BenchmarkTools

n = 10
triA = rand(Int, n, n)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
triB = rand(Int, n, n)
symB = [triB[sort([i, j])...] for i = 1:n, j = 1:n]

A = Tensor(Dense(Dense(Element(0))), symA)    
B = Tensor(Dense(Dense(Element(0))), symB) 
C = Tensor(Dense(Dense(Element(0))), zeros(n, n))
ref = Tensor(Dense(Dense(Element(0))), zeros(n, n))

eval(@finch_kernel mode=fastfinch function symsym_ref(C, A, B)
    C .= 0
    for j=_, k=_, i=_
        C[i, j] += A[i, k] * B[k, j]
    end
end)

eval(@finch_kernel mode=fastfinch function symsym_opt1(C, A, B)
    C .= 0
    for j=_, k=_, i=_
        if i <= k
            C[i, j] += A[i, k] * B[k, j]
        end
        if i < k 
            C[k, j] += A[i, k] * B[i, j]
        end
    end
end)

eval(@finch_kernel mode=fastfinch function symsym_opt2(C, A, B)
    C .= 0
    for j=_, k=_, i=_
        if i <= k && k <= j
            if i != k && k != j
                C[i, j] += A[i, k] * B[k, j]
                C[k, j] += A[k, i] * B[i, j]
                C[i, k] += A[i, j] * B[j, k]
                # what is the logic as to why we don't include the below case??
                # C[k, k] += A[k, i] * B[j, k] 
            end
            if i == k && k != j
                C[i, j] += A[i, k] * B[k, j]
                C[i, k] += A[i, j] * B[j, k]
            end
            if i != k && k == j
                C[i, j] += A[i, k] * B[k, j]
                C[k, j] += A[k, i] * B[i, j]
            end
            if i == k && k == j
                C[i, j] += A[i, k] * B[k, j]
            end
        end
    end
end)


function main()
    @btime symsym_ref($ref, $A, $B)

    @btime symsym_opt1($C, $A, $B)
    @info "check opt1" C == ref

    check = Scalar(true)
    @btime symsym_opt2($C, $A, $B)
    @finch for j=_, i=_
        if i <= j
            check[] &= C[i, j] == ref[i, j]
        end
    end
    @info "check opt2" check[]
end

main()