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

eval(@finch_kernel mode=:fast function symsym_ref(C, A, B)
    C .= 0
    for j=_, k=_, i=_
        C[i, j] += A[i, k] * B[k, j]
    end
end)

eval(@finch_kernel mode=:fast function symsym_opt1(C, A, B)
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

eval(@finch_kernel mode=:fast function symsym_opt2(C, A, B)
    C .= 0
    for j=_, k=_, i=_
        if i <= k && k <= j
            if i != k && k != j
                C[i, j] += A[i, k] * B[k, j] 
                C[k, j] += A[k, i] * B[i, j] 
                C[i, k] += A[i, j] * B[j, k] 
                C[k, i] += A[k, j] * B[j, i]
                C[j, k] += A[j, i] * B[i, k] 
                C[j, i] += A[j, k] * B[k, i]
            end
            if i == k && k != j
                C[i, j] += A[i, k] * B[k, j]
                C[i, k] += A[i, j] * B[j, k]
                C[j, k] += A[j, i] * B[i, k] 
            end
            if i != k && k == j
                C[i, j] += A[i, k] * B[k, j]
                C[k, j] += A[k, i] * B[i, j]
                C[k, i] += A[k, j] * B[j, i]
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

    @btime symsym_opt2($C, $A, $B)
    @info "check opt2" C == ref
end

main()