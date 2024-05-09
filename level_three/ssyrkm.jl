using Finch
using BenchmarkTools
using SparseArrays

n = 100

triA = rand(Int, n, n)
# triA = fill(1, (n, n))
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]

A = Fiber!(Dense(Dense(Element(0))), symA)    
B = Fiber!(Dense(Dense(Element(0))), rand(Int, n, n)) 
C = Fiber!(Dense(Dense(Element(0))), zeros(n, n)) 
T = Fiber!(Dense(Dense(Element(0))), zeros(n, n))    

eval(@finch_kernel mode=:fast function ssyrkm_ref1(C, A, B, T)
    T .= 0
    for j = _, k = _, i = _
        T[i, j] += A[i, k] * A[k, j]
    end
    C .= 0
    for j = _, k = _, i = _
        C[i, j] += T[i, k] * B[k, j]
    end
end)

eval(@finch_kernel mode=:fast function ssyrkm_ref2(C, A, B)
    C .= 0
    for j = _, k = _, l = _, i = _
        C[i, j] += A[i, k] * A[k, l] * B[l, j]
    end
end)

eval(@finch_kernel mode=:fast function ssyrkm_opt1(C, A, B)
    C .= 0
    for j = _, k = _, l = _, i = _
        let temp = A[i, k] * A[k, l]
            if i >= l 
                C[i, j] += temp * B[l, j]
            end
            if i > l
                C[l, j] += temp * B[i, j]
            end
        end
    end
end)

function main()
    ref = Fiber!(Dense(Dense(Element(0))), zeros(n, n))
    @btime ssyrkm_ref1($ref, $A, $B, $T)

    @btime ssyrkm_ref2($C, $A, $B)
    @info "check ref2" C == ref

    @btime ssyrkm_opt1($C, $A, $B)
    @info "check opt1" C == ref
end

main()