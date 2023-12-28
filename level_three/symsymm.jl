using Finch
using BenchmarkTools
using SparseArrays

n = 10

triA = rand(Int, n, n)
# triA = fill(1, (n, n))
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
# triB = rand(Int, n, n)
triB = fill(1, (n, n))
symB = [triB[sort([i, j])...] for i = 1:n, j = 1:n]

A = Fiber!(Dense(Dense(Element(0))), symA)    
B = Fiber!(Dense(Dense(Element(0))), symB) 
# C = Fiber!(Dense(Dense(Element(0))), rand(Int, n, n))    
C = Fiber!(Dense(Dense(Element(0))), fill(1, (n, n)))    
D = Fiber!(Dense(Dense(Element(0))), zeros(n, n))
T = Fiber!(Dense(Dense(Element(0))), zeros(n, n))

eval(@finch_kernel mode=fastfinch function symsymm_ref1(D, A, B, C, T)
    T .= 0
    for j = _, k = _, i = _
        T[i, j] += A[i, k] * B[k, j]
    end
    D .= 0
    for j = _, k = _, i = _
        D[i, j] += T[i, k] * C[k, j]
    end
end)

eval(@finch_kernel mode=fastfinch function symsymm_ref2(D, A, B, C)
    D .= 0
    for j = _, k = _, l = _, i = _
        D[i, j] += A[i, k] * B[k, l] * C[l, j]
    end
end)

eval(@finch_kernel mode=fastfinch function symsymm_opt1(D, A, B, C)
    D .= 0
    for j = _, k = _, l = _, i = _
        if i >= l 
            D[i, j] += A[i, k] * B[k, l] * C[l, j]
        end
        if i > l
            D[l, j] += A[l, k] * B[k, i] * C[i, j]
        end
    end
end)

function main()
    ref = Fiber!(Dense(Dense(Element(0))), zeros(n, n))
    @btime symsymm_ref1($ref, $A, $B, $C, $T)

    @btime symsymm_ref2($D, $A, $B, $C)
    @info "check ref2" D == ref

    @btime symsymm_opt1($D, $A, $B, $C)
    @info "check opt1" D == ref
end

main()