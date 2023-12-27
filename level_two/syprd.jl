using Finch
using BenchmarkTools
using SparseArrays

n = 10000
triA = rand(Int, n, n)
symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
x = rand(Int, n)

A = Fiber!(Dense(SparseList(Element(0))), symA)  
X = Fiber!(Dense(Element(0)), x)
C = Scalar(0)
D = Scalar(0)

eval(@finch_kernel mode=fastfinch function syprd_ref(C, A, X)
    C .= 0
    for j=_, i=_ 
        C[] += X[i] * A[i, j] * X[j] 
    end
end)

eval(@finch_kernel mode=fastfinch function syprd_opt1(C, D, A, X)
    C .= 0
    D .= 0
    for j=_ 
        let X_j = X[j]
            for i=_ 
                let a = A[i, j]
                    if i < j 
                        C[] += X[i] * a * X_j
                    end
                    if i == j
                        D[] += X[i] * a * X_j
                    end
                end
            end
        end
    end
end)



function main()
    ref = Scalar(0)
    @btime syprd_ref($ref, $A, $X)

    C = Scalar(0)
    D = Scalar(0)
    @btime syprd_opt1($C, $D, $A, $X)

    @info ref[] == 2 * C[] + D[]
end

main()