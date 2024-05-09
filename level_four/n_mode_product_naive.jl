using Finch
using BenchmarkTools


n = 100
triA = rand(Int, n, n, n)
symA = [triA[sort([i, j, k])...] for i = 1:n, j = 1:n, k = 1:n]
A = Fiber!(Dense(SparseList(SparseList(Element(0)))), symA)

# n = 3
# A = zeros(n, n, n)
# # A[:,:,1] = [1 2 3; -1 0 -2; 2 1 3]
# # A[:,:,2] = [2 1 3; 1 2 3; -1 -2 -3]
# # A[:,:,3] = [0 1 2; 3 -1 0; 2 2 1]
# A[:,:,1] = [0 0 0; 0 4 -2; 0 -2 2]
# A[:,:,2] = [0 4 -2; 4 0 0; -2 0 2]
# A[:,:,3] = [0 -2 2; -2 0 4; 2 4 0]

# A = Fiber!(Dense(Dense(SparseList(Element(0.0)))), fsprand((n, n, n), 0.1))
C = Fiber!(Dense(Dense(Dense(Element(0)))), zeros((n, n, n)))
X = Fiber!(Dense(Dense(Element(0))), rand(Int, n, n))
# X = [1 2 3; -1 2 3; -3 -2 1]
# TIME: ~35ms
eval(@finch_kernel mode=:fast function mode2_product_naive(C, A, X)
    C .= 0
    for l=_, j=_, k=_, i=_ 
        C[i, j, l] += A[i, k, l] * X[k, j] 
    end
end)

eval(@finch_kernel mode=:fast function mode1_product_naive(C, A, X)
    C .= 0
    for l=_, j=_, k=_, i=_ 
        C[i, j, l] += A[k, j, l] * X[k, i] 
    end
end)

eval(@finch_kernel mode=:fast function mode3_product_naive(C, A, X)
    C .= 0
    for l=_, k=_, j=_, i=_ 
        C[i, j, l] += A[i, j, k] * X[k, l] 
    end
end)

mode1_C = Fiber!(Dense(Dense(Dense(Element(0), n), n), n))
@btime mode1_product_naive($mode1_C, $A, $X)

mode2_C = Fiber!(Dense(Dense(Dense(Element(0), n), n), n))
@btime mode2_product_naive($mode2_C, $A, $X)

mode3_C = Fiber!(Dense(Dense(Dense(Element(0), n), n), n))
@btime mode3_product_naive($mode3_C, $A, $X)

check = Scalar(true)
check = (mode1_C == mode2_C) && (mode2_C == mode3_C)
@info "check" check[]