# @finch begin
#     for k=_, j=_, l=_, i=_ 
#         C[i, j, k] = A[i, l, k] * X[j, l] 
#     end
# end

@finch begin
    for k=_, l=_, j=_, i=_ 
        C[i, j, l] = A[i, k, l] * X[k, j] 
    end
end

A = Fiber!(Dense(Dense(SparseList(Element(0.0)))), fsprand((10, 10, 10), 0.1))
C = Fiber!(Dense(Dense(Dense(Element(0.0)))), zeros((10, 10, 10)))
X = Fiber!(Dense(Dense(Element(0.0))), rand(10, 10))
temp2 = Scalar(0.0)
@finch begin
    for l=_, j=_, i=_
        temp2 .= 0
        temp1 = X[i, j]
        for k=_ 
            # temp3 = A[i, l, k]
            if uptrimask[k+1, i]
                C[k, j, l] += temp1 * A[k, i, l]
            end
            if uptrimask[k, l]
                temp2[] += X[k, j] * A[k, i, l]
            end
        end
        C[i, j, l] += temp2[]
    end
end