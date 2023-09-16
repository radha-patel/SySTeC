using Finch
using BenchmarkTools

A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1));
x = Fiber!(Dense(Element(0.0)), rand(10000));
y = Fiber!(Dense(Element(0.0)), zeros(10000));
temp2 = Scalar(0.0)

eval(@finch_kernel function ssymv(y, A, x, temp2)
    for j = _
        temp1 = x[j]
        temp2 .= 0
        for i = _
            temp3 = A[i, j]
            if uptrimask[i, j]
                y[i] += temp1 * temp3
            end
            if uptrimask[i, j - 1]
                temp2[] += temp3 * x[i]
            end
        end
        y[j] += temp2[]
    end
end)

@btime ssymv($y, $A, $x, $temp2)