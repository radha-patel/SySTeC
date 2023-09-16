using Finch
using BenchmarkTools

A = Fiber!(Dense(SparseList(Element(0.0))), fsprand((10000, 10000), 0.1));
x = Fiber!(Dense(Element(0.0)), rand(10000));
y = Fiber!(Dense(Element(0.0)), zeros(10000));

eval(@finch_kernel function symv(y, A, x) 
    for j = _, i = _
        y[i] += A[i, j] * x[j]
    end
end)

@btime symv($y, $A, $x)