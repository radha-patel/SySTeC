using Finch
using BenchmarkTools

for n in (10, 100, 1000)
    println("N = ", n)
    print("Copy dense to dense: ")
    c = Fiber!(Dense(Dense(Element(0.0))), fsprand((n, n), 0.1))
    t = Fiber!(Dense(Dense(Element(0.0))), zeros(n, n))
    eval(@finch_kernel mode=:fast function copy_transpose1(c, t) 
        for j = _, i = _ 
            t[j, i] = c[i, j] 
        end
    end)
    @btime copy_transpose1($c, $t)

    print("Copy sparse to dense: ")
    c = Fiber!(Dense(SparseList(Element(0.0))), fsprand((n, n), 0.1))
    t = Fiber!(Dense(Dense(Element(0.0))), zeros(n, n))
    eval(@finch_kernel mode=:fast function copy_transpose2(c, t)
        for j = _, i = _ 
            t[j, i] = c[i, j] 
        end
    end)
    @btime copy_transpose2($c, $t)

    print("Copy sparse to sparse: ")
    c = Fiber!(Dense(SparseList(Element(0.0))), fsprand((n, n), 0.1))
    t = Fiber!(Dense(SparseList(Element(0.0))), zeros(n, n))
    eval(@finch_kernel mode=:fast function copy_transpose3(c, t)
        t .= 0
        for j = _, i = _ 
            t[j, i] = c[i, j] 
        end
    end)
    @btime copy_transpose3($c, $t)

    print("Copy hash to hash: ")
    c = Fiber!(SparseHash{2}(Element(0.0)), fsprand((n, n), 0.1))
    t = Fiber!(SparseHash{2}(Element(0.0)), zeros(n, n))
    eval(@finch_kernel mode=:fast function copy_transpose4(c, t)
        for j = _, i = _ 
            t[j, i] = c[i, j] 
        end
    end)
    @btime copy_transpose4($c, $t)
    println("")
end