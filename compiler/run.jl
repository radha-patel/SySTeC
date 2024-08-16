include("symmetrize.jl")

function generate_code()
    y = :y
    x = :x
    A = :A
    B = :B 
    C = :C

    i = index(:i)
    j = index(:j)
    k = index(:k)
    l = index(:l)
    m = index(:m)
    n = index(:n)

    ex = @finch_program y[i] += A[i, j] * x[j]
    func_name = "ssymv_finch_opt_helper"
    symmetric_tns = [A]
    loop_order = [i, j]
    filename = "ssymv.jl"
    execute(ex, func_name, symmetric_tns, loop_order, filename)

    ex = @finch_program y[] += x[i] * A[i, j] * x[j]
    func_name = "syprd_finch_opt_helper"
    symmetric_tns = [A]
    loop_order = [i, j]
    filename = "syprd.jl"
    execute(ex, func_name, symmetric_tns, loop_order, filename)

    ex = @finch_program C[i, j] += A[i, k] * A[j, k]
    func_name = "ssyrk_finch_opt_helper"
    symmetric_tns = []
    loop_order = [i, j, k]
    filename = "ssyrk.jl"
    execute(ex, func_name, symmetric_tns, loop_order, filename)

    ex = @finch_program C[i, j] += A[i, k] * B[k, j]
    func_name = "ssymm_finch_opt_helper"
    symmetric_tns = [A]
    loop_order = [j, i, k]
    filename = "ssymm.jl"
    execute(ex, func_name, symmetric_tns, loop_order, filename)

    ex = @finch_program C[i, j, l] += A[k, j, l] * B[k, i]
    func_name = "ttm_finch_opt_helper"
    symmetric_tns = [A]
    loop_order = [i, j, k, l]
    filename = "ttm.jl"
    execute(ex, func_name, symmetric_tns, loop_order, filename)

    ex = @finch_program C[i, j] += A[i, k, l] * B[l, j] * B[k, j]
    func_name = "mttkrp_dim3_finch_opt_helper"
    symmetric_tns = [A]
    loop_order = [j, i, k, l]
    filename = "mttkrp_dim3.jl"
    execute(ex, func_name, symmetric_tns, loop_order, filename)

    ex = @finch_program C[i, j] += A[i, k, l, m] * B[l, j] * B[k, j] * B[m, j]
    func_name = "mttkrp_dim4_finch_opt_helper"
    symmetric_tns = [A]
    loop_order = [j, i, k, l, m]
    filename = "mttkrp_dim4.jl"
    execute(ex, func_name, symmetric_tns, loop_order, filename)

    ex = @finch_program C[i, j] += A[i, k, l, m, n] * B[l, j] * B[k, j] * B[m, j] * B[n, j]
    func_name = "mttkrp_dim5_finch_opt_helper"
    symmetric_tns = [A]
    loop_order = [j, i, k, l, m, n]
    filename = "mttkrp_dim5.jl"
    execute(ex, func_name, symmetric_tns, loop_order, filename)
end

generate_code()