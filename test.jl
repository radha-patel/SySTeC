using Test

function files_are_equal(file1::String, file2::String)::Bool
    content1 = open(read, file1)
    content2 = open(read, file2)
    return content1 == content2
end

if files_are_equal("tests/ssymv_correct.jl", "ssymv.jl")
    println("[SUCCESS] SSYMV generated correctly")
else
    println("[FAIL] SSYMV is incorrect")
end

if files_are_equal("tests/ssyrk_correct.jl", "ssyrk.jl")
    println("[SUCCESS] SSYRK generated correctly")
else
    println("[FAIL] SSYRK is incorrect")
end

if files_are_equal("tests/ssymm_correct.jl", "ssymm.jl")
    println("[SUCCESS] SSYMM generated correctly")
else
    println("[FAIL] SSYMM is incorrect")
end

if files_are_equal("tests/ttm_correct.jl", "ttm.jl")
    println("[SUCCESS] TTM generated correctly")
else
    println("[FAIL] TTM is incorrect")
end

if files_are_equal("tests/mttkrp_dim3_correct.jl", "mttkrp_dim3.jl")
    println("[SUCCESS] MTTKRP Dim 3 generated correctly")
else
    println("[FAIL] MTTKRP Dim 3 is incorrect")
end

if files_are_equal("tests/mttkrp_dim4_correct.jl", "mttkrp_dim4.jl")
    println("[SUCCESS] MTTKRP Dim 4 generated correctly")
else
    println("[FAIL] MTTKRP Dim 4 is incorrect")
end

if files_are_equal("tests/mttkrp_dim5_correct.jl", "mttkrp_dim5.jl")
    println("[SUCCESS] MTTKRP Dim 5 generated correctly")
else
    println("[FAIL] MTTKRP Dim 5 is incorrect")
end