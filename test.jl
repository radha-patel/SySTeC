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