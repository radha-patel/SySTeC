#!/usr/bin/env julia
if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
end

using MatrixDepot
using BenchmarkTools
using ArgParse
using DataStructures
using JSON
using SparseArrays
using Printf
using LinearAlgebra

s = ArgParseSettings("Benchmark SSYMV experiments.")

@add_arg_table! s begin
    "--output", "-o"
        arg_type = String
        help = "output file path"
        default = "spmv_results.json"
    "--dataset", "-d"
        arg_type = String
        help = "dataset keyword"
        default = "oski"
end

parsed_args = parse_args(ARGS, s)

datasets = Dict(
    "oski" => [
        "Boeing/ct20stif",
        "Simon/olafu",
        "Boeing/bcsstk35",
        "Boeing/crystk02",
        "Boeing/crystk03",
        "Nasa/nasasrb",
        "Rothberg/3dtube",
        "Simon/raefsky4",
        "Mulvey/finan512",
        "Pothen/pwt",
        "Cote/vibrobox",
        "HB/saylr4",
        "Gupta/gupta1"
    ],
)

sizes = [10, 100, 1000, 10000, 20000]
sparsity = [0.1]

include("../level_two/ssymv.jl")

results = []

for n in sizes
    triA = fsprand(Int, (n, n), 0.1)
    symA = [triA[sort([i, j])...] for i = 1:n, j = 1:n]
    A = Fiber!(Dense(SparseList(Element(0))), symA);
    x = Fiber!(Dense(Element(0)), rand(Int, n));
    y = Fiber!(Dense(Element(0)), zeros(n));
    ref = Fiber!(Dense(Element(0)), zeros(n));

# for mtx in datasets[parsed_args["dataset"]]
#     A = Fiber!(Dense(SparseList(Element(0.0))), matrixdepot(mtx));
#     (n, n) = size(A)
#     x = Fiber!(Dense(Element(0.0)), rand(n));
#     y = Fiber!(Dense(Element(0.0)), zeros(n));
#     ref = Fiber!(Dense(Element(0.0)), zeros(n))

    for (key, method) in [
        "ref" => ssymv_ref_helper,
        "opt" => ssymv_opt_helper
    ] 
        @info "testing" key n
        if key == "ref"
            res = method(ref, A, x)
        else
            res = method(y, A, x)
        end

        time = res.time
        
        check = Scalar(true)
        @finch for i=_
            check[] &= ref[i] == y[i]
        end
        # check = (norm(ref - y)/norm(ref) < 0.1)
        @info "check" key == "ref" || check[]

        @info "results" time
        push!(results, OrderedDict(
            "time" => time,
            "method" => key,
            "kernel" => "spmv",
            "size" => n,
            # "mtx" => mtx,
        ))
        write(parsed_args["output"], JSON.json(results, 4))
    end
end