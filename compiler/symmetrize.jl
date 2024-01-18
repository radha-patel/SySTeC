using Finch
using Finch: freshen, JuliaContext
using Finch.FinchNotation
using RewriteTools
using RewriteTools.Rewriters
using Combinatorics
using AbstractTrees

import Finch.FinchNotation: and, or

isindex(ex::FinchNode) = ex.kind === index
order_canonically(idxs) = sort(idxs, by = i->i.val)
iscommutative(op) = typeof(op.val) == typeof(*) ? true : false
get_var_name(tn, idxs) = Symbol(tn.val, "_", join([idx.val for idx in idxs]))

function triangularize(idxs)
    idxs = order_canonically(idxs)
    conditions = []
    for i in 1:length(idxs)-1
        push!(conditions, call(<=, idxs[i], idxs[i+1]))
    end
    return length(conditions) > 1 ? call(and, conditions...) : conditions[1]
end

# TODO: add condition
function permute_symmetries(ex, idxs)
    permuted_exs = []
    for perm in permutations(idxs)
        ex_2 = Rewrite(Postwalk(@rule ~idx::isindex => begin
            idx in idxs ? perm[findfirst(i -> i == idx, idxs)] : idx
        end))(ex)
        push!(permuted_exs, ex_2)
    end
    return sieve(triangularize(idxs), block(permuted_exs...))
end

function get_permutable_idxs(rhs, issymmetric)
    permutable = []
    Postwalk(@rule access(~tn::issymmetric, ~mode, ~idxs...) => push!(permutable, idxs))(rhs)
    return permutable
end

function normalize(ex, issymmetric)
    _normalize = Rewrite(Postwalk(Chain([
        # Sort indices of symmetric matrices in canonical order (alphabetical)
        (@rule access(~tn::issymmetric, ~mode, ~idxs...) => access(tn, mode, order_canonically(idxs)...)),
        # Sort operands in canonical order
        (@rule call(~op::iscommutative, ~tns...) => call(op, sort(tns, by = tn->hash(tn))...))
    ])))
    _normalize(ex)
end

function transform(ex)
    transformed = false
    prev = ex 
    while prev != ex || !transformed 
        _transform = Rewrite(Postwalk(Chain([
            # Group equivalent assignments
            (@rule block(~s1..., assign(~lhs, +, ~rhs), ~s2..., assign(~lhs, +, ~rhs), ~s3...) =>
                block(s1..., assign(lhs, +, call(*, 2, rhs)), s2..., s3...)),
            # Consolidate identical reads
            (@rule block(~s1...) => begin
                counts = Dict()
                ex = block(~s1...)
                for node in PostOrderDFS(block(~s1...)) 
                    counts[node] = get(counts, node, 0) + 1
                end
                for (node, count) in counts
                    if @capture(node, access(~tn, reader, ~idxs...)) && count > 1
                        ctx = JuliaContext()
                        var = freshen(ctx, get_var_name(tn, idxs))
                        ex = Postwalk(@rule node => var)(ex)
                        ex = define(var, access(tn, reader, idxs...), ex)
                    end
                end
                ex
            end)
        ])))
        transformed = true
        prev = ex
        ex = _transform(ex)
    end
    ex
end

function get_intermediate_output(tn, count)
    ctx = JuliaContext()
    var_name = Symbol("_", tn.val, count)
    var = freshen(ctx, var_name)
end

function find_swaps(_A, B)
    A = copy(_A)
    swaps = []
    for i in 1:length(A)
        if A[i] != B[i]
            match_idx = findfirst(x -> x == B[i], A)
            push!(swaps, (i, match_idx))
            A[i], A[match_idx] = A[match_idx], A[i]
        end
    end
    return swaps
end

function get_indices_to_replicate(lhs1, lhs2)
    @capture lhs1 access(~tn1, updater, ~idxs1...)
    @capture lhs2 access(~tn2, updater, ~idxs2...)
    find_swaps(idxs1, idxs2)
end


function transform_with_post_processing(ex)
    replicate = Dict()
    count = 1
    transformed = false
    prev = ex 
    while prev != ex || !transformed 
        _transform = Rewrite(Postwalk(Chain([
            # Group equivalent assignments
            (@rule block(~s1..., assign(~lhs1, +, ~rhs), ~s2..., assign(~lhs2, +, ~rhs), ~s3...) => begin
                @capture lhs1 access(~tn1, updater, ~idxs1...)
                @capture lhs2 access(~tn2, updater, ~idxs2...)
                swaps = find_swaps(idxs1, idxs2)
                if haskey(replicate, swaps)
                    _tn = replicate[swaps]
                else
                    _tn = get_intermediate_output(tn1, count)
                    replicate[swaps] = _tn
                    count += 1
                end
                _lhs = access(_tn, updater, idxs1...)
                block(s1..., assign(_lhs, +, rhs), s2..., s3...)
            end)
        ])))
        transformed = true
        prev = ex
        ex = _transform(ex)
    end
    ex, replicate
end

# TODO: account for partial symmetry
function find_symmetry(ex, symmetric_tns)
    @capture ex assign(access(~lhs, updater, ~idxs...), ~op, ~rhs)

    # Helper methods
    issymmetric(tn) = tn.val in symmetric_tns

    permutable_idxs = get_permutable_idxs(rhs, issymmetric)
    # TODO: assumption that there is one symmetric matrix
    permuted = permute_symmetries(ex, permutable_idxs[1])
    normalized = normalize(permuted, issymmetric)
    transformed = transform(normalized)
    ex, replicate = transform_with_post_processing(transformed)
    display(ex)
    print(replicate)
end
