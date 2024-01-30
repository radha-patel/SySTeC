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


function triangularize(idxs, strict)
    idxs = order_canonically(idxs)
    conditions = []
    for i in 1:length(idxs)-1
        op = strict[i] == 1 ? (<) : (<=)
        push!(conditions, call(op, idxs[i], idxs[i+1]))
    end
    return length(conditions) > 1 ? call(and, conditions...) : conditions[1]
end

function make_strict(cond)
    Rewrite(Postwalk(Chain([
        (@rule call(<=, ~a, ~b) => call(<, a, b))
    ])))(cond)
end

# TODO: add condition
function permute_symmetries(ex, idxs, diagonals)
    permuted_exs = []
    for perm in permutations(idxs)
        ex_2 = Rewrite(Postwalk(@rule ~idx::isindex => begin
            idx in idxs ? perm[findfirst(i -> i == idx, idxs)] : idx
        end))(ex)
        push!(permuted_exs, ex_2)
    end
    return diagonals ? block(permuted_exs...) : sieve(triangularize(idxs), block(permuted_exs...))
end

function permute_indices(ex, idxs, permutations)
    permuted_exs = []
    for perm in permutations
        ex_2 = Rewrite(Postwalk(@rule ~idx::isindex => begin
            idx in idxs ? perm[findfirst(i -> i == idx, idxs)] : idx
        end))(ex)
        push!(permuted_exs, ex_2)
    end
    return block(permuted_exs...)
end

function get_permutable_idxs(rhs, issymmetric)
    permutable = Dict()
    Postwalk(@rule access(~tn::issymmetric, ~mode, ~idxs...) => begin 
        permutable_idxs = get(permutable, tn.val, Set())
        permutable[tn.val] = union(permutable_idxs, Set(idxs))
    end)(rhs)
    return [collect(idxs) for idxs in values(permutable)]
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

function transform(ex, diagonals)
    _transform = Rewrite(Postwalk(Chain([
        # Group equivalent assignments
        (@rule block(~s1..., assign(~lhs, +, ~rhs), ~s2..., assign(~lhs, +, ~rhs), ~s3...) =>
            block(s1..., assign(lhs, +, call(*, 2, rhs)), s2..., s3...)),
        # TODO: move this to separate step
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
    Fixpoint(_transform)(ex)
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
    Fixpoint(_transform)(ex), replicate
end

function add_conditions(ex, idxs)
    Rewrite(Postwalk(@rule block(~s...) => begin
        sieves = []
        @assert length(s) % length(idxs) == 0 || length(s) == 1
        if length(s) == 1
            clause = sieve(triangularize(idxs), block(s[1]))
            push!(sieves, clause)
        else
            n = Int(length(s) / length(idxs))
            clause = sieve(triangularize(idxs), block(s[1:n]...))
            push!(sieves, clause)

            range = length(idxs) - 1
            for i in 1:range
                strict = zeros(range)
                strict[i] = 1
                clause = sieve(triangularize(idxs, strict), block(s[n*i + 1 : n*i + n]...)) 
                push!(sieves, clause)
            end
        end

        # Account for invisible output symmetry
        sieves_2 = []
        for i in 1:length(sieves)
            # TODO: generalize for multiple updates in sieve
            if @capture sieves[i] sieve(~cond, block(assign(~lhs, +, call(*, 2, ~tn))))
                s = assign(lhs, +, tn)
                push!(sieves_2, sieve(cond, s))
                push!(sieves_2, sieve(make_strict(cond), s))
            else
                push!(sieves_2, sieves[i])
            end
        end

        block(sieves_2...)
    end))(ex)
end

function recursively_generate_conditions(pairs)
    idx_1, idx_2 = pairs[1][1], pairs[1][2]
    if length(pairs) == 1
        return [call(!=, idx_1, idx_2), call(==, idx_1, idx_2)]
    elseif length(pairs) == 2
        all_conds = []
        for partial_cond in recursively_generate_conditions(pairs[2:end])
            push!(all_conds, call(and, call(!=, idx_1, idx_2), partial_cond))
            push!(all_conds, call(and, call(==, idx_1, idx_2), partial_cond))
        end
        return all_conds
    else
        all_conds = []
        for partial_cond in recursively_generate_conditions(pairs[2:end])
            @capture partial_cond call(and, ~partial_conds...)
            push!(all_conds, call(and, call(!=, idx_1, idx_2), partial_conds...))
            push!(all_conds, call(and, call(==, idx_1, idx_2), partial_conds...))
        end
        return all_conds
    end
end

function get_conditions(idxs)
    idxs = order_canonically(idxs)
    pairs = [[idxs[i], idxs[i+1]] for i in 1:length(idxs)-1]
    recursively_generate_conditions(pairs)
end

function get_subsymmetry(cond)
    subsymmetry = []
    @capture cond call(and, ~conds...)
    for cond in conds
        if @capture cond call(==, ~idx_1, ~idx_2)
            if length(subsymmetry) > 0 && idx_1 in subsymmetry[1]
                push!(subsymmetry[1], idx_2)
            else
                push!(subsymmetry, [idx_1, idx_2])
            end
        end
    end
    subsymmetry
end

sort_permutations = perms -> sort(collect(perms), by=idxs -> [i.val for i in idxs])
function get_permutations(idxs, subsymmetry)
    perms = Dict()
    for perm in permutations(idxs)
        translated_perm = deepcopy(perm)
        for group in subsymmetry
            for idx in group[2:end]
                translated_perm[findfirst(i -> i == idx, translated_perm)] = group[1]
            end
        end
        perms_2 = get(perms, translated_perm, Set())
        perms[translated_perm] = union(perms_2, Set([perm]))
    end
    return [sort_permutations(perms_2)[1] for perms_2 in values(perms)]
end

function add_updates(ex, conds, permutable_idxs)
    for cond in conds
        subsymmetry = get_subsymmetry(cond)
        perms = get_permutations(permutable_idxs, subsymmetry)
        updates = permute_indices(ex, permutable_idxs, perms)
        display(cond)
        display(updates)
    end
end

# function consolidate(ex)
#     _consolidate = Rewrite(Postwalk(Chain([
#         (@rule block(~s1..., sieve(~cond1, ~body), ~s2..., sieve(~cond2, ~body), ~s3...) => begin
#             if @capture cond1 call(and, ~c1...) && @capture cond2 call(and, ~c2...)
#                 replicated_conds = []
#                 single_conds = []
#                 for (_c1, _c2) in zip(c1, c2)
#                     @capture _c1 call(~op1, ~i1...)
#                     @capture _c2 call(~op2, ~i2...)
#                     if op1 == literal(<=) || op2 == literal(<=)
#                         push!(single_conds, call(==, i1...))
#                     else
#                         push!(single_conds, call(!=, i1...))
#                     end
#                     push!(replicated_conds, call(<, i1...))
#                 end
#             end

#             @capture body assign(~lhs, +, ~rhs)
#             clause_1 = sieve(call(and, replicated_conds...), assign(lhs, +, call(*, 2, rhs)))
#             clause_2 = sieve(call(and, single_conds...), body)
#             block(clause_1, clause_2, s1..., s2..., s3...)
#         end)
#     ])))
#     Fixpoint(_consolidate)(ex)
#     # _consolidate(ex)
# end

# TODO: account for partial symmetry
function symmetrize(ex, symmetric_tns, include_diagonals)
    @capture ex assign(access(~lhs, updater, ~idxs...), ~op, ~rhs)

    # Helper methods
    issymmetric(tn) = tn.val in symmetric_tns

    permutable_idxs = get_permutable_idxs(rhs, issymmetric)
    # TODO: assumption that there is one symmetric matrix
    original_ex = ex
    permuted = permute_symmetries(ex, permutable_idxs[1], include_diagonals)
    normalized = normalize(permuted, issymmetric)
    transformed = transform(normalized, include_diagonals)
    ex, replicate = transform_with_post_processing(transformed)
    if include_diagonals
        ex = add_conditions(ex, permutable_idxs[1])
    end
    conds = get_conditions(permutable_idxs[1])
    add_updates(original_ex, conds, [k, j, l])
    # ex = consolidate(ex)
    println("REPLICATE: ", replicate)
    display(ex)
end
