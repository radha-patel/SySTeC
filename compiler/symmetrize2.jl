using Finch
using Finch: freshen, JuliaContext
using Finch.FinchNotation
using RewriteTools
using RewriteTools.Rewriters
using Combinatorics
using AbstractTrees

import Finch.FinchNotation: and, or

#TODO: figure out how to handle having repeats of the same symmetric matrix, but with different indices 
"""
    get_permutable_idxs(rhs, issymmetric)

Returns list of lists, where each list consists of indices that are in the
same symmetry group and can be permuted with each other.
"""
function get_permutable_idxs(rhs, issymmetric)
    permutable = Dict()
    Postwalk(@rule access(~tn::issymmetric, ~mode, ~idxs...) => begin 
        permutable_idxs = get(permutable, tn.val, Set())
        permutable[tn.val] = union(permutable_idxs, Set(idxs))
    end)(rhs)
    return [collect(idxs) for idxs in values(permutable)]
end


# TODO: "canonical" order should maybe be loop order, instead of alphabetical order?
"""
    order_canonically(idxs)

Return list of indices sorted in alphabetical order. 
"""
order_canonically(idxs) = sort(idxs, by = i->i.val)


"""
    recursively_generate_conditions(pairs)

Given list of lists of pairs of indices, return list of all combinations of == and !=
relationships between indices in each pair.
"""
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


"""
    get_conditions(idxs)

Return list of all possible combinations of matching/non-matching indices.   
"""
function get_conditions(idxs)
    idxs = order_canonically(idxs)
    pairs = [[idxs[i], idxs[i+1]] for i in 1:length(idxs)-1]
    recursively_generate_conditions(pairs)
end


"""
    get_subsymmetry(cond)

Given an expression of equivalent and non-equivalent indices, return list of lists grouping
together indices that are equivalent (and thus are swappable in an expression).
"""
function get_subsymmetry(cond)
    subsymmetry = []
    if !(@capture cond call(and, ~conds...))
        conds = [cond]
    end
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


"""
    sort_permutations(perms)

    Sort a list of lists of indices in "canonical" (i.e. alphabetical) order.
"""
sort_permutations = perms -> sort(collect(perms), by=idxs -> [i.val for i in idxs])


"""
    get_permutations(idxs, subsymmetry)

    Returns list of the caconical (i.e. alphabetical) form of all unique permutations 
    of `idxs` based on `subsymmetry` grouping (groups of equivalent indices)
"""
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


"""
    permute_indices(ex, idxs, permutations)

    Returns a block with length(permutations) expressions where expression i has its
    indices permutated such that each index idxs[j] is replaced with index 
    permutations[i][j]
"""
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


"""
    normalize(ex, issymmetric)

    Returns normalized form of `ex` 
"""
function normalize(ex, issymmetric)
    _normalize = Rewrite(Postwalk(Chain([
        # Sort indices of symmetric matrices in canonical order (alphabetical)
        (@rule access(~tn::issymmetric, ~mode, ~idxs...) => access(tn, mode, order_canonically(idxs)...)),
        # Sort operands in canonical order
        (@rule call(~op::iscommutative, ~tns...) => call(op, sort(tns, by = tn->hash(tn))...))
    ])))
    _normalize(ex)
end


"""
    add_updates(ex, conds, permutable_idxs, issymmetric)

    Returns a block containing sieves with each possible combination of equivalent
    and non-equivalent indices as the condition and all the updates that need to be 
    applied to the output given that combination of indices.   
"""
function add_updates(ex, conds, permutable_idxs, issymmetric)
    sieves = []
    for cond in conds
        subsymmetry = get_subsymmetry(cond)
        perms = get_permutations(permutable_idxs, subsymmetry)
        updates = permute_indices(ex, permutable_idxs, perms)
        normalized_updates = normalize(updates, issymmetric)
        push!(sieves, sieve(cond, normalized_updates))
    end
    return block(sieves...)
end


"""
    symmetrize2(ex, symmetric_tns)

Rewrite ex to exploit symmetry in the tensors marked as symmetric in symmetric_tns
"""
function symmetrize2(ex, symmetric_tns)
    # helper methods
    issymmetric(tn) = tn.val in symmetric_tns

    @capture ex assign(access(~lhs, updater, ~idxs...), ~op, ~rhs)

    permutable_idxs = get_permutable_idxs(rhs, issymmetric)
    conditions = get_conditions(permutable_idxs[1])
    add_updates(ex, conditions, permutable_idxs[1], issymmetric)
end