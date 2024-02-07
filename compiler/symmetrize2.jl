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
function recursively_generate_conditions(pairs, diagonals=true)
    idx_1, idx_2 = pairs[1][1], pairs[1][2]
    if length(pairs) == 1
        return diagonals ? [call(!=, idx_1, idx_2), call(==, idx_1, idx_2)] : [call(!=, idx_1, idx_2)]
    elseif length(pairs) == 2
        all_conds = []
        for partial_cond in recursively_generate_conditions(pairs[2:end], diagonals)
            push!(all_conds, call(and, call(!=, idx_1, idx_2), partial_cond))
            if diagonals
                push!(all_conds, call(and, call(==, idx_1, idx_2), partial_cond))
            end
        end
        return all_conds
    else
        all_conds = []
        for partial_cond in recursively_generate_conditions(pairs[2:end], diagonals)
            @capture partial_cond call(and, ~partial_conds...)
            push!(all_conds, call(and, call(!=, idx_1, idx_2), partial_conds...))
            if diagonals
                push!(all_conds, call(and, call(==, idx_1, idx_2), partial_conds...))
            end
        end
        return all_conds
    end
end


"""
    get_conditions(idxs)

Return list of all possible combinations of matching/non-matching indices.   
"""
function get_conditions(idxs, diagonals=true)
    idxs = order_canonically(idxs)
    pairs = [[idxs[i], idxs[i+1]] for i in 1:length(idxs)-1]
    recursively_generate_conditions(pairs, diagonals)
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

Returns list of the canonical (i.e. alphabetical) form of all unique permutations 
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


#TODO: this step can probably be combined/happen within add_updates 
#TODO: enumerate all permutations of b1_s and b2_s and check if pairs can be made the same
function group_sieves(ex)
    Rewrite(Postwalk(@rule block(~s1..., sieve(~c1, ~b1), ~s2..., sieve(~c2, ~b2), ~s3...) => begin
        @capture b1 block(~b1_s...)
        @capture b2 block(~b2_s...)

        if length(b1_s) != length(b2_s)
            return nothing
        end

        # rewrite b1 and b2 using only canonical index of equivalent indices
        b1_subsymmetry = get_subsymmetry(c1)
        b2_subsymmetry = get_subsymmetry(c2)

        b1_2 = b1
        for group in b1_subsymmetry
            for idx in group[2:end]
                b1_2 = Rewrite(Postwalk(@rule ~idx_2::isindex => idx_2 == idx ? group[1] : idx_2))(b1_2)
            end
        end

        b2_2 = b2
        for group in b2_subsymmetry
            for idx in group[2:end]
                b2_2 = Rewrite(Postwalk(@rule ~idx_2::isindex => idx_2 == idx ? group[1] : idx_2))(b2_2)
            end
        end
        
        # iterate through equivalent indices in b1 and rewrite b1 to include only one of index
        # from each set of equivalent indices
        for group in b1_subsymmetry
            for idx in group[2:end]
                b1_3 = Rewrite(Postwalk(@rule ~idx_2::isindex => idx_2 == group[1] ? idx : idx_2))(b1_2)
                display(b1_3)
                display(b2_2)
                if b1_3 == b2_2 
                    println("EQUIVALENT SIEVES FOUND")
                end
            end
        end

        #TODO: fix

    end))(ex)
    return nothing
end


"""
    group_assignments(ex)
    
Return rewritten ex with equivalent assignments in a block grouped together.
"""
function group_assignments(ex)
    Fixpoint(Rewrite(Postwalk(@rule block(~s1..., assign(~lhs, +, ~rhs), ~s2..., assign(~lhs, +, ~rhs), ~s3...) =>
                        block(s1..., assign(lhs, +, call(*, 2, rhs)), s2..., s3...))))(ex)
end


"""
    find_swaps(_A, B)

Return list of lists of length two consisting of pair of indices that need
to be swapped in `_A` to make `_A` equivalent to `B`.
"""
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


"""
    get_intermediate_output(tn, count)

Return new symbol with name of form "tn_count".
"""
function get_intermediate_output(tn, count)
    ctx = JuliaContext()
    var_name = Symbol("_", tn.val, count)
    var = freshen(ctx, var_name)
end


# TODO: what if multiple indices need to be swapped in intermediate output? what is a kernel in which this would be the case?
"""
    exploit_output_replication_base(ex)

Given block `ex` consisting of updates that need to be performed for all nondiagonal 
coordinates, determine where output symmetry exists and introduce an intermediate 
output tensor that includes only canonical coordinates of otuput and can be replicated 
in post-processing.

Returns new expression, dictionary mapping lists of indices that need to be swapped
to the intermediate tensor in which they should be swapped, and a boolean representing
whether or not the output is fully symmetry across an axis. 
"""
function exploit_output_replication_base(ex)
    @capture ex block(~s...)
    update_count = length(s)

    replicate = Dict()
    count = 1
    ex = Fixpoint(Rewrite(Postwalk(@rule block(~s1..., assign(~lhs1, +, ~rhs), ~s2..., assign(~lhs2, +, ~rhs), ~s3...) => begin
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
    ))(ex)

    fully_replicable = false
    @capture ex block(~s...)
    if length(s) == update_count / 2 
        fully_replicable = true
    end

    ex, replicate, fully_replicable
end


"""
    exploit_output_replication_edge(ex, fully_replicable, replicate)

Given block `ex` consisting of updates that need to be performed to handle a set
of diagonal coordinates, boolean `fully_replicable` representing whether the output
is fully symmetric across an axis in the output, and `replicate` mapping lists of indices
that need to be swapped to intermediate output tensors, return new expression with 
intermediate outputs used wherever possible.
"""
function exploit_output_replication_edge(ex, fully_replicable, replicate)
    @capture ex block(~s...)
    update_count = length(s)

    count = 1
    # TODO: standardize variable names
    ex = Fixpoint(Rewrite(Postwalk(@rule block(~s1..., assign(~lhs1, +, ~rhs), ~s2..., assign(~lhs2, +, ~rhs), ~s3...) => begin
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
    ))(ex)

    # if fully replicable replace all instances of original output tensor
    # with tensor that will later be replicated 
    if fully_replicable
        #TODO: should confirm that replicate dict only has one value if fully_replicable
        tn_2 = first(values(replicate))
        ex = Fixpoint(Rewrite(Postwalk(@rule assign(~lhs, +, ~rhs) => begin
                if @capture lhs access(~tn, updater, ~idxs...)
                    if tn != tn_2
                        lhs = access(tn_2, updater, idxs...)
                    end
                end
                assign(lhs, +, rhs)
            end)
        ))(ex)
    end
    ex, replicate
end


"""
    is_base(cond)

Returns boolean representing whether `cond` is an expression with only != relations.
"""
function is_base(cond)
    if !(@capture cond call(and, ~conds...))
        conds = [cond]
    end
    for cond in conds
        if @capture cond call(==, ~a, ~b)
            return false
        end
    end
    return true
end


# TODO: make sure that we keep canonical coordinates of output (replicate for non-canonical)
# TODO: how to handle/mark fully_replicable case if there are multiple axes of output symmetry 
# TODO: does being fully_replicable mean that there is only one entry in replicate dict
"""
    exploit_output_replication(ex) 

Return an expression that exploits output symmetry and dictionary that maps lists of indices
across which values need to be replicated to the intermediate output in which these values need
to be replicated.
"""
function exploit_output_replication(ex) 
    replicate = Dict()
    fully_replicable = false

    ex = Rewrite(Postwalk(@rule ex sieve(~cond::is_base, ~body) => begin
        body, replicate, fully_replicable = exploit_output_replication_base(body)
        return sieve(cond, body)
    end))(ex)

    Rewrite(Postwalk(@rule ex sieve(~cond::((c) -> !is_base(c)), ~body) => begin
        body, replicate = exploit_output_replication_edge(body, fully_replicable, replicate)
        return sieve(cond, body)
    end))(ex)
end


"""
    consolidate_reads(ex)

Replace multiple access to a tensor with a let statement outside block.
"""
function consolidate_reads(ex)
    Rewrite(Prewalk(@rule block(~s1...) => begin
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
    end))(ex)
end


"""
    triangularize(ex, permutable_idxs)

Wrap `ex` in clause that restricts accesses to canonical triangle of `permutable_idxs`.
"""
function triangularize(ex, permutable_idxs, diagonals=true)
    # should only have one sieve in ex if no diagonals
    if !diagonals && @capture ex block(sieve(~cond, ~body))
        ex = body
    end

    idxs = order_canonically(permutable_idxs)
    conditions = []
    for i in 1:length(idxs)-1
        push!(conditions, diagonals ? call(<=, idxs[i], idxs[i+1]) : call(<, idxs[i], idxs[i+1]))
    end
    condition = length(conditions) > 1 ? call(and, conditions...) : conditions[1]
    return sieve(condition, ex)
end

# TODO: there should be a better/safer way to do this (e.g. somehow indicate nesting of ops) - something for future
function consolidate_operators(ops)
    if all(op -> in(op, ops), [literal(<=), literal(!=)]) && length(ops) == 2
        return literal(<)
    elseif all(op -> in(op, ops), [literal(<=), literal(==)]) && length(ops) == 2
        return literal(==)
    elseif all(op -> in(op, ops), [literal(<=), literal(==), literal(!=)]) && length(ops) == 3
        return literal(<=)
    else
        throw(ArgumentError("Unexpected combination of operators"))
    end
end

function consolidate_comparisons(conditions)
    idx_comparisons = Dict()

    for condition in conditions
        if !(@capture condition call(and, ~pairwise_conditions...))
            pairwise_conditions = [condition]
        end
        for pair in pairwise_conditions
            @capture pair call(~op, ~idxs...)
            ops = get(idx_comparisons, idxs, Set())
            push!(ops, op)
            idx_comparisons[idxs] = ops 
        end
    end

    conditions_2 = []
    for (idx, comparisons) in pairs(idx_comparisons)
        op = consolidate_operators(comparisons)
        push!(conditions_2, call(op, idx[1], idx[2]))
    end

    return length(conditions_2) == 1 ? conditions_2[1] : call(and, conditions_2...)
end

# TODO: still need to determine/evaluate when this is necessary/not necessary to do
# TODO: docstrings for this & helpers
function consolidate_conditions(ex)
    conditions_map = Dict()
    @capture ex sieve(~base_condition, ~body)

    Postwalk(@rule sieve(~body_condition, ~body_2) => begin
        @capture body_2 block(~updates...)
        for update in updates
            update_conditions = get(conditions_map, update, Set())
            push!(update_conditions, body_condition)
            push!(update_conditions, base_condition)
            conditions_map[update] = update_conditions 
        end
    end)(body)
    
    updates_map = Dict()
    for (update, conditions) in pairs(conditions_map)
        consolidated_condition = consolidate_comparisons(conditions)
        updates = get(updates_map, consolidated_condition, [])
        push!(updates, update)
        updates_map[consolidated_condition] = updates
    end

    sieves = []
    for (condition, updates) in pairs(updates_map)
        push!(sieves, sieve(condition, block(updates...)))
    end

    return block(sieves...)
end


# TODO: implement
function nest_conditions(ex)
    # determine which conditions are nestable
    # nest accordingly 
end


"""
    symmetrize2(ex, symmetric_tns)

    Rewrite ex to exploit symmetry in the tensors marked as symmetric in symmetric_tns
"""
function symmetrize3(ex, symmetric_tns, diagonals=true)
    # helper methods
    issymmetric(tn) = tn.val in symmetric_tns

    @capture ex assign(access(~lhs, updater, ~idxs...), ~op, ~rhs)

    permutable_idxs = get_permutable_idxs(rhs, issymmetric)
    conditions = get_conditions(permutable_idxs[1], diagonals)
    ex = add_updates(ex, conditions, permutable_idxs[1], issymmetric)
    # group_sieves(ex) # TODO: fix group_sieves
    ex = group_assignments(ex)
    ex = exploit_output_replication(ex)
    ex = triangularize(ex, permutable_idxs[1], diagonals)
    @info "before consolidating conditions"
    display(ex)
    ex = consolidate_conditions(ex)
    ex = consolidate_reads(ex) # TODO: need to figure out best place to do this (and how - prewalk or postwalk?)
    # TODO: compare # of conditions from before/after consolidate_conditions and keep version with less
    @info "after consolidating conditions"
    display(ex)
end