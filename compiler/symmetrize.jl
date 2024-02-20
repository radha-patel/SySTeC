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


# TODO: docstrings / reorganize
is_index(ex::FinchNode) = ex.kind === index
is_commutative(op) = typeof(op.val) == typeof(*) ? true : false

"""
    permute_indices(ex, idxs, permutations)

Returns a block with length(permutations) expressions where expression i has its
indices permutated such that each index idxs[j] is replaced with index 
permutations[i][j]
"""
function permute_indices(ex, idxs, permutations)
    permuted_exs = []
    for perm in permutations
        ex_2 = Rewrite(Postwalk(@rule ~idx::is_index => begin
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
        (@rule call(~op::is_commutative, ~tns...) => call(op, sort(tns, by = tn->hash(tn))...))
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


# TODO: may need an additional layer of nesting here
"""
    get_possible_swaps(subsymmetry)

Given the subsymmetry of a tensor, return a list of lists consisting of the pairs
of indices (i.e. equivalent indices) to access the same value of the tensor.
"""
function get_possible_swaps(subsymmetry)
    swaps = []
    for group in subsymmetry
        append!(swaps, collect(combinations(group, 2)))
    end
    return swaps
end


"""
    updates_count(body)

Returns the number of updates that are performed in the expression `body`
where an "update" consists of adding to the output at a particular coordinate
the product of the input matrices for some set index values.  
"""
function updates_count(body)
    count = 0
    Postwalk(@rule assign(~lhs, +, ~rhs) => begin
        if @capture rhs call(*, 2, ~rhs_2)
            count += 2
        else
            count += 1
        end
    end)(body)
    return count
end


"""
    swap_indices(ex, idx_1, idx_2, issymmetric)

Returns expression `ex` with all instances of indices `idx_1` and `idx_2` swapped.
"""
function swap_indices(ex, idx_1, idx_2, issymmetric) 
    Rewrite(Postwalk(@rule access(~tn::((t) -> (!issymmetric(t))), ~mode, ~idxs...) => begin
        swapped_idxs = deepcopy(idxs)
        for i_1 in 1:length(idxs)
            if idxs[i_1] == idx_1
                swapped_idxs[i_1] = idx_2
            elseif idxs[i_1] == idx_2
                swapped_idxs[i_1] = idx_1
            end
        end
        access(tn, mode, swapped_idxs...)
    end))(ex)
end


"""
    are_updates_identical(updates_1, udpates_2) 

Given lists of updates `updates_1` and `updates_2`, return whether the updates in 
`updates_1` are some permutation of `updates_2``.
"""
function are_updates_identical(updates_1, udpates_2) 
    for updates_1_perm in permutations(updates_1)
        if updates_1_perm == udpates_2
            return true
        end
    end
    return false
end


"""
    decouple_assignments(ex, subsymmetry)

Upgroup groups of two assignments and rewrite one with a pair of equivalent indices
swapped in all nonsymmetric tensors. Return list of resulting assignment expressions.
"""
function decouple_grouped_updates(ex, subsymmetry, issymmetric)
    assignments = []
    Postwalk(@rule assign(~lhs, +, ~rhs) => begin
        if @capture rhs call(*, 2, ~rhs_2)
            idx_1, idx_2 = subsymmetry[1] # TODO: what if we have a more complicated subsymmetry?
            lhs_swapped = swap_indices(lhs, idx_1, idx_2, issymmetric)
            rhs_3 = swap_indices(rhs_2, idx_1, idx_2, issymmetric)
            push!(assignments, assign(lhs, +, rhs_2))
            push!(assignments, assign(lhs_swapped, +, normalize(rhs_3, issymmetric)))
        else
            push!(assignments, assign(lhs, +, rhs))
        end
    end)(ex)
    return assignments
end


#TODO: this step can probably be combined/happen within add_updates 
"""
    group_sieves(ex, issymmetric)

Identifies sieves that perform the same updates and groups them together. Returns
a boolean representing whether any equivalent sieves were found and the resulting
expression as a result of performing this transform. 
"""
function group_sieves(ex, issymmetric)
    grouped = false
    ex = Rewrite(Postwalk(@rule block(~s1..., sieve(~c1, ~b1), ~s2..., sieve(~c2, ~b2), ~s3...) => begin
        if updates_count(b1) != updates_count(b2)
            return nothing
        end

        b1_subsymmetry = get_subsymmetry(c1)
        b2_subsymmetry = get_subsymmetry(c2)

        b1_s = decouple_grouped_updates(b1, b1_subsymmetry, issymmetric)
        b2_s = decouple_grouped_updates(b2, b2_subsymmetry, issymmetric)

        if are_updates_identical(b1_s, b2_s)
            grouped = true
            return block(s1..., sieve(call(or, c1, c2), block(b1_s...)), s2..., s3...)
        end
        # TODO: combine conditions here -> 1. maybe just throw in or and combine in another transform, or 2. "merge" here
    end))(ex)
    return grouped, ex
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
            idxs_3 = sort_permutations([idxs1, idxs2])[1]
            _lhs = access(_tn, updater, idxs_3...)
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
            idxs_3 = sort_permutations([idxs1, idxs2])[1]
            _lhs = access(_tn, updater, idxs_3...)
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

Returns boolean representing whether `cond` is an expression that restricts updates to only
non-diagonal coordinates of output.
"""
function is_base(cond)
    conds = []
    Postwalk(@rule call(~op::is_comparison, ~idx_1, ~idx_2) => begin
        push!(conds, call(op, idx_1, idx_2))
    end)(cond)
    for cond in conds
        if (@capture cond call(==, ~a, ~b)) || (@capture cond call(<=, ~a, ~b))
            return false
        end
    end
    return true
end


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
    get_tn_var_name(tn, idxs)

Returns name of variable for value of `tn` at `idxs` coordinate. 
"""
get_tn_var_name(tn, idxs) = Symbol(tn.val, "_", join([idx.val for idx in idxs]))


"""
    op_str(op)

Returns string representation of `op`.
"""
function op_str(op)
    if op == literal(==)
        return "eq"
    elseif op == literal(!=)
        return "neq"
    elseif op == literal(<=)
        return "leq"
    elseif op == literal(<)
        return "lt"
    else
        throw(ArgumentError("Unexpected operator"))
    end
end
    

"""
    get_comp_var_name(tn, idx_1, idx_2)

Returns name of variable to represent comparison between `idx_1` and `idx_2` 
"""
function get_comp_var_name(op, idx_1, idx_2)
    @capture idx_1 call(identity, ~idx_1)
    @capture idx_2 call(identity, ~idx_2)
    
    Symbol(idx_1.val, idx_2.val, "_", op_str(op))
end


"""
    consolidate_reads(ex)

Replace the following with a let statement outside block:
    - multiple accesses to a tensor
    - repetitions of a comparison (including the negation)
"""
# TODO: would be nice to put let statements on separate lines for readability
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
                var = freshen(ctx, get_tn_var_name(tn, idxs))
                ex = Postwalk(@rule node => var)(ex)
                ex = define(var, access(tn, reader, idxs...), ex)
            end
        end
        for (node, count) in counts
            if @capture(node, call(~op::is_comparison, ~idx_1, ~idx_2)) && count > 1
                ctx = JuliaContext()
                var = freshen(ctx, get_comp_var_name(op, idx_1, idx_2))
                ex = Postwalk(@rule node => var)(ex)
                if is_eq_op(op) 
                    negation = call(literal(!=), idx_1, idx_2)
                    ex = Postwalk(@rule negation => call(literal(!=), var))(ex)
                end
                ex = define(var, call(op, idx_1, idx_2), ex)
            end
        end
        # TODO: how to best consolidate the other nodes?
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
"""
    consolidate_operators(ops)

Given list of operators `ops`, each representing an expression of the form (j ops[i] k) for 
1 <= i < length(ops), return operator for a singular expression (j `operator` k) that would
satisfy the combination of expressions.

Here, the combination of expressions is of the form (j ops[1] k) && ((j ops[2] k) || (j ops[3] k) ...)
"""
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


"""
    sort_conditions(conds)

Sort conditions comparing indices such that the first indices are in canonical order.
"""
sort_conditions = conds -> sort(conds, by=cond -> begin 
    @capture cond call(~op, ~idx_1, ~idx_2)
    idx_1.val
end)


"""
    make_strict(cond)

For a condition expression of the form `idx_1 <= idx_2 <= idx_3 ...` returns equivalent
expression but with all unstrict inequalities made strict.
"""
function make_strict(cond)
    if !(@capture cond call(and, ~conds...))
        conds = [cond]
    end
    strict_conds = []
    for cond in conds 
        @capture cond call(<=, ~idx_1, ~idx_2)
        push!(strict_conds, call(<, idx_1, idx_2))
    end
    return length(strict_conds) == 1 ? strict_conds[1] : call(and, strict_conds...)
end


"""
    consolidate_comparisons(conditions)

Given a list of Finch boolean expressions `conditions`, return a single
Finch boolean expression that satisfies all the conditions.
"""
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

    # TODO: sort conditions here
    conditions_2 = sort_conditions(conditions_2)
    return length(conditions_2) == 1 ? conditions_2[1] : call(and, conditions_2...)
end


# TODO: still need to determine/evaluate when this is necessary/not necessary to do
"""
    consolidate_conditions(ex)

Reorganize sieve conditions and bodies in `ex` such that each update
operate happens only once (no repeats).
"""
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


"""
    insert_loops(ex, idxs, loop_order)

Wrap `ex` with loop with all indices.
"""
# TODO: how to optimally order indices? should we be interspersing non_permutable_idxs in permutable_idxs
function insert_loops(ex, permutable_idxs, loop_order = [])
    if isempty(loop_order)
        permutable_idxs = order_canonically(permutable_idxs)
        all_idxs = get_idxs(ex)
        non_permutable_idxs = idxs_not_in_set(all_idxs, permutable_idxs)
        loop_order = vcat(non_permutable_idxs, permutable_idxs)
    end

    ex_with_loops = ex
    for idx in loop_order
        ex_with_loops = loop(idx, virtual(Dimensionless()), ex_with_loops)
    end
    return ex_with_loops
end


"""
    is_comparison(op)

Returns whether `op` is a comparison operator (not including negations -- e.g. !=).
"""
is_comparison(op) = op in [literal(<), literal(<=), literal(==)]


"""
    is_eq_op(op)

Returns whether `op` is a an == operator.
"""
is_eq_op(op) = op == literal(==)


"""
    insert_identity(ex)

Inserts identity() wrapper around all indices in `ex` that correspond to a comparison done
for a diagonal-related update.
"""
function insert_identity(ex)
    if @capture ex sieve(~triangular_condition, ~ex_2)
        ex_2 = Rewrite(Postwalk(@rule call(~op::is_comparison, ~idx_1, ~idx_2) => begin
            call(op, call(identity, idx_1), call(identity, idx_2))
        end))(ex_2)
        ex = sieve(triangular_condition, ex_2)
    end
    ex 
end


"""
    separate_loop_nests(ex)

Returns two expressions, where first expression uses the values on the nondiagonals
of the symmetric tensors and the second expression uses the values on the diagonals.
"""
function separate_loop_nests(ex)
    base = []
    edge = []
    if @capture ex sieve(~triangle_cond, ~ex_2)
        Postwalk(@rule sieve(~cond, ~body) => begin 
            if is_base(cond)
                strict_triangle_cond = make_strict(triangle_cond)
                push!(base, sieve(strict_triangle_cond, body))
            else
                push!(edge, sieve(cond, body))
            end
        end)(ex_2)
    end
    base[1], sieve(triangle_cond, block(edge...))
end


"""
    conditions_count(ex)

Returns the number of sieves in `ex`.
"""
function conditions_count(ex)
    count = 0
    Postwalk(@rule sieve(~cond, ~body) => begin count += 1 end)(ex)
    return count
end


"""
    get_idxs(ex)

Returns set of all indices in `ex`.
"""
function get_idxs(ex)
    s = Set()
    Postwalk(@rule ~idx::is_index => push!(s, idx))(ex)
    return s
end


"""
    idxs_not_in_set(s, idxs)

Returns list of indices in `s` that are not in `idxs`
"""
function idxs_not_in_set(s, idxs)
    return [idx for idx in s if !(idx in idxs)]
end


"""
    get_transposed_tn_name(tn, count)

Return new name for the `count`-th transposed form of tensor `tn`.
"""
get_transposed_tn_name(tn, count) = count == 1 ? Symbol(tn.val, "_T") : Symbol(tn.val, "_T", "_", count)


"""
    transpose_operands(ex, issymmetric, loop_order)

Returns `ex` with nonsymmetric operands transposed such that accesses are concordant with
loop order and dictionary mapping tuples of the form `(tn, swaps)` to the new name of the 
transposed tensor where `tn` is the name of the original tensor and `swaps` is a list of 
lists of length two corresponding to the indices that need to be transposed in `tn` to get
the new tensor.
"""
function transpose_operands(ex, issymmetric, loop_order)
    @capture ex assign(access(~lhs, ~updater, ~idxs...), ~op, ~rhs)

    transposed = Dict()
    count = 0
    rhs = Rewrite(Postwalk(@rule access(~tn::((t) -> !issymmetric(t)), ~mode, ~idxs...) => begin
        idxs_depths = [(findfirst(x -> x == idx, loop_order), idx) for idx in idxs]
        transposed_idxs = [tup[2] for tup in sort(idxs_depths)]
        swaps = Set(find_swaps(idxs, transposed_idxs))

        if isempty(swaps)
            return nothing
        end

        if haskey(transposed, (tn, swaps))
            new_tn = transposed[(tn, swaps)]
        else
            count += 1
            new_tn = get_transposed_tn_name(tn, count)
            transposed[(tn, swaps)] = new_tn
        end
        return access(new_tn, mode, transposed_idxs...)
    end))(rhs)

    ex = assign(access(lhs, updater, idxs...), op, rhs)
    return ex, transposed
end

"""
    symmetrize2(ex, symmetric_tns)

    Rewrite ex to exploit symmetry in the tensors marked as symmetric in symmetric_tns
"""
# TODO: transpose some tensors?
# TODO: option to generate code to initialize tensors 
# TODO: stress testing - for what sizes/complexities of tensors/kernels does this work?
function symmetrize(ex, symmetric_tns, loop_order=[], diagonals=true)
    # helper methods
    issymmetric(tn) = tn.val in symmetric_tns

    @capture ex assign(access(~lhs, updater, ~idxs...), ~op, ~rhs)
    if !isempty(loop_order)
        all_idxs = get_idxs(ex)
        @assert all_idxs == Set(loop_order) "loop order does not include all indices"
        ex, transposed = transpose_operands(ex, issymmetric, loop_order)
    end

    # TODO: given multiple symmetric matrices, how many is it worth optimizing for
    permutable_idxs = get_permutable_idxs(rhs, issymmetric)
    permutable_idxs = collect(Set(Iterators.flatten(permutable_idxs)))

    conditions = get_conditions(permutable_idxs, diagonals)
    ex = add_updates(ex, conditions, permutable_idxs, issymmetric)
    ex = group_assignments(ex)
    ex = exploit_output_replication(ex)
    # TODO: grouping takes a LONG time - constrain when we actually do this?
    grouped, ex = group_sieves(ex, issymmetric)
    # grouped = false
    ex = triangularize(ex, permutable_idxs, diagonals)
    if !grouped 
        ex_2 = consolidate_conditions(ex)
        # TODO: maybe there is a better metric to determine which expression to keep?
        ex = conditions_count(ex_2) < conditions_count(ex) ? ex_2 : ex
        ex = consolidate_reads(ex) # TODO: need to figure out best place to do this (and how - prewalk or postwalk?)
        ex = insert_loops(ex, permutable_idxs, loop_order)
        display(ex) 
    else
        ex_base, ex_edge = separate_loop_nests(ex)
        ex_base = consolidate_reads(ex_base)
        ex_edge = consolidate_reads(ex_edge)
        ex_edge = insert_identity(ex_edge)
        ex_base = insert_loops(ex_base, permutable_idxs, loop_order)
        ex_edge = insert_loops(ex_edge, permutable_idxs, loop_order)
        display(ex_base)
        display(ex_edge)
    end
end