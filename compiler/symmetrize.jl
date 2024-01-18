using Finch
using Finch.FinchNotation
using RewriteTools
using RewriteTools.Rewriters
using Combinatorics

import Finch.FinchNotation: and, or

isindex(ex::FinchNode) = ex.kind === index
order_canonically(idxs) = sort(idxs, by = i->i.val)

iscommutative(::typeof(or)) = true
iscommutative(::typeof(and)) = true
iscommutative(::typeof(|)) = true
iscommutative(::typeof(&)) = true
iscommutative(::typeof(+)) = true
iscommutative(::typeof(*)) = true
iscommutative(::typeof(min)) = true
iscommutative(::typeof(max)) = true
iscommutative(alg) = false

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
        (@rule call(~op, ~tns...) => call(op, sort(tns, by = tn->hash(tn))...))
    ])))
    _normalize(ex)
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
    display(normalized)

    permuted_exs = []
    map = Dict{Vector{FinchNode}, Set{FinchNode}}([])
    for idxs in permutable_idxs
        for permutation in permutations(idxs)
            new_ex = deepcopy(ex)
            operands = []
            Postwalk((x) -> if ((@capture x access(~tn, reader, ~original_idxs...)) || (@capture x access(~tn, updater, ~original_idxs...)))
                _original_idxs = deepcopy(original_idxs)
                if !(tn.val in symmetric_tns)
                    for i in 1:length(idxs)
                        idxs_to_replace = findall(y -> y == idxs[i], _original_idxs)
                        for j in idxs_to_replace
                            original_idxs[j] = permutation[i]
                        end
                    end
                end
                if @capture x access(~tn, reader, ~original_idxs...)
                    push!(operands, x)
                end
            end)(new_ex)
            display(new_ex)
            push!(permuted_exs, new_ex)
            # Sort operands in canonical order (ensure that operator is commutative before doing this)
            sort!(operands, by = operands -> hash(operands))
        
            if @capture new_ex assign(~lhs, ~op, ~rhs)
                lhs_set = get(map, operands, Set([]))
                push!(lhs_set, lhs)
                map[operands] = lhs_set
            end
        end
    end

    for update in keys(map)
        if length(map[update]) > 1 
            println("The following positions are the same (output symmetry):")
            for lhs in map[update]
                println(lhs)
            end
        end
    end
end
