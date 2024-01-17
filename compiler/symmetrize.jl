using Finch
using Finch.FinchNotation
using RewriteTools
using RewriteTools.Rewriters
using Combinatorics

function find_symmetry(ex, symmetric_tns)
    @capture ex assign(access(~lhs, updater, ~idxs...), ~op, ~rhs)

    permutable_idxs = []
    # Sort indices of symmetric matrices in canonical order (alphabetical)
    Postwalk((x) -> if @capture x access(~tn, reader, ~idxs...)
        if tn.val in symmetric_tns
            sort!(idxs, by = i -> i.val)
            push!(permutable_idxs, idxs)
        end
    end)(rhs)

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
