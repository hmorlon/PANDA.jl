"""
    sample_tips(tree::Tree, f::Float64 ; root_age = true)

Sample tips from a tree. Each tip is kept with probability f.

# Arguments
- `tree::Tree`: the phylogeny.
- `f::Float64`: the sampling probability.

# Keyword arguments
- `root_age::Bool`: Should extinct lineages be removed from the output tree? Default to false.
"""
function sample_tips(tree::Tree, f::Float64 ; root_age = true)
    function combine(left, right, bl, att)
        extant = (left.extant || right.extant)
        return Tree([left,right], bl, att, extant, "")
    end
    function aux(subtree)
        u = 0.5
        if subtree.extant
            if length(subtree.offsprings) == 0
                if f < 1.
                    u = rand()
                end
                if u < f
                    return subtree
                else
                    return Tree()
                end
            else
                combine(aux(subtree.offsprings[1]), aux(subtree.offsprings[2]), subtree.branch_length, subtree.attributes)
            end
        else
            return subtree
        end
    end

    if root_age
        left_tree = Tree()
        while (!left_tree.extant)
            left_tree = aux(tree.offsprings[1])
        end
        right_tree = Tree()
        while (!right_tree.extant)
            right_tree = aux(tree.offsprings[2])
        end
        new_tree = Tree([left_tree,right_tree], tree.branch_length, tree.attributes, tree.label)
        tip_rates = extract_tip_rates_light(new_tree, return_extinct = false)
        return (prune_extinct_lineages(new_tree), tip_rates)
    else
        new_tree = aux(tree)
        tip_rates = extract_tip_rates_light(new_tree, return_extinct = false)
        return (prune_extinct_lineages(new_tree), tip_rates)
    end
end
