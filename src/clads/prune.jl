"""
    sample_tips(tree::Tree, f::Float64 ; root_age = true)

Sample tips from a tree. Each tip is kept with probability `f`. The fuction returns a tupple `(phylo, rates)`, where phylo is a `Tree` object containing the sampled tree, and rates a vector of `Floats` with its tip rates.

# Arguments
- `tree::Tree`: the phylogeny.
- `f::Float64`: the sampling probability.

# Keyword arguments
- `root_age::Bool`: If true, the sampling is constrained so that the tree keeps the same root age : at least one tip is sampled in each of the two tree subtrees.
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

function extract_tip_rates_light(tree::Tree ; id = 1, return_extinct = true)

    function aux(subtree, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            if subtree.extant || return_extinct
                pushfirst!(x,subtree.attributes[id])
            end
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
        end
    end

    x = Array{Float64,1}(undef,0)
    aux(tree,x)
    return x
end

function prune_extinct_lineages(tree::Tree)
    function combine(left, right, bl, att)
        if ! left.extant
            if ! right.extant
                return Tree()
            else
                return Tree(right.offsprings, right.branch_length +bl,
                    att, true, right.label)
            end
        elseif ! right.extant
            return Tree(left.offsprings, left.branch_length + bl,
                att, true, left.label)
        else
            return Tree([left, right], bl, att, true, "")
        end
    end

    function aux(subtree)
        if length(subtree.offsprings) == 0
            return subtree
        else
            combine(aux(subtree.offsprings[1]), aux(subtree.offsprings[2]), subtree.branch_length, subtree.attributes)
        end
    end

    aux(tree)
end
