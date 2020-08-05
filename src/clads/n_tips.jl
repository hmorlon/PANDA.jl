function n_extant_tips(tree::Tree)
    function aux(subtree)
        if subtree.n_nodes < 2
            if subtree.extant
                return 1
            else
                return  0
            end
        else
            return aux(subtree.offsprings[1]) + aux(subtree.offsprings[2])
        end
    end

    aux(tree)
end


function n_extant_tips(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1})

    n::Int64 = 0

    for et in edge_trees
        n += et.tip_number
    end

    return n
end


function n_left(tree::Tree)
    function aux(subtree, x)
        if length(subtree.offsprings) == 0
            pushfirst!(x,0)
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
            pushfirst!(x,subtree.offsprings[1].n_nodes)
        end
    end

    x = Array{Int64,1}(undef,0)
    aux(tree, x)
    return x
end
