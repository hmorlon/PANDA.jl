
function extract_branch_lengths(tree::Tree)
    function aux(subtree, bl)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(bl,subtree.branch_length)
        else
            aux(subtree.offsprings[2],bl)
            aux(subtree.offsprings[1],bl)
            pushfirst!(bl,subtree.branch_length)
        end
    end

    bl =  Array{Float64,1}(undef,0)
    aux(tree, bl)
    return bl
end

function get_node_depth(tree, edge_id)
    function aux(subtree, edge)
        if subtree.n_nodes < 2
            return 0.
        else
            if edge == 1
                return aux(subtree.offsprings[1], edge) + subtree.offsprings[1].branch_length
            elseif (edge - 1) > subtree.offsprings[1].n_nodes
                return aux(subtree.offsprings[2], edge - 1 - subtree.offsprings[1].n_nodes)
            else
                return aux(subtree.offsprings[1], edge - 1)
            end
        end
    end

    aux(tree, edge_id + 1)
end
