function tip_labels(tree::Tree)
    function aux(subtree, x)
        if length(subtree.offsprings) == 0
            pushfirst!(x,subtree.label)
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
        end
    end

    x = Array{String,1}(undef,0)
    aux(tree, x)
    return x
end

function tips(tree::Tree)
    function aux(subtree, x)
        if length(subtree.offsprings) == 0
            pushfirst!(x,true)
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
            pushfirst!(x,false)
        end
    end

    x = Array{Bool,1}(undef,0)
    aux(tree, x)
    return x
end



function update_rates(tree::Tree, rates::Array{T,1} ; id = 1) where {T<:Number}
    function aux(sub_tree, sub_rates)
        if sub_tree.n_nodes < 2
            root_rate = popfirst!(sub_rates)
            current_rates = sub_tree.attributes
            current_rates[id] = root_rate
            return Tree(sub_tree.offsprings, sub_tree.branch_length, current_rates, sub_tree.n_nodes, sub_tree.label)
        else
            root_rate = popfirst!(sub_rates)
            t1 = aux(sub_tree.offsprings[1], sub_rates)
            t2 = aux(sub_tree.offsprings[2], sub_rates)
            current_rates = sub_tree.attributes
            current_rates[id] = root_rate
            return Tree([t1,t2], sub_tree.branch_length, current_rates, sub_tree.n_nodes, sub_tree.label)
        end
    end

    new_rates = deepcopy(rates)
    aux(tree, new_rates)
end

function update_rates!(tree::Tree, rates::Array{T,1} ; id = 1) where {T<:Number}
    function aux(sub_tree, sub_rates)
        if sub_tree.n_nodes < 2
            root_rate = popfirst!(sub_rates)
            sub_tree.attributes[id] = root_rate
        else
            root_rate = popfirst!(sub_rates)
            sub_tree.attributes[id] = root_rate
            t1 = sub_tree.offsprings[1]
            aux(t1, sub_rates)
            t2 = sub_tree.offsprings[2]
            aux(t2, sub_rates)
            #return Tree([t1,t2], sub_tree.branch_length, root_rate, sub_tree.n_nodes)
        end
    end

    new_rates = deepcopy(rates)
    aux(tree, new_rates)
end

function update_rates!(tree::Tree, rates::Array{T,1}, edge_trees ; id = 1) where {T<:Number}
    update_rates!(tree,rates, id = id)

    for i in 1:length(edge_trees)
        edge_trees[i].tree.attributes[id] = rates[i+1]
    end

    return tree, edge_trees
end


function update_rates!(tree::Tree, rate::T ; id = 1) where {T<:Number}
    function aux(sub_tree)
        if sub_tree.n_nodes < 2
            sub_tree.attributes[id] *= rate
        else
            sub_tree.attributes[id] *= rate
            t1 = sub_tree.offsprings[1]
            aux(t1)
            t2 = sub_tree.offsprings[2]
            aux(t2)
        end
    end

    aux(tree)
end

function get_parent_edge(tree::Tree, edge_id::Int64)
    if edge_id == 0
        return NaN
    end

    function aux(subtree, edge)
        n_edges_left = subtree.offsprings[1].n_nodes + 1
        if edge == 2 || edge == (n_edges_left + 1)
            return 1
        elseif edge > (n_edges_left + 1)
            return n_edges_left + aux(subtree.offsprings[2], edge - n_edges_left)
        else
            return 1 + aux(subtree.offsprings[1], edge - 1)
        end
    end

    aux(tree, edge_id + 1) - 1
end
