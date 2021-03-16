
function extract_rates(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1}; id = 1)
    function aux(subtree, rates,ri)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(rates,subtree.attributes[id]*ri)
        else
            aux(subtree.offsprings[2],rates,ri)
            aux(subtree.offsprings[1],rates,ri)
            pushfirst!(rates,subtree.attributes[id]*ri)
        end
    end

    rates =  Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, rates, edge_trees[i].stem_rate[1])
    end
    pushfirst!(rates, tree.attributes[1])
    return rates
end


function extract_relative_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, rates::Array{Float64,1}; id = 1)

    function aux(subtree, parent_rate, x)
        if subtree.n_nodes > 0
            if length(subtree.offsprings) == 0
                pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
            else
                log_rate = deepcopy(log(subtree.attributes[id]))
                aux(subtree.offsprings[2], log_rate, x)
                aux(subtree.offsprings[1], log_rate, x)
                pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
            end
        end
    end

    relative_rates = Array{Float64,1}(undef,0)

    i = 0
    for et in edge_trees
        i+= 1
        parent_edge = et.parent_edge

        parent_rate = rates[1]
        if parent_edge > 0
            parent_rate = edge_trees[parent_edge].tip_rate * rates[parent_edge + 1]
        end

        r = log(rates[i + 1]) - log(parent_rate)
        pushfirst!(relative_rates,r)
        if length(et.tree.offsprings) > 0
            aux(et.tree.offsprings[1], 0., relative_rates)
            aux(et.tree.offsprings[2], 0., relative_rates)
        end
    end
    return relative_rates
end

function extract_relative_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}; id = 1)
    function aux(subtree, parent_rate, x)
        if subtree.n_nodes > 0
            if length(subtree.offsprings) == 0
                pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
            else
                log_rate = deepcopy(log(subtree.attributes[id]))
                aux(subtree.offsprings[2], log_rate, x)
                aux(subtree.offsprings[1], log_rate, x)
                pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
            end
        end
    end

    relative_rates = Array{Float64,1}(undef,0)

    for et in edge_trees
        parent_edge = et.parent_edge

        parent_rate = tree.attributes[1]
        if parent_edge > 0
            parent_rate = edge_trees[parent_edge].tip_rate * edge_trees[parent_edge].stem_rate[1]
        end


        r = log(et.stem_rate[1]) - log(parent_rate)
        pushfirst!(relative_rates,r)

        if length(et.tree.offsprings) > 1#0000000
            aux(et.tree.offsprings[1], 0., relative_rates)
            aux(et.tree.offsprings[2], 0., relative_rates)
        end
    end

    return relative_rates
end

function extract_rates(tree::Tree; id=1)
    function aux(subtree, rates)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(rates,subtree.attributes[id])
        else
            aux(subtree.offsprings[2],rates)
            aux(subtree.offsprings[1],rates)
            pushfirst!(rates,subtree.attributes[id])
        end
    end

    rates = Array{Float64,1}(undef,0)
    aux(tree, rates)
    return(rates)
end

#=
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
=#

function extract_tip_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1},
    tips_id::Array{Bool,1}, rates::Array{Float64,1} ; id = 1, return_extinct = true)

    tip_rates = Array{Float64,1}(undef,0)
    for i in 1:length(edge_trees)
        if tips_id[i]
            r = sample(edge_trees[i].tip_rates) * edge_trees[i].stem_rate[1]
            push!(tip_rates,r)
        end
    end
    return tip_rates
end
