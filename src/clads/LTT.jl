
function LTT(tree::Tree)
    function aux(subtree, node_times, events, current_depth)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth + subtree.branch_length)
            if subtree.extant
                pushfirst!(events, 0)
            else
                pushfirst!(events, -1)
            end
        else
            aux(subtree.offsprings[2],node_times,events, current_depth + subtree.branch_length)
            aux(subtree.offsprings[1],node_times,events, current_depth + subtree.branch_length)
            pushfirst!(node_times,current_depth+ subtree.branch_length)
            pushfirst!(events,1)
        end
    end

    node_times = Array{Float64,1}(undef,0)
    events = Array{Int64,1}(undef,0)
    aux(tree, node_times, events, 0.)
    indices = sortperm(node_times)
    events = events[indices]
    node_times = node_times[indices]
    if length(events) > 1
        while events[end-1] == 0
            pop!(events)
            pop!(node_times)
        end
    end
    events[1] += 1
    events = cumsum(events, dims = 1)


    return node_times, events
end

function LTT(tree::Tree, times::Array{Float64,1})
    n = length(times)

    function aux(subtree::Tree, from_root::Float64, id_time::Int64, ltt::Array{Int64,1})
        if id_time <= n
            if (subtree.branch_length + from_root) <= times[id_time]
                if subtree.n_nodes < 2
                    if !subtree.extant
                        ltt[id_time] += -1
                    end
                else
                    ltt[id_time] += 1
                    aux(subtree.offsprings[2],from_root+ subtree.branch_length, id_time, ltt)
                    aux(subtree.offsprings[1],from_root+ subtree.branch_length, id_time, ltt)
                end
            else
                aux(subtree,from_root, id_time+1, ltt)
            end
        end

    end

    ltt = fill(0, length(times))
    ltt[1] = 1
    aux(tree, 0., 1, ltt)
    ltt = cumsum(ltt, dims = 1)


    return times, ltt

end

function LTT(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, times::Array{Float64,1})
    n = length(times)

    function aux(subtree::Tree, from_root::Float64, id_time::Int64, ltt::Array{Int64,1}, scale::Float64)
        if id_time <= n
            if (subtree.branch_length/scale + from_root) <= times[id_time]
                if subtree.n_nodes < 2
                    if !subtree.extant
                        ltt[id_time] += -1
                    end
                else
                    ltt[id_time] += 1
                    aux(subtree.offsprings[2],from_root+ subtree.branch_length/scale, id_time, ltt, scale)
                    aux(subtree.offsprings[1],from_root+ subtree.branch_length/scale, id_time, ltt, scale)
                end
            else
                aux(subtree,from_root, id_time+1, ltt, scale)
            end
        end

    end

    ltt = fill(0, length(times))
    ltt[1] = 1
    aux(tree, 0., 1, ltt, 1.)
    #println("1 $ltt ")
    i=0
    md = edge_trees[1].stem_depth
    for et in edge_trees# in 1:length(edge_trees)
        if et.tree.n_nodes > 1
            aux(et.tree, md - et.stem_depth, 1, ltt,et.scale)
        end
    end
    #println("2 $ltt ")
    ltt = cumsum(ltt, dims = 1)
    #println("3 $ltt ")

    return times, ltt
end
