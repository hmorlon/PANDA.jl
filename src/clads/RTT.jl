function branch_time_id(tree::Tree, times::Array{Float64,1})
    nb = node_depths_base(tree)
    ne = node_depths(tree)

    rd= maximum(ne)
    ids = Array{Array{Int64,1},1}(undef,0)

    for t in times
        if t <= rd
            idt = Array{Int64,1}(undef,0)

            for i in 1:length(nb)
                if nb[i] <= t < ne[i]
                    push!(idt,i)
                end
            end
            push!(ids,idt)
        end
    end

    return ids
end

function branch_time_id(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, times::Array{Float64,1})
    nb = node_depths_base(tree, edge_trees)

    ne = node_depths(tree,edge_trees)

    rd= maximum(ne)
    ids = Array{Array{Int64,1},1}(undef,0)

    for t in times
        if t <= rd
            idt = Array{Int64,1}(undef,0)

            for i in 1:length(nb)
                if nb[i] <= t < ne[i]
                    push!(idt,i+1)
                end
            end
            push!(ids,idt)
        end
    end

    return ids
end

function time_rates(ids::Array{Array{Int64,1},1}, row::Array{Float64,1})
    mean_rates = Array{Float64,1}(undef,0)

    for id in ids
        mr = 0.
        n = 0.
        for i in id
            n += 1.
            mr += row[i]
        end
        push!(mean_rates, mr/n)
    end
    return mean_rates
end

function time_rates(tree::Tree, row::Array{Float64,1}, times::Array{Float64,1})
    ids = branch_time_id(tree, times)

    time_rates(ids, row)
end

function time_rates(tree::Tree, times::Array{Float64,1})
    ids = branch_time_id(tree, times)
    row = extract_rates(tree)
    time_rates(ids, row)
end


function time_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, times::Array{Float64,1})
    ids = branch_time_id(tree, edge_trees, times)
    row = extract_rates(tree, edge_trees)
    time_rates(ids, row)
end
