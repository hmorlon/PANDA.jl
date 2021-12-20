#=
New DA structure that allows minimal manipulation
=#

struct EdgeTreeRates2
    tree::Tree                          # the grafted tree
    stem_rate::Array{Float64,1}         # λ_i : every other rates are given divided by λ_i
    tip_rate::Float64                   # the tip rate
    tip_id::Int64                       # tip_id
    tip_rates::Array{Float64,1}         # rates at t_i (for internal branches)
    effective_bl::Float64               # t_j * λ_j / λ_i
    tip_number::Int64                   # extant tips
    scale::Float64
    #relative_rates::Array{Float64,1}
    depth::Float64
    stem_depth::Float64
    parent_edge::Int64
end

function EdgeTreeRates2(tree::Tree, stem_rate::Float64, tip_rate::Float64,
    tip_id::Int64, tip_rates::Array{Float64,1}, tip_number::Int64, depth::Float64,stem_depth::Float64, parent_edge::Int64)
    effective_bl = sum(extract_branch_lengths(tree) .* extract_rates(tree)) / stem_rate
    return EdgeTreeRates2(tree, [stem_rate], tip_rate, tip_id, tip_rates,
        effective_bl, tip_number, depth, stem_depth,parent_edge)
end

function EdgeTreeRates2(tree::Tree, stem_rate::Float64, tip_rate::Float64,
    tip_id::Int64, tip_rates::Array{Float64,1},
    effective_bl::Float64, tip_number::Int64,depth::Float64,stem_depth::Float64,
    parent_edge::Int64)
    return EdgeTreeRates2(tree, [stem_rate], tip_rate, tip_id, tip_rates, effective_bl, tip_number, depth,stem_depth, parent_edge)
end

function EdgeTreeRates2(tree::Tree, stem_rate::Array{Float64,1}, tip_rate::Float64,
    tip_id::Int64, tip_rates::Array{Float64,1},
    effective_bl::Float64, tip_number::Int64,
    depth::Float64,stem_depth::Float64,parent_edge::Int64)

    return EdgeTreeRates2(tree, stem_rate, tip_rate, tip_id, tip_rates, effective_bl,
        tip_number, stem_rate[1], depth,stem_depth, parent_edge)
end

#=
A few functions to manipulate it
=#

function change_bl(tree::Tree, scale::Float64)
    function aux(sub_tree)
        if sub_tree.n_nodes < 2
            return Tree(sub_tree.offsprings,
                sub_tree.branch_length / scale,
                sub_tree.attributes.*scale,
                sub_tree.n_nodes,
                sub_tree.extant)
        else
            left = aux(sub_tree.offsprings[1])
            right = aux(sub_tree.offsprings[2])
            return Tree([left,right],
                sub_tree.branch_length / scale,
                sub_tree.attributes.*scale,
                sub_tree.n_nodes,
                sub_tree.extant)
        end
    end

    aux(tree)
end

function change_bl(tree::Tree, scale::Float64, stem_rate::Float64)
    function aux(sub_tree)
        if sub_tree.n_nodes < 2
            return Tree(sub_tree.offsprings,
                sub_tree.branch_length / scale,
                sub_tree.attributes.*stem_rate,
                sub_tree.n_nodes,
                sub_tree.extant)
        else
            left = aux(sub_tree.offsprings[1])
            right = aux(sub_tree.offsprings[2])
            return Tree([left,right],
                sub_tree.branch_length / scale,
                sub_tree.attributes.*stem_rate,
                sub_tree.n_nodes,
                sub_tree.extant)
        end
    end

    aux(tree)
end


function graft_edge_trees(tree, edge_trees::Array{EdgeTreeRates2,1})

    function aux(subtree, sub_edge_trees)
        if subtree.n_nodes < 2
            return change_bl(sub_edge_trees[1].tree, sub_edge_trees[1].scale, sub_edge_trees[1].stem_rate[1])
        else
            n_edges_left = subtree.offsprings[1].n_nodes + 1
            sub_left = deepcopy(sub_edge_trees[2:n_edges_left])
            sub_right = deepcopy(sub_edge_trees[(n_edges_left+1):subtree.n_nodes])
            subtree_left = aux(subtree.offsprings[1], sub_left)
            subtree_right = aux(subtree.offsprings[2], sub_right)
            return graft_to_edge(Tree([subtree_left, subtree_right], 0., subtree.attributes),
                change_bl(sub_edge_trees[1].tree, sub_edge_trees[1].scale, sub_edge_trees[1].stem_rate[1]), 0, sub_edge_trees[1].tip_id)
        end
    end

    n_edges_left = tree.offsprings[1].n_nodes
    left = deepcopy(edge_trees[1:n_edges_left])
    right = deepcopy(edge_trees[(n_edges_left+1):(tree.n_nodes - 1)])
    tree_left = aux(tree.offsprings[1], left)
    tree_right = aux(tree.offsprings[2], right)
    return Tree([tree_left, tree_right], tree.branch_length, tree.attributes)
end

function node_depths_base(tree::Tree, edge_trees::Array{EdgeTreeRates2,1})
    function aux(subtree, node_times, current_depth, att)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length/att, att)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length/att, att)
            pushfirst!(node_times,current_depth)
        end
    end

    node_times = node_depths_base(tree)[2:end]
    node_times_full = Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, node_times_full, node_times[i], edge_trees[i].scale)
    end
    #pushfirst!(node_times,0.)
    return node_times_full
end


function node_depths(tree::Tree, edge_trees::Array{EdgeTreeRates2,1})
    function aux(subtree, node_times, current_depth, sa)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth + subtree.branch_length/sa)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length/sa, sa)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length/sa, sa)
            pushfirst!(node_times,current_depth+ subtree.branch_length/sa)
        end
    end

    node_times = node_depths_base(tree)[2:end]
    node_times_full = Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, node_times_full, node_times[i], edge_trees[i].scale)
    end
    #pushfirst!(node_times,0.)
    return node_times_full
end

function node_depths_base(tree::Tree)
    function aux(subtree, node_times, current_depth)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length)
            pushfirst!(node_times,current_depth)
        end
    end

    node_times = Array{Float64,1}(undef,0)
    aux(tree, node_times, 0.)
    #pushfirst!(node_times,0.)
    return node_times
end

function node_depths(tree::Tree)
    function aux(subtree, node_times, current_depth)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth + subtree.branch_length)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length)
            pushfirst!(node_times,current_depth+ subtree.branch_length)
        end
    end

    node_times = Array{Float64,1}(undef,0)
    aux(tree, node_times, 0.)
    #pushfirst!(node_times,0.)
    return node_times
end


function graft_to_tip(tree,subtree, tip_id)

    function aux(t1, t2, tip)
        if t1.n_nodes == 1
            return Tree(t2.offsprings, t2.branch_length + t1.branch_length, t1.attributes, t2.n_nodes, t2.extant)
        else
            t11 = t1.offsprings[1]
            t12 = t1.offsprings[2]
            t11_tips = (t11.n_nodes + 1) / 2
            if t11_tips < tip
                t12 = aux(t12, t2, tip - t11_tips)
                return Tree([t11, t12], t1.branch_length, t1.attributes)
            else
                t11 = aux(t11, t2, tip)
                return Tree([t11, t12], t1.branch_length, t1.attributes)
            end
        end
    end

    aux(tree,subtree, tip_id)
end

function graft_to_edge(tree, sub_tree, edge, tip)

    if sub_tree.n_nodes < 2
        return Tree(tree.offsprings, tree.branch_length + sub_tree.branch_length, tree.attributes, tree.n_nodes, sub_tree.extant)
    end

    function aux(t1, t2, edge_id)
        if t1.n_nodes == 1
            return t2
        elseif edge_id == 1
            new_t1 = Tree(t1.offsprings, 0., t1.attributes, t1.n_nodes, t1.extant)
            return graft_to_tip(t2, new_t1, tip)
        else
            t11 = t1.offsprings[1]
            t12 = t1.offsprings[2]
            if t11.n_nodes < (edge_id-1)
                t12 = aux(t12, t2, edge_id-1-t11.n_nodes)
                return Tree([t11, t12], t1.branch_length, t1.attributes)
            else
                t11 = aux(t11, t2, edge_id-1)
                return Tree([t11, t12], t1.branch_length, t1.attributes)
            end
        end
    end

    return aux(tree, sub_tree, edge+1)
end

function graft_to_tips(tree, subtree_list, tips)
    if length(tips) == 0
        return tree
    end

    ntips = Int64((tree.n_nodes + 1) / 2)
    subtrees = fill(Tree(),ntips)
    #println("$tips, $ntips")

    for i in 1:ntips
        for j in 1:length(tips)
            if i==tips[j]
                subtrees[i] = subtree_list[j]
                break
            end
        end
    end


    function aux(t1, t2)
        if t1.n_nodes <= 1
            #=if ! t2[1].extant
                println("is it ok ? $(t2[1]) ; $(t1.extant)")
            end=#
            return Tree(t2[1].offsprings, t2[1].branch_length + t1.branch_length, t1.attributes, t2[1].extant)
        else
            t11 = t1.offsprings[1]
            t12 = t1.offsprings[2]
            t11_tips = Int64((t11.n_nodes + 1) / 2)

            t11 = aux(t11, t2[1:t11_tips])
            t12 = aux(t12, t2[(t11_tips+1):end])

            return Tree([t11, t12], t1.branch_length, t1.attributes)
        end
    end

    new_tree = aux(tree,subtrees)
    return new_tree
end

function init_edge_tree_rates2(tree, rates)
    bl = extract_branch_lengths(tree)
    edge_trees = Array{EdgeTreeRates2,1}(undef,0)
    lefts = n_left(tree)
    for i in 2:tree.n_nodes
        nd = get_node_depth(tree, i-1)
        parent_edge = get_parent_edge(tree, i-1)

        ne = 0
        if lefts[i]==0
            ne = 1
        end

        push!(edge_trees,
            EdgeTreeRates2(Tree(Array{Tree,1}(undef,0), bl[i] * rates[i], 1.),       # tree
                rates[i],                                                           # stem_rate
                1.,                                                                 # tip_rate
                1,                                                                  # tip_id
                [1.],                                                               # tip_rates
                bl[i],                                                              # effective_bl
                ne,                                                                  # tip_number
                nd,
                nd + bl[i],
                parent_edge))
    end
    return edge_trees
end

function sim_ClaDS2_time(root_age,σ,α,ε,λ0,u,sf,lf ; return_if_extinct = true, max_node_number = Inf, make_tree = true, prune_extinct = false, return_if_max = true)
    # accesory functions that will be called latter
    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ)
        λ * rand(lambda_law, 2)
    end
    #=function new_rates(λ::Float64)
        X = (rand(lambda_law, 2))
        λ * X#* sort2!(X)
    end=#
    #println()
    function aux_ext(time, rate, n_max, s)
        #println(n_max)
        if n_max == 0
            #println("max")
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0, -Inf
        end
        node_time = randexp()/(rate*(1+ε))
        if node_time <= 0
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0, -Inf
        end
        if node_time > time
            s += lf
            #print(" $s $u ;")
            if u<s
                return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max + 1, s
            else
                #print(" !")
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0, s
            end
        else
                is_dead = dead()
                left_time = time - node_time
                if time == left_time
                    return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
                end
                if is_dead
                    #print("dead")
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1, s
                else
                    offspring_rates = new_rates(rate)
                    if offspring_rates[1] > offspring_rates[2]
                        offspring_rates = [offspring_rates[2],offspring_rates[1]]
                    end
                    #print("$depth ;")
                    left_tree = aux_ext(left_time, offspring_rates[1], n_max - 1, s)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_ext(left_time, offspring_rates[2], left_tree[2], left_tree[3])
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2], right_tree[3]
                end
        end
    end
    tree = aux_ext(root_age, λ0, max_node_number, sf)


    if prune_extinct
        if (tree[1].n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), -Inf
        else
            return prune_extinct_lineages(tree[1]), tree[3]
        end
    else
        if (tree[1].n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), -Inf
        else
            return tree[1], tree[3]
        end
    end
end

function sim_ClaDS2_time_rates(root_age,σ,α,ε,λ0 ; return_if_extinct = true, max_node_number = 1_000,
    make_tree = true, prune_extinct = false, return_if_max = true, return_na = false)

    # accesory functions that will be called latter
    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ)
        λ * rand(lambda_law, 2)
    end
    #=function new_rates(λ::Float64)
        X = (rand(lambda_law, 2))
        λ * X#* sort2!(X)
    end=#
    if return_na
        function aux_na(time, rate, n_max, rates, live)
            #println(n_max)
            if n_max == 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
            end
            node_time = randexp()/(rate*(1+ε))
            if node_time <= 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
            end
            if node_time > time
                #offsprings, branch_length, attributes, n_nodes
                push!(rates, rate)
                push!(live, true)
                return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max
            else
                is_dead = dead()
                if is_dead
                    #print("dead")
                    push!(live, false)
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    if time == left_time
                        return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
                    end
                    #print("$depth ;")
                    left_tree = aux_na(left_time, offspring_rates[1], n_max - 1, rates, live)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_na(left_time, offspring_rates[2], left_tree[2], rates, live)
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2]
                end
            end
        end

        rates = Array{Float64,1}(undef,0)
        alive = Array{Bool,1}(undef,0)

        tree = aux_na(root_age, λ0, max_node_number, rates, alive)[1]
        if (tree.n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), Array{Float64,1}(undef,0), Array{Bool,1}(undef,0)
        else
            return tree, rates, alive
        end
    else
        function aux_ext(time, rate, n_max, rates)
            #println(n_max)
            if n_max == 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
            end
            node_time = randexp()/(rate*(1+ε))
            if node_time > time
                #offsprings, branch_length, attributes, n_nodes
                push!(rates, rate)
                return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max
            else
                is_dead = dead()
                if is_dead
                    #print("dead")
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    #print("$depth ;")
                    left_tree = aux_ext(left_time, offspring_rates[1], n_max - 1, rates)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_ext(left_time, offspring_rates[2], left_tree[2], rates)
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2]
                end
            end
        end

        rates = Array{Float64,1}(undef,0)

        tree = aux_ext(root_age, λ0, max_node_number, rates)[1]
        if (tree.n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), Array{Float64,1}(undef,0)
        else
            return tree, rates
        end
    end

end


function sim_ClaDS2_time_unsampled_rates(root_age,σ,α,ε,λ0, u,fl ;
    sampling_proba = 1., return_if_extinct = true, max_node_number = 1_000,
    make_tree = false, prune_extinct = false)

    while true
        tree, rates, s = sim_ClaDS2_time_rates_tip(root_age,σ,α,ε,λ0,u,sampling_proba,fl,
            return_if_extinct = true, max_node_number = max_node_number,
            prune_extinct = false)
        if tree.n_nodes == -1
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), false, Inf, rates
        else
            n_alive = length(rates)
            if s < u
                return tree, false, n_alive, rates
            else
                return tree, true, n_alive, rates
            end
        end
    end
end


function sim_ClaDS2_time_rates_tip(root_age,σ,α,ε,λ0,u,sf,fl ; return_if_extinct = true, max_node_number = 1_000,
    make_tree = true, prune_extinct = false, return_if_max = true, return_na = false)

    # accesory functions that will be called latter
    #print(" $(ε/(1+ε))")
    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ)
        λ * rand(lambda_law, 2)
    end
    #=function new_rates(λ::Float64)
        X = (rand(lambda_law, 2))
        λ * X#* sort2!(X)
    end=#
    if return_na
    else
        function aux_ext(time, rate, n_max, rates, s ,n)
            if n_max == 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0,0,0
            end
            node_time = randexp()/(rate*(1+ε))
            #println(node_time)
            if node_time <= 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0,0,0
            end#print(" $s $n $u;")# $time $node_time $rate;")
            if node_time > time
                #print("zzzz")
                new_s = s
                #offsprings, branch_length, attributes, n_nodes
                push!(rates, rate)
                if n>0
                    new_s = s + log(1-sf)+log(n+1)-log(n)
                else
                    new_s = s + log(sf)
                end
                if new_s >= s
                    #print("a")
                    return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max, new_s , n+1
                elseif u < new_s
                    #print("b")
                    return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max, new_s , n+1
                else
                    #print("c")
                    return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0,new_s ,n+1
                end
            else
                is_dead = dead()
                if is_dead
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1, s, n
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    if time == left_time
                        return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0,0,0
                    end
                    #print("$depth ;")
                    left_tree = aux_ext(left_time, offspring_rates[1], n_max - 1, rates, s, n)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_ext(left_time, offspring_rates[2], left_tree[2], rates, left_tree[3], left_tree[4])
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2], right_tree[3], right_tree[4]
                end
            end
        end

        rates = Array{Float64,1}(undef,0)

        tree = aux_ext(root_age, λ0, max_node_number, rates,fl,0)
        #if ((tree[1].n_nodes == -1) & return_if_max) | (tree[4] == 0) | (u > tree[3])
        #    return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), Array{Float64,1}(undef,0)
        #else
        return tree[1], rates, tree[3]
        #end
    end

end

function update_edges_ETR2!(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1},
    lefts::Array{Int64,1}; it_rates = 1, prior_ε = "lognormal", logε0 = 0., sd = 0.5)

    if prior_ε == "uniform"
        ε = draw_ε_crown_priorUnif(tree, edge_trees, lefts)
    elseif prior_ε == "uniformInf"
        ε = draw_ε_crown(tree, edge_trees, lefts)
    elseif prior_ε == "ClaDS0"
        ε = 0.
    else
        ε = draw_ε_crown_priorln(tree, edge_trees, lefts, logε0 = logε0, sd = sd)
    end

    for j in 1:it_rates

        draw_λi_quad!(rates, edge_trees, σ, α, ε, tree, lefts)
        relative_rates = extract_relative_rates(tree, edge_trees, rates)

        σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
        α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)

        if prior_ε == "uniform"
            ε = draw_ε_crown_priorUnif(tree, edge_trees, lefts)
        elseif prior_ε == "uniformInf"
            ε = draw_ε_crown(tree, edge_trees, lefts)
        elseif prior_ε == "ClaDS0"
            ε = 0.
        else
            ε = draw_ε_crown_priorln(tree, edge_trees, lefts, logε0 = logε0, sd = sdε)
        end
    end

    return ε, σ, α
end
