
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

function sim_ClaDS2_time_unsampled(root_age,σ,α,ε,λ0 ;
    sampling_proba = 1., return_if_sampled = false, return_if_extinct = true, max_node_number = 1_000, not_sampled = true,
    make_tree = false, prune_extinct = false)

    while true
        tree = sim_ClaDS2_time(root_age,σ,α,ε,λ0 ,
            return_if_extinct = true, max_node_number = max_node_number,
            prune_extinct = false)
        if tree.n_nodes == -1
            if return_if_sampled
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), false, Inf
            end
        else
            n_alive = n_extant_tips(tree)
            if not_sampled
                if n_alive == 0
                    log_sampling_proba = 0
                else
                    log_sampling_proba = n_alive * log(1-sampling_proba)
                end
            else
                if n_alive == 0
                    log_sampling_proba = -Inf
                elseif n_alive == 1
                    log_sampling_proba = log(sampling_proba)
                else
                    log_sampling_proba = log(n_alive) + log(sampling_proba) + (n_alive - 1) * log(1-sampling_proba)
                end
            end

            u = log(rand())
            unsampled = u < log_sampling_proba
            if unsampled | return_if_sampled
                return tree, unsampled, n_alive
            end
        end
    end
end

function sim_ClaDS2_time_unsampled_rates(root_age,σ,α,ε,λ0 ;
    sampling_proba = 1., return_if_sampled = false, return_if_extinct = true, max_node_number = 1_000, not_sampled = true,
    make_tree = false, prune_extinct = false)

    while true
        tree, rates = sim_ClaDS2_time_rates(root_age,σ,α,ε,λ0 ,
            return_if_extinct = true, max_node_number = max_node_number,
            prune_extinct = false)
        if tree.n_nodes == -1
            if return_if_sampled
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), false, Inf, rates
            end
        else
            n_alive = length(rates)
            if not_sampled
                if n_alive == 0
                    log_sampling_proba = 0
                else
                    log_sampling_proba = n_alive * log(1-sampling_proba)
                end
            else
                if n_alive == 0
                    log_sampling_proba = -Inf
                elseif n_alive == 1
                    log_sampling_proba = log(sampling_proba)
                else
                    log_sampling_proba = log(n_alive) + log(sampling_proba) + (n_alive - 1) * log(1-sampling_proba)
                end
            end

            u = log(rand())
            unsampled = u < log_sampling_proba
            if unsampled | return_if_sampled
                return tree, unsampled, n_alive, rates
            end
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


function sim_ClaDS2_ntips_aux(n,σ,α,ε,λ0 ; return_if_extinct = false, make_tree = true, prune_extinct = true, tree_only = true)

    # accesory functions that will be called latter

    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ)
        λ * rand(lambda_law, 2)
    end

    while true

        # initialise internal variables
        current_node = 3
        alive = [1,2]
        rates = new_rates(λ0)
        branches = [[[0;1]]; [[0;2]]]
        branch_lengths = [0., 0.]
        dead_branches = [false,false]
        node_time = 0
        times = [0.]
        time_diffs = [0.]
        n_lineages = [1]

        time_int = randexp()/(sum(rates[alive])*(1+ε))
        node_time += time_int

        while 0 < length(alive) < n
            individual= splice!(alive, sample(eachindex(alive),Weights(rates[alive])))

            is_dead=dead()
            push!(times, node_time)
            push!(time_diffs, time_int)

            for node in alive
                branch_lengths[node] += time_int
            end
            if is_dead
                branch_lengths[individual] += time_int

                dead_branches[individual]=true
                push!(n_lineages,length(alive))
            else
                sample_rates = new_rates(rates[individual])
                push!(dead_branches,false,false)
                push!(rates, sample_rates[1], sample_rates[2])
                push!(branch_lengths, 0., 0.);
                push!(alive,current_node,current_node+1)
                push!(n_lineages,length(alive))
                push!(branches,[individual;current_node], [individual;current_node+1])

                current_node = current_node+2
                branch_lengths[individual] += time_int

            end

            time_int = randexp()/(sum(rates[alive])*(1+ε))
            node_time += time_int
            if time_int==0
            end
        end

        push!(times, node_time)
        push!(time_diffs, time_int)
        push!(n_lineages,length(alive))

        for i in alive
            branch_lengths[i] += time_int
        end
        if length(alive) > 0 || return_if_extinct
            if make_tree
                if tree_only
                    return build_tree(branches, branch_lengths, rates, extinct = dead_branches, prune_extinct = prune_extinct, root_attributes=[λ0])
                else
                    t = build_tree(branches, branch_lengths, rates, extinct = dead_branches, prune_extinct = prune_extinct, root_attributes=[λ0])
                    return t, times, n_lineages, time_diffs
                end
            else
                return branches, branch_lengths, rates, dead_branches, times, n_lineages
            end
        end
    end
end


#=

=#

function build_tree(branches::Array{Int64,2}, branch_lengths::Array{Float64,2}, attributes::Array{Array{T2,1},1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T2,1}(undef,0)) where {T2<:Number}
    n_attributes = length(attributes)
    if length(extinct) == 0
        dead_branches = fill(false, size(branches)[1])
    else
        dead_branches = extinct
    end

    if size(branches)[1] == 0
        if return_void
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),stem_age,root_attributes)
        end
    end
    if size(branches)[1] == 1
        if ! dead_branches[1]
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, true)
        elseif prune_extinct
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, false)
        end
    end

    root = -1
    for i in branches[:,1]
        has_no_parent = true
        for j in branches[:,2]
            if i == j
                has_no_parent = false
                break
            end
        end
        if has_no_parent
            root = i
            break
        end
    end

    max_node = maximum(branches[:,1:2])+1
    offsprings = fill(-1,max_node,2)
    parent_edges = fill(-1,max_node)
    daughter_edges = fill(-1,max_node)

    for i in 1:size(branches)[1]
        individual = branches[i,1]+1
        parent_edges[individual] = i
        daughter_edges[branches[i,2]+1] = i
        if offsprings[individual,1] < 0
            offsprings[individual,1] = branches[i,2]
        else
            offsprings[individual,2] = branches[i,2]
        end

        if dead_branches[i]
            individual = branches[i,2]+1
            offsprings[individual,:]= [-2 -2]
        end
    end

    if prune_extinct
        function aux(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes])
            elseif offsprings[node+1,1]==-2
                return Tree()
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes #repeat([Array{Float64}(undef,0)],n_attributes)
                end

                tree_left = aux(offsprings[node+1,1])
                tree_right = aux(offsprings[node+1,2])

                if tree_left.n_nodes + tree_right.n_nodes == 0
                    return Tree()
                elseif tree_left.n_nodes == 0
                    tree_right = Tree(tree_right.offsprings, tree_right.branch_length + add_time,
                        add_attribute, tree_right.n_nodes)

                    return tree_right
                elseif tree_right.n_nodes == 0
                    tree_left = Tree(tree_left.offsprings, tree_left.branch_length + add_time,
                        add_attribute, tree_left.n_nodes)
                    return tree_left
                else
                    return Tree([tree_left, tree_right], add_time, add_attribute)
                end
            end
        end

        return aux(root)

    else
        function aux_all(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], true)
            elseif offsprings[node+1,1]==-2
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], false)
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes#repeat([Array{Float64}(undef,0)],n_attributes)
                end
                tree_left = aux_all(offsprings[node+1,1])
                tree_right = aux_all(offsprings[node+1,2])
                return Tree([tree_left, tree_right], add_time, add_attribute)
            end
        end

        return aux_all(root)

    end
end

function build_tree(branches::Array{Array{T1,1}}, branch_lengths::Array{Float64,2}, attributes::Array{Array{T2,1},1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T2,1}(undef,0)) where {T1<:Number, T2<:Number}

    n_attributes = length(attributes)
    if length(extinct) == 0
        dead_branches = fill(false, size(branches)[1])
    else
        dead_branches = extinct
    end

    if size(branches)[1] == 0
        if return_void
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),stem_age,root_attributes)
        end
    end
    if size(branches)[1] == 1
        if ! dead_branches[1]
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, true)
        elseif prune_extinct
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, false)
        end
    end

    root = -1
    for b1 in branches
        has_no_parent = true
        for b2 in branches
            if b1[1] == b2[2]
                has_no_parent = false
                break
            end
        end
        if has_no_parent
            root = b1[1]
            break
        end
    end

    max_node = -1
    for b in branches
        max_node = max(b[1], b[2], max_node)
    end
    max_node += 1
    offsprings = fill(-1,max_node,2)
    parent_edges = fill(-1,max_node)
    daughter_edges = fill(-1,max_node)

    for i in 1:length(branches)
        individual = branches[i][1]+1
        parent_edges[individual] = i
        daughter_edges[branches[i][2]+1] = i
        if offsprings[individual,1] < 0
            offsprings[individual,1] = branches[i][2]
        else
            offsprings[individual,2] = branches[i][2]
        end

        if dead_branches[i]
            individual = branches[i][2]+1
            offsprings[individual,:]= [-2 -2]
        end
    end

    if prune_extinct
        function aux(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes])
            elseif offsprings[node+1,1]==-2
                return Tree()
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes #repeat([Array{Float64}(undef,0)],n_attributes)
                end

                tree_left = aux(offsprings[node+1,1])
                tree_right = aux(offsprings[node+1,2])

                if tree_left.n_nodes + tree_right.n_nodes == 0
                    return Tree()
                elseif tree_left.n_nodes == 0
                    tree_right = Tree(tree_right.offsprings, tree_right.branch_length + add_time,
                        add_attribute, tree_right.n_nodes)

                    return tree_right
                elseif tree_right.n_nodes == 0
                    tree_left = Tree(tree_left.offsprings, tree_left.branch_length + add_time,
                        add_attribute, tree_left.n_nodes)
                    return tree_left
                else
                    return Tree([tree_left, tree_right], add_time, add_attribute)
                end
            end
        end

        return aux(root)

    else
        function aux_all(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], true)
            elseif offsprings[node+1,1]==-2
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], false)
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes#repeat([Array{Float64}(undef,0)],n_attributes)
                end
                tree_left = aux_all(offsprings[node+1,1])
                tree_right = aux_all(offsprings[node+1,2])
                return Tree([tree_left, tree_right], add_time, add_attribute)
            end
        end

        return aux_all(root)

    end

end

function build_tree(branches::Array{Int64,2}, branch_lengths::Array{Float64,2}, attributes::Array{Array{T2,1},1}, tip_labels::Array{String,1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T2,1}(undef,0)) where {T2<:Number}

    n_attributes = length(attributes)
    if length(extinct) == 0
        dead_branches = fill(false, size(branches)[1])
    else
        dead_branches = extinct
    end

    if size(branches)[1] == 0
        if return_void
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),stem_age,root_attributes)
        end
    end
    if size(branches)[1] == 1
        if ! dead_branches[1]
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, true)
        elseif prune_extinct
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, false)
        end
    end

    root = -1
    for i in branches[:,1]
        has_no_parent = true
        for j in branches[:,2]
            if i == j
                has_no_parent = false
                break
            end
        end
        if has_no_parent
            root = i
            break
        end
    end

    max_node = maximum(branches[:,1:2])+1
    offsprings = fill(-1,max_node,2)
    parent_edges = fill(-1,max_node)
    daughter_edges = fill(-1,max_node)

    for i in 1:size(branches)[1]
        individual = branches[i,1]+1
        parent_edges[individual] = i
        daughter_edges[branches[i,2]+1] = i
        if offsprings[individual,1] < 0
            offsprings[individual,1] = branches[i,2]
        else
            offsprings[individual,2] = branches[i,2]
        end

        if dead_branches[i]
            individual = branches[i,2]+1
            offsprings[individual,:]= [-2 -2]
        end
    end

    if prune_extinct
        function aux(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], tip_labels[node])
            elseif offsprings[node+1,1]==-2
                return Tree()
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes #repeat([Array{Float64}(undef,0)],n_attributes)
                end

                tree_left = aux(offsprings[node+1,1])
                tree_right = aux(offsprings[node+1,2])

                if tree_left.n_nodes + tree_right.n_nodes == 0
                    return Tree()
                elseif tree_left.n_nodes == 0
                    tree_right = Tree(tree_right.offsprings, tree_right.branch_length + add_time,
                        add_attribute, tree_right.n_nodes, tree_right.label)

                    return tree_right
                elseif tree_right.n_nodes == 0
                    tree_left = Tree(tree_left.offsprings, tree_left.branch_length + add_time,
                        add_attribute, tree_left.n_nodes, tree_left.label)
                    return tree_left
                else
                    return Tree([tree_left, tree_right], add_time, add_attribute)
                end
            end
        end

        return aux(root)

    else
        function aux_all(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], true, tip_labels[node])
            elseif offsprings[node+1,1]==-2
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], false)
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes#repeat([Array{Float64}(undef,0)],n_attributes)
                end
                tree_left = aux_all(offsprings[node+1,1])
                tree_right = aux_all(offsprings[node+1,2])
                return Tree([tree_left, tree_right], add_time, add_attribute)
            end
        end

        return aux_all(root)

    end
end

function build_tree(branches::Array{Int64,2}, branch_lengths::Array{Float64,2}, attribute::Array{T,1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T,1}(undef,0)) where {T<:Number}
    build_tree(branches, branch_lengths, [attribute], extinct=extinct, prune_extinct=prune_extinct,
        return_void = return_void, stem_age = stem_age, root_attributes= root_attributes)
end

function build_tree(branches::Array{Int64,2}, branch_lengths::Array{Float64,2}; extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T,1}(undef,0))
    attributes = [1:length(branch_lengths)...]
    build_tree(branches, branch_lengths, [attributes], extinct=extinct,
        prune_extinct=prune_extinct, return_void = return_void, stem_age = stem_age, root_attributes= root_attributes)
end





struct EdgeTree
    tree::Tree
    tip_id::Int64
    n::Int64
end

struct EdgeTreeRates
    tree::Tree
    tip_id::Int64 #for tip edges it is the number of alive species
    rate::Float64
    rates::Array{Float64,1}
end

function rasterize(tree::Tree)
    function aux(subtree)
        if length(subtree.offsprings) == 0
            return subtree
        else
            left = aux(subtree.offsprings[1])
            right = aux(subtree.offsprings[2])

            if left.n_nodes < right.n_nodes
                return Tree([left,right], subtree.branch_length,subtree.attributes,subtree.n_nodes,subtree.extant,subtree.label)
            else
                return Tree([right,left], subtree.branch_length,subtree.attributes,subtree.n_nodes,subtree.extant,subtree.label)
            end
        end
    end

    return aux(tree)
end

function n_tip(tree::Tree)
    return Int64((tree.n_nodes+1)/2)
end

function n_tip(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    n = n_tip(tree)
    for i in 1:length(edge_trees)
        n += n_tip(edge_trees[i].tree) - 1
    end
    return n
end

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

function n_extant_tips(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    n = n_extant_tips(tree)
    for i in 1:length(edge_trees)
        n += n_extant_tips(edge_trees[i].tree) - 1
    end
    return n
end

function n_extinct(tree::Tree)
    n_tip(tree) - n_extant_tips(tree)
end

function n_extinct(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    n = 0
    for i in 1:length(edge_trees)
        n += n_extinct(edge_trees[i].tree)
    end
    return n
end

function is_tip(tree::Tree)
    length(tree.offsprings) == 0
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


function add_extinct(tree::Tree ; tol = 0.0001)
    function max_depth(tr)
        if length(tr.offsprings) == 0
            return tr.branch_length
        else
            return tr.branch_length + max(max_depth(tr.offsprings[2]),max_depth(tr.offsprings[1]))
        end
    end

    rd = max_depth(tree)
    function aux(tr,depth)
        if length(tr.offsprings) == 0
            if (depth + tr.branch_length) >= (rd * (1-tol))
                return Tree(tr.offsprings, tr.branch_length, tr.attributes, tr.n_nodes,
                    true, tr.label)
            else
                return Tree(tr.offsprings, tr.branch_length, tr.attributes, tr.n_nodes,
                    false, tr.label)
            end
        else
            t1 = aux(tr.offsprings[1], depth + tr.branch_length)
            t2 = aux(tr.offsprings[2], depth + tr.branch_length)
            return Tree([t1,t2], tr.branch_length, tr.attributes, tr.n_nodes,
                tr.label)
        end
    end

    aux(tree,0.)
end

function extant_tips(tree::Tree)
    function aux(subtree, x)
        if length(subtree.offsprings) == 0 && subtree.extant
            pushfirst!(x,true)
        elseif length(subtree.offsprings) == 0
            pushfirst!(x,false)
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

function extract_rates(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}; id = 1)
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

    rates =  Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, rates)
    end
    pushfirst!(rates, tree.attributes[1])
    return rates
end

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

function extract_branch_lengths(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})

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
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree,bl)
    end
    pushfirst!(bl, tree.branch_length)
    return bl
end

function which_extant(tree::Tree)
    function aux(subtree, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            if subtree.extant
                pushfirst!(x, 1)
            else
                pushfirst!(x, 0)
            end
        else
            aux(subtree.offsprings[2],x)
            aux(subtree.offsprings[1],x)
            if subtree.extant
                pushfirst!(extant, 1)
            else
                pushfirst!(extant, 0)
            end
        end
    end

    x = Array{Int64,1}(undef,0)
    aux(tree, x)
    return x
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

function get_parent_edges(tree::Tree)
    parents = Array{Int64,1}(undef,0)
    for i in 1:(tree.n_nodes-1)
        push!(parents,get_parent_edge(tree, i))
    end
    return parents
end

function get_parent_rate(tree::Tree, edge_id::Int64, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    parent_edge = get_parent_edge(tree, edge_id)
    if parent_edge == 0
        return tree.attributes[1]
    else
        if edge_trees[parent_edge].tree.n_nodes < 2
            return edge_trees[parent_edge].tree.attributes[1]
        else
            return extract_tip_rates(edge_trees[parent_edge].tree, edge_trees[parent_edge].tip_id)
        end
    end
end

function get_parent_rate(tree::Tree, edge_id::Int64, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, rates::Array{Float64,1})
    parent_edge = get_parent_edge(tree, edge_id)
    if parent_edge == 0
        return rates[1]
    else
        if edge_trees[parent_edge].tree.n_nodes < 2
            return rates[parent_edge+1]#
        else
            return extract_tip_rates(edge_trees[parent_edge].tree, edge_trees[parent_edge].tip_id)
        end
    end
end

function get_parent_rate(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, rates::Array{Float64,1}, parent_edge::Int64)
    if parent_edge == 0
        return rates[1]
    else
        if edge_trees[parent_edge].tree.n_nodes < 2
            return rates[parent_edge+1]#
        else
            return extract_tip_rates(edge_trees[parent_edge].tree, edge_trees[parent_edge].tip_id)
            #tip_rates = extract_tip_rates(edge_trees[parent_edge].tree)
            #return tip_rates[edge_trees[parent_edge].tip_id]
        end
    end
end

function get_parent_rates(tree::Tree; id = 1)
    function aux(subtree, root_rate, rates)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(rates,root_rate)
        else
            new_root = deepcopy(subtree.attributes[id])
            aux(subtree.offsprings[2],new_root,rates)
            aux(subtree.offsprings[1],new_root,rates)
            pushfirst!(rates,root_rate)
        end
    end

    root_rate = tree.attributes[id]
    x = Array{Float64, 1}(undef,0)
    aux(tree,root_rate,x)
    return(x)
end

function extract_relative_rates(tree::Tree; id = 1)
    function aux(subtree, parent_rate, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
        else
            log_rate = deepcopy(log(subtree.attributes[id]))
            aux(subtree.offsprings[2], log_rate, x)
            aux(subtree.offsprings[1], log_rate, x)
            pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
        end
    end

    x = Array{Float64,1}(undef,0)
    aux(tree, 0., x)
    return x
end

function extract_relative_rates(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, rates; id = 1)
    function aux(subtree, parent_rate, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
        else
            log_rate = deepcopy(log(subtree.attributes[id]))
            aux(subtree.offsprings[2], log_rate, x)
            aux(subtree.offsprings[1], log_rate, x)
            pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
        end
    end

    relative_rates = []
    for edge_id in 1:length(edge_trees)
        parent_rate = log(get_parent_rate(tree, edge_id, edge_trees))
        r = log(rates[edge_id+1]) - parent_rate
        if isnan(r) || r == -Inf
            println("$edge_id, $r, $(rates[edge_id+1]), $(get_parent_rate(tree, edge_id, edge_trees))")
        end
        if edge_trees[edge_id].tree.n_nodes == 1
            pushfirst!(relative_rates,r)
        else
            aux(edge_trees[edge_id].tree,parent_rate, relative_rates)
        end
    end
    return relative_rates
end

function extract_relative_rates(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}; id = 1)
    rates = extract_rates(tree)
    extract_relative_rates(tree, edge_trees, rates, id = id)
end

function extract_partial_relative_rates(tree::Tree, edge_trees, rates; id = 1)

    relative_rates = []
    for edge_id in 1:length(edge_trees)
        parent_rate = log(get_parent_rate(tree, edge_id, edge_trees))
        r = log(rates[edge_id+1]) - parent_rate
        if isnan(r) || r == -Inf
            println("$edge_id, $r, $(rates[edge_id+1]), $(get_parent_rate(tree, edge_id, edge_trees))")
        end
        pushfirst!(relative_rates,r)
    end
    return relative_rates
end

function extract_partial_relative_rates(tree::Tree, edge_trees; id = 1)
    rates = extract_rates(tree)
    extract_partial_relative_rates(tree, edge_trees, rates, id = id)
end

function make_ape(tree::Tree ; id = 1)
    n_nodes = tree.n_nodes
    node_id = (n_nodes + 3)/2

    function aux(subtree, root, next_node, tip)
        if subtree.n_nodes == 0
            return Array{Int64,2}(undef,0,2), Array{Float64,1}(undef,0), Array{Float64,1}(undef,0), next_node, tip, Array{String,1}(undef,0)
        elseif length(subtree.offsprings) == 0
            return [root tip], subtree.branch_length, subtree.attributes[id], next_node, tip+1, subtree.label
        else
            tupple1 = aux(subtree.offsprings[1], next_node, next_node+1, tip)
            tupple2 = aux(subtree.offsprings[2], next_node, tupple1[4], tupple1[5])
            if subtree.branch_length == 0.
                return [tupple1[1]; tupple2[1]],
                    [tupple1[2]; tupple2[2]],
                    [tupple1[3]; tupple2[3]],
                    tupple2[4], tupple2[5],
                    [tupple1[6]; tupple2[6]]
            else
                return [[root next_node]; tupple1[1]; tupple2[1]],      # edges
                    [subtree.branch_length; tupple1[2]; tupple2[2]],    # branch lengths
                    [subtree.attributes[id]; tupple1[3]; tupple2[3]],   # rates
                    tupple2[4], tupple2[5],                             # auxiliary variables
                    [tupple1[6]; tupple2[6]]                            # tip labels
            end
        end
    end

    tupple = aux(tree, -1, node_id, 1)
    return tupple[1], tupple[2], tupple[3], tupple[6]
end

function Rsave_tree(tree::Tree ; file = "tree.Rdata")
    edges, branch_lengths, rates, tip_labels = make_ape(tree)
    ntip = (tree.n_nodes + 1)/2

    @rput edges
    @rput branch_lengths
    @rput rates
    @rput ntip
    @rput file
    @rput tip_labels

    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
        save(tree, rates, file = file)
    """)
end

function get_daughter_edges(tree, edge_id, lefts)

    function aux(subtree, edge, add, x)
        if subtree.n_nodes <= 1
            pushfirst!(x, edge + add)
        elseif edge == 0
            aux(subtree.offsprings[2], 0, add + 1 + lefts[add + 1], x)
            aux(subtree.offsprings[1], 0, add + 1, x)
            pushfirst!(x, add)
        elseif edge < lefts[add + 1]
            aux(subtree.offsprings[1], edge - 1, add + 1, x)
        else
            aux(subtree.offsprings[2], edge - 1 - lefts[add + 1], add + lefts[add + 1] + 1, x)
        end
    end

    x = Array{Int64,1}(undef,0)
    aux(tree, edge_id, 0,x)
    return x
end

function get_daughter_edges(tree, edge_id)
    lefts = n_left(tree)

    get_daughter_edges(tree, edge_id, lefts)
end

function daughters(tree, lefts)
    daughter_edges = Array{Array{Int64,1},1}(undef,tree.n_nodes-1)
    for i in 2:tree.n_nodes
        daughter_edges[i-1] = get_daughter_edges(tree, i-1, lefts)
    end

    return daughter_edges
end

function daughters(tree)
    lefts = n_left(tree)
    daughter_edges = Array{Array{Int64,1},1}(undef,tree.n_nodes)
    for i in 2:tree.n_nodes
        daughter_edges[i-1] = get_daughter_edges(tree, i-1, lefts)
    end

    return daughter_edges
end

function get_daughter_rates(tree, edge_id ; rate_id = 1)

    function aux(subtree, edge)
        if edge == 1
            if subtree.n_nodes < 2
                return [NaN,NaN]
            else
                return [subtree.offsprings[1].attributes[rate_id] ; subtree.offsprings[2].attributes[rate_id]]
            end
        else
            t11 = subtree.offsprings[1]
            t12 = subtree.offsprings[2]
            if t11.n_nodes < (edge-1)
                return aux(t12, edge-1-t11.n_nodes)
            else
                return aux(t11, edge-1)
            end
        end
    end

    aux(tree, edge_id + 1)
end

function get_daughter_rates(tree, edge_id, edge_trees ; rate_id = 1)

    if edge_id == 0
        return [tree.offsprings[1].attributes[rate_id] ; tree.offsprings[2].attributes[rate_id]]
    elseif edge_trees[edge_id].tree.n_nodes > 1
        return [edge_trees[edge_id].tree.offsprings[1].attributes[rate_id] ; edge_trees[edge_id].tree.offsprings[2].attributes[rate_id]]
    end

    function aux(subtree, edge)
        if edge == 1
            if subtree.n_nodes < 2
                return [NaN,NaN]
            else
                return [subtree.offsprings[1].attributes[rate_id] ; subtree.offsprings[2].attributes[rate_id]]
            end
        else
            t11 = subtree.offsprings[1]
            t12 = subtree.offsprings[2]
            if t11.n_nodes < (edge-1)
                return aux(t12, edge-1-t11.n_nodes)
            else
                return aux(t11, edge-1)
            end
        end
    end

    aux(tree, edge_id + 1)
end

function get_daughter_rates(tree, edge_id, edge_trees, rates, lefts ; rate_id = 1, with_edge_tree = true)

    edge = edge_id +1
    if edge_id == 0
        return [rates[edge+1]; rates[edge+1+lefts[edge]]]
    elseif edge_trees[edge_id].tree.n_nodes > 1 && with_edge_tree
        return [edge_trees[edge_id].tree.offsprings[1].attributes[rate_id] ; edge_trees[edge_id].tree.offsprings[2].attributes[rate_id]]
    elseif lefts[edge]==0
        return [NaN,NaN]
    else
        return [rates[edge+1]; rates[edge+1+lefts[edge]]]
    end
end

function get_daughter_rates(tree, edge_id, edge_trees, rates ; rate_id = 1)

    lefts = n_left(tree)
    get_daughter_rates(tree, edge_id, edge_trees, rates, lefts, rate_id = rate_id)
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

function get_base_node_depth(tree, edge_id)
    function aux(subtree, edge)
        if subtree.n_nodes < 2
            return subtree.branch_length
        else
            if edge == 1
                return aux(subtree.offsprings[1], edge) + subtree.branch_length
            elseif (edge - 1) > subtree.offsprings[1].n_nodes
                return aux(subtree.offsprings[2], edge - 1 - subtree.offsprings[1].n_nodes)
            else
                return aux(subtree.offsprings[1], edge - 1)
            end
        end
    end

    aux(tree, edge_id + 1)
end

function extract_tip_rates(tree::Tree ; id = 1, return_extinct = true)

    function aux(subtree, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            if subtree.extant || return_extinct
                pushfirst!(x,subtree.attributes[id])
            else
                pushfirst!(x,NaN)
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

function extract_live_tip_rates(tree::Tree ; id = 1)

    function aux(subtree, x)
        if length(subtree.offsprings) == 0
            if subtree.extant
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

function extract_tip_rates_bu1(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, tips_id::Array{Bool,1} ; id = 1, return_extinct = true)
    function aux(subtree,x)
        if length(subtree.offsprings)==0
            push!(x,subtree.attributes[id])
        else
            if subtree.offsprings[1].extant
                aux(subtree.offsprings[1],x)
            else
                aux(subtree.offsprings[2],x)
            end
        end
    end
    tip_rates = Array{Float64,1}(undef,0)
    for i in 1:length(edge_trees)
        if tips_id[i]
            aux(edge_trees[i].tree, tip_rates)
        end
    end
    return tip_rates
end



function extract_tip_rates(tree::Tree, edge_trees::Array{EdgeTree,1}, tips_id::Array{Bool,1}, rates::Array{Float64,1} ; id = 1, return_extinct = true)

    tip_rates = Array{Float64,1}(undef,0)
    for i in 1:length(edge_trees)
        if tips_id[i]
            push!(tip_rates,sample(extract_live_tip_rates(edge_trees[i].tree)))
        end
    end
    return tip_rates
end

function extract_tip_rates(tree::Tree, edge_trees::Array{EdgeTreeRates,1},
    tips_id::Array{Bool,1}, rates::Array{Float64,1} ; id = 1, return_extinct = true)
    tip_rates = Array{Float64,1}(undef,0)
    for i in 1:length(edge_trees)
        if tips_id[i]
            r = sample(edge_trees[i].rates)*rates[i+1]/edge_trees[i].rate
            push!(tip_rates,r)
        end
    end
    return tip_rates
end

function extract_tip_rates(tree::Tree, tip_id::Int64 ; id = 1, return_extinct = true)

    function aux(subtree, tip)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            if subtree.extant || return_extinct
                return subtree.attributes[id]
            else
                return NaN
            end
        elseif ((subtree.offsprings[1].n_nodes+1)/1) >= tip
            aux(subtree.offsprings[1], tip)
        else
            aux(subtree.offsprings[2], tip - (subtree.offsprings[1].n_nodes+1)/1)
        end
    end

    return aux(tree, tip_id * 2)
end

function cut_tree(tree::Tree, time; prune_extinct = false) #time from the root

    function aux(subtree, sub_time)
        if subtree.n_nodes < 2
            if subtree.branch_length < sub_time
                if prune_extinct && ! subtree.extant
                    return Tree()
                else
                    return subtree
                end
            else
                return Tree(subtree.offsprings, sub_time, subtree.attributes, subtree.n_nodes, true)
            end
        else
            if subtree.branch_length > sub_time
                return Tree(Array{Tree,1}(undef,0), sub_time, subtree.attributes)
            else
                left_tree = aux(subtree.offsprings[1], sub_time - subtree.branch_length)
                right_tree = aux(subtree.offsprings[2], sub_time - subtree.branch_length)
                if ! left_tree.extant && prune_extinct
                    if ! right_tree.extant
                        return Tree()
                    else
                        return Tree(right_tree.offsprings, right_tree.branch_length + subtree.branch_length,
                            subtree.attributes)
                    end
                elseif ! right_tree.extant && prune_extinct
                    return Tree(left_tree.offsprings, left_tree.branch_length + subtree.branch_length,
                        subtree.attributes)
                else
                    return Tree([left_tree, right_tree], subtree.branch_length, subtree.attributes)
                end
            end

        end
    end

    aux(tree, time)
end
