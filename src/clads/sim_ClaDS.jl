#=
Simulate a ClaDS2 tree conditionned on tip number
=#

"""
    sim_ClaDS2_ntips(n::Int64,σ::Float64,α::Float64,ε::Float64,λ0::Float64 ; ...)

Simulate a tree from the ClaDS model conditionned on the number of tips at present

# Arguments
- `n::Int64`: the goal number of tips at present.
- `σ::Float64`: the stochasticity parameter.
- `α::Float64`: the trend parameter.
- `ε::Float64`: turnover rate (extinction / speciation).
- `λ0::Float64`: initial speciation rate.

# Keyword arguments
- `prune_extinct::Bool`: Should extinct lineages be removed from the output tree? Default to false.
"""
function sim_ClaDS2_ntips(n::Int64,σ::Float64,α::Float64,ε::Float64,λ0::Float64 ;
    prune_extinct = true, sed = 0.001,
    max_time = 5, max_simulation_try = 100)

    # accesory functions that will be called latter

    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ)
        λ * rand(lambda_law, 2)
    end

    N = n*max_time                  # maximal number of species before killing the simulation
    s = λ0*sed
    breaked = 0

    while breaked < max_simulation_try

        # initialise the tree to one alive lineage with rate λ0
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
        scale = 1.

        while 0 < length(alive) < N
            time_int = randexp()/(sum(rates[alive])*(1+ε))          # time of the new event

            if time_int==0                                          # to avoid numerical issues
                breaked += 1
                break
            end

            if length(alive) == n
                samp = s * time_int                     # probability to sample the tree in this interval with n species
                node_time += time_int
                u = rand()

                if u < samp*scale                       # meaning the tree is sampled
                    time_int *= rand()                  # the end time is uniformly sampled along the branch

                    push!(times, node_time)
                    push!(time_diffs, time_int)
                    push!(n_lineages,length(alive))

                    for i in alive                      # set all extant branch end to hte present time
                        branch_lengths[i] += time_int
                    end

                    t = build_tree(branches, branch_lengths, rates, extinct = dead_branches, prune_extinct = prune_extinct, root_attributes=[λ0])
                    return t
                else
                    scale /= (1. - samp)                # so that sampling is uniform on all time points wit n species
                end
            else
                node_time += time_int
            end

            individual= splice!(alive, sample(eachindex(alive),Weights(rates[alive])))

            is_dead=dead()                                          # check if the new event is a death or a speciation
            push!(times, node_time)
            push!(time_diffs, time_int)

            for node in alive
                branch_lengths[node] += time_int
            end

            # update internal variables
            if is_dead
                branch_lengths[individual] += time_int
                dead_branches[individual] = true
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
        end

        if length(alive) > 0
            s *=2
        end
    end

    # if the function finished before a good tree was found, returns an empty tree
    return Tree()
end

#=
Simulate a ClaDS2 tree conditionned on time
=#

function sim_ClaDS2_time(root_age,σ,α,ε,λ0 ; return_if_extinct = true, max_node_number = Inf, make_tree = true, prune_extinct = false, return_if_max = true)

    # accesory functions that will be called latter
    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ::Float64)
        λ * rand(lambda_law, 2)
    end

    function aux_ext(time, rate, n_max)
        #println(n_max)
        if n_max == 0
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
        end
        node_time = randexp()/(rate*(1+ε))
        if node_time > time
                #offsprings, branch_length, attributes, n_nodes
                return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max
        else
                is_dead = dead()
                if is_dead
                    print("dead")
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    #print("$depth ;")
                    left_tree = aux_ext(left_time, offspring_rates[1], n_max - 1)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_ext(left_time, offspring_rates[2], left_tree[2])
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2]
                end
        end
    end
    tree = aux_ext(root_age, λ0, max_node_number)[1]


    if prune_extinct
        if (tree.n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false)
        else
            return prune_extinct_lineages(tree)
        end
    else
        if (tree.n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false)
        else
            return tree
        end
    end
end


#=
Simulate a ClaDS2 tree
=#
