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

function R_tree(tree::Tree ; file = "tree.Rdata")
    edges, branch_lengths, temp_rates, tip_labels = make_ape(tree)
    ntip = (tree.n_nodes + 1)/2

    @rput edges
    @rput branch_lengths
    @rput temp_rates
    @rput ntip
    @rput file
    @rput tip_labels

    reval("""
        temp_tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
    """)
end

function save_ClaDS_in_R(co::CladsOutput, path::String ; maxit = Inf, burn = 0.25)

    chains = co.chains
    rtt_chains = co.rtt_chains
    sig_map = co.σ_map
    al_map = co.α_map
    eps_map = co.ε_map
    lambda0_map = co.λ0_map
    lambdai_map = co.λi_map
    lambdatip_map = co.λtip_map
    DTT_mean = co.DTT_mean
    RTT_map = co.RTT_map
    time_points = co.time_points
    enhanced_trees = co.enhanced_trees
    gelm_id, gelm_stat  = co.gelm

    R_tree(co.tree)
    reval("""tree = temp_tree""")

    reval("""enhanced_trees = list()""")
    for i in 1:length(co.enhanced_trees)
        @rput i
        R_tree(co.enhanced_trees[i])
        reval("""enhanced_trees[[i]] = list(tree = temp_tree, rates = temp_rates)""")
    end

    @rput chains
    @rput rtt_chains
    @rput sig_map
    @rput al_map
    @rput eps_map
    @rput lambda0_map
    @rput lambdai_map
    @rput lambdatip_map
    @rput DTT_mean
    @rput RTT_map
    @rput time_points
    @rput gelm_id
    @rput gelm_stat

    @rput path
    reval("""
        CladsOutput = list(
            tree = tree,
            rtt_chains = rtt_chains,
            sig_map = sig_map,
            al_map = al_map,
            eps_map = eps_map,
            lambda0_map = lambda0_map,
            lambdai_map = lambdai_map,
            lambdatip_map = lambdatip_map,
            DTT_mean = DTT_mean,
            RTT_map = RTT_map,
            time_points = time_points,
            enhanced_trees = enhanced_trees,
            gelm = c(gelm_id, gelm_stat)
        )

        save(CladsOutput, file = path)
    """)
end
