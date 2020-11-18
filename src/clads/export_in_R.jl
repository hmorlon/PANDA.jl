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
        class(temp_tree) = "phylo"
    """)
end

"""
    save_ClaDS_in_R(co::CladsOutput, path::String ; maxit = Inf, ...)

Save the output of a ClaDS run as a Rdata file.

# Arguments
- `co::CladsOutput`: the result of a run of `infer_ClaDS`.
- `path::String`: the path the file should be saved to.

# .Rdata file
The .Rdata file contains an object called `CladsOutput` that is a list with the following fields :
- `tree`: the phylogeny, saved as a `phylo` object
- `chains`: the MCMC chains
- `rtt_chains`: the MCMCs with the rate through time information
- `sig_map`: the MAP estimate for the σ parameter
- `al_map`: the MAP estimate for the α parameter
- `eps_map`: the MAP estimate for the ε parameter
- `lambda0_map`: the MAP estimate for the λ0 parameter
- `lambdai_map`: the MAP estimate for the branch-specific speciation rates
- `lambdatip_map`: the MAP estimate for the tip speciation rates
- `DTT_mean`: the point estimate for the diversity through time
- `RTT_map` : the rate through time estimates
- `time_points` : the times at which `DTT_mean` and `RTT_map` are computed
- `enhanced_trees` : a sample of the complete tree distribution, as a list of `phylo` ojects
- `gelm` : the gelman parameter
"""
function save_ClaDS_in_R(co::CladsOutput, path::String)

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
            chains = chains,
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
