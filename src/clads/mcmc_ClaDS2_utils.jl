struct CladsOutput
    tree::Tree
    chains::Array{Array{Array{Float64,1},1},1}
    rtt_chains::Array{Array{Array{Float64,1},1},1}
    σ_map::Float64
    α_map::Float64
    ε_map::Float64
    λ0_map::Float64
    λi_map::Array{Float64,1}
    λtip_map::Array{Float64,1}
    DTT_mean::Array{Float64,1}
    RTT_map::Array{Float64,1}
    time_points::Array{Float64,1}
    enhanced_trees::Array{Tree,1}
    gelm::Tuple{Int64,Float64}
    current_state # all the rest that is mising to continue the run
end

function CladsOutput()
    return CladsOutput(
        Tree(),
        Array{Array{Array{Float64,1},1},1}(undef,0),
        Array{Array{Array{Float64,1},1},1}(undef,0),
        0.,0.,0.,0.,[0.],[0.],[0.],[0.],
        [0.],[Tree()],(0,10.),0.
    )
end

function plot_coda(sampler ; burn = 0, thin = 1, id_par = [1:4...], smooth = false)
    chains = sampler[1]
    @rput chains
    @rput thin
    @rput burn
    @rput id_par
    @rput smooth

    reval("""
        require(coda)
        N = length(chains)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
            it = seq(ini, n_row, thin)

        plot_chains = mcmc.list(lapply(1:N, function(i){
        mcmc(sapply((id_par), function(j){
        if(j<4) {
            chains[[i]][[j]][it]
        }else{
            log(chains[[i]][[j]][it])
        }
        }))}))

        plot(plot_chains, smooth = smooth)
    """)
end

function add_to_chain!(chain, param)
    for i in 1:length(param)
        push!(chain[i], param[i])
    end
end

function initialize_ClaDS2(tree::Tree ; ini_par = [], initialize_rates = 0, ltt_steps = 10, enhance_method = "reject", n_chains = 3, n_trees = 10) where {N}

    new_tree = Tree(deepcopy(tree.offsprings), 0., deepcopy(tree.attributes))

    root_depth = maximum(node_depths(new_tree)) * 1.01
    live_nd = node_depths_base(new_tree)

    ltt_times = [0:ltt_steps...] * root_depth/ltt_steps
    ltt_extant = LTT(tree,ltt_times)
    if length(ini_par) == 0
        ini_par = [[[10. ^(1-i),1,0.5];fill(0.000001, tree.n_nodes); fill(0.000001, Int64((tree.n_nodes+1)/2));[0, n_extant_tips(tree)] ; ltt_extant[2]] for i in 1:n_chains]
    end

    edge_trees_s = [init_edge_tree_rates2(update_rates(new_tree, ini_par[i][4:(tree.n_nodes + 3)]), ini_par[i][4:(end-2)]) for i in 1:n_chains]


    for i in 1:n_chains
        new_tree = Tree(deepcopy(tree.offsprings), 0., deepcopy(tree.attributes))
        update_rates!(new_tree, ini_par[i][4:(tree.n_nodes + 3)])

        rates = ini_par[i][4:(tree.n_nodes + 3)]
        ε = ini_par[i][3]
        σ = ini_par[i][1]
        α = ini_par[i][2]
        edge_trees = edge_trees_s[i]
        for j in 1:initialize_rates
            rates, ε, σ, α = update_edges_ini!(new_tree, edge_trees, σ, α, ε, rates, it_rates = 1, with_ε = false)
            update_rates!(new_tree, rates)
        end
        ini_par[i][4:(tree.n_nodes + 3)] = rates
        ini_par[i][3] = ε
        edge_trees_s[i] = edge_trees
    end

    extant_branch_lengths = extract_branch_lengths(tree)
    trees = [deepcopy(update_rates(new_tree, ini_par[i][4:end])) for i in 1:n_chains]
    chains = [[[deepcopy(ini_par[i][j])] for j in 1:length(ini_par[1])] for i in 1:n_chains]
    param = deepcopy(ini_par)
    mean_rates_chains = [[time_rates(tree,ini_par[i][4:(tree.n_nodes+3)],ltt_times)] for i in 1:n_chains]
    println(" ")
    return chains, param, edge_trees_s, trees, extant_branch_lengths, ltt_times, live_nd, mean_rates_chains, Array{Tree,1}(undef,0)
end

function add_iter_ClaDS2(sampler, n_reccord::Int64; thin = 1, fs = 1., plot_tree = 0, print_state = 0, quad = 1,
    max_node_number = 1_000, max_try = 100_000, it_edge_tree = 1, print_all = false, it_rates = 1, enhance_method = "reject", n_trees = 10)

    chain_s, param_s, edge_trees_s, tree_s, extant_branch_lengths, ltt_times, live_nd, mean_rates_chains, enhanced_trees = sampler

    n_chains = length(chain_s)
    tips_id = tips(tree_s[1])[2:end]
    lefts = n_left(tree_s[1])
    ntips = Int64((tree_s[1].n_nodes+1)/2)

    if plot_tree>0
        reval("par(mfrow=c(3,4), mar=c(5,2,2,2))")
    end

    n_par = tree_s[1].n_nodes + 3

    for k in 1:n_chains
        chain = chain_s[k]
        tree = tree_s[k]

        σ = param_s[k][1]
        α = param_s[k][2]
        ε = param_s[k][3]
        rates = param_s[k][4:(tree.n_nodes + 3)]

        edge_trees = edge_trees_s[k]

        relative_rates = Array{Float64,1}(undef,0)

        for i in 1:n_reccord
            for j in 1:thin
                for l in 1:it_edge_tree
                    update_edge_trees!(edge_trees, tree, σ, α, ε, fs, rates, extant_branch_lengths, enhance_method= enhance_method,
                        max_node_number = max_node_number, max_try = max_try, keep_if_any = true,
                        do_tips =  true)


                    relative_rates = extract_relative_rates(tree, edge_trees, rates)
                    draw_λ0_slicing!(rates, edge_trees, σ, α, ε, lefts)
                    σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
                    α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)
                    ε = draw_ε_crown(tree, edge_trees, lefts)
                    draw_λ0!(tree::Tree, ε::Float64, rates::Array{Float64,1},
                        edge_trees::Array{EdgeTreeRates2,1})
                    tree.attributes[1] = rates[1]
                end

                if it_rates > 0
                    ε, σ, α = update_edges_ETR2!(tree, edge_trees, σ, α, ε, rates,lefts,it_rates = it_rates )
                end

                if plot_tree>0
                    if i%plot_tree == 0 && j == thin
                        enhanced_tree = graft_edge_trees(tree, edge_trees)
                        plot_ClaDS(enhanced_tree)#, extant_edges[2:end], ln = false)
                    end
                end
                if print_state>0
                    if i%print_state == 0 && (j == thin || print_all)
                        relative_rates = extract_relative_rates(tree, edge_trees)
                        println("$i : sd $(sqrt(var(relative_rates))), mean $(mean(relative_rates)), σ $σ, α $α, ε $ε, λ0 $(rates[1]), n_ext $(extract_nextinct(tree,edge_trees)) $(n_tip(tree,edge_trees)) $(n_extant_tips(tree,edge_trees))")
                    end
                end

                if length(enhanced_trees) < n_trees
                    enhanced_tree = graft_edge_trees(tree, edge_trees)
                    push!(enhanced_trees, enhanced_tree)
                else
                    u = rand()
                    if u < 0.1*n_trees/i
                        enhanced_tree = graft_edge_trees(tree, edge_trees)
                        itree = sample(1:n_trees)
                        enhanced_trees[itree] = enhanced_tree
                    end
                end

                #println(" 3 nnn")

            end
            #println("$σ ")

            tip_rates = extract_tip_rates(tree, edge_trees, tips_id, rates)
            param_s[k][1] = σ
            param_s[k][2] = α
            param_s[k][3] = ε
            param_s[k][4:n_par] = rates
            param_s[k][(n_par+1):(n_par+ntips)] = tip_rates
            param_s[k][n_par + ntips + 1] = extract_nextinct(tree,edge_trees)
            param_s[k][n_par + ntips + 2] = n_tip(tree,edge_trees) - extract_nextinct(tree,edge_trees)
            #enhanced_tree = graft_edge_trees(tree, edge_trees)
            #ltt = LTT(enhanced_tree, ltt_times)[2]
            ltt = LTT(tree, edge_trees, ltt_times)[2]
            param_s[k][(n_par + ntips + 3):end] = ltt
            edge_trees_s[k] = edge_trees
            add_to_chain!(chain_s[k], param_s[k])
            #println(time_rates(tree,edge_trees,ltt_times))
            #println(mean_rates_chains)
            for imr in 1
                #add_to_chain!(mean_rates_chains[k], time_rates(tree,edge_trees,ltt_times))
                push!(mean_rates_chains[k], time_rates(tree,edge_trees,ltt_times))
            end
            tree_s[k] = tree
        end
    end

    println(" ")
    return chain_s, param_s, edge_trees_s, tree_s, extant_branch_lengths, ltt_times, live_nd, mean_rates_chains, enhanced_trees
end


function gelman_est(mcmc::Array{Array{Array{Float64,1},1},1}, Npar::Int64 ; thin = 1, burn = 0)
    nChains = length(mcmc)
    Nit = length(mcmc[1][1])
    ini = Int64(floor(1 + Nit * burn))

    it = [ini]
    i = ini + thin
    while i<= Nit
        push!(it,i)
        i += thin
    end
    Nit = length(it)

    # compute means
    M = [zeros(Npar) for i in 1:nChains]

    for n in 1:nChains
        for i in 1:Npar
            if ((i == 1 )|| (i == 3))
                for j in it
                    M[n][i] += sum(mcmc[n][i][j])
                end
            else
                for j in it
                    M[n][i] += sum(log(mcmc[n][i][j]))
                end
            end
        end
    end

    M ./= Nit

    m = zeros(Npar)
    for n in 1:nChains
        m .+= M[n]
    end

    m ./= nChains

    # compute vars
    Σ = [zeros(Npar) for i in 1:nChains]

    for n in 1:nChains
        for i in 1:Npar
            if ((i == 1 )|| (i == 3))
                for j in it
                    Σ[n][i] += (mcmc[n][i][j] - M[n][i])^2
                end
            else
                for j in it
                    Σ[n][i] += (log(mcmc[n][i][j]) - M[n][i])^2
                end
            end

        end
    end

    Σ ./= (Nit-1)

    # between chains variance
    B = zeros(Npar)

    for n in 1:nChains
        for i in 1:Npar
            B[i] += (M[n][i] - m[i])^2
        end
    end

    B .* (nChains + 1)/(nChains * (nChains - 1))


    # within chains variance

    W = zeros(Npar)
    for n in 1:nChains
        W .+= Σ[n]
    end
    W ./= nChains

    Wd = W .* (Nit-1)/Nit


    # pooled variance
    V = Wd .+ B

    # PSRF
    PSRF = V ./ W

    #return PSRF
    id = 0
    maxi = 0.

    for i in 1:Npar
        if PSRF[i] > maxi
            id = i
            maxi = PSRF[i]
        end
    end

    return id, sqrt(maxi)
end
