function infer_ClaDS(tree::Tree, n_reccord::Int64; ini_par = [], initialize_rates = 0, goal_gelman = 1.05,
    thin = 1, burn = 1/4, f = 1., plot_tree = 0, print_state = 0, max_node_number = 100, plot_chain = false,quad = 1,
    max_try = 10_000, it_edge_tree = 30, print_all = false, it_rates = 3, former_run = CladsOutput(), plot_burn = NaN, ltt_steps = 50,
    save_as_R = false, Rfile = "coda.Rdata", max_it_number = Inf, enhance_method = "MHrr", end_it = Inf, n_chains = 3, n_trees = 10)

    ntips = Int64((tree.n_nodes + 1)/2)
    if isnan(plot_burn)
        plot_burn = burn
    end

    if length(former_run.chains) == 0
        sampler = initialize_ClaDS2(tree , ini_par = ini_par, initialize_rates = initialize_rates, ltt_steps = ltt_steps, enhance_method = enhance_method, n_chains = n_chains, n_trees = n_trees)
    else
        sampler = (former_run.chains,
            former_run.current_state[1],
            former_run.current_state[2],
            former_run.current_state[3],
            former_run.current_state[4],
            former_run.time_points,
            former_run.current_state[5],
            former_run.rtt_chains,
            former_run.enhanced_trees)

            former_run = CladsOutput()
    end

    if length(f) == 1
        fs = fill(f,tree.n_nodes-1)
    elseif length(f) == ntips
        fs = sample_fractions(tree,f)[2:end]
    else
        fs = f
    end

    gelman = 10.
    MAPS=[]
    nit = 0

    while (gelman > goal_gelman) & (nit < end_it)

        sampler = add_iter_ClaDS2(sampler, n_reccord, thin = thin, fs = fs, plot_tree = plot_tree, print_state = print_state,
            max_node_number = max_node_number, max_try = max_try, it_edge_tree = it_edge_tree,
            print_all = print_all, it_rates = it_rates, enhance_method = enhance_method, quad = quad, n_trees = n_trees)

        nit += n_reccord



        chains = sampler[1]
        npar = sampler[4][1].n_nodes + 3
        @rput chains
        @rput npar

        if goal_gelman == 0.
            id_gelman = 0
            gelman = 10
        else
            gelman = gelman_est(sampler[1], npar, burn = burn)
            id_gelman = gelman[1]
            gelman = gelman[2]
        end

        if plot_chain
            plot_coda(sampler, burn = plot_burn, id_par=[1,2,3,max(4,id_gelman)])
        else
            chains = sampler[1]
            @rput chains
            @rput thin
            @rput burn
            id_par = 1
            @rput id_par
            reval("""
                require(coda)
                n_row = length(chains[[1]][[1]])
                ini = floor(burn * n_row)
                if(ini == 0) ini = 1
                    it = seq(ini, n_row, 1)
            """)
        end

        @rput ntips
        @rput n_chains

        reval("""
            npar2 = length(chains[[1]])
            n_row = length(chains[[1]][[1]])
            ini = floor(burn * n_row)
            if(ini == 0) ini = 1
            it = c(unique(ceiling(seq(ini, n_row-1,length.out = 100))),n_row)

            map_chains = mcmc.list(lapply(1:n_chains, function(i){
                mcmc(sapply(1:npar2, function(j){
                    if(j<=3 || j>(npar+ntips)) {
                        chains[[i]][[j]][it]
                    }else{
                        log(chains[[i]][[j]][it])
                    }
            }))}))

            unlist_chains = (sapply(1:npar2, function(k){
                x = c()
                for (n in 1:n_chains){
                    x = c(x,map_chains[[n]][,k] )
                    }
                return(x)
                }))

            MAPS=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
                return(D[[1]][which.max(D[[2]])])})
            #MAPS[2] = exp(MAPS[2])
        """)
        @rget MAPS
        if save_as_R
            chains = sampler[1]
            sampler_to_Rdata(tree, (sampler, MAPS), Rfile , sample_fraction = f, max_it_number = max_it_number)
        end
        println("iteration $(length(sampler[1][1][1])) ; gelman = $gelman for param $id_gelman")
        println("   σ = $(MAPS[1]), α = $(MAPS[2]), ε = $(MAPS[3]), λ0 = $(MAPS[4])")
        println("")
    end

    gelm = (0,0)
    npar = sampler[4][1].n_nodes + 3

    if goal_gelman > 0.
        gelm = gelman_est(sampler[1], npar, burn = burn)
    end

    @rput ntips
    @rput n_chains
    @rput npar

    reval("""
        npar2 = length(chains[[1]])
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, 1)

        map_chains = mcmc.list(lapply(1:n_chains, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j==1 || j==3 || j>(npar+ntips)) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))

        unlist_chains = (sapply(1:npar2, function(k){
            x = c()
            for (n in 1:n_chains){
                x = c(x,map_chains[[n]][,k] )
                }
            return(x)
            }))

        MAPS=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
            return(D[[1]][which.max(D[[2]])])})
        MAPS[2] = exp(MAPS[2])

        mean_ltt = sapply((npar + ntips + 3):npar2, function(i){mean(unlist_chains[,i])})
        unlist_chains = list()
        map_chains = list()
    """)
    @rget MAPS
    @rget mean_ltt

    mean_rates_chains = sampler[8]
    @rput mean_rates_chains
    reval("""
        nmr = length(mean_rates_chains[[1]][[1]])
        map_chains = mcmc.list(lapply(1:n_chains, function(i){
            mcmc(sapply(1:nmr, function(j){
                res = c()
                for(t in it){
                    res = c(res,log(mean_rates_chains[[i]][[t]][j]))
                    }
                return(res)

                }))}))

        unlist_chains = (sapply(1:nmr, function(k){
            x = c()
            for (n in 1:n_chains){
                x = c(x,map_chains[[n]][,k] )
                }
            return(x)
            }))
        RTT_map = sapply(1:nmr, function(i){D=density(unlist_chains[,i]);
            return(D[[1]][which.max(D[[2]])])})
        RTT_map = exp(RTT_map)
    """)
    @rget RTT_map

    output = CladsOutput(
        tree,
        sampler[1],                         #chains
        mean_rates_chains,                  #rtt_chains
        MAPS[1], MAPS[2], MAPS[3], MAPS[4], #σ_map,α_map,ε_map,λ0_map
        MAPS[5:npar],                       #λi_map
        MAPS[(npar+1):(npar+ntips)],      #λtip_map
        mean_ltt,                           #DTT_mean
        RTT_map,                            #RTT_map
        sampler[6],                         #time_points
        sampler[9],                         #enhanced_trees
        gelm,                               #gelman statistics
        (sampler[2], sampler[3], sampler[4], sampler[5], sampler[7])
    #=
    current_state=#
    )
    return output


end
