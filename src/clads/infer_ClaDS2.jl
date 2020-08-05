function infer_ClaDS(tree::Tree, n_reccord::Int64; ini_par = [], initialize_rates = 0, goal_gelman = 1.05,
    thin = 10, burn = 1/4, f = 1., plot_tree = 0, print_state = 0, max_node_number = 100, plot_chain = false,quad = 1,
    max_try = 10_000, it_edge_tree = 10, print_all = false, it_rates = 1, sampler = [], plot_burn = NaN, ltt_steps = 50,
    save_as_R = false, Rfile = "coda.Rdata", max_it_number = Inf, enhance_method = "MHrr", end_it = Inf, n_chains = 3)

    ntips = Int64((tree.n_nodes + 1)/2)
    if isnan(plot_burn)
        plot_burn = burn
    end

    if length(sampler) == 0
        sampler = initialize_ClaDS2(tree , ini_par = ini_par, initialize_rates = initialize_rates, ltt_steps = ltt_steps, enhance_method = enhance_method, n_chains = n_chains)
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
            print_all = print_all, it_rates = it_rates, enhance_method = enhance_method, quad = quad)

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
    """)
    @rget MAPS
    return sampler, MAPS, gelm
end
