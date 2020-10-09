#=
R ploting function
=#

function create_plot_ClaDS_ape()
    reval("""
    require(fields)
    require(ape)
    require(RColorBrewer)

    plot_ClaDS=function(phylo,rate1,rate2=NULL,same.scale=T,main=NULL,lwd=3,log=T, show_labels=F, minr = Inf, maxr = 0, show_legend = T,...){
        Colors = colorRampPalette(rev(c('darkred',brewer.pal(n = 8, name = "Spectral"),'darkblue')))(100)

     if(nrow(phylo[[1]]) <=1){
        plot(1000,xlim = c(0,1), ylim = c(0,2), axes = F, xlab = "", ylab = "")
        lines(c(0,1), c(1,1), lwd=lwd, col="steelblue2")
      }else{
        if(is.null(rate2)){
          if(log) {
              rate1=log(rate1)
              minr = log(minr)
              maxr = log(maxr)
              minr = min(c(minr, log(0.95)+max(rate1),min(rate1)))
              maxr = max(c(maxr, log(1.05)+min(rate1), max(rate1)))
              }else{
                minr = min(c(minr, 0.99*max(rate1),min(rate1)))
                maxr = max(c(maxr, 1.01*min(rate1), max(rate1)))
              }

            col = round( (rate1 - minr) / (maxr - minr)*99   )+1
            plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,main=main,edge.width =lwd,...)

            if(log){
              m10=floor(minr/log(10))
              M10=ceiling(maxr/log(10))
              if((M10-m10)<4){
                ticks=c(1,2,5)
              }else{
                ticks=1
              }

              ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
              lt=length(ticks[ticks>exp(minr) & ticks<exp(maxr)])
              if(lt<4){
                  ticks=c(round(exp(minr),digit=2),max(0,-1*m10+(lt<2)),ticks,round(exp(maxr),digit = 2))
                  }
              if(show_legend) image.plot(z = c(minr,maxr),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
            }else{
              if(show_legend) image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
            }

        }else{
          if(log){
            rate1=log(rate1)
            rate2=log(rate2)
          }
          if(same.scale){
            min=min(min(rate1),min(rate2))
            max=max(max(rate1),max(rate2))
            par(mfrow=c(1,2))
            col = round(( (rate1 - min) / (max-min))*99   )+1
            plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
            col = round(( (rate2 - min) / (max-min))*99   )+1
            plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
            par(mfrow=c(1,1))
            if(log){
              m10=floor(min/log(10))
              M10=ceiling(max/log(10))
              if((M10-m10)<4){
                ticks=c(1,2,5)
              }else{
                ticks=1
              }
              ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
              lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
              if(lt<4){

                  ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))
                  }
              # ticks=seq(min,max,length.out = 5)
              if(show_legend) image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
            }else{
              if(show_legend) image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T)
            }
          }else{
            par(mfrow=c(1,2))
            if(isTRUE(all.equal(rep(rate1[1],length(rate1)),rate1))){
              col=rep(1,length(rate1))
              plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
              if(log){

                if(show_legend) image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T)
              }else{
                if(show_legend) image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
              }
            }else{
              col = round(( (rate1 - min(rate1)) / (max(rate1)-min(rate1)))*99   )+1
              plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
              if(log){
                min=min(rate1)
                max=max(rate1)
                m10=floor(min/log(10))
                M10=ceiling(max/log(10))
                if((M10-m10)<4){
                  ticks=c(1,2,5)
                }else{
                  ticks=1
                }
                ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
                lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
                if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
                if(show_legend) image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
              }else{
                if(show_legend) image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
              }
            }
            if(isTRUE(all.equal(rep(rate2[1],length(rate2)),rate2))){
              col=rep(1,length(rate2))
              plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
              if(log){
                if(show_legend) image.plot(z = c(exp(rate2[1]),2*exp(rate2[1])),col = Colors, horizontal=T,legend.only = T)
              }else{
                if(show_legend) image.plot(z = c(rate2[1],2*rate2[1]),col = Colors, horizontal=T,legend.only = T)
              }
            }else{
              col = round(( (rate2 - min(rate2)) / (max(rate2)-min(rate2)))*99   )+1
              plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
              if(log){
                minr=min(rate2)
                maxr=max(rate2)
                m10=floor(minr/log(10))
                M10=ceiling(maxr/log(10))
                if((M10-m10)<4){
                  ticks=c(1,2,5)
                }else{
                  ticks=1
                }
                ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
                lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
                if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
                if(show_legend) image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
              }else{
                if(show_legend) image.plot(z = as.matrix(rate2),col = Colors, horizontal=T,legend.only = T)
              }
            }
          }
          par(mfrow=c(1,1))
        }
      }
      if(!show_legend) return(list(Colors = Colors, col = col, minr = minr, maxr = maxr))
    }
    """)
end


#=
Plot a tree with branch rates from Julia
3 methods
=#

# specifying only the tree, .attribute[id] used as rates
function plot_ClaDS(tree::Tree ; id = 1, ln=true, show_labels = false, options="")
    plot_tree = Tree(tree.offsprings, 0., tree.attributes, tree.n_nodes)

    opt = options
    if length(options) > 0
        if options[1] != ","[1]
            opt = ", $options"
        end
    end

    edges, branch_lengths, rates, tip_labels = Tree2ape(plot_tree, id = 1)
    ntip = (tree.n_nodes + 1)/2
    @rput edges
    @rput branch_lengths
    @rput rates
    @rput ntip
    @rput ln
    @rput show_labels
    @rput tip_labels

    create_plot_ClaDS_ape()
    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
        colors = plot_ClaDS(tree, rates, log=ln, show_labels = show_labels $opt);
    """)
end

# specifying the tree and a vector of rates
function plot_ClaDS(tree::Tree, rates ; id = 1, ln=true, round = false, options = "", show_labels=false)
    plot_tree = Tree(tree.offsprings, 0., tree.attributes, tree.n_nodes)

    opt = options
    if length(options) > 0
        if options[1] != ","[1]
            opt = ", $options"
        end
    end

    edges, branch_lengths, tree_rates, tip_labels = Tree2ape(plot_tree, id = 1)
    ntip = (tree.n_nodes + 1)/2
    @rput edges
    @rput branch_lengths
    @rput rates
    @rput ntip
    @rput ln
    @rput show_labels
    @rput tip_labels

    create_plot_ClaDS_ape()
    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
        plot_ClaDS(tree, rates, log=ln, show_labels = show_labels $opt)
    """)
end

# specifying the tree and two vectors of rates, both kind of rates are plotted on the same color scale
function plot_ClaDS(tree::Tree, rates1, rates2 ; id = 1, ln=true)
    plot_tree = Tree(tree.offsprings, 0., tree.attributes, tree.n_nodes)

    reval("""
    source("/Users/maliet/ownCloud/My_folder/ClaDS_Julia_list_of_lists/plot_ClaDS_2trees.R")
    """)
    edges, branch_lengths, new_rates = Tree2ape(plot_tree, id = 1)
    ntip = (tree.n_nodes + 1)/2

    @rput edges
    @rput branch_lengths
    @rput rates1
    @rput rates2
    @rput ntip
    @rput ln

    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = 1:ntip)
        class(tree) = "phylo"
        leg = plot_ClaDS_noLeg(tree, rates1, rates2, log=ln $opt)
        leg = plot_ClaDS_noLeg(tree, rates2, rates1, log=ln $opt)
        image.plot(z = leg[[1]],col = leg[[2]], horizontal=F,legend.only = T,axis.args=leg[[5]], legend.mar=4.5)
    """)
end


# diversity through time
function plot_DTT(co::CladsOutput; nyl = 1000, n_ltt = 100, alpha_col = 0.05, burn = 0.25,
	thin = 1, axes_cex = 1., lwd = 1., lab_cex = 1.) where {T <: Number}

	extant_tree = co.tree
	chains = co.chains
    npar = extant_tree.n_nodes + 3 + (extant_tree.n_nodes+1)/2
    ltt_times = co.time_points
    maps = co.DTT_mean

    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput alpha_col
    @rput burn
    @rput thin
    @rput live_ltt
	@rput axes_cex
	@rput lab_cex
    @rput nyl
    reval("""
        require(coda)
        id_ltt = (npar+3):length(chains[[1]])
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(id_ltt, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:length(id_ltt), function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))
		y_max = max(unlist_chains)
        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:length(id_ltt), function(i){if(T){log(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = log(c(2,y_max)), xlab = "time",
        ylab = "nb lineages", cex.lab = lab_cex)
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        if(nyl < length(y_lab)){
            y_lab = y_lab[unique(floor(seq(1,length(y_lab), length.out = nyl)))]
            }
        lines(ltt_times, log(live_ltt), col = "black", lwd = 6, lty = 1)

        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))

        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){log(chains[[i]][[k]][j])})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
		quant = sapply(1:ncol(Ys), function(i){quantile(log(unlist_chains[,i]), probs = c(0.05,0.95))})

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[1:length(id_ltt)], col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[1:length(id_ltt)], col = "darkseagreen4", lwd = 3, lty = 6)

        axis(1, cex.axis = axes_cex, lwd = lwd, lwd.ticks = lwd)
        axis(2, at = log(y_lab), lab = y_lab, cex.axis = axes_cex, lwd = lwd, lwd.ticks = lwd)

    """)
end


function plot_RTT(co::CladsOutput ;nplot = 50, miny = -1, maxy = 1, alpha_col = 0.05, burn = 0., axes_cex = 1., lwd = 1., lab_cex = 1., lay = 5)

	times = co.time_points[2:end]
	RTT_map = co.RTT_map
    N = length(co.rtt_chains[1])
    t = times#[1:length(tr)]
	ini_plot = Int64(2+ floor(burn * N))
    plot_id = Array{Int64,1}(floor.(range(ini_plot, length=nplot, stop=N)))
    mean_mr = zeros(length(times))
    a = Array{Float64,1}(undef,0)
    map_mr = [deepcopy(a) for i in 1:length(times)]
	range_RTT = [minimum(minimum(minimum(co.rtt_chains)[2:end]));maximum(maximum(maximum(co.rtt_chains)[2:end]))]
    @rput t
    @rput miny
    @rput maxy
    @rput alpha_col
	@rput lab_cex
	@rput axes_cex
	@rput lay
	@rput RTT_map
	@rput range_RTT

    reval("""
		ylim = log(range_RTT)#+c(miny,maxy)
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = ylim,
        xlab = "time", ylab = "mean rate", cex.lab = lab_cex)
    """)
    n = 0
    colors = ["deepskyblue2" ; "deepskyblue2" ; "deepskyblue2"]

    #colors = ["deepskyblue2" ; "orange2" ; "orangered3"]
    for i in 1:3
        color = colors[i]
        @rput color
        for k in ini_plot:N
            y = co.rtt_chains[i][k]
            for j in 1:length(y)
                push!(map_mr[j],log(y[j]))
            end
			if k ∈ plot_id
            	@rput y
            	reval("""lines(t, log(y), col = alpha(color, alpha = alpha_col), lwd = 2)""")
			end
            mean_mr .+= y
            n += 1
        end
    end

    mean_mr ./= n
    @rput mean_mr
    @rput map_mr
    reval("""
        maps = log(RTT_map) #sapply(map_mr, function(x){D=density(x); return(D[[1]][which.max(D[[2]])])})
		#quant = sapply(map_mr, function(x){quantile(x, probs = c(0.05,0.95))})
		quant = sapply(map_mr, function(x){quantile(x, probs = c(0.05,0.95))})

    """)
    reval("""

        lines(t, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(t, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(t, maps, col = "darkseagreen1", lwd = 5, lty = 1)
        lines(t, maps, col = "darkseagreen4", lwd = 3, lty = 6)
        axis(1, cex.axis = axes_cex, lwd = lwd, lwd.ticks = lwd)
		labels_a2 = unique(signif(exp(seq(ylim[1], ylim[2], length.out = lay)), digits = 1))
		axis(2, at = log(labels_a2), labels = labels_a2, cex.axis = axes_cex, lwd = lwd, lwd.ticks = lwd)
    """)

end

function plot_density(co::CladsOutput, id_par::Int64; burn = 0.25)
	n_par = length(co.chains[1])
	n_tips = n_extant_tips(co.tree)
	n_edges = (n_tips - 1)*2
	nit = length(co.chains[1][1])
	init = Int64(ceil(burn * nit))
	color = "orange"

	if id_par >= 5 + n_edges + n_tips
		id_par += 2
	end

	if id_par <= n_par
		if id_par == 1				#σ
			name = "sigma"
			map = co.σ_map
		elseif id_par == 2			#α
			name = "alpha"
			map = co.α_map
		elseif id_par == 3			#ε
			name = "epsilon"
			map = co.ε_map
		elseif id_par == 4			#λ_0
			name = "lambda_0"
			map = co.λ0_map
		elseif id_par < 5 + n_edges #λ_i
			i = id_par - 4
			name = "lambda_$i"
			map = co.λi_map[i]
		elseif id_par < 5 + n_edges + n_tips						#λ_tip_i
			i = id_par - 4 - n_edges
			tip_name = tip_labels(co.tree)[i]
			name = "lambda_tip $i"
			if length(tip_name)>0
				name = "lambda_tip $tip_name"
			end
			map = co.λtip_map[i]
		else
			color = "#3498db"
			i = id_par - (6 + n_edges + n_tips)
			name = "number of lineages (time $(co.time_points[i]))"
			map = co.DTT_mean[i]
		end

		chain = [co.chains[i][id_par][j] for i in 1:length(co.chains) for j in init:length(co.chains[1][1])]
	else

		i = id_par - n_par
		name = "mean rate (time $(co.time_points[i]))"
		map = co.RTT_map[i]

		chain = [co.rtt_chains[i][j][id_par - n_par] for i in 1:length(co.rtt_chains) for j in init:length(co.rtt_chains[1])]
	end
	@rput chain
	@rput name
	@rput map
	@rput id_par
	@rput n_par
	@rput color
	reval("""
		if (id_par < 5 ){
			d = density(chain)
			plot(d\$x, d\$y, type = 'l', lwd = 3, xlab = name, ylab = "density")
		}else{
			d = density(log(chain))
			plot(exp(d\$x), d\$y, type = 'l', lwd = 3, xlab = name, ylab = "density", log = 'x')
		}
		abline(v = map, lwd = 3, col = color)

		color_q = adjustcolor("coral3", alpha.f = 0.1)
		quant = quantile(chain, c(0.025,0.975))
		polygon(quant[c(1,2,2,1)], c(-1,-1,100000,100000), col = color_q, border = NA)

	""")
end

function plot_chain(co::CladsOutput, id_par::Int64; burn = 0.25)
	n_par = length(co.chains[1])
	n_tips = n_extant_tips(co.tree)
	n_edges = (n_tips - 1)*2
	nit = length(co.chains[1][1])
	init = Int64(ceil(burn * nit))
	if id_par >= 5 + n_edges + n_tips
		id_par += 2
	end

	if id_par <= n_par
		if id_par == 1				#σ
			name = "sigma"
			map = co.σ_map
		elseif id_par == 2			#α
			name = "alpha"
			map = co.α_map
		elseif id_par == 3			#ε
			name = "epsilon"
			map = co.ε_map
		elseif id_par == 4			#λ_0
			name = "lambda_0"
			map = co.λ0_map
		elseif id_par < 5 + n_edges #λ_i
			i = id_par - 4
			name = "lambda_$i"
			map = co.λi_map[i]
		elseif id_par < 5 + n_edges + n_tips						#λ_tip_i
			i = id_par - 4 - n_edges
			tip_name = tip_labels(co.tree)[i]
			name = "lambda_tip $i"
			if length(tip_name)>0
				name = "lambda_tip $tip_name"
			end
			map = co.λtip_map[i]
		else
			i = id_par - (6 + n_edges + n_tips)
			name = "number of lineages (time $(co.time_points[i]))"
			map = co.DTT_mean[i]
		end

		chain = [co.chains[i][id_par][init:length(co.chains[1][1])] for i in 1:length(co.chains)]
	else

		i = id_par - n_par
		name = "mean rate (time $(co.time_points[i]))"
		map = co.RTT_map[i]

		chain = Array{Array{Float64,1},1}(undef,0)
		for i in 1:length(co.rtt_chains)
			push!(chain, [co.rtt_chains[i][j][id_par - n_par] for j in init:length(co.rtt_chains[1])])
		end
	end
	@rput chain
	@rput name
	@rput map
	@rput id_par
	@rput n_par

	reval("""
		require("coda")
		traceplot(mcmc.list(lapply(chain, mcmc)), ylab = name)

	""")
end

function convert_id(co::CladsOutput, id::String)
	n_par = length(co.chains[1])
	n_tips = n_extant_tips(co.tree)
	n_edges = (n_tips - 1)*2

	if occursin("σ",id) | occursin("sig",id)
		par_id = 1
	elseif occursin("α",id) | occursin("al",id)
		par_id = 2
	elseif occursin("ε",id) | occursin("ϵ",id) | occursin("eps",id) | occursin("turn",id)
		par_id = 3
	elseif occursin("λ0",id) | occursin("lambda0",id) | occursin("lambda_0",id) | occursin("lamb_0",id)
		par_id = 4
	elseif occursin(r"λ*tip",id) | occursin(r"lamb*tip",id)
		number = 0
		tip_id=split(id,"tip_")[end]
		re = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";
		if occursin(re, tip_id)
			number=parse(Int, tip_id)
		else
			dist = Inf
			tip_labs = tip_labels(co.tree)
			for itl in 1:length(tip_labs)
				d = Levenshtein()(tip_id, tip_labs[itl])
				if d < dist
					dist = d
					number = itl
					if d == 0
						break
					end
				end
			end
		end
		par_id = 4 + n_edges + number
	elseif occursin("λ",id) | occursin("lamb",id)
		number_as_string=split(id,"_")[end]
		number=parse(Int, number_as_string)
		par_id = 4 + number
	elseif occursin("lin",id) | occursin("div",id)
		number_as_string=split(id,"_")[end]
		number=parse(Int, number_as_string)
		par_id = 4 + n_edges + n_tips + number
	elseif occursin("rat",id)
		number_as_string=split(id,"_")[end]
		number=parse(Int, number_as_string)
		par_id = n_par - 2 + number
	else
		par_id = 1
	end

	return par_id
end

function plot_density(co::CladsOutput, id_par::String; burn = 0.25)
	id = convert_id(co, id_par)
	plot_density(co, id, burn = burn)
end

function plot_chain(co::CladsOutput, id_par::String; burn = 0.25)
	id = convert_id(co, id_par)
	plot_chain(co, id, burn = burn)
end
