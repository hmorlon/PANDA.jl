#=
R ploting function
=#

function create_plot_ClaDS_ape()
    reval("""
    require(fields)
    require(ape)
    require(RColorBrewer)

    plot_ClaDS=function(phylo,rate1,rate2=NULL,same.scale=T,main=NULL,lwd=1,log=T, show_labels=F, minr = Inf, maxr = 0, show_legend = T,...){
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
function plot_ClaDS(tree::Tree ; id = 1, ln=true, lwd=3, show_labels = false, options="")
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
    @rput lwd
    @rput show_labels
    @rput tip_labels

    create_plot_ClaDS_ape()
    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
        colors = plot_ClaDS(tree, rates, log=ln, lwd=lwd, show_labels = show_labels $opt);
    """)
end

# specifying the tree and a vector of rates
function plot_ClaDS(tree::Tree, rates ; id = 1, ln=true, lwd=3, round = false, options = "", show_labels=false)
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
    @rput lwd
    @rput show_labels
    @rput tip_labels

    create_plot_ClaDS_ape()
    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
        plot_ClaDS(tree, rates, log=ln, lwd=lwd, show_labels = show_labels $opt)
    """)
end

# specifying the tree and two vectors of rates, both kind of rates are plotted on the same color scale
function plot_ClaDS(tree::Tree, rates1, rates2 ; id = 1, ln=true, lwd=3)
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
    @rput lwd

    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = 1:ntip)
        class(tree) = "phylo"
        leg = plot_ClaDS_noLeg(tree, rates1, rates2, log=ln, lwd=lwd)
        leg = plot_ClaDS_noLeg(tree, rates2, rates1, log=ln, lwd=lwd)
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
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})

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
    plot_id = Array{Int64,1}(floor.(range(2+ floor(burn * N), length=nplot, stop=N)))
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
        for k in plot_id
            y = co.rtt_chains[i][k]
            for j in 1:length(y)
                push!(map_mr[j],log(y[j]))
            end
            @rput y
            reval("""lines(t, log(y), col = alpha(color, alpha = alpha_col), lwd = 2)""")
            mean_mr .+= y
            n += 1
        end
    end

    mean_mr ./= n
    @rput mean_mr
    @rput map_mr
    reval("""
        maps = log(RTT_map) #sapply(map_mr, function(x){D=density(x); return(D[[1]][which.max(D[[2]])])})
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
