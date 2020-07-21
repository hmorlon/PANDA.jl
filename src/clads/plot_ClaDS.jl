#=
R ploting function
=#

function create_plot_ClaDS_ape()
    reval("""
    require(fields)
    require(ape)
    require(RColorBrewer)

    plot_ClaDS=function(phylo,rate1,rate2=NULL,same.scale=T,main=NULL,lwd=1,log=F, show_labels=F,...){
        Colors = colorRampPalette(rev(c('darkred',brewer.pal(n = 8, name = "Spectral"),'darkblue')))(100)

     if(nrow(phylo[[1]]) <=1){
        plot(1000,xlim = c(0,1), ylim = c(0,2), axes = F, xlab = "", ylab = "")
        lines(c(0,1), c(1,1), lwd=lwd, col="steelblue2")
      }else{
        if(is.null(rate2)){
          if(log) rate1=log(rate1)
          if(isTRUE(all.equal(rep(as.numeric(rate1[1]),length(rate1)),as.numeric(rate1)))){
            col=rep(1,length(rate1))
            plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,main=main,edge.width =lwd,...)
            if(log){
              image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T)
            }else{
              image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
            }
          }else{
            col = round( (rate1 - min(rate1)) / diff(range(rate1))*99   )+1
            plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,main=main,edge.width =lwd,...)
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
              image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
            }else{
              image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
            }
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
              if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
              # ticks=seq(min,max,length.out = 5)
              image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
            }else{
              image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T)
            }
          }else{
            par(mfrow=c(1,2))
            if(isTRUE(all.equal(rep(rate1[1],length(rate1)),rate1))){
              col=rep(1,length(rate1))
              plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
              if(log){

                image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T)
              }else{
                image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
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
                image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
              }else{
                image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
              }
            }
            if(isTRUE(all.equal(rep(rate2[1],length(rate2)),rate2))){
              col=rep(1,length(rate2))
              plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
              if(log){
                image.plot(z = c(exp(rate2[1]),2*exp(rate2[1])),col = Colors, horizontal=T,legend.only = T)
              }else{
                image.plot(z = c(rate2[1],2*rate2[1]),col = Colors, horizontal=T,legend.only = T)
              }
            }else{
              col = round(( (rate2 - min(rate2)) / (max(rate2)-min(rate2)))*99   )+1
              plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
              if(log){
                min=min(rate2)
                max=max(rate2)
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
                image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
              }else{
                image.plot(z = as.matrix(rate2),col = Colors, horizontal=T,legend.only = T)
              }
            }
          }
          par(mfrow=c(1,1))
        }
      }

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
        plot_ClaDS(tree, rates, log=ln, lwd=lwd, show_labels = show_labels $options)
    """)
end

# specifying the tree and a vector of rates
function plot_ClaDS(tree::Tree, rates ; id = 1, ln=true, lwd=3, round = false)
    plot_tree = Tree(tree.offsprings, 0., tree.attributes, tree.n_nodes)

    create_plot_ClaDS_ape()
    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
        plot_ClaDS(tree, rates, log=ln, lwd=lwd, show_labels = show_labels $options)
    """)
    if round
        reval("""
        source("/Users/maliet/ownCloud/My_folder/ClaDS_Julia/plot_ClaDS_round.R")
        """)
    else
        reval("""
        source("/Users/maliet/ownCloud/My_folder/ClaDS_Julia/plot_ClaDS.R")
        """)
    end

    edges, branch_lengths, new_rates = Tree2ape(plot_tree, id = 1)
    ntip = (tree.n_nodes + 1)/2

    @rput edges
    @rput branch_lengths
    @rput rates
    @rput ntip
    @rput ln
    @rput lwd

    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = 1:ntip)
        class(tree) = "phylo"
        plot_ClaDS(tree, rates, log=ln, lwd=lwd)
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
