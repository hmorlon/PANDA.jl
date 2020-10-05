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

function print_CladsOutput(co::CladsOutput)
    n = n_tip(co.tree)
    nit = length(co.chains[1][1])
    to_print = string("ClaDS run for a $n tipped tree","\n",
        "$nit iterations, the gelman statistics in $(co.gelm[2])","\n",
        "inferred parameters : σ = $(co.σ_map), α = $(co.α_map), ε = $(co.ε_map), λ0 = $(co.λ0_map)"
        )
    println(to_print)
end


function plot_CladsOutput(co::CladsOutput ; method = "tree",
    nplot = 50, alpha_col = 0.05, options = "", id_par = 1)

    if method == "tree"
        plot_ClaDS(co.tree,co.λi_map, options = options)
    elseif method == "DTT"
        plot_DTT(co, n_ltt = nplot, alpha_col = alpha_col)
    elseif method == "RTT"
        plot_RTT(co, nplot = nplot, alpha_col = alpha_col)
    elseif method == "density"
        plot_density(co, id_par)
    elseif method == "chain"
        plot_chain(co, id_par)
    end
end
