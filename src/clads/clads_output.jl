"""
A structure containig the informations about the resulot of a ClaDS run. Contains the following fields :

- `tree::Tree`: the phylogeny on which the analysis was performed.
- `chains::Array{Array{Array{Float64,1},1},1}` : the mcmc chains
- `rtt_chains::Array{Array{Array{Float64,1},1},1}` : the mcmc chains with the rate through time
- `σ_map::Float64` : the σ parameter estimate
- `α_map::Float64` : the α parameter estimate
- `ε_map::Float64` : the ε parameter estimate
- `λ0_map::Float64` : the initial speciation rate estimate
- `λi_map::Array{Float64,1}` : the estimates of the branh specific speciation rates
- `λtip_map::Array{Float64,1}` : the estimates of the tip speciation rates
- `DTT_mean::Array{Float64,1}` : the diversity through time estimates
- `RTT_map::Array{Float64,1}` : the rate through time estimates
- `time_points::Array{Float64,1}` : the times at which `DTT_mean` and `RTT_map` are computed
- `enhanced_trees::Array{Tree,1}` : a sample of the complete tree distribution
- `gelm::Tuple{Int64,Float64}` : the gelman parameter
- `current_state` : other variables, by `infer_ClaDS` used to continue the run
"""
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
    reason_for_stop::String
    current_state # all the rest that is mising to continue the run
end

# function CladsOutput()
#     return CladsOutput(
#         Tree(),
#         Array{Array{Array{Float64,1},1},1}(undef,0),
#         Array{Array{Array{Float64,1},1},1}(undef,0),
#         0.,
#         0.,
#         0.,
#         0.,
#         [0.],
#         [0.],
#         [0.],
#         [0.],
#         [0.],
#         [Tree()],
#         (0,10.),
#         0.,
#         "",
#     )
# end

function CladsOutput()
    return CladsOutput(
        Tree(),
        Array{Array{Array{Float64,1},1},1}(undef, 0),
        Array{Array{Array{Float64,1},1},1}(undef, 0),
        0.0,
        0.0,
        0.0,
        0.0,
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [Tree()],
        (0, 10.0),
        "",
        ([], [], [], [], [])  # This initializes current_state
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


"""
    plot_CladsOutput(co::CladsOutput ; method = "tree", ...)

Plots various aspects of the output of ClaDS

# Arguments
- `co::CladsOutput` : A `CladsOutput` object, the output of a ClaDS run.

# Keyword arguments
- `method::String` : A `String` indicating what aspect of the output should be plotted, see details.
"""
function plot_CladsOutput(co::CladsOutput ; method = "tree",
    nplot = 50, alpha_col = 0.05, options = "", id_par = 1, force_ylim = Array{Float64,1}(undef,0))

    if method == "tree"
        plot_ClaDS(co.tree,co.λi_map, options = options)
    elseif method == "DTT"
        plot_DTT(co, n_ltt = nplot, alpha_col = alpha_col, force_ylim = force_ylim)
    elseif method == "RTT"
        plot_RTT(co, nplot = nplot, alpha_col = alpha_col, force_ylim = force_ylim)
    elseif method == "density"
        plot_density(co, id_par)
    elseif method == "chain"
        plot_chain(co, id_par)
    end
end


function tip_rate(co::CladsOutput, sp_name::String)
    number = 0
    dist = Inf
    tip_labs = tip_labels(co.tree)
    for itl in 1:length(tip_labs)
        d = Levenshtein()(sp_name, tip_labs[itl])
        if d < dist
            dist = d
            number = itl
            if d == 0
                break
            end
        end
    end
    return co.λtip_map[number]
end
