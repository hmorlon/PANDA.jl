#=

  "ClaDS: Cladogenetic XX" submodule package

=#

module ClaDS

using Random: randexp
using Distributions: Bernoulli, LogNormal, Weights, sample, Normal, InverseGamma
using RCall: @rput, @rget, reval
using Statistics: mean, var

# other submodules dependencies
using ..Utils

# files
include("clads/tree_class.jl")
include("clads/clads_output.jl")
include("clads/enhance_tree_utils.jl")
include("clads/sim_ClaDS_utils.jl")
include("clads/sim_ClaDS.jl")
include("clads/load_tree.jl")
include("clads/plot_ClaDS.jl")
include("clads/infer_ClaDS2.jl")
include("clads/mcmc_ClaDS2_utils.jl")
include("clads/sample_fractions.jl")
include("clads/LTT.jl")
include("clads/n_tips.jl")
include("clads/Tree_utils.jl")
include("clads/branch_lengths.jl")
include("clads/RTT.jl")
include("clads/enhance_tree.jl")
include("clads/rates.jl")
include("clads/gibbs.jl")

reval("""
        list.of.packages <- c("ape", "coda")   
        new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
        if(length(new.packages)) install.packages(new.packages)
    """)

end # module ClaDS
