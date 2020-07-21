#=

  "ClaDS: Cladogenetic XX" submodule package

=#

module ClaDS

using Random: randexp
using Distributions: Bernoulli, LogNormal, Weights, sample
using RCall: @rput, reval

# other submodules dependencies
using ..Utils

# files
include("clads/tree_class.jl")
include("clads/sim_ClaDS_utils.jl")
include("clads/sim_ClaDS.jl")
include("clads/plot_ClaDS.jl")


export sim_ClaDS2_ntips, plot_ClaDS

end # module ClaDS
