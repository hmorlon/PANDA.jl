#=

  "PANDA: Phylogenetic ANalyses of DiversificAtion" module package

=#

module PANDA

__precompile__(true) #put back to true when checking done

#=
 Submodules
=#

# Utilities
include("Utils.jl")

# ESSE: Environmental and State Dependent Diversification
include("ESSE.jl")

# TRIBE: Trait and Range Interspecific Biogeographic Evolution
include("TRIBE.jl")

# ClaDS: Cladogenetic XX
include("ClaDS.jl")


#=
 Exported functions
=#

using .ESSE: esse, simulate_sse
export esse, simulate_sse

using .TRIBE: tribe, simulate_tribe
export tribe, simulate_tribe

using .ClaDS: Tree, tip_labels, n_tips,
    sim_ClaDS2_ntips, plot_ClaDS, infer_ClaDS,
    CladsOutput, print_CladsOutput, plot_CladsOutput, load_tree,
    save_ClaDS_in_R
export Tree, tip_labels, n_tips, 
    sim_ClaDS2_ntips, plot_ClaDS, infer_ClaDS,
    CladsOutput, print_CladsOutput, plot_CladsOutput, load_tree,
    save_ClaDS_in_R

end
