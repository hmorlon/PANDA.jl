#=

  "PANDA: Phylogenetic ANalyses of DiversificAtion" module package

=#

module PANDA

__precompile__(true) #put back to true when checking done

#=
 Submodules
=#

# ClaDS: Cladogenetic XX
include("ClaDS.jl")


#=
 Exported functions
=#
using .ClaDS: Tree, tip_labels, n_tips, sample_tips,
    sim_ClaDS2_ntips, plot_ClaDS, infer_ClaDS,
    CladsOutput, print_CladsOutput, plot_CladsOutput, tip_rate,
    load_tree, save_ClaDS_in_R
export Tree, tip_labels, n_tips, sample_tips,
    sim_ClaDS2_ntips, plot_ClaDS, infer_ClaDS,
    CladsOutput, print_CladsOutput, plot_CladsOutput, tip_rate,
    load_tree, save_ClaDS_in_R
end
