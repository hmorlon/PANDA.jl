using Revise
using PANDA
using RCall
using Random
#include("PANDA.jl")

Random.seed!(813)
t20 = sim_ClaDS2_ntips(100,0.2,1.,0.01,10.);
plot_ClaDS(t20, options = "type = 'fan'")

infered = infer_ClaDS(t20, 100, plot_chain = true,
    plot_tree=0,print_state=50, thin = 1)

infered[1][9]

reval("par(mfrow=c(3,4))")
for i in 1:10
    plot_ClaDS(infered[1][9][i])
end
