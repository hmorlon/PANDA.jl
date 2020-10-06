# ClaDS tutorial

The `ClaDS` module implements the Data Augmentation inference method for the ClaDS model, that allows estimating branch specific speciation rates on a reconstructed phylogeny.

## Loading a tree

You can import a phylogeny to the environment using the `load_tree` function. Currently supported extensions include `.tre` and `.nex`.

```julia
my_tree = load_tree(tree_path)
```

## Running ClaDS

The parameter inference is ran with the function [`infer_ClaDS`](@ref)

```julia
output = infer_ClaDS(my_tree)
```

### Incomplete sampling

By default, the function considers that the clade was perfectly sampled, i.e. that all the species alive at present time are included in the phylogeny. If it is not the case, the sampling fraction can be specified through the keyword argument `f`. `f` can be a `Float`, in which case the sampling fraction is taken as homogeneous on the whole phylogeny.

```julia
output = infer_ClaDS(my_tree, f = 0.94)
```

Alternatively, different sampling fractions can be specified for different subclades. To do so, `f` should be passed as a `Array{Float64}` of length `n`, where `n` is the number of tip in the phylogeny. `f[i]` is the sampling fraction of the subclade that contains tip `i`. If the `Tree` object has tip labels (which can be accessed using `tip_labels(my_tree)`, the sampling fractions in `f` are in the same order as the tip labels, and `f[i]` is the sampling fraction of the subclade that contains the tip with label `tip_labels(my_tree)[i]`.

In the following example, the left subtree of `my_tree`is assign the sampling fraction `0.3` and its right subtree the sampling fraction `0.8`.

```julia
#=
 create a vector of size n_tips(my_tree) where the first
 n_tips(my_tree.offsprings[1]) elements are equal to 0.3
 and the rest to 0.8
=#
f = [ i < n_tips(my_tree.offsprings[1]) ? 0.3 : 0.8 for i in 1:n_tips(my_tree)]
output = infer_ClaDS(my_tree, f = f)
```

### Result

The result is a `CladsOutput` object, that contains the following fields:
- `tree`: the `Tree` object on which the inference was performed.
- `chains`: the resulting mcmc chains.
- `rtt_chains`: the mcmc chains of  the mean rate through time.
- `σ_map`, `α_map`, `ε_map`, `λ0_map`: the MAP estimates of the model's parameters.
- `λi_map`, `λtip_map`: the MAP estimates of the branch-specific and present rates.
- `time_points`: Time points at which the number of lineages through time and rate through time are computed. The number of time points can be specified using the keyword argument `ltt_steps`.
- `DTT_mean`: Estimate of the number of lineages through time.
- `RTT_map`: Estimate of the mean rate through time.
- `enhanced_tree`: Sample from the complete phylogeny distribution. Their number can be specified through the keyword argument `n_trees`.
- `gelm`: Evaluation of the gelman statistics.

### Plot the branch specific rates

It can be plotted using the `plot_CladsOutput` function. By default, this function plots the reconstructed phylogeny painted with the inferred branch-specific speciation rates, but other methods are available.

```julia
plot_CladsOutput(output)
```

### Diversity through time plot

Using the keyword argument `method = "DTT"`, the function plots the estimate of the number of lineages through time.

```julia
plot_CladsOutput(output, method = "DTT")
```

On this plot, we have:
- black line: the LTT plot (number of lineages through time in the reconstructed phylogeny)
- thin blue lines: individual MCMC iterations
- thick blue lines: the $95\%$ confidence interval
- dotted green line: the point estimates

### Mean rate through time plot

Using the keyword argument `method = "RTT"`, the function plots the estimate of the mean speciation rate through time.

```julia
plot_CladsOutput(output, method = "RTT")
```

Similarly to the diversity through time plot, we have here:
- thin blue lines: individual MCMC iterations
- thick blue lines: the $95\%$ confidence interval
- dotted green line: the point estimates

### Marginal posterior densities
