# PANDA.jl    Phylogenetic ANalyses of DiversificAtion

```@meta
CurrentModule = PANDA
using PANDA
```

Implements macroevolutionary analyses on phylogenetic trees.

## Installation

You can install the package by typing

```julia
julia> using Pkg
julia> Pkg.add("PANDA")
```

PANDA uses R functions and packages for plotting. If you want to be able to use the plotting functions, [the R language](https://www.r-project.org/) needs to be installed on your computer. You will also need a few R packages to be installed, including : [ape](https://cran.r-project.org/web/packages/ape/index.html), [coda](https://cran.r-project.org/web/packages/coda/index.html), [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [fields](https://cran.r-project.org/web/packages/fields/index.html). You can install them from a R session by typing

```r
install.packages("ape", "coda", "RColorBrewer", "fields")
```

You will then be able to load PANDA to Julia by typing

```julia
julia> using PANDA
```

## ClaDS

The ClaDS module implements the inference of ClaDS parameters on a phylogeny using data augmentation. A step by step presentation of how to perform the inference is available in the [manual](clads/tutorial.md)

```@autodocs
Modules = [PANDA.ClaDS]
Pages   = ["load_tree.jl", "infer_ClaDS2.jl", "clads_output.jl", "plot_ClaDS.jl","tree_class.jl", "Tree_utils.jl","clads_output.jl"]
```
