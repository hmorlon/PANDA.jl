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

Once installed, it can be loaded to your environment with

```julia
julia> using PANDA
```

## ClaDS

The ClaDS module implements the inference of ClaDS parameters on a phylogeny using data augmentation.


```@autodocs
Module = [PANDA.ClaDS]
```

```@docs
PANDA.ClaDS.load_tree
```

```@docs
PANDA.ClaDS.infer_ClaDS
```

```@docs
PANDA.ClaDS.plot_ClaDS
```

```@docs
PANDA.ClaDS.sim_ClaDS2_ntips
```
