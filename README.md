# PANDA: Phylogenetic ANalyses of DiversificAtion
Implements macroevolutionary analyses on phylogenetic trees.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hmorlon.github.io/PANDA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hmorlon.github.io/PANDA.jl/dev)
[![Build Status](https://travis-ci.com/hmorlon/PANDA.jl.svg?branch=master)](https://travis-ci.com/hmorlon/PANDA.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/hmorlon/PANDA.jl?svg=true)](https://ci.appveyor.com/project/hmorlon/PANDA-jl)
[![Codecov](https://codecov.io/gh/hmorlon/PANDA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/hmorlon/PANDA.jl)



## Installation

First, Julia needs to be installed on your computer. If this is not the case, you can install it by following the instructions provided on [this page](https://julialang.org/downloads/)

Within a Julia session, you can install the package by typing

```julia
julia> using Pkg
julia> Pkg.add("PANDA")
```

PANDA uses R functions and packages for plotting. If you want to be able to use the plotting functions, [the R language](https://www.r-project.org/) needs to be installed on your computer. You will also need a few R packages to be installed, including : [ape](https://CRAN.R-project.org/package=ape), [coda](https://CRAN.R-project.org/package=coda), [RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer), [fields](https://CRAN.R-project.org/package=fields). You can install them from a R session by typing

```r
> install.packages("ape", "coda", "RColorBrewer", "fields")
```

You will then be able to load PANDA to Julia by typing

```julia
julia> using PANDA
```

## Help

The help of the latest released version of the package is available using [this link](https://hmorlon.github.io/PANDA.jl/stable)
