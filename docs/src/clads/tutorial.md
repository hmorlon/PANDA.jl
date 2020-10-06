# ClaDS tutorial

## Loading a tree

```julia
my_tree = load_tree(tree_path)
```

## Running ClaDS

The parameter inference is ran with the function [`infer_ClaDS`](@ref)

## Tree simulation

```julia
t20 = sim_ClaDS2_ntips(200,0.2,1.,0.01,10.);
plot_ClaDS(t20, options = "type = 'fan'")
```
