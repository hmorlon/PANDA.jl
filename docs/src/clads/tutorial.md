# ClaDS tutorial

## Loading a tree

```julia
my_tree = load_tree(tree_path)
```

## Running ClaDS

The parameter inference is ran with the function [`infer_ClaDS`](@ref)

```julia
output = infer_ClaDS(my_tree, 200)
```

The result is a `CladsOutput` object, that contains the following fields:
- `tree::Tree`: the `Tree` object on which the inference was performed
