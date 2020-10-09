#=
The tree structure used for all the ClaDS functions
=#

"""
A phylogeny object for the module ClaDS. It is represented as a branch with a given length and optional attributes, and its daughter trees.

- `offsprings::Array{Tree,1}`: the two daughter trees.
- `branch_length::Float64`: the length of the branch
- `attributes::Array{T,1} where {T<:Number}`: soem attributes of the branch. Used to store the speciation rates.
- `n_nodes::Int64`: the number of nodes in the tree (internal nodes + tips)
- `extant::Bool`: does the tree have any extant species?
- `label::String`: if the tree is a tip (ie n_nodes == 1), contains the name of the corresponding species.
"""
struct Tree
    offsprings::Array{Tree,1}
    branch_length::Float64
    attributes::Array{T,1} where {T<:Number}
    n_nodes::Int64
    extant::Bool
    label::String
end

#=
Constructors (12 methods with various default options)
=#

# specify all
function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::T, n_nodes::Int64, extant::Bool, label::String) where {T<:Number}
    Tree(offsprings, branch_length, [attributes], n_nodes, extant, label)
end

# all except tree labels, set to empty strings
function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::T, n_nodes::Int64, extant::Bool) where {T<:Number}
    Tree(offsprings, branch_length, attributes, n_nodes, extant, "")
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::Array{T,1}, n_nodes::Int64, extant::Bool) where {T<:Number}
    Tree(offsprings, branch_length, attributes, n_nodes, extant, "")
end

# without extant information (set to true)
function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::Array{T,1}, n_nodes::Int64, label::String) where {T<:Number}
    Tree(offsprings, branch_length, attributes, n_nodes, true, label)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::T, n_nodes::Int64, label::String) where {T<:Number}
    Tree(offsprings, branch_length, [attributes], n_nodes, true, label)
end

# without tip labels (set to empty strings) nor extant information (set to true)
function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::Array{T,1}, n_nodes::Int64) where {T<:Number}
    Tree(offsprings, branch_length, attributes, n_nodes, true)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::T, n_nodes::Int64) where {T<:Number}
    Tree(offsprings, branch_length, [attributes], n_nodes, true)
end


# whithout specification of the node number (set to the sum of the node number of the offspring +1)
function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes)
    n_nodes = length(branch_length)
    extant = false
    for i in offsprings
        n_nodes += i.n_nodes
        extant = extant || i.extant
    end
    extant = extant || n_nodes == 1
    Tree(offsprings, branch_length, attributes, n_nodes, extant)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes, label::String)
    n_nodes = length(branch_length)
    extant = false
    for i in offsprings
        n_nodes += i.n_nodes
        extant = extant || i.extant
    end
    extant = extant || n_nodes == 1
    Tree(offsprings, branch_length, attributes, n_nodes, extant, label)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes, extant::Bool)
    n_nodes = length(branch_length)
    for i in offsprings
        n_nodes += i.n_nodes
    end
    Tree(offsprings, branch_length, attributes, n_nodes, extant)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes, extant::Bool, label::String)
    n_nodes = length(branch_length)
    for i in offsprings
        n_nodes += i.n_nodes
    end
    Tree(offsprings, branch_length, attributes, n_nodes, extant, label)
end

# empty tree
function Tree()
    Tree(Array{Tree,1}(undef,0), 0., Array{Float64,1}(undef,0), 0, false)
end
