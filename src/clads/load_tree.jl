function make_ape(tree::Tree ; id = 1)
    n_nodes = tree.n_nodes
    node_id = (n_nodes + 3)/2

    function aux(subtree, root, next_node, tip)
        if subtree.n_nodes == 0
            return Array{Int64,2}(undef,0,2), Array{Float64,1}(undef,0), Array{Float64,1}(undef,0), next_node, tip, Array{String,1}(undef,0)
        elseif length(subtree.offsprings) == 0
            return [root tip], subtree.branch_length, subtree.attributes[id], next_node, tip+1, subtree.label
        else
            tupple1 = aux(subtree.offsprings[1], next_node, next_node+1, tip)
            tupple2 = aux(subtree.offsprings[2], next_node, tupple1[4], tupple1[5])
            if subtree.branch_length == 0.
                return [tupple1[1]; tupple2[1]],
                    [tupple1[2]; tupple2[2]],
                    [tupple1[3]; tupple2[3]],
                    tupple2[4], tupple2[5],
                    [tupple1[6]; tupple2[6]]
            else
                return [[root next_node]; tupple1[1]; tupple2[1]],      # edges
                    [subtree.branch_length; tupple1[2]; tupple2[2]],    # branch lengths
                    [subtree.attributes[id]; tupple1[3]; tupple2[3]],   # rates
                    tupple2[4], tupple2[5],                             # auxiliary variables
                    [tupple1[6]; tupple2[6]]                            # tip labels
            end
        end
    end

    tupple = aux(tree, -1, node_id, 1)
    return tupple[1], tupple[2], tupple[3], tupple[6]
end

function build_tree(branches::Array{Int64,2}, branch_lengths, attributes::Array{Array{T2,1},1}, tip_labels::Array{String,1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T2,1}(undef,0)) where {T2<:Number}

    n_attributes = length(attributes)
    if length(extinct) == 0
        dead_branches = fill(false, size(branches)[1])
    else
        dead_branches = extinct
    end

    if size(branches)[1] == 0
        if return_void
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),stem_age,root_attributes)
        end
    end
    if size(branches)[1] == 1
        if ! dead_branches[1]
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, true)
        elseif prune_extinct
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, false)
        end
    end

    root = -1
    for i in branches[:,1]
        has_no_parent = true
        for j in branches[:,2]
            if i == j
                has_no_parent = false
                break
            end
        end
        if has_no_parent
            root = i
            break
        end
    end

    max_node = maximum(branches[:,1:2])+1
    offsprings = fill(-1,max_node,2)
    parent_edges = fill(-1,max_node)
    daughter_edges = fill(-1,max_node)

    for i in 1:size(branches)[1]
        individual = branches[i,1]+1
        parent_edges[individual] = i
        daughter_edges[branches[i,2]+1] = i
        if offsprings[individual,1] < 0
            offsprings[individual,1] = branches[i,2]
        else
            offsprings[individual,2] = branches[i,2]
        end

        if dead_branches[i]
            individual = branches[i,2]+1
            offsprings[individual,:]= [-2 -2]
        end
    end

    if prune_extinct
        function aux(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], tip_labels[node])
            elseif offsprings[node+1,1]==-2
                return Tree()
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes #repeat([Array{Float64}(undef,0)],n_attributes)
                end

                tree_left = aux(offsprings[node+1,1])
                tree_right = aux(offsprings[node+1,2])

                if tree_left.n_nodes + tree_right.n_nodes == 0
                    return Tree()
                elseif tree_left.n_nodes == 0
                    tree_right = Tree(tree_right.offsprings, tree_right.branch_length + add_time,
                        add_attribute, tree_right.n_nodes, tree_right.label)

                    return tree_right
                elseif tree_right.n_nodes == 0
                    tree_left = Tree(tree_left.offsprings, tree_left.branch_length + add_time,
                        add_attribute, tree_left.n_nodes, tree_left.label)
                    return tree_left
                else
                    return Tree([tree_left, tree_right], add_time, add_attribute)
                end
            end
        end

        return aux(root)

    else
        function aux_all(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], true, tip_labels[node])
            elseif offsprings[node+1,1]==-2
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], false)
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes#repeat([Array{Float64}(undef,0)],n_attributes)
                end
                tree_left = aux_all(offsprings[node+1,1])
                tree_right = aux_all(offsprings[node+1,2])
                return Tree([tree_left, tree_right], add_time, add_attribute)
            end
        end

        return aux_all(root)

    end
end

function load_tree(file)
    if file[(end-3):end] == ".tre"
        @rput file
        reval("""
            require(ape)
            temp_tree = read.tree(file)
        """)
    elseif file[(end-3):end] == ".nex"
        @rput file
        reval("""
            require(ape)
            temp_tree = read.nexus(file)
        """)
    end
    @rget temp_tree

    reval(""" temp_tree = list() """)
    return ape2Tree(temp_tree)
end
