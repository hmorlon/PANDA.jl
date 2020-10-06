function sample_fractions(tree, sampling_at_tips)
    function aux(subtree, sub_sampling, x)
        if subtree.n_nodes < 2
            s = pop!(sub_sampling)
            pushfirst!(x,s)
            return s
        else
            s_left = aux(subtree.offsprings[2], sub_sampling,x)
            n_left = (subtree.offsprings[2].n_nodes+1)/2
            s_right = aux(subtree.offsprings[1], sub_sampling,x)
            n_right = (subtree.offsprings[1].n_nodes+1)/2
            s = (n_left + n_right)/(n_left/s_left + n_right/s_right)
            pushfirst!(x,s)

            return s
        end
    end

    fs = Array{Float64,1}(undef,0)
    aux(tree, deepcopy(sampling_at_tips), fs)
    return fs
end
