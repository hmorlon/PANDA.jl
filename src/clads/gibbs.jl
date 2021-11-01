function newton_LambertW(x ; rel_tol=0.0001)
    w0 = 0
    if x<10
        w1 = 1
    else
        w1 = max(1,log(x)-log(log(x))+ log(log(x))/log(x))
    end

    while abs((w0-w1)/min(abs(w0),abs(w1))) > rel_tol
        w0 = w1
        w1 = w1 - (w1 - x * exp(-w1))/(1+w1)
    end
    return w1
end

function slicing(x0_eff, σ_eff, bl_eff, n, former_λ ; n_it = 5, bounds = 5)
    a = x0_eff + n * σ_eff^2
    t_eff = bl_eff * σ_eff^2

    if (exp(a) * t_eff) < 10000
        best_λi = a - newton_LambertW(t_eff * exp(a), rel_tol=0.000001)
    else
        lt = log(t_eff) + a
        best_λi = a - lt + log(lt) - (log(lt)/lt)
    end
    max_λi = best_λi + bounds * σ_eff
    min_λi = best_λi - bounds * σ_eff
    extant_λi = max_λi - min_λi

    max_f = n * best_λi - bl_eff * exp(best_λi) - (best_λi - x0_eff)^2/(2 * σ_eff^2)
    j = 1.
    function fx(x ; bl = bl_eff, ne = n, x0 = x0_eff, s = 2*σ_eff^2)
         exp(x * ne - bl * exp(x) - (x-x0)^2/s - max_f)
     end

    λ = best_λi
    fλ = fx(λ)

    if -1e10 < max_f < 1e10
        for j in 1:n_it
            min_eff = min_λi
            max_eff = max_λi

            u = fλ * rand()
            reject = true
            if u < 1
                n_iter = 0
                while reject
                    n_iter += 1
                    λ = rand() * (max_eff - min_eff) + min_eff
                    fλ = fx(λ)
                    reject = (fλ < u)
                    if reject
                        if λ < best_λi
                            min_eff = λ
                        else
                            max_eff = λ
                        end
                    end

                end
            else
                λ = best_λi
            end
        end
    else
        λ = log(former_λ)
        #new_λ = former_λ
    end

    new_λ = exp(λ)
    if isnan(new_λ) || new_λ <= 0.
        new_λ = former_λ#rates[i + 1]
    end

    return new_λ
end

function slicing0(S, n_events, former_λ0 ; n_it = 5, bounds = 5)

    best_λ0 = log(n_events/S)
    max_λ0 = best_λ0+20
    min_λ0 = best_λ0-20
    extant_λ0 = max_λ0 - min_λ0

    max_f = n_events * (best_λ0 - 1)

    function fx(x ; ne = n_events, s = S)
         exp(x * ne - s * exp(x) - max_f)
     end

    λ = best_λ0
    fλ = fx(λ)

    if -1e10 < max_f < 1e10
        for j in 1:n_it
            min_eff = min_λ0
            max_eff = max_λ0

            u = fλ * rand()
            reject = true
            if u < Inf
                while reject
                    λ = rand() * (max_eff - min_eff) + min_eff
                    if max_eff - min_eff <= 0
                        println(" oh oh $j, $(max_eff - min_eff), $fλ, $u, $max_f, $min_eff, $max_eff, $best_λ0, $min_λ0, $max_λ0")
                        break
                    end
                    fλ = fx(λ)
                    reject = fλ < u
                    if reject
                        if λ < best_λ0
                            min_eff = λ
                        else
                            max_eff = λ
                        end
                    end

                end
            else
                λ = log(former_λ0)
            end
        end
    else
        new_λ = former_λ0
    end

    new_λ0 = exp(λ)

    if isnan(new_λ0) || new_λ0 <= 0.
        new_λ0 = former_λ0#rates[i + 1]
    end

    return new_λ0
end

function extract_S(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1})

    S::Float64 = 0.
    for et in edge_trees
        S += et.stem_rate[1] * et.effective_bl
    end

    return S

end

function extract_nextinct(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1})

    n::Int64 = 2 - n_tip(tree)
    i=0
    for et in edge_trees
        n += n_tip(et.tree) - et.tip_number
    end

    return n
end


function n_tip(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1})
    n = n_tip(tree)
    for i in 1:length(edge_trees)
        n += n_tip(edge_trees[i].tree) - 1
    end
    return n
end

function n_tip(tree::Tree)
    return Int64((tree.n_nodes+1)/2)
end


function draw_λ0!(tree::Tree, ε::Float64, rates::Array{Float64,1},
    edge_trees::Array{EdgeTreeRates2,1})

    S = extract_S(tree, edge_trees)
    S *= 1+ε
    n = n_tip(tree, edge_trees)
    n_events = 2 * n - 2 - n_extant_tips(tree,edge_trees)
    former_λ0 = rates[1]
    S /= former_λ0
    λ0 = slicing0(S, n_events, former_λ0)

    #rates[1] = λ0
    ratio = λ0 / former_λ0
    rates .*= ratio

    update_rates!(tree, ratio)
    for i in 1:length(edge_trees)
        edge_trees[i].stem_rate[1] *= ratio
    end
end

function draw_σ(relative_rates, α ; α0 = 1., β0 = 0.01)
    n = length(relative_rates)
    m = mean(relative_rates)
    α_n = α0 + n/2
    log_α = log(α)#m #log(m)
    β_n = β0

    for rate in relative_rates
        β_n += (rate - log_α)^2 /2
    end

    if isnan(β_n)
        println("inv gam $α0 $β0; $α_n , $β_n ;    ")
        print(relative_rates)
    end
    sqrt(rand(InverseGamma(α_n,β_n)))
end

function draw_α(relative_rates, σ ; α_0 = 0., σ_0 = 1)
    # know posterior for the mean of a normal with normal prior, see for example
    # www.ams.sunysb.edu/~zhu/ams570/Bayesian_Normal.pdf

    # so here we have a normal prior on log(α) (with default parameters (α_0,σ_0)=(0,1))
    n = length(relative_rates)
    #α_n = (α_0 * σ^2 + σ_0 * sum(map(x -> x^2, relative_rates)))/(σ^2 + n * σ_0^2)
    σ_n = sqrt((σ^2 * σ_0^2)/(σ^2 + n * σ_0^2))
    α_n = (α_0 * σ^2 + σ_0^2 * sum(relative_rates))/(σ^2 + n * σ_0^2)


    a = exp(rand(Normal(α_n, σ_n)))
    if a==0
        println("$α_n , $σ_n, $α_0 ; $σ ; $relative_rates")
    end
    return a
end

function draw_ε_crown(S::Float64, n_extinct::Int64, n_cond::Int64 ; n_it = 10) where {T<:Number}

    n_extinct_eff = n_extinct + 1  # because we want a flat prior on the nonloged value
    n_cond_eff = n_cond - 1
    δ = sqrt((n_extinct_eff + n_cond_eff - S)^2 + 4*n_extinct_eff*S)
    best_ε = log(max((n_extinct_eff + n_cond_eff - S + δ)/(2*S),(n_extinct_eff + n_cond_eff - S - δ)/(2*S)))
    max_f = best_ε * n_extinct_eff +
        log(exp(best_ε)+1) * n_cond_eff -
        S * exp(best_ε)

    function fx(x ; ne = n_extinct_eff, nc = n_cond_eff, s = S)
         exp(x * ne  + log(exp(x)+1) * nc - s * exp(x) - max_f)
     end


    λ = best_ε
    fλ = fx(λ)

    for j in 1:n_it
        u = fλ * rand()
        reject = true

        min_eff = best_ε - 10
        max_eff = best_ε + 10

        while reject
            λ = rand() * (max_eff - min_eff) + min_eff
            fλ = fx(λ)

            #println("$fλ , $u , $λ, $best_ε ")

            reject = fλ < u

            if reject
                if λ < best_ε
                    min_eff = λ
                else
                    max_eff = λ
                end
            end
        end
    end

    if !isnan(λ)
        return exp(λ)
    else
        return 0.001
    end
end

function draw_ε_crown(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, lefts::Array{Int64,1})
    S = extract_S(tree, edge_trees)
    n_extinct = extract_nextinct(tree, edge_trees)
    crown_edges = [1, 1 + lefts[1]]

    n_cond = 0
    if tree.offsprings[1].n_nodes > 1
        n_cond += 1
    end
    if tree.offsprings[2].n_nodes > 1
        n_cond += 1
    end
    bl = [tree.offsprings[1].branch_length, tree.offsprings[2].branch_length]
    for edge in 1:2
        if edge_trees[crown_edges[edge]].tree.n_nodes > 1
            n_cond -= 1
            ltt = LTT(edge_trees[crown_edges[edge]].tree)
            s = edge_trees[crown_edges[edge]].scale
            b = bl[edge]*s
            for e in 1:length(ltt[2])
                if ltt[2][e] == 1
                    n_cond += 1
                end
                if ltt[1][e] >= b
                    break
                end

            end
        end
    end

    return draw_ε_crown(S, n_extinct, n_cond)
end

function draw_ε_crown_priorUnif(S::Float64, n_extinct::Int64, n_cond::Int64 ; n_it = 10) where {T<:Number}

    n_extinct_eff = n_extinct + 1  # because we want a flat prior on the nonloged value
    n_cond_eff = n_cond - 1
    δ = sqrt((n_extinct_eff + n_cond_eff - S)^2 + 4*n_extinct_eff*S)
    best_ε = log(max((n_extinct_eff + n_cond_eff - S + δ)/(2*S),(n_extinct_eff + n_cond_eff - S - δ)/(2*S)))
    best_ε = min(best_ε,log(990))
    max_f = best_ε * n_extinct_eff +
        log(exp(best_ε)+1) * n_cond_eff -
        S * exp(best_ε)

    xmax = log(1000)

    function fx(x ; ne = n_extinct_eff, nc = n_cond_eff, s = S)
        if x > xmax
            0
        else
            exp(x * ne  + log(exp(x)+1) * nc - s * exp(x) - max_f)
        end
     end


    λ = best_ε
    fλ = fx(λ)
    print(fλ)
    for j in 1:n_it
        u = fλ * rand()
        reject = true

        min_eff = best_ε - 10
        max_eff = min(best_ε + 10, xmax)

        while reject
            λ = rand() * (max_eff - min_eff) + min_eff
            fλ = fx(λ)

            #println("$fλ , $u , $λ, $best_ε ")

            reject = fλ < u

            if reject
                if λ < best_ε
                    min_eff = λ
                else
                    max_eff = λ
                end
            end
        end
    end

    if !isnan(λ)
        return exp(λ)
    else
        return 0.001
    end
end

function draw_ε_crown_priorUnif(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, lefts::Array{Int64,1})
    S = extract_S(tree, edge_trees)
    n_extinct = extract_nextinct(tree, edge_trees)
    crown_edges = [1, 1 + lefts[1]]

    n_cond = 0
    if tree.offsprings[1].n_nodes > 1
        n_cond += 1
    end
    if tree.offsprings[2].n_nodes > 1
        n_cond += 1
    end
    bl = [tree.offsprings[1].branch_length, tree.offsprings[2].branch_length]
    for edge in 1:2
        if edge_trees[crown_edges[edge]].tree.n_nodes > 1
            n_cond -= 1
            ltt = LTT(edge_trees[crown_edges[edge]].tree)
            s = edge_trees[crown_edges[edge]].scale
            b = bl[edge]*s
            for e in 1:length(ltt[2])
                if ltt[2][e] == 1
                    n_cond += 1
                end
                if ltt[1][e] >= b
                    break
                end

            end
        end
    end

    return draw_ε_crown_priorUnif(S, n_extinct, n_cond)
end

function draw_ε_crown_priorln(S::Float64, n_extinct::Int64, n_cond::Int64 ; n_it = 10, logε0 = 0., sd = 0.5) where {T<:Number}

    n_extinct_eff = n_extinct
    n_cond_eff = n_cond - 1
    δ = sqrt((n_extinct_eff + n_cond_eff - S)^2 + 4*n_extinct_eff*S)

    best_ε = log(max((n_extinct_eff + n_cond_eff - S + δ)/(2*S),(n_extinct_eff + n_cond_eff - S - δ)/(2*S)))
    max_f = best_ε * n_extinct_eff +
        log(exp(best_ε)+1) * n_cond_eff -
        S * exp(best_ε)

    function fx(x ; ne = n_extinct_eff, nc = n_cond_eff, s = S)
         exp(x * ne  + log(exp(x)+1) * nc - s * exp(x) - max_f - (x - logε0)^2/(2*sd^2))
     end
     function df(x;  ne = n_extinct_eff, nc = n_cond_eff, s = S)
         ne + exp(x)/(exp(x) + 1) * nc - s * exp(x) - (x - logε0)/(sd^2)
    end


    λ = best_ε
    fλ = fx(λ)

    for j in 1:n_it
        u = fλ * rand()
        reject = true

        min_eff = best_ε - 10
        max_eff = best_ε + 10

        while reject
            λ = rand() * (max_eff - min_eff) + min_eff
            fλ = fx(λ)

            #println("$fλ , $u , $λ, $best_ε ")

            reject = fλ < u

            dfλ = df(λ)
            if reject
                if dfλ > 0
                    min_eff = λ
                else
                    max_eff = λ
                end
            end
        end
    end

    if !isnan(λ)
        return exp(λ)
    else
        return 0.001
    end
end

function draw_ε_crown_priorln(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, lefts::Array{Int64,1} ; logε0 = 1., sd = 0.5)
    S = extract_S(tree, edge_trees)
    n_extinct = extract_nextinct(tree, edge_trees)
    crown_edges = [1, 1 + lefts[1]]

    n_cond = 0
    if tree.offsprings[1].n_nodes > 1
        n_cond += 1
    end
    if tree.offsprings[2].n_nodes > 1
        n_cond += 1
    end
    bl = [tree.offsprings[1].branch_length, tree.offsprings[2].branch_length]
    for edge in 1:2
        if edge_trees[crown_edges[edge]].tree.n_nodes > 1
            n_cond -= 1
            ltt = LTT(edge_trees[crown_edges[edge]].tree)
            s = edge_trees[crown_edges[edge]].scale
            b = bl[edge]*s
            for e in 1:length(ltt[2])
                if ltt[2][e] == 1
                    n_cond += 1
                end
                if ltt[1][e] >= b
                    break
                end

            end
        end
    end

    return draw_ε_crown_priorln(S, n_extinct, n_cond,  logε0 = logε0, sd = sd)
end

function draw_λi_quad!(rates::Array{Float64,1}, edge_trees::Array{EdgeTreeRates2,1},
    σ::Float64, α::Float64, ε::Float64, tree::Tree, lefts::Array{Int64,1}; n_it = 5, bounds = 10)

    log_α = log(α)
    function aux(subtree, id, parent_rate, sampled_rates)

        if subtree.n_nodes < 2
            former_λ = subtree.attributes[1]
            x0_eff = parent_rate + log_α

            n = 0
            bl_eff = 0.
            if id>0
                bl_eff = edge_trees[id].effective_bl * (1 + ε)
                #s = sum(extract_branch_lengths(edge_trees[id].tree) .*extract_rates(edge_trees[id].tree) / edge_trees[id].scale)
                n = edge_trees[id].tree.n_nodes - edge_trees[id].tip_number
            end

            λ = slicing(x0_eff, σ, bl_eff, n, former_λ)
            pushfirst!(sampled_rates, λ)


            return  n, bl_eff*λ

        else

            tip_λ = subtree.attributes[1]
            former_λ = subtree.attributes[1]
            if id>0
                tip_λ *= edge_trees[id].tip_rate
            end
            tip_λ = log(tip_λ)

            n_right, bl_right = aux(subtree.offsprings[2], id + 1 + lefts[id+1], tip_λ, sampled_rates)
            n_left, bl_left = aux(subtree.offsprings[1], id + 1, tip_λ, sampled_rates)

            x0_eff = parent_rate + log_α
            n = n_left + n_right
            bl_eff = (bl_right + bl_left) / former_λ

            λ = 0.

            if id>0
                bl_eff += edge_trees[id].effective_bl * (1 + ε)
                n += edge_trees[id].tree.n_nodes - edge_trees[id].tip_number

                λ = slicing(x0_eff, σ, bl_eff, n, former_λ)
            else
                λ = slicing0(bl_eff, n, former_λ)
            end

            ratio = λ/former_λ

            for k in 1:(subtree.n_nodes - 1)
                sampled_rates[k] *= ratio
            end

            pushfirst!(sampled_rates, λ)
            return  n, bl_eff*λ
        end
    end

    r = Array{Float64,1}(undef,0)
    aux(tree, 0, 0., r);
    rates[1] = r[1]

    for i in 1:length(edge_trees)
        rates[i+1] = r[i+1]
        edge_trees[i].stem_rate[1] = rates[i+1]
    end

    update_rates!(tree,rates)
end

function draw_λi_quad!(rates::Array{Float64,1}, edge_trees::Array{EdgeTreeRates2,1}, σ::Float64, α::Float64, ε::Float64, tree::Tree; n_it = 5, bounds = 10)
    lefts = n_left(tree)
    draw_λi_quad!(rates, edge_trees, σ, α, ε, tree, lefts, n_it = n_it, bounds = bounds)
end

function draw_λ0_slicing!(rates::Array{Float64,1}, edge_trees::Array{EdgeTreeRates2,1}, σ::Float64, α::Float64, ε::Float64, lefts; n_it = 5, bounds = 20)

    parent_λ = NaN

    daughter_λ = map(log,[rates[2],rates[2+lefts[1]]])
    x0_eff = (-2*log(α)+daughter_λ[1]+daughter_λ[2])/2
    σ_eff = σ / sqrt(2)

    law = Normal(x0_eff, σ_eff)
    new_λ = exp(rand(law,1)[1])
    rates[1] = new_λ

    return new_λ
end
