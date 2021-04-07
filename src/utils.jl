#####
##### Statistics helper functions
#####

"""
Returns the log probability of the data under a uniform multinomial distribution.
"""
function ll_uniform_multinomial(data)
    k = length(data)
    n = sum(data)
    n == 0 && return 0.
    logpdf(Multinomial(n, ones(k)./k), data)
end
"""
Returns the maximum likelihood estimate of the Multinomial rates and the log
probability of the counts under the Multinomial with these rates
"""
function MLE(counts)
    Σ = sum(counts)
    Σ == 0 && return zeros(length(counts)), 0.
    mle = counts ./ Σ
    mle, logpdf(Multinomial(Σ, mle), counts)
end
"""
Returns the normalization constant of a Dirichlet distribution with parameters
counts .+ 1
"""
function logmultinomialnormalization(counts)
    N = sum(counts)
    loggamma(N + 1) - +(map(x -> loggamma(x + 1), counts)...)
end
"""
Returns the log probability of the data under a Dirichlet Multinomial with
parameters `alpha` (`= ones(length(counts))` by default).
"""
function logpdf_dirichletmultinomial(counts, alpha = ones(length(counts)))
    n = sum(counts)
    n == 0 && return 0.
    logpdf(DirichletMultinomial(n, alpha), counts)
end

function _f(y, j)
    j == 1 && return prod(y)
    d = length(y)
    j == d + 1 && return 1 - y[1]
    prod(y[i] for i in 1:d-j+1) * (1 - y[d-j+2])
end
function get_constraints(args...)
    y -> begin
        (&)([a[2](_f(y, a[1]), _f(y, a[3])) for a in args]...)
    end
end
function _integrand(x, constraints)
    n = cumsum(x)
    K = length(n)
    y -> constraints(y) ? prod(y[i]^(n[K - i] - 1)*(1 - y[i])^(x[K - i + 1] - 1) for i in eachindex(y)) : 0.
end
function logpdf_dirichlet_int(x; constraints = y -> true)
    f = _integrand(x, constraints)
    K = length(x)
    res = hcubature(f, SVector{K-1}(zeros(K-1)), SVector{K-1}(ones(K-1)), maxevals = 10^5)
    res[2]/res[1] > .05 && @warn "Relative error in numerical integration = $(res[2]/res[1])"
    log(res[1])
end
logpdf_dirichletmultinomial_int(counts, c; alpha = ones(length(counts))) =
    logpdf_dirichletmultinomial_int(counts, get_constraints(c...), alpha = alpha)
function logpdf_dirichletmultinomial_int(counts, constraints::Function = y -> true;
                                         alpha = ones(length(counts)))
    x = counts .+ alpha
    logpdf_dirichlet_int(x, constraints = constraints) + loggamma(sum(counts) + 1) - sum(loggamma.(counts .+ 1)) - logpdf_dirichlet_int(alpha, constraints = constraints)
end


"""
Akaike information criterion for parameters `mle` with log likelihood `ll`
"""
AIC(n, ll) = 2n - 2ll


#####
##### helper functions
#####

"""
    partialsumsperbird(data, compartement_indep, food_indep, N)

Sums those counts that are assumed to have the same rate and returns additionally
the log probability that the summed counts are produced by a uniform Multinomial,
(see Methods, Eq. 2).


## Example
In the following we compute the sums and the log probability of uniform
distribution under the assumption that 30 decisions are taken in total and the
rates are independent of the compartment, i.e. they are the same for each
compartment, but they are different for each type of food. Data for different
compartments is in the rows, for different foodtypes in the columns.

```julia-repl
julia> partialsumsperbird([1 2;
                    3 4;
                    5 6], true, false, 30)
([9, 9, 12], -7.311519521215035)
```
"""
function partialsumsperbird(data, compartement_indep, food_indep, N)
    if compartement_indep && food_indep
        v = [sum(data)]
        ll = ll_uniform_multinomial(data[:])
    elseif compartement_indep
        v = sum(data, dims = 1)[:]
        ll = +([ll_uniform_multinomial(data[:, i]) for i in 1:length(v)]...)
    elseif food_indep
        v = sum(data, dims = 2)[:]
        ll = +([ll_uniform_multinomial(data[i, :]) for i in 1:length(v)]...)
    else
        v = data[:]
        ll = 0
    end
    if N > 0
        n = N - sum(data)
        @assert n >= 0
        push!(v, n), ll
    else
        v, ll
    end
end

"""
    partialsums(data, bird_indep, compartement_indep, food_indep, N)

Performs the partial sums per bird and, if the rates are assumed to be
independent of the bird identity, performs the partial sums over birds and adds
the correction factors (see Methods, Eq. 3).
"""
function partialsums(data, bird_indep, compartement_indep, food_indep, N)
    countsll = [partialsumsperbird(data[b], compartement_indep, food_indep, N)
                for b in 1:length(data)]
    counts = [c[1] for c in countsll]
    ll = +([c[2] for c in countsll]...)
    if bird_indep
        sum(counts),
        ll + +(map(logmultinomialnormalization, counts)...) -
            logmultinomialnormalization(sum(counts))
    else
        counts, ll
    end
end

"""
    logcategoricalposterior(data, bird_indep, compartement_indep, food_indep, N = 40)

Computes the log probability of the posterior given assumptions about the
independence of the different parameters. The strategy is to sum all counts that
are assumed to have the same rate and add the log probability that the individual
counts are generated with a uniform Multinomial.
"""
function logcategoricalposterior(data; bird_indep, compartement_indep, food_indep, N = 40, group = nothing, constraints = nothing, hyp = constraints)
    counts, ll = partialsums(data, bird_indep, compartement_indep, food_indep, N)
    if hyp == :cch
        if group === nothing
            constraints = ((1, <, 2), (3, <, 2), (4, >, 5), (6, >, 5))
        elseif group == :FN
            constraints = ((1, <, 2), (2, >, 3))
        elseif group == :NF
            constraints = ((1, >, 2), (2, <, 3))
        end
    elseif hyp == :fph2
        if group === nothing
            constraints = ((1, >, 2), (3, >, 2), (4, <, 5), (6, <, 5))
        elseif group == :FN
            constraints = ((1, >, 2), (2, <, 3))
        elseif group == :NF
            constraints = ((2, >, 1), (2, >, 3))
        end
    elseif hyp == :fph1
        if group === nothing
            constraints = ((1, >, 4), (1, >, 2), (1, >, 5), (1, >, 3), (1, >, 6))
        elseif group == :FN
            constraints = ((1, >, 2), (1, >, 3))
        elseif group == :NF
            constraints = ((2, >, 1), (2, >, 3))
        end
    elseif constraints !== nothing
        constraints = constraints
    else
        constraints = y -> true
    end
    if bird_indep
        if hyp === nothing
            logpdf_dirichletmultinomial(counts) + ll
        else
            logpdf_dirichletmultinomial_int(counts, constraints) + ll
        end
    else
        if hyp === nothing
            +(map(logpdf_dirichletmultinomial, counts)...) + ll
        else
            +(map(c -> logpdf_dirichletmultinomial_int(c, constraints), counts)...) + ll
        end
    end
end

"""
    compute_posterior(data, models; choices = [0], group = nothing)

This is the main entry point to compute the posteriors.
"""
function compute_posterior(data, models; choices = [0], group = nothing, unnormalizedlog = false)
    logposteriors = zeros(length(models))
    for k in choices
        logposteriors .+= [logcategoricalposterior(data; N = k,
                                                   group = group, c...) for c in models]
    end
    logposteriors ./= length(choices)                # uniform prior over choices
    unnormalizedlog && return logposteriors
    logposteriors .-= minimum(logposteriors)         # shift for numerical stability
    posteriors = exp.(logposteriors)/sum(exp.(logposteriors))  # exp and normalize
end

"""
    AIC(data, bird_indep, compartement_indep, food_indep, N = 40)

Computes the AIC given assumptions about the independence of the different
parameters.
"""
function AIC(data, bird_indep, compartement_indep, food_indep, N = 40)
    counts, ll = partialsums(data, bird_indep, compartement_indep, food_indep, N)
    if bird_indep
        mle, llmle = MLE(counts)
        AIC(length(mle), llmle + ll), mle
    else
        n = 0
        llmle = 0
        totalmle = []
        for c in counts
            mle, llmle_tmp = MLE(c)
            append!(totalmle, mle)
            n += length(mle)
            llmle += llmle_tmp
        end
        AIC(n, llmle + ll), totalmle
    end
end

"""
Returns an array where the numbers in row 1 correspond to the cached items in
the the first compartment etc. and numbers in column 1 corresponds to the
cached items of the first type the bird was prefed with etc.
"""
function sort_semantically(data)
    comp = copy(data.comp_order)
    prefed = copy(data.prefed_order)
    res = copy(data.cached)
    while comp[1] != :A
        comp = circshift(comp, 1)
        if size(res, 2) > 1
            res = circshift(res, (-1, 0))
        else
            res = circshift(res, -1)
        end
    end
    while prefed[1] != :P && size(res, 2) > 1
        prefed = circshift(prefed, 1)
        res = circshift(res, (0, 1))
    end
    (comp_order = comp, prefed_order = prefed, cached = res)
end

function load_data(file_name)
    data = YAML.load_file(file_name)
    result = Dict()
    for (key, content) in data
        cached = content["cached"]
        # Data is in row first form in the yaml file.
        length(cached) > 3 && (cached = Array(reshape(cached, :, 3)'))
        result[Symbol(key)] = (prefed_order = Symbol.(content["prefed_order"]),
                               comp_order = Symbol.(content["comp_order"]),
                               cached = cached)
    end
    data_sorted = Dict(name => sort_semantically(content) for (name, content) in result)
    (data = result,
     data_sorted = data_sorted,
     cached = [b.cached for b in values(data_sorted)])
end

