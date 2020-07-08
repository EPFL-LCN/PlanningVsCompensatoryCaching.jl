import PlanningVsCompensatoryCaching: load_data, sort_semantically,
datatotextable, plotallindividualdata, logcategoricalposterior, plotcomparison,
pgfsave, compute_posterior

@isdefined(RESULT_DIR) || (RESULT_DIR = joinpath(@__DIR__, "results/"))

@info "Analysing Experiment 1..."

data1, data1_sorted, cached1 = load_data(joinpath(@__DIR__, "data", "experiment1.yaml"))

# Export raw and sorted data.
write(joinpath(RESULT_DIR, "datatable.tex"), datatotextable(data1, data1_sorted))

# Define qualitatively some hypotheses and plot them together with the data.
data_sorted_with_hypothesis = deepcopy(data1_sorted)
data_sorted_with_hypothesis[:FPH1] = (comp_order = [:A, :B, :C],
                                      prefed_order = [:P, :M],
                                      cached = [5 1; 1 1; 1 1])
data_sorted_with_hypothesis[:FPH2] = (comp_order = [:A, :B, :C],
                                      prefed_order = [:P, :M],
                                      cached = [3 1; 1 3; 3 2])
data_sorted_with_hypothesis[:CCH] = (comp_order = [:A, :B, :C],
                                     prefed_order = [:P, :M],
                                     cached = [2 3; 4 1; 1 3])
plotallindividualdata(data_sorted_with_hypothesis, dir = RESULT_DIR,
                      xlabelat = "FPH1", titles = Dict(:FPH1 => "FPH 1",
                                                       :FPH2 => "FPH 2"));

models = [[(bird_indep = b, compartement_indep = c, food_indep = f, hyp = nothing)
           for b in (true, false), f in (true, false), c in (true, false) if !(b && c && f)]...;
          [(bird_indep = b, compartement_indep = false, food_indep = false, hyp = h)
           for b in (true, false), h in (:cch, :fph1, :fph2)]...]

# Compute Posterior for each condition.
posteriors = compute_posterior(cached1, models)
@show posteriors

# Generate comparison plots.
p1 = plotcomparison(posteriors, models, title = "\$P(M|D)\$", ymax = 5);
pgfsave(joinpath(RESULT_DIR, "posteriors.tex"), p1, include_preamble = false)


"""
Conclusion

1 bird clear FPH 1 (Lisbon)
2 birds unclear (Jerusalem, Quito)
3 birds rather CCH (Caracas, Washington, Wellington)

BUT: I don't want to overinterpret these results since the data is quite well
explained with a model that distributes uniformly among compartments
"""


#####
##### Data Analysis: Experiment 2
#####

import PlanningVsCompensatoryCaching: load_data, sort_semantically,
datatotextable2, plotallindividualdata, logcategoricalposterior, plotcomparison,
pgfsave, compute_posterior

@info "Analysing Experiment 2..."

data2, data2_sorted, cached2 = load_data(joinpath(@__DIR__, "data", "experiment2.yaml"))

# Export raw and sorted data.
write(joinpath(RESULT_DIR, "datatable2.tex"), datatotextable2(data2, data2_sorted))

# Define qualitatively some hypotheses and plot them together with the data.
data_sorted_with_hypothesis = deepcopy(data2_sorted)
data_sorted_with_hypothesis[:FPH1] = (comp_order = [:A, :B, :C],
                                      prefed_order = [:F, :N],
                                      cached = [4; 1; 1])
data_sorted_with_hypothesis[:FPH2] = (comp_order = [:A, :B, :C],
                                      prefed_order = [:F, :N],
                                      cached = [2; 1; 2])
data_sorted_with_hypothesis[:CCH] = (comp_order = [:A, :B, :C],
                                     prefed_order = [:F, :N],
                                     cached = [1; 2; 1])
data_sorted_with_hypothesis[:FPH1NF] = (comp_order = [:A, :B, :C],
                                        prefed_order = [:N, :F],
                                        cached = [1; 3; 1])
data_sorted_with_hypothesis[:FPH2NF] = (comp_order = [:A, :B, :C],
                                        prefed_order = [:N, :F],
                                        cached = [1; 3; 1])
data_sorted_with_hypothesis[:CCHNF] = (comp_order = [:A, :B, :C],
                                       prefed_order = [:N, :F],
                                       cached = [2; 1; 2])
plotallindividualdata(data_sorted_with_hypothesis, exp = 2, xlabelat = "FPH1",
                      y_axis_line_style = "{->}", dir = RESULT_DIR, ylabel = "cached",
                      titles = Dict(:FPH1 => "FPH 1", :FPH2 => "FPH 2",
                                    :CCHNF => "CCH", :FPH1NF => "FPH 1",
                                    :FPH2NF => "FPH 2"));

models = [[(bird_indep = b, compartement_indep = c, food_indep = true, hyp = nothing)
           for b in (true, false), c in (true, false) if !(b && c)]...;
          [(bird_indep = b, compartement_indep = false, food_indep = true, hyp = h)
           for b in (true, false), h in (:cch, :fph1, :fph2)]...]

# Compute Posterior for each condition and 40 caching trials.
idx_fn = findall(x -> x ∈ (:Caracas, :Washington, :Quito), collect(keys(data2)))
idx_nf = findall(x -> x ∈ (:Wellington, :Rome, :Lisbon), collect(keys(data2)))
lp_fn = compute_posterior(cached2[idx_fn], models, group = :FN, unnormalizedlog = true)
lp_nf = compute_posterior(cached2[idx_nf], models, group = :NF, unnormalizedlog = true)
posteriors = (x -> x ./ sum(x))(exp.(lp_fn .+ lp_nf))

@show posteriors

p1 = plotcomparison(posteriors, models, title = "\$P(M|D)\$", altpattern = "north east lines", ymax = 5, exp2 = true, legendtoname = "modelcomparisonlegend2", legendentries = ["bird independent", "bird dependent"]);
pgfsave(joinpath(RESULT_DIR, "posteriors2.tex"), p1, include_preamble = false)

"""
Conclusion

3 birds rather FPH 1 (Caracas, Washington, Quito)
3 birds rather CCH (Rome, Lisbon Wellington)

BUT: I don't want to overinterpret these results since the data is quite well
explained with a model that distributes uniformly among compartments
"""

