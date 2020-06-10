# Testing two competing hypotheses for Eurasian jays' caching for the future: planning versus compensatory caching

This repository contains the data and the code of the data analysis of

Piero Amodio, Johanni Brea,Benjamin Farrar, Ljerka Ostojic and Nicola S. Clayton, (2020),
"Testing two competing hypotheses for Eurasian jays' caching for the future: planning versus compensatory caching".

To run the data analysis open a [julia](https://julialang.org) REPL and run
```julia
# load the code
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/jbrea/PlanningVsCompensatoryCaching.jl"))
import PlanningVsCompensatoryCaching
# set-up results directory
RESULT_DIR = "results"
mkdir(RESULT_DIR)
# run the script
DIR = joinpath(dirname(pathof(PlanningVsCompensatoryCaching)), "..")
Pkg.activate(DIR)
include(joinpath(DIR, "script.jl"))
```
