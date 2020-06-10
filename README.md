# Testing two competing hypotheses for Eurasian jays' caching for the future: planning versus compensatory caching

This repository contains the data and the code of the data analysis of

Piero Amodio, Johanni Brea,Benjamin Farrar, Ljerka Ostojic and Nicola S. Clayton, (2020),
"Testing two competing hypotheses for Eurasian jays' caching for the future: planning versus compensatory caching".

To run the data analysis open a [julia](https://julialang.org) REPL and run
```julia
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/jbrea/PlanningVsCompensatoryCaching.jl"))
Pkg.activate("PlanningVsCompensatoryCaching")
import PlanningVsCompensatoryCaching
include(joinpath(dirname(pathof(PlanningVsCompensatoryCaching))))
```
