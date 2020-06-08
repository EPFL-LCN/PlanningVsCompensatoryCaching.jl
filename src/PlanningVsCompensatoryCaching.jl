module PlanningVsCompensatoryCaching
using SpecialFunctions, PGFPlotsX, LinearAlgebra, Distributions, DataStructures, YAML, HCubature, StaticArrays
push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usetikzlibrary{patterns}")

include("plotting.jl")
include("utils.jl")

end # module
