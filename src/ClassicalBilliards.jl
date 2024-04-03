module ClassicalBilliards

using StaticArrays 
using LinearAlgebra
using CoordinateTransformations, Rotations
using ForwardDiff
using Roots
using Makie
using SavingPlotting

abstract type AbsCurve end
abstract type AbsDomain end
abstract type AbsBilliard end
abstract type AbsBoundaryCondition end

export AbsCurve, AbsDomain, AbsBilliard, AbsBoundaryCondition

include("geometry/geometry.jl")
include("dynamics/dynamics.jl")

include("plottingmakie.jl")

end