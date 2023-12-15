using StaticArrays 
using LinearAlgebra
using CoordinateTransformations, Rotations

abstract type AbsCurve end
abstract type AbsDomain end
abstract type AbsBilliard end
abstract type AbsBoundaryCondition end

include("geometry.jl")


