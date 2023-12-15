include("linesegment.jl")
include("circlesegment.jl")


struct Domain{T} <: AbsDomain where T<:Real
    boundary::Vector{AbsCurve}
    id::Int64
end

function is_inside(domain::D, pt) where {D<:AbsDomain}
    d = [is_inside(crv, pt) for crv in domain.boundary]
    return d
end

function is_inside(domain::D, pts::AbstractArray) where {D<:AbsDomain}
    d = [is_inside(crv, pts) for crv in domain.boundary]
    return  reduce(hcat,d)
end


#check if points inside for general curves
function is_inside(curve::C, pt) where {C<:AbsCurve}
    return domain_fun(curve, pt) .< zero(eltype(pt)) 
end

function is_inside(curve::C, pts::AbstractArray) where {C<:AbsCurve}
    let
    d = domain_fun(curve, pts)
    return d .< zero(eltype(pts[1])) 
    end
end

struct Billiard{T} <: AbsBilliard where T<:Real
    domains::Vector{AbsDomain}
    #symmetries::
end
