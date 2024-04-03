using Revise
using ClassicalBilliards
#include("plottingmakie.jl")
using CairoMakie
using StaticArrays

using Roots
using BenchmarkTools



f = Figure(resolution=(1000,1000))
ax = Axis(f[1,1])
plot_curve!(ax, circle)
plot_curve!(ax, x_segment)
plot_curve!(ax, y_segment)
display(f)

f = Figure(resolution=(1000,1000))
plot_domain_fun!(f, y_segment)
display(f)


f = Figure(resolution=(1000,1000))
ax = Axis(f[1,1])
for crv in rectangle_dom.boundary
    plot_curve!(ax, crv)
end
display(f)

f = Figure(resolution=(1000,1000))
plot_domain!(f, circle_dom; xlim=(-1.0,2.0))
plot_domain!(f, rectangle_dom; xlim=(-1.0,2.0))
display(f)



p = Particle(1.2,0.6,-1.0,1.2)
pts = [p.r]
ns = []
f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
for i in 1:1000
    iterate_bounce!(p, circle_dom; dt = 1.0)
    #println(p.r)
    push!(pts, p.r)
end

plot_curve!(ax, circle)
plot_curve!(ax, x_segment)
plot_curve!(ax, y_segment)
lines!(ax, getindex.(pts,1), getindex.(pts,2) )
#arrows!(ax,getindex.(pts,1),getindex.(pts,2), getindex.(ns,1),getindex.(ns,2), color = :black, lengthscale = 0.1)
display(f)



f = Figure(resolution=(1000,1000))
ax, hmap = plot_domain_fun!(f, circle)

pts = curve(circle, range(0.0,1.0,100))
g = domain_gradient_vector(circle, pts)
arrows!(ax,getindex.(pts,1),getindex.(pts,2), getindex.(g,1),getindex.(g,2), color = :black, lengthscale = 0.1)
display(f)


#origin = SVector{2}([0.0,0.0])
normal = SVector{2}([0.4,1.0])
n = -normal./norm(normal)
incident = SVector{2}([1.0,-0.5])
out = incident - 2*dot(incident,n)*n 
2*dot(incident,n)*n 


f = Figure(resolution=(1000,1000))
ax = Axis(f[1,1])
arrows!(ax,[0.0],[0.0],[incident[1]],[incident[2]], color = :blue, lengthscale = 1.0)
arrows!(ax,[0.0],[0.0],[n[1]],[n[2]], color = :black, lengthscale = 1.0)
arrows!(ax,[0.0],[0.0],[-n[2]],[n[1]], color = :black, lengthscale = 1.0)
arrows!(ax,[0.0],[0.0],[out[1]],[out[2]], color = :red, lengthscale = 1.0)
#arrows!(ax,getindex.(data1,1),getindex.(data1,2), getindex.(data2,1),getindex.(data2,2), color = :black, lengthscale = 0.1)

display(f)


exit_time, idx = find_intersection(p, circle_dom; dt = 0.01)

n = domain_gradient_vector(circle_dom.boundary[idx], propagate(p,exit_time))

########################333
p = Particle(0.0,0.0,1.0,1.0)
ts = range(0.0,1.0,200)

pts  =propagate(p, ts)


f = Figure(resolution=(1000,1000))
ax, hmap = plot_domain_fun!(f, y_segment)
lines!(ax, getindex.(pts,1), getindex.(pts,2))
display(f)

fun = domain_fun(circle,pts)
f = Figure(resolution=(1000,1000))
ax = Axis(f[1,1])
lines!(ax, ts, fun)
r(t) = domain_fun(circle, propagate(p, t))
z = find_zero(t -> domain_fun(circle,propagate(p, t)), (0.0, 1.0))
scatter!(ax,z,r(z))
display(f)

p = Particle(1.0,0.5,-1.0,1.1)
ts = range(0.0,1.0,200)
pts  = propagate(p, ts)
f = Figure(resolution=(1000,1000))
ax, hmap = plot_domain!(f, circle_dom; xlim=(-1.0,2.0))
lines!(ax, getindex.(pts,1), getindex.(pts,2))
#plot_domain!(f, rectangle_dom; xlim=(-1.0,2.0))
display(f)

d, tim = find_domain_exit_curves(p, circle_dom, 0.1)

using LinearAlgebra
function collision!(ax, particle::P, curve::C, time) where {P<:Particle, C<:AbsCurve}
    let v = particle.v #particle initial velocity
    r = propagate(particle, time) #colision point
    
    g = domain_gradient_vector(curve, r)
    n =  g./norm(g)
    
    v_new = v .- 2*dot(v,n)*n 
    arrows!(ax,[r[1]],[r[2]],[-n[1]],[-n[2]], color = :black, lengthscale = 0.1)
    arrows!(ax,[r[1]],[r[2]],[v[1]],[v[2]], color = :blue, lengthscale = 0.1)
    arrows!(ax,[r[1]],[r[2]],[v_new[1]],[v_new[2]], color = :red, lengthscale = 0.1)
    println("collision")
    println(particle.v)
    particle.v = v_new
    println(particle.v)
    particle.r = r
    particle.time += time
    end
end

function iterate_bounce!(ax, particle::P, domain::D; dt = 1.0 ) where {P<:AbsParticle, D<:AbsDomain}
    collision_time, idx = find_intersection(particle, domain; dt)
    #println(collision_time)
    println(domain.boundary[idx])
    collision!(ax, particle, domain.boundary[idx], collision_time)
end

p = Particle(1.0,0.6,-1.0,1.2)
pts = [p.r]
ns = []
f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
for i in 1:4
    println(i)
    iterate_bounce!(ax, p, circle_dom; dt = 1.0)
    #println(p.r)
    push!(pts, p.r)
end

plot_curve!(ax, circle)
plot_curve!(ax, x_segment)
plot_curve!(ax, y_segment)
lines!(ax, getindex.(pts,1), getindex.(pts,2) )
#arrows!(ax,getindex.(pts,1),getindex.(pts,2), getindex.(ns,1),getindex.(ns,2), color = :black, lengthscale = 0.1)
display(f)