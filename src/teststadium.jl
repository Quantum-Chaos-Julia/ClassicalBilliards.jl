using Revise
using ClassicalBilliards
#include("plottingmakie.jl")
using GLMakie
using StaticArrays
using LinearAlgebra
using Roots
using BenchmarkTools


stad = Stadium(0.5)

#trajectory(p,stad,100)


#scatter!(ax, p.r[1], p.r[2])
#arrows!(ax, [p.r[1]], [p.r[2]], [p.v[1]], [p.v[2]]; lengthscale = 0.05)
p = Particle(1.0,0.6,-1.0,1.2;domain_idx=1)
pts,vel,ts = trajectory(p,stad,100; dt = 0.1)

pts

f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_billiard!(ax,stad)
lines!(ax, getindex.(pts,1), getindex.(pts,2))
scatter!(ax, getindex.(pts,1), getindex.(pts,2), color=:black)
arrows!(ax, getindex.(pts,1), getindex.(pts,2), getindex.(vel,1), getindex.(vel,2);color=:black, lengthscale = 0.05)

display(f)
p.time

p = Particle(1.0,0.6,-1.0,1.2;domain_idx=1)
@btime trajectory(p,stad,100000; dt = 0.01)
time

f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_billiard!(ax,stad;dens = 100.0)
scatter!(ax, p.r[1], p.r[2])
arrows!(ax, [p.r[1]], [p.r[2]], [p.v[1]], [p.v[2]]; lengthscale = 0.05)
display(f)


#p = Particle(1.0,0.6,-1.0,1.2;domain_idx=1)
#pts1,vel1,ts = trajectory(p,stad,100000; dt = 0.01)

#ts

#append!(pts,pts1)