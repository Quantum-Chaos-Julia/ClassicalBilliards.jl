using Revise
using ClassicalBilliards
using BenchmarkTools
using GLMakie

limacon = Limacon(1.0)

f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_billiard!(ax,limacon; plot_virtual = true)
display(f)
#trajectory(p,stad,100)


#scatter!(ax, p.r[1], p.r[2])
#arrows!(ax, [p.r[1]], [p.r[2]], [p.v[1]], [p.v[2]]; lengthscale = 0.05)
p = Particle(0.5,0.5,1.1,1.2;subdomain=1)
T = 10000

f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_trajectory!(ax, p, limacon, T; dt = 0.0001, traj_args = Dict(:alpha=>0.5), full_domain=false, plot_virtual = true)
display(f)


p = Particle(1.0,0.6,-1.0,1.2;subdomain=1)
@btime trajectory(p,limacon,100000; dt = 0.001)


f = Figure(size=(1000,1000))
plot_domain_fun!(f, limacon.subdomains[1].boundary[1]; xlim=(-1.0,1.0),ylim=(-1.0,1.0), dens=100.0, hmargs=Dict(),cmap=:binary) where {C<:AbsCurve}


using StaticArrays
f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
curve = limacon.subdomains[1].boundary[1]
xlim=(-1.5,1.5)
ylim=(-1.5,1.5)
dens=1000.0
d = one(dens)/dens
x_grid = range(xlim... ; step=d)
y_grid = range(ylim... ; step=d)
pts = [SVector(x,y) for y in y_grid for x in x_grid]

Z = reshape(domain_fun(curve,pts),length(x_grid),length(y_grid))
m = findmax(Z)[1]
vmax = 1.0
range_val = (-m*vmax,m*vmax)
hmap = heatmap!(ax,x_grid, y_grid, Z; colormap = Reverse(:balance),colorrange=range_val)
plot_billiard!(ax,limacon; plot_virtual = true)
display(f)