using Revise
using ClassicalBilliards
#include("plottingmakie.jl")
using GLMakie
using StaticArrays
using LinearAlgebra
using Roots
using BenchmarkTools

using LaTeXStrings
using MathTeXEngine # required for texfont

textheme = Theme(fontsize = 8,
            fonts=(; regular=texfont(:text),
                        bold=texfont(:bold),
                        italic=texfont(:italic),
                        bold_italic=texfont(:bolditalic)),
            Axis = (;xgridvisible = false,
            ygridvisible = false, yticksize=3,xticksize=3,yticklabelpad=2.0,xticklabelpad=2.0, xlabelpadding = 0.0),
            #lines = (;grid=false)
                        
            )
set_theme!(textheme)
#=
function check_subdomain!(particle,billiard)
    idx = 1
    for dom in billiard.subdomains
        if all(is_inside(dom, particle.r))
            particle.subdomain = idx
        end
    idx +=1
    end
end
=#
#check_subdomain

mushroom = Mushroom(0.75)

f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_billiard!(ax,mushroom; plot_virtual = true)
display(f)
#trajectory(p,stad,100)


#scatter!(ax, p.r[1], p.r[2])
#arrows!(ax, [p.r[1]], [p.r[2]], [p.v[1]], [p.v[2]]; lengthscale = 0.05)
p = Particle(0.5,0.5,-1.0,-1.3;subdomain=1)

T = 10000
#pts, vel, ts = trajectory(p,mushroom,T; dt = 0.01)

f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_trajectory!(ax, p, mushroom, T; dt = 0.01,traj_args = Dict(:alpha=>0.25), full_domain=true, plot_virtual = false, plot_velocity=false)
display(f)


f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_trajectory!(ax, p, mushroom, T; dt = 0.01,traj_args = Dict(:alpha=>0.25), full_domain=true, plot_virtual = true, plot_velocity=true)
display(f)


