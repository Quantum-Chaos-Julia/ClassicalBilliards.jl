include("ClassicalBilliards.jl")
include("plottingmakie.jl")
using CairoMakie

#curve and billiard ploting
function plot_curve!(ax, crv::AbsCurve; dens = 20.0)
    L = crv.length
    grid = max(round(Int, L*dens),3)
    t = range(0.0,1.0, grid)
    pts = curve(crv,t)
    lines!(ax, pts, color = :black )
    ax.aspect=DataAspect()
end

function plot_domain_fun!(f, curve::C; xlim=(-1.0,1.0),ylim=(-1.0,1.0), dens=100.0, hmargs=Dict(),cmap=:binary) where {C<:AbsCurve}
    d = one(dens)/dens
    x_grid = range(xlim... ; step=d)
    y_grid = range(ylim... ; step=d)
    pts = [SVector(x,y) for y in y_grid for x in x_grid]
    Z = reshape(domain_fun(curve,pts),length(x_grid),length(y_grid))
    hmap, ax = plot_heatmap_balaced!(f,x_grid,y_grid,Z) 
    ax.aspect=DataAspect()
    return ax, hmap
end


function plot_domain!(f, domain::D; xlim=(-1.0,1.0),ylim=(-1.0,1.0), dens=100.0, hmargs=Dict(),cmap=:binary) where {D<:AbsDomain}
    d = one(dens)/dens
    x_grid = range(xlim... ; step=d)
    y_grid = range(ylim... ; step=d)
    pts = [SVector(x,y) for y in y_grid for x in x_grid]
    #Z = is_inside(domain,pts)
    Z = reshape([all(r) for r in eachrow(is_inside(domain,pts))],length(x_grid),length(y_grid))
    hmap, ax = plot_heatmap!(f, x_grid, y_grid, Z, cmap = cmap, vmax=1.0 ,hmargs...)
    ax.aspect=DataAspect()
end


function plot_domain!(f, billiard::D; xlim=(-1.0,1.0),ylim=(-1.0,1.0), dens=100.0, hmargs=Dict(),cmap=:binary) where {B<:AbsClassicalBilliard}
    d = one(dens)/dens
    x_grid = range(xlim... ; step=d)
    y_grid = range(ylim... ; step=d)
    pts = [SVector(x,y) for y in y_grid for x in x_grid]
    #Z = is_inside(domain,pts)
    domains = [[all(r) for r in eachrow(is_inside(domain,pts))] for domain in billiard.domains]
    for dom in domains
        Z = reshape(dom,length(x_grid),length(y_grid))
        hmap, ax = plot_heatmap!(f, x_grid, y_grid, Z, cmap = cmap, vmax=1.0 ,hmargs...)
    ax.aspect=DataAspect()
end


half_width = 0.5
circle = CircleSegment(1.0, pi/2, 0.0, half_width, 0.0)
x_segment = SymLineSegment(SVector(half_width,0.0),SVector(half_width+1.0,0.0))
y_segment = VirtLineSegment(SVector(half_width, 1.0), SVector(half_width,0.0))




f = Figure(resolution=(1000,1000))
ax = Axis(f[1,1])
plot_curve!(ax, circle)
plot_curve!(ax, x_segment)
plot_curve!(ax, y_segment)
display(f)

f = Figure(resolution=(1000,1000))
plot_domain_fun!(f, y_segment)
display(f)

dom1 =  = Domain{Float64}([circle,x_segment,y_segment],1)


f = Figure(resolution=(1000,1000))
plot_domain!(f, stadium)
display(f)

