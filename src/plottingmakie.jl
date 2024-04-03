#include("../abstracttypes.jl")
#include("../utils/gridutils.jl")

#helper functions
export plot_domain!, plot_domain_fun!, plot_curve!, plot_billiard!
#curve and billiard ploting
function plot_billiard!(ax, billiard::AbsBilliard; dens = 20.0)
    for dom in billiard.domains
        for crv in dom.boundary
            plot_curve!(ax, crv; dens)
        end
    end
end

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
    return ax, hmap
end


function plot_domain!(f, billiard::B; xlim=(-1.0,1.0),ylim=(-1.0,1.0), dens=100.0, hmargs=Dict(),cmap=:binary) where {B<:AbsBilliard}
    d = one(dens)/dens
    x_grid = range(xlim... ; step=d)
    y_grid = range(ylim... ; step=d)
    pts = [SVector(x,y) for y in y_grid for x in x_grid]
    #Z = is_inside(domain,pts)
    domains = [[all(r) for r in eachrow(is_inside(domain,pts))] for domain in billiard.domains]
    for dom in domains
        Z = reshape(dom,length(x_grid),length(y_grid))
        hmap, ax = plot_heatmap!(f, x_grid, y_grid, Z, cmap = cmap, vmax=1.0 ,hmargs...)
    end
    ax.aspect=DataAspect()
    return ax, hmap
end
