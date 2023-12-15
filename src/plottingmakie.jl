#include("../abstracttypes.jl")
#include("../utils/gridutils.jl")
using Makie
using StaticArrays
#helper functions
function plot_heatmap!(f,x,y,Z ;vmax = 1.0,log=(false,-5.0), cmap=Reverse(:gist_heat),hmargs=Dict(),axargs=Dict())
    if log[1]
        X = log10.(Z)
        ax = Axis(f[1,1],axargs...)        
        m = findmax(X)[1]
        range_val = (log[2],m*vmax)
        hmap = heatmap!(ax,x, y, X, colormap = cmap, colorrange=range_val, hmargs...)
        ax.aspect=DataAspect()
        Colorbar(f[1,2], colormap = cmap, limits = Float64.(range_val),tellheight=true)
        rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
    else
        ax = Axis(f[1,1],axargs...)        
        m = findmax(Z)[1]
        range_val = (0,m*vmax)
        hmap = heatmap!(ax,x, y, Z, colormap = cmap, colorrange=range_val, hmargs...)
        ax.aspect=DataAspect()
        Colorbar(f[1,2], colormap = cmap, limits = Float64.(range_val),tellheight=true)
        rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
    end
    return hmap, ax
end

function plot_heatmap_balaced!(f,x,y,Z ;vmax = 1.0, cmap=Reverse(:balance),hmargs=Dict(),axargs=Dict())
    ax = Axis(f[1,1],axargs...)        
    m = findmax(abs, Z)[1]
    range_val = (-m*vmax,m*vmax)
    hmap = heatmap!(ax,x, y, Z, colormap = cmap, colorrange=range_val, hmargs...)
    ax.aspect=DataAspect()
    Colorbar(f[1,2], colormap = cmap, limits = Float64.(range_val),tellheight=true)
    rowsize!(f.layout, 1, ax.scene.px_area[].widths[2])
    return hmap, ax
end
