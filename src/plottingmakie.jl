#plotting_parameters
default_linewidth = 0.75
dens = 25.0
boundary_args = Dict(:color=>:black, :linewidth=>1.5, :linestyle=>:solid)
virt_args = Dict(:color=>:grey, :linewidth=>default_linewidth, :linestyle=>:dot)
sym_args = Dict(:color=>:black, :linewidth=>default_linewidth, :linestyle=>:dash)
traj_args = Dict(:linewidth=>0.75, :linestyle=>:solid, :transparency => true, :alpha => 0.1)

#helper functions
export plot_domain!, plot_domain_fun!, plot_curve!, plot_billiard!, plot_trajectory!
#curve and billiard ploting
function plot_billiard!(ax, billiard::B; plot_full_domain = true, plot_virtual = true, dens = dens)  where B<:AbsBilliard
    for dom in billiard.subdomains
        for crv in dom.boundary
            if typeof(crv) <: AbsVirtLine
                if plot_virtual
                    plot_curve!(ax, crv; dens)
                end
            else
                if plot_full_domain
                    L = crv.length
                    grid = max(round(Int, L*dens),3)
                    t = range(0.0,1.0, grid)
                    pts = curve(crv,t)
                    lines!(ax, pts; boundary_args... )
                    for sym in billiard.symmetries[2:end]
                        sym_pts = [sym(pt) for pt in pts]
                        lines!(ax, sym_pts; boundary_args... )
                    end 
                else
                    plot_curve!(ax, crv; dens)
                end
            end
        end
    end
    ax.aspect=DataAspect()
end

function plot_domain!(ax, domain::D; dens = dens)  where D<:AbsDomain
    for crv in domain.boundary
        plot_curve!(ax, crv; dens)
    end
    ax.aspect=DataAspect()
end

function plot_curve!(ax, crv::AbsCurve;  dens = dens)
    L = crv.length
    grid = max(round(Int, L*dens),3)
    t = range(0.0,1.0, grid)
    pts = curve(crv,t)
    lines!(ax, pts; boundary_args... )
    ax.aspect=DataAspect()
end

function plot_curve!(ax, crv::C; dens = dens)  where C<:VirtLineSegment
    L = crv.length
    grid = max(round(Int, L*dens),3)
    t = range(0.0,1.0, grid)
    pts = curve(crv,t)
    lines!(ax, pts; virt_args... )
    ax.aspect=DataAspect()
end

function plot_curve!(ax, crv::C; dens = dens)  where C<:SymLineSegment
    L = crv.length
    grid = max(round(Int, L*dens),3)
    t = range(0.0,1.0, grid)
    pts = curve(crv,t)
    lines!(ax, pts; sym_args... )
    ax.aspect=DataAspect()
end

#trajectory plotting
function plot_trajectory!(ax, particle, billiard, T; dt = 0.01, full_domain=true, traj_args = traj_args, plot_virtual = false, plot_velocity = false)
    pts, vel, ts = trajectory(particle, billiard, T; dt, full_domain)
    plot_billiard!(ax,billiard; plot_virtual, plot_full_domain = full_domain)
    lines!(ax, getindex.(pts,1), getindex.(pts,2); traj_args...)
    if plot_velocity
        arrows!(ax, getindex.(pts,1), getindex.(pts,2), getindex.(vel,1), getindex.(vel,2); lengthscale = 0.05)
    end
end

#domain plotting

function plot_domain_fun!(f, curve::C; xlim=(-1.0,1.0),ylim=(-1.0,1.0), dens=100.0, hmargs=Dict(),cmap=:binary) where {C<:AbsCurve}
    d = one(dens)/dens
    x_grid = range(xlim... ; step=d)
    y_grid = range(ylim... ; step=d)
    pts = [SVector(x,y) for y in y_grid for x in x_grid]
    Z = reshape(domain_fun(curve,pts),length(x_grid),length(y_grid))
    hmap, ax = plot_heatmap_balaced!(f,x_grid,y_grid,Z; hmargs...) 
    #ax.aspect=DataAspect()
    #return hmap, ax
end
#=
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
    subdomains = [[all(r) for r in eachrow(is_inside(domain,pts))] for domain in billiard.subdomains]
    for dom in subdomains
        Z = reshape(dom,length(x_grid),length(y_grid))
        hmap, ax = plot_heatmap!(f, x_grid, y_grid, Z, cmap = cmap, vmax=1.0 ,hmargs...)
    end
    ax.aspect=DataAspect()
    return ax, hmap
end
=#
