abstract type AbsParticle end

mutable struct Particle{T}  <: AbsParticle where T<:Real
    r::SVector{2,T}
    v::SVector{2,T}
    m::T #mass of particle
    time::T #current time of the particle
    subdomain::Int64 #id of current subdomain
    sym_sector::Int64 #id of current symmetry sector
    v_b::SVector{2,T} # velocity vector in local coordinate system on the boundary (v_t, v_n) v_t is tangential and v_n normal component.
end

Particle(x,y,vx,vy; m = 1.0, time = 0.0, subdomain=1, sym_sector=1, v_b = SVector(vx,vy)) = Particle{typeof(x)}(SVector(x,y),SVector(vx,vy),m,time,subdomain,sym_sector, v_b)
Particle(point,velocity; m = 1.0, time = 0.0, subdomain=1, sym_sector=1, v_b = SVector(vx,vy)) = Particle{eltype(point)}(point,velocity,m,time,subdomain,sym_sector, v_b)

export Particle, propagate

function propagate(part::P, t) where P<:Particle
    return @. part.r + part.v*t
end

function propagate(part::P, ts::AbstractArray) where P<:Particle
    return [@. part.r .+ part.v.*t for t in ts] 
end


function find_domain_exit_curves(particle::P, domain::D, dt) where {P<:AbsParticle, D<:AbsDomain}
    local_time = 0.0 #time since last bounce
    pt = propagate(particle, local_time + dt) #particle position after dt
    local_time += dt #increase time by dt
    d = is_inside(domain, pt) #check if particle is inside all curves
    while all(d) #while inside all curves
        pt = propagate(particle, local_time + dt)
        d = is_inside(domain, pt)
        local_time += dt
    end
    return d, local_time
end


function find_intersection_time(particle::P, curve::C, exit_time) where {P<:AbsParticle, C<:AbsCurve}
    #time since last bounce
    let
        r(t) = domain_fun(curve, propagate(particle, t))
        t = find_zeros(r, (0.0+1e-14, exit_time)) #10.0*eps(exit_time)
        return t[1]  #take shortest time
    end
end
#=
function find_intersection_time(particle::P, line::C, exit_time) where {P<:Particle, C<:AbsLine}
    #time since last bounce
    let
        r(t) = domain_fun(curve, propagate(particle, t))
        t = find_zeros(r, (0.0+1e-14, exit_time)) #10.0*eps(exit_time)
        return t[1]  #take shortest time
    end
end
=#

function find_intersection(particle::P, domain::D; dt = 1.0) where {P<:AbsParticle, D<:AbsDomain}
    #time since last bounce
    boundary = domain.boundary
    N =  length(boundary)
    inside, approx_exit_time = find_domain_exit_curves(particle, domain, dt) #d gives the curves that were exited
    
    exit_curves = boundary[.~inside]
    crv_idx = collect(1:N)[.~inside]
    idx = crv_idx[1]
    exit_time = approx_exit_time
    for (i,crv) in enumerate(exit_curves)
        t = find_intersection_time(particle, crv, approx_exit_time)
        if t < exit_time
            exit_time = t
            idx = crv_idx[i]
        end    
    end 
    return exit_time, idx
end

function collision_rule!(particle::P, curve::C, n) where {P<:AbsParticle, C<:AbsCurve}
    let v = particle.v #particle initial velocity
        v_n = dot(v,n)*n #normal velocity
        v_t = v .- v_n #tangetial velocity
        v_new = v .- 2*v_n 
        particle.v = v_new
        particle.v_b = SVector(norm(v_t),norm(v_n))
    end
end

function collision_rule!(particle::P, sym_line::C, n) where {P<:AbsParticle, C<:SymLineSegment}
    #collisions with symetry line result in change of symmetry sector
    let v = particle.v , id = particle.sym_sector, sym = sym_line.symmetry #particle initial velocity
        v_n = dot(v,n)*n #normal velocity
        v_t = v .- v_n #tangetial velocity
        v_new = v .- 2*v_n 
        particle.v = v_new
        particle.v_b = SVector(norm(v_t),norm(v_n))
        new_id = new_sector(sym, id) #change symmetry sector
        particle.sym_sector = new_id
    end
end

function collision_rule!(particle::P, virt_line::C, n) where {P<:AbsParticle, C<:VirtLineSegment}
    let v = particle.v
        v_n = dot(v,n)*n #normal velocity
        v_t = v .- v_n #tangetial velocity
        particle.v_b = SVector(norm(v_t),norm(v_n))
        particle.subdomain = virt_line.domain_exit_idx
    end
end

function collision!(particle::P, curve::C, time) where {P<:Particle, C<:AbsCurve}
    r = propagate(particle, time) #colision point
    g = domain_gradient_vector(curve, r)
    n =  g./norm(g) 
    collision_rule!(particle, curve, n) #modify velocity using rule corresponding to boundary
    particle.r = r
    particle.time += time
    
end


function iterate_bounce!(particle::P, billiard::B; dt = 1.0 ) where {P<:AbsParticle, B<:AbsBilliard}
    domain = billiard.subdomains[particle.subdomain]
        
    collision_time, idx = find_intersection(particle, domain; dt)
    crv = domain.boundary[idx]
    #particle.curve_idx = idx
    collision!(particle, crv, collision_time)
    #=
    if typeof(crv) <: AbsVirtLine #repeat iteration on internal domain lines and symmetry lines
        domain = billiard.subdomains[particle.subdomain]
        
        collision_time, idx = find_intersection(particle, domain; dt)
        crv = domain.boundary[idx]
        collision!(particle, crv, collision_time)
    end
    =#
end

function trajectory(particle::P, billiard::B, T::Int; dt = 1.0, full_domain=true) where {P<:AbsParticle, B<:AbsBilliard}
    let p = particle
        pts = Vector{typeof(p.r)}(undef,T+1)
        vel = Vector{typeof(p.v)}(undef,T+1)
        ts = Vector{typeof(p.time)}(undef,T+1)
        #sym_ids = Vector{Int64}(undef,T+1)
        pts[1] = p.r
        vel[1] = p.v
        ts[1] = p.time
        #sym_ids[1] = p.sym_sector
        #println("0, r=$(p.r)")
        for i in 1:T
            iterate_bounce!(p, billiard; dt)
            ##println("$i, r=$(p.r), dom_idx=$(p.subdomain)")
            if full_domain
                id = p.sym_sector
                pts[i+1] = billiard.symmetries[id](p.r)
                vel[i+1] = billiard.symmetries[id](p.v)
            else
                pts[i+1] = p.r
                vel[i+1] = p.v
            end
            ts[i+1] = p.time
            #sym_ids[i+1] = p.sym_sector
        end
        return pts, vel, ts
    end
end

export trajectory, iterate_bounce!