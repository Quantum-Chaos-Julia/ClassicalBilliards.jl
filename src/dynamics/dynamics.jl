include("intersections.jl")
export find_domain_exit_curves, find_intersection_time, find_intersection
include("colissionrules.jl")
export collision_rule!


function collision!(particle::P, curve::C, time) where {P<:PointParticle, C<:AbsCurve}
    r = propagate(particle, time) #colision point
    g = domain_gradient_vector(curve, r)
    n =  g./norm(g) 
    collision_rule!(particle, curve.bc, n) #modify velocity using rule corresponding to boundary
    particle.r = r
    particle.time += time
end


function iterate_bounce!(particle::P, billiard::B; dt = 1.0 ) where {P<:AbsParticle, B<:AbsBilliard}
    domain = billiard.fundamental_domain.subdomains[particle.subdomain]    
    collision_time, idx = find_intersection(particle, domain; dt)
    crv = domain.boundary[idx]
    #particle.curve_idx = idx
    collision!(particle, crv, collision_time)
    if typeof(crv.bc) <: Transparent
        iterate_bounce!(particle, billiard; dt)
    end
end

function trajectory(particle::P, billiard::B, T::Int; dt = 1.0, full_domain=false) where {P<:AbsParticle, B<:AbsBilliard}
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

export trajectory, iterate_bounce!, colission!