function find_domain_exit_curves(particle::P, domain::D, dt) where {P<:AbsParticle, D<:AbsSimpleDomain}
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

function find_intersection(particle::P, domain::D; dt = 1.0) where {P<:AbsParticle, D<:AbsSimpleDomain}
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