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


function find_intersection_times(particle, curve, exit_time)
    #time since last bounce
    println("generic_times")
    let
        r(t) = domain_fun(curve, propagate(particle, t))
        t = find_zeros(r, (0.0+1e-14, exit_time)) #10.0*eps(exit_time)
        return t  #returnes all intersections
    end
end 


function find_intersection_times(particle::P, line::C, exit_time) where {P<:PointParticle, C<:LineSegment}
    #time since last bounce
    let
        pt1 = line.pt1
        pt0 = line.pt0
        r = particle.r
        v = particle.v
        m = pt1 .- pt0
        return (pt1[1]*pt0[2] - pt0[1]*pt1[2] + m[2]*r[1] - m[1]*r[2])/(m[1]*v[2]-m[2]*v[1])
    end
end



function find_intersection(particle::P, domain::D; dt = 0.1) where {P<:AbsParticle, D<:SimpleDomain}
    #time since last bounce
    boundary = domain.boundary
    N =  length(boundary)
    inside, approx_exit_time = find_domain_exit_curves(particle, domain, dt) #d gives the curves that were exited
    
    exit_curves = boundary[.~inside]
    crv_idx = collect(1:N)[.~inside]
    idx = crv_idx[1]
    exit_time = approx_exit_time
    for (i,crv) in enumerate(exit_curves)
        times = find_intersection_times(particle, crv, approx_exit_time)
        if  isempty(times)
            #println("Trajectory is broken at time $(particle.time)")
            continue
        else
            t = times[1]
            if t < exit_time
                exit_time = t
                idx = crv_idx[i]
            end
        end    
    end 
    return exit_time, idx
end


function line_polar(r,v,theta)
    type = eltype(v)
    if v[1] == zero(type)
        return @. r[1] / cos(theta)
    end
    if v[2] == zero(type)
        return @. r[2] / sin(theta)
    end
    return @. (v[2]*r[1] - v[1]*r[2])/(v[2]*cos(theta) - v[1]*sin(theta))
end

function determine_brackets(r,v; eps=1e-12) 
    dir = cross(r,v)
    theta0 = atan(r[2],r[1])
    pole = atan(v[2],v[1])
    poles = rem2pi.([-pi, pole-pi, pole, pole+pi,  pi], RoundNearest)
    poles = sort(poles)
    unique = [true for p in poles]
    for i in 1:(length(poles)-1)
        if isapprox(poles[i],poles[i+1])
            unique[i] = false
        end
    end
    poles = poles[unique]
    brackets = [(poles[1],poles[2]),(poles[2],poles[3]),(poles[3],poles[4])]
    if poles[1] < theta0 < poles[2]
        if dir > zero(dir)
            brackets =  [(theta0 + eps, poles[2])]
        else
            brackets =  [(-1.0*pi, theta0 - eps), (poles[3],1.0*pi)]
        end
    elseif poles[2] < theta0 < poles[3]
        if dir > zero(dir)
            brackets =  [(theta0 + eps, poles[3])]
        else
            brackets =  [(poles[2], theta0 - eps)]
        end
    elseif poles[3] < theta0 < poles[4]
        if dir > zero(dir)
            brackets = [(theta0 + eps,1.0*pi),(-1.0*pi, poles[2])]
        else
            brackets = [(poles[3], theta0 - eps)]
        end
    end
    return brackets
end


function find_intersection_angles(particle::PointParticle{T}, curve::PolarSegment{N,T,BC} ; eps=1e-12) where {N, T, BC}
    #time since last bounce
    let  r = particle.r, v = particle.v, R = curve.R, coef=curve.coef
        brackets = determine_brackets(r,v; eps)
        fun(theta) = polar_radius(R,coef,theta) - line_polar(r,v,theta)
        angles = Vector{eltype(r)}()
        for b in brackets
            angles1 = find_zeros(fun, b) #10.0*eps(exit_time)
            append!(angles, angles1)
        end
        return angles  #returns all angles
    end
end 

function find_intersection_times(particle::PointParticle{T}, curve::PolarSegment{N,T,BC}, exit_time) where {N, T, BC}
    let  r = particle.r, v = particle.v, #R = curve.R, coef=curve.coef
        angles = find_intersection_angles(particle, curve)
        pt0 = particle.r
        pts_fun(phi) =  line_polar(r,v,phi) .* SVector{2,Float64}([cos(phi),sin(phi)])
        pts = pts_fun.(angles)
        #pts = line_polar.(r,v,angles)
        speed = norm(particle.v)
        return sort([hypot((pt-pt0)...)/speed for pt in pts])
    end
end

function find_intersection(particle::P, domain::D; dt = 0.1) where {P<:PointParticle, D<:PolarDomain}
    #time since last bounce
    boundary = domain.boundary
    times = find_intersection_times(particle, boundary[1], 100.0)
    if isempty(times)
        return 0.0, 1
    end
    bounce_time = times[1]
    idx  = 1
    for (i,crv) in enumerate(boundary[2:end])
        times = find_intersection_times(particle, crv, 100.0)
        t = times[1]
        if t < bounce_time
            bounce_time = t
            idx = i
        end   
    end 
    return bounce_time, idx
end