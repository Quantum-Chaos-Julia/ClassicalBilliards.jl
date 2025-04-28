function collision_rule!(particle::P, bc::BC, n) where {P<:AbsParticle, BC<:SpecularReflection}
    let v = particle.v #particle initial velocity
        v_n = dot(v,n)*n #normal velocity
        v_t = v .- v_n #tangetial velocity
        v_new = v .- 2*v_n 
        particle.v = v_new
        particle.v_b = SVector(norm(v_t),norm(v_n))
    end
end

function collision_rule!(particle::P, bc::BC, n) where {P<:AbsParticle, BC<:ReflectionSymmetry}
    #collisions with symetry line result in change of symmetry sector
    let v = particle.v , id = particle.sym_sector, sym = bc.symmetry #particle initial velocity
        v_n = dot(v,n)*n #normal velocity
        v_t = v .- v_n #tangetial velocity
        v_new = v .- 2*v_n 
        particle.v = v_new
        particle.v_b = SVector(norm(v_t),norm(v_n))
        new_id = change_sector(bc, id) #change symmetry sector
        particle.sym_sector = new_id
    end
end

#reflections
#sector change when crossing y axis means x coordinate is reflected
# N_sectors is number of symmetry sectors 2 - for only x symmetry, 4 - for x and y symmetry 
function change_sector(bc::ReflectionSymmetry{S}, id) where S<:YAxisReflection
    if bc.N_sectors == 2
        new_id = (id == 1) ? 2 : 1
    else # N_sectors = 4 case
        if id == 1
            new_id = 2
        elseif id == 2
            new_id = 1
        elseif id == 3
            new_id = 4
        elseif id == 4
            new_id = 3
        end
    end
    return new_id
end

#sector change when crossing x axis means y coordinate is reflected
# N_sectors is number of symmetry sectors 2 - for only x symmetry, 4 - for x and y symmetry 
function change_sector(bc::ReflectionSymmetry{S}, id) where S<:XAxisReflection
    if bc.N_sectors == 2
        new_id = (id == 1) ? 2 : 1
    else # N_sectors = 4 case
        if id == 1
            new_id = 4
        elseif id == 2
            new_id = 3
        elseif id == 3
            new_id = 2
        elseif id == 4
            new_id = 1
        end
    end
    return new_id
end

function collision_rule!(particle::P, bc::BC, n) where {P<:AbsParticle, BC<:Transparent}
    #collisions with transparent lines result in change of subdomain
    let v = particle.v
        v_n = dot(v,n)*n #normal velocity
        v_t = v .- v_n #tangetial velocity
        particle.v_b = SVector(norm(v_t),norm(v_n))
        particle.subdomain = bc.next_id
    end
end
