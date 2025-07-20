
function invert_point(pt, crv::C) where {C<:AbsCurve}
    inv_x(theta) = curve(crv,theta)[1] - pt[1] 
    inv_y(theta) = curve(crv,theta)[2] - pt[2]
    roots_y = find_zeros(inv_y, (0.0, 1.0))
    roots_x = find_zeros(inv_x, (0.0, 1.0))
    return roots_x, roots_y
    #if roots_y.length>1
    #    roots_x = find_zeros(inv_x, (0.0, 1.0))
    #    check = [isapprox(th_x, th_y)]
    #else
    #    return roots_X[1]
    #end
end

