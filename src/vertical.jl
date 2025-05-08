using StaticArrays

include("dubinsmaneuver2d.jl")

function Vertical(qi::SVector{3,Float64}, qf::SVector{3,Float64}, rhomin::Float64, pitchmax::SVector{2,Float64})::DubinsManeuver2D
    maneuver = DubinsManeuver2D(qi, qf, rhomin, DubinsStruct(0.0, 0.0, 0.0, Inf, ""))

    dx = maneuver.qf[1] - maneuver.qi[1]
    dy = maneuver.qf[2] - maneuver.qi[2]
    D = sqrt(dx^2 + dy^2)

    # Distance normalization
    d = D / maneuver.rhomin

    # Normalize the problem using rotation
    rotationAngle = mod2pi(atan(dy, dx))
    a = mod2pi(maneuver.qi[3] - rotationAngle)
    b = mod2pi(maneuver.qf[3] - rotationAngle)

    # CSC
    pathLSL = _LSL(maneuver)
    pathRSR = _RSR(maneuver)
    pathLSR = _LSR(maneuver, pitchmax)
    pathRSL = _RSL(maneuver, pitchmax)
    _paths = [pathLSL, pathRSR, pathLSR, pathRSL]

    # for such a short vector, its actually faster not to sort it
    best_path = nothing
    best_length = Inf

    for p in _paths
        # Check if the turns are too long (do not meet pitch constraint)
        if abs(p.t) < π && abs(p.q) < π
            # Check the inclination based on pitch constraint
            center_angle = maneuver.qi[3] + ((p.case[1] == 'L') ? p.t : -p.t)
            if center_angle < pitchmax[1] || center_angle > pitchmax[2]
                continue
            end
            if p.length < best_length
                best_length = p.length
                best_path = p
            end
        end
    end

    if best_path === nothing
        maneuver.maneuver = DubinsStruct(Inf, Inf, Inf, Inf, "XXX")
    else
        maneuver.maneuver = best_path
    end

    return maneuver
end

########## LSL ##########
function _LSL(self::DubinsManeuver2D)::DubinsStruct
    theta1 = self.qi[3]
    theta2 = self.qf[3]

    if theta1 <= theta2
        # start/end points
        p1 = self.qi[1:2]
        p2 = self.qf[1:2]

        radius = self.rhomin

        c1, s1 = radius * cos(theta1), radius * sin(theta1)
        c2, s2 = radius * cos(theta2), radius * sin(theta2)

        # origins of the turns
        o1 = SVector(p1[1] - s1, p1[2] + c1)
        o2 = SVector(p2[1] - s2, p2[2] + c2)

        diff = o2 - o1
        center_distance = sqrt(diff[1]^2 + diff[2]^2)
        centerAngle = atan(diff[2], diff[1])

        t = mod2pi(-theta1 + centerAngle)
        p = center_distance / radius
        q = mod2pi(theta2 - centerAngle)

        if t > pi
            t = 0.0
            q = theta2 - theta1
            turn_end_y = o2[2] - radius * cos(theta1)
            diff_y = turn_end_y - p1[2]
            if abs(theta1) > 1e-5 && (diff_y < 0 == theta1 < 0)
                p = diff_y / sin(theta1) / radius
            else
                t = p = q = Inf
            end
        end
        if q > pi
            t = theta2 - theta1
            q = 0.0
            turn_end_y = o1[2] - radius * cos(theta2)
            diff_y = p2[2] - turn_end_y
            if abs(theta2) > 1e-5 && (diff_y < 0 == theta2 < 0)
                p = diff_y / sin(theta2) / radius
            else
                t = p = q = Inf
            end
        end
    else
        t = p = q = Inf
    end

    length = (t + p + q) * self.rhomin
    case = "LSL"

    return DubinsStruct(t, p, q, length, case)
end

########## RSR ##########
function _RSR(self::DubinsManeuver2D)::DubinsStruct
    theta1 = self.qi[3]
    theta2 = self.qf[3]

    if theta2 <= theta1
        # start/end points
        p1 = self.qi[1:2]
        p2 = self.qf[1:2]

        radius = self.rhomin

        c1, s1 = radius * cos(theta1), radius * sin(theta1)
        c2, s2 = radius * cos(theta2), radius * sin(theta2)

        # origins of the turns
        o1 = SVector(p1[1] + s1, p1[2] - c1)
        o2 = SVector(p2[1] + s2, p2[2] - c2)

        diff = o2 - o1
        center_distance = sqrt(diff[1]^2 + diff[2]^2)
        centerAngle = atan(diff[2], diff[1])

        t = mod2pi(theta1 - centerAngle)
        p = center_distance / radius
        q = mod2pi(-theta2 + centerAngle)

        if t > pi
            t = 0.0
            q = -theta2 + theta1
            turn_end_y = o2[2] + radius * cos(theta1)
            diff_y = turn_end_y - p1[2]
            if abs(theta1) > 1e-5 && (diff_y < 0 == theta1 < 0)
                p = diff_y / sin(theta1) / radius
            else
                t = p = q = Inf
            end
        end
        if q > pi
            t = -theta2 + theta1
            q = 0.0
            turn_end_y = o1[2] + radius * cos(theta2)
            diff_y = p2[2] - turn_end_y
            if abs(theta2) > 1e-5 && (diff_y < 0 == theta2 < 0)
                p = diff_y / sin(theta2) / radius
            else
                t = p = q = Inf
            end
        end
    else
        t = p = q = Inf
    end

    length = (t + p + q) * self.rhomin
    case = "RSR"

    return DubinsStruct(t, p, q, length, case)
end

########## LSR ##########
function _LSR(self::DubinsManeuver2D, pitchmax::SVector{2,Float64})::DubinsStruct
    theta1 = self.qi[3]
    theta2 = self.qf[3]

    # start/end points
    p1 = self.qi[1:2]
    p2 = self.qf[1:2]

    radius = self.rhomin

    c1, s1 = radius * cos(theta1), radius * sin(theta1)
    c2, s2 = radius * cos(theta2), radius * sin(theta2)

    # origins of the turns
    o1 = SVector(p1[1] - s1, p1[2] + c1)
    o2 = SVector(p2[1] + s2, p2[2] - c2)

    diff = o2 - o1
    center_distance = sqrt(diff[1]^2 + diff[2]^2)

    # not constructible
    if center_distance < 2 * radius
        diff = SVector(sqrt(4.0 * radius * radius - diff[2] * diff[2]), diff[2])
        alpha = π / 2.0
    else
        alpha = asin(2.0 * radius / center_distance)
    end

    centerAngle = atan(diff[2], diff[1]) + alpha

    if centerAngle < pitchmax[2]
        t = mod2pi(-theta1 + centerAngle)
        p = sqrt(max(0.0, center_distance * center_distance - 4.0 * radius * radius)) / radius
        q = mod2pi(-theta2 + centerAngle)
    else
        centerAngle = pitchmax[2]
        t = mod2pi(-theta1 + centerAngle)
        q = mod2pi(-theta2 + centerAngle)

        # points on boundary between C and S segments
        c1, s1 = radius * cos(centerAngle), radius * sin(centerAngle)
        c2, s2 = radius * cos(centerAngle), radius * sin(centerAngle)
        w1 = o1 - SVector(-s1, c1)
        w2 = o2 - SVector(s2, -c2)

        p = (w2[2] - w1[2]) / sin(centerAngle) / radius
    end

    length = (t + p + q) * self.rhomin
    case = "LSR"

    return DubinsStruct(t, p, q, length, case)
end


########## RSL ##########
function _RSL(self::DubinsManeuver2D, pitchmax::SVector{2,Float64})::DubinsStruct
    theta1 = self.qi[3]
    theta2 = self.qf[3]

    # start/end points
    p1 = self.qi[1:2]
    p2 = self.qf[1:2]

    radius = self.rhomin

    c1, s1 = radius * cos(theta1), radius * sin(theta1)
    c2, s2 = radius * cos(theta2), radius * sin(theta2)

    # origins of the turns
    o1 = SVector(p1[1] + s1, p1[2] - c1)
    o2 = SVector(p2[1] - s2, p2[2] + c2)

    diff = o2 - o1
    center_distance = sqrt(diff[1]^2 + diff[2]^2)

    # not constructible
    if center_distance < 2 * radius
        diff = SVector(sqrt(4.0 * radius * radius - diff[2] * diff[2]), diff[2])
        alpha = π / 2.0
    else
        alpha = asin(2.0 * radius / center_distance)
    end

    centerAngle = atan(diff[2], diff[1]) - alpha

    if centerAngle > pitchmax[1]
        t = mod2pi(theta1 - centerAngle)
        p = sqrt(max(0.0, center_distance * center_distance - 4.0 * radius * radius)) / radius
        q = mod2pi(theta2 - centerAngle)
    else
        centerAngle = pitchmax[1]
        t = mod2pi(theta1 - centerAngle)
        q = mod2pi(theta2 - centerAngle)

        # points on boundary between C and S segments
        c1, s1 = radius * cos(centerAngle), radius * sin(centerAngle)
        c2, s2 = radius * cos(centerAngle), radius * sin(centerAngle)
        w1 = o1 - SVector(s1, -c1)
        w2 = o2 - SVector(-s2, c2)

        p = (w2[2] - w1[2]) / sin(centerAngle) / radius
    end

    length = (t + p + q) * self.rhomin
    case = "RSL"

    return DubinsStruct(t, p, q, length, case)
end

