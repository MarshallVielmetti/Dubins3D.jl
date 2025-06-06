"""
Classical 2D Dubins Curve
"""
struct DubinsStruct{F<:Real}
    t::F
    p::F
    q::F
    length::F
    case::String
end

"""
Classical 2D Dubins Curve
"""
mutable struct DubinsManeuver2D{F<:Real}
    qi::SVector{3,F}
    qf::SVector{3,F}
    rhomin::F
    maneuver::DubinsStruct{F}
end

"""
Classical 2D Dubins Curve
"""
function DubinsManeuver2D(qi::SVector{3,F}, qf::SVector{3,F}; rhomin::F=F(1.0), minLength::Union{Nothing,F}=nothing, disable_CCC::Bool=false)::DubinsManeuver2D{F} where {F<:Real}
    maneuver = DubinsManeuver2D{F}(qi, qf, rhomin, DubinsStruct{F}(F(0.0), F(0.0), F(0.0), F(Inf), ""))

    dx = maneuver.qf[1] - maneuver.qi[1]
    dy = maneuver.qf[2] - maneuver.qi[2]
    D = sqrt(dx^2 + dy^2)

    # Distance normalization
    d = D / maneuver.rhomin

    # Normalize the problem using rotation
    rotationAngle = mod2pi(atan(dy, dx))
    a = mod2pi(maneuver.qi[3] - rotationAngle)
    b = mod2pi(maneuver.qf[3] - rotationAngle)

    sa, ca = sincos(a)
    sb, cb = sincos(b)

    # CSC
    pathLSL = _LSL(maneuver, a, b, d, sa, ca, sb, cb)
    pathRSR = _RSR(maneuver, a, b, d, sa, ca, sb, cb)
    pathLSR = _LSR(maneuver, a, b, d, sa, ca, sb, cb)
    pathRSL = _RSL(maneuver, a, b, d, sa, ca, sb, cb)

    if disable_CCC
        _paths = [pathLSL, pathRSR, pathLSR, pathRSL]
    else
        # CCC
        pathRLR = _RLR(maneuver, a, b, d, sa, ca, sb, cb)
        pathLRL = _LRL(maneuver, a, b, d, sa, ca, sb, cb)
        _paths = [pathLSL, pathRSR, pathLSR, pathRSL, pathRLR, pathLRL]
    end

    if (abs(d) < maneuver.rhomin * 1e-5 && abs(a) < maneuver.rhomin * 1e-5 && abs(b) < maneuver.rhomin * 1e-5)
        dist_2D = maximum(abs.(maneuver.qi[1:2] - maneuver.qf[1:2]))
        if dist_2D < maneuver.rhomin * 1e-5
            pathC = _C(maneuver)
            _paths = [pathC]
        end
    end

    a(x) = x.length
    sort!(_paths, by=a)

    if (minLength === nothing)
        maneuver.maneuver = _paths[1]
    else
        for p in _paths
            if p.length >= minLength
                maneuver.maneuver = p
                break
            end
        end

        if (maneuver.maneuver === nothing)
            inf = Inf
            maneuver.maneuver = DubinsManeuver2D(inf, inf, inf, inf, "XXX")
        end
    end

    return maneuver
end

########## LSL ##########
@inline function _LSL(self::DubinsManeuver2D{F}, a::F, b::F, d::F, sa::F, ca::F, sb::F, cb::F)::DubinsStruct{F} where {F<:Real}
    aux = atan(cb - ca, d + sa - sb)
    t = mod2pi(-a + aux)
    p = sqrt(2 + d^2 - 2 * cos(a - b) + 2 * d * (sa - sb))
    q = mod2pi(b - aux)
    length = (t + p + q) * self.rhomin
    case = "LSL"
    return DubinsStruct{F}(t, p, q, length, case)
end

########## RSR ##########
@inline function _RSR(self::DubinsManeuver2D{F}, a::F, b::F, d::F, sa::F, ca::F, sb::F, cb::F)::DubinsStruct{F} where {F<:Real}
    aux = atan(ca - cb, d - sa + sb)
    t = mod2pi(a - aux)
    p = sqrt(2 + d^2 - 2 * cos(a - b) + 2 * d * (sb - sa))
    q = mod2pi(mod2pi(-b) + aux)
    length = (t + p + q) * self.rhomin
    case = "RSR"
    return DubinsStruct{F}(t, p, q, length, case)
end

########## LSR ##########
@inline function _LSR(self::DubinsManeuver2D{F}, a::F, b::F, d::F, sa::F, ca::F, sb::F, cb::F)::DubinsStruct{F} where {F<:Real}
    aux1 = -2 + d^2 + 2 * cos(a - b) + 2 * d * (sa + sb)
    if (aux1 > 0)
        p = sqrt(aux1)
        aux2 = atan(-ca - cb, d + sa + sb) - atan(-2 / p)
        t = mod2pi(-a + aux2)
        q = mod2pi(-mod2pi(b) + aux2)
    else
        t = p = q = F(Inf)
    end
    length = (t + p + q) * self.rhomin
    case = "LSR"
    return DubinsStruct{F}(t, p, q, length, case)
end

########## RSL ##########
@inline function _RSL(self::DubinsManeuver2D{F}, a::F, b::F, d::F, sa::F, ca::F, sb::F, cb::F)::DubinsStruct{F} where {F<:Real}
    aux1 = d^2 - 2 + 2 * cos(a - b) - 2 * d * (sa + sb)
    if (aux1 > 0)
        p = sqrt(aux1)
        aux2 = atan(ca + cb, d - sa - sb) - atan(2 / p)
        t = mod2pi(a - aux2)
        q = mod2pi(mod2pi(b) - aux2)
    else
        t = p = q = F(Inf)
    end
    length = (t + p + q) * self.rhomin
    case = "RSL"
    return DubinsStruct{F}(t, p, q, length, case)
end

########## RLR ##########
@inline function _RLR(self::DubinsManeuver2D{F}, a::F, b::F, d::F, sa::F, ca::F, sb::F, cb::F)::DubinsStruct{F} where {F<:Real}
    aux = (6 - d^2 + 2 * cos(a - b) + 2 * d * (sa - sb)) / 8
    if (abs(aux) <= 1)
        p = mod2pi(-acos(aux))
        t = mod2pi(a - atan(ca - cb, d - sa + sb) + p / 2)
        q = mod2pi(a - b - t + p)
    else
        t = p = q = F(Inf)
    end
    length = (t + p + q) * self.rhomin
    case = "RLR"
    return DubinsStruct{F}(t, p, q, length, case)
end


########## LRL ##########
@inline function _LRL(self::DubinsManeuver2D{F}, a::F, b::F, d::F, sa::F, ca::F, sb::F, cb::F)::DubinsStruct{F} where {F<:Real}
    aux = (6 - d^2 + 2 * cos(a - b) + 2 * d * (-sa + sb)) / 8
    if (abs(aux) <= 1)
        p = mod2pi(-acos(aux))
        t = mod2pi(-a + atan(-ca + cb, d + sa - sb) + p / 2)
        q = mod2pi(b - a - t + p)
    else
        t = p = q = F(Inf)
    end
    length = (t + p + q) * self.rhomin
    case = "LRL"
    return DubinsStruct{F}(t, p, q, length, case)
end

########## C ##########
@inline function _C(self::DubinsManeuver2D{F})::DubinsStruct{F} where {F<:Real}
    t = F(0.0)
    p = F(2 * pi)
    q = F(0.0)
    length = (t + p + q) * self.rhomin
    case = "RRR"
    return DubinsStruct{F}(t, p, q, length, case)
end

function getCoordinatesAt(self::DubinsManeuver2D{F}, offset::F)::SVector{3,F} where {F<:Real}
    # Offset normalizado
    noffset = offset / self.rhomin

    qi = SVector{3,F}(F(0.0), F(0.0), self.qi[3])

    # Gerando as configurações intermediárias            
    l1 = self.maneuver.t
    l2 = self.maneuver.p
    q1 = getPositionInSegment(self, l1, qi, self.maneuver.case[1]) # Final do segmento 1
    q2 = getPositionInSegment(self, l2, q1, self.maneuver.case[2]) # Final do segmento 2

    # Obtendo o restante das configurações
    if (noffset < l1)
        q = getPositionInSegment(self, noffset, qi, self.maneuver.case[1])
    elseif (noffset < (l1 + l2))
        q = getPositionInSegment(self, noffset - l1, q1, self.maneuver.case[2])
    else
        q = getPositionInSegment(self, noffset - l1 - l2, q2, self.maneuver.case[3])
    end

    return SVector{3,F}(q[1] * self.rhomin + self.qi[1],
        q[2] * self.rhomin + self.qi[2],
        mod2pi(q[3]))
end

function getPositionInSegment(self::DubinsManeuver2D{F}, offset::F, qi::SVector{3,F}, case::Char)::SVector{3,F} where {F<:Real}
    if (case == 'L')
        return SVector{3,F}(
            qi[1] + sin(qi[3] + offset) - sin(qi[3]),
            qi[2] - cos(qi[3] + offset) + cos(qi[3]),
            qi[3] + offset
        )
    elseif (case == 'R')
        return SVector{3,F}(
            qi[1] - sin(qi[3] - offset) + sin(qi[3]),
            qi[2] + cos(qi[3] - offset) - cos(qi[3]),
            qi[3] - offset
        )
    elseif (case == 'S')
        return SVector{3,F}(
            qi[1] + cos(qi[3]) * offset,
            qi[2] + sin(qi[3]) * offset,
            qi[3]
        )
    end
    return SVector{3,F}(F(0.0), F(0.0), F(0.0))  # Should never reach here
end

function getSamplingPoints(self::DubinsManeuver2D{F}, res::F=F(0.1))::Vector{SVector{3,F}} where {F<:Real}
    points = Vector{SVector{3,F}}()
    range = F(0.0):res:self.maneuver.length
    for offset in range
        push!(points, getCoordinatesAt(self, offset))
    end
    return points
end

