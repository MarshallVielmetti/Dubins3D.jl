using StaticArrays

include("vertical.jl")

"""
    DubinsManeuver3D struct

This struct contains all necessary information about the maneuver.
  * qi - initial configuration (x, y, z, heading, pitch)
  * qf - final configuration (x, y, z, heading, pitch)
  * rhomin - minimum turning radius
  * pitchlims - limits of the pitch angle [pitch_min, pitch_max] 
    where pitch_min < 0.0   
  * path - array containing horizontal and vertical Dubins paths
  * length - total length of the 3D maneuver 
"""
mutable struct DubinsManeuver3D{F<:Real}
    qi::SVector{5,F}  # Fixed-size array for initial configuration
    qf::SVector{5,F}  # Fixed-size array for final configuration

    rhomin::F
    pitchlims::SVector{2,F}  # Fixed-size array for pitch limits

    path::Vector{DubinsManeuver2D}  # Keeping this dynamic as it depends on the number of paths
    length::F
end

function DubinsManeuver3D(qi::Vector, qf::Vector, rhomin, pitchlims::Vector)
    @assert length(qi) >= 5
    @assert length(qf) >= 5
    @assert length(pitchlims) >= 2
    F = promote_type(eltype(qi), eltype(qf), typeof(rhomin), eltype(pitchlims))
    return DubinsManeuver3D(
        SVector{5,F}(qi[1:5]),
        SVector{5,F}(qf[1:5]),
        F(rhomin),
        SVector{2,F}(pitchlims[1:2])
    )
end

function wrapToPi(x::F) where {F<:Real}
    return atan(sin(x), cos(x))
end

"""
    DubinsManeuver3D(qi, qf, rhomin, pitchlims)

Create 3D Dubins path between two configurations qi, qf
    * qi - initial configuration (x, y, z, heading, pitch)
    * qf - final configuration (x, y, z, heading, pitch)
    * rhomin - minimum turning radius
    * pitchlims - limits of the pitch angle [pitch_min, pitch_max] 
"""
function DubinsManeuver3D(qi::SVector{5,F}, qf::SVector{5,F},
    rhomin::F, pitchlims::SVector{2,F}) where {F<:Real}
    maneuver = DubinsManeuver3D(qi, qf, rhomin, pitchlims, Vector{DubinsManeuver2D}(undef, 2), -1.0)

    # Pre-allocate paths for reuse
    dlat_best = DubinsManeuver2D(SVector{3,F}(0, 0, 0), SVector{3,F}(0, 0, 0), 1.0,
        DubinsStruct(0.0, 0.0, 0.0, Inf, ""))
    dlon_best = DubinsManeuver2D(SVector{3,F}(0, 0, 0), SVector{3,F}(0, 0, 0), 1.0,
        DubinsStruct(0.0, 0.0, 0.0, Inf, ""))

    result = Vector{DubinsManeuver2D}(undef, 2)

    # Multiplication factor of rhomin in [1, 1000]
    a = 1.0
    b = 1.0

    found_valid_a = try_to_construct!(result, maneuver, maneuver.rhomin * a)

    found_valid_b = try_to_construct!(result, maneuver, maneuver.rhomin * b)
    while !found_valid_b
        b *= 2.0
        found_valid_b = try_to_construct!(result, maneuver, maneuver.rhomin * b)
        if b > 1000.0
            error("No Dubins path exists")
        end
    end

    best_length = result[2].maneuver.length
    maneuver.length = best_length

    # Binary search with less allocations
    step::F = F(0.1)
    while abs(step) > F(1e-10)
        c = b + step
        c = max(c, F(1.0))
        valid = try_to_construct!(result, maneuver, maneuver.rhomin * c)

        if valid && result[2].maneuver.length < best_length
            b = c
            best_length = result[2].maneuver.length
            maneuver.length = best_length
            step *= F(2.0)
            continue
        end

        step *= F(-0.1)
    end

    # Assign the best paths to maneuver.path
    maneuver.path[1] = result[1]
    maneuver.path[2] = result[2]

    return maneuver
end

function compute_sampling(self::DubinsManeuver3D{F}; numberOfSamples::Integer=1000) where {F<:Real}
    Dlat, Dlon = self.path
    points = Vector{SVector{5,F}}(undef, numberOfSamples)
    lena = Dlon.maneuver.length

    step = lena / F(numberOfSamples - 1)
    for i in 1:numberOfSamples
        offsetLon = F(i - 1) * step
        qSZ = getCoordinatesAt(Dlon, offsetLon)
        qXY = getCoordinatesAt(Dlat, qSZ[1])
        points[i] = SVector{5,F}(qXY[1], qXY[2], qSZ[2], wrapToPi(qXY[3]), wrapToPi(qSZ[3]))
    end

    return points
end

"""
    compute_at_len(self::DubinsManeuver3D; offsetLon::F=0.0)

Compute the coordinates a specific distance along the path
"""
function compute_at_len(self::DubinsManeuver3D, offset::F=0.0) where {F<:Real}
    @assert offset >= 0
    @assert offset <= self.length

    Dlat, Dlon = self.path
    offsetLon = (Dlon.maneuver.length / self.length) * offset

    qSZ = getCoordinatesAt(Dlon, offsetLon)
    qXY = getCoordinatesAt(Dlat, qSZ[1])
    return SVector{5,F}(qXY[1], qXY[2], qSZ[2], wrapToPi(qXY[3]), wrapToPi(qSZ[3]))
end

function try_to_construct!(result::Vector{DubinsManeuver2D}, self::DubinsManeuver3D{F}, horizontal_radius::F) where {F<:Real}
    qi2D = SVector{3,F}(self.qi[1], self.qi[2], self.qi[4])
    qf2D = SVector{3,F}(self.qf[1], self.qf[2], self.qf[4])

    Dlat = DubinsManeuver2D(qi2D, qf2D; rhomin=horizontal_radius)

    vertical_curvature = sqrt(F(1.0) / self.rhomin^2 - F(1.0) / horizontal_radius^2)
    if vertical_curvature < F(1e-5)
        return false
    end
    vertical_radius = F(1.0) / vertical_curvature

    qi3D = SVector{3,F}(0.0, self.qi[3], self.qi[5])
    qf3D = SVector{3,F}(Dlat.maneuver.length, self.qf[3], self.qf[5])

    Dlon = DubinsManeuver2D(qi3D, qf3D; rhomin=vertical_radius)

    if Dlon.maneuver.case == "RLR" || Dlon.maneuver.case == "LRL"
        return false
    end

    if Dlon.maneuver.case[1] == 'R'
        if self.qi[5] - Dlon.maneuver.t < self.pitchlims[1]
            return false
        end
    else
        if self.qi[5] + Dlon.maneuver.t > self.pitchlims[2]
            return false
        end
    end

    result[1] = Dlat
    result[2] = Dlon
    return true
end

function getLowerBound(qi, qf, rhomin=1.0, pitchlims=[-pi / 4, pi / 2])
    F = promote_type(eltype(qi), eltype(qf), typeof(rhomin), eltype(pitchlims))
    qi_sv = SVector{5,F}(qi...)
    qf_sv = SVector{5,F}(qf...)
    return getLowerBound(qi_sv, qf_sv, F(rhomin), SVector{2,F}(pitchlims[1], pitchlims[2]))
end

function getLowerBound(qi::SVector{5,F}, qf::SVector{5,F}, rhomin::F=F(1.0), pitchlims::SVector{2,F}=SVector{2,F}(-pi / 4, pi / 2)) where {F<:Real}
    maneuver = DubinsManeuver3D(qi, qf, rhomin, pitchlims, Vector{DubinsManeuver2D}(), F(-1.0))

    spiral_radius = rhomin * ((cos(max(-pitchlims[1], pitchlims[2])))^2)

    qi2D = @SVector F[maneuver.qi[1], maneuver.qi[2], maneuver.qi[4]]
    qf2D = @SVector F[maneuver.qf[1], maneuver.qf[2], maneuver.qf[4]]

    Dlat = DubinsManeuver2D(qi2D, qf2D; rhomin=spiral_radius)

    qi3D = @SVector F[0.0, maneuver.qi[3], maneuver.qi[5]]
    qf3D = @SVector F[Dlat.maneuver.length, maneuver.qf[3], maneuver.qf[5]]

    Dlon = Vertical(qi3D, qf3D, maneuver.rhomin, maneuver.pitchlims)

    if Dlon.maneuver.case == "XXX"
        # TODO - update Vertical1D such that it compute the shortest prolongation
        maneuver.length = 0.0
        return maneuver
    end

    maneuver.path = [Dlat, Dlon]
    maneuver.length = Dlon.maneuver.length
    return maneuver
end

function getUpperBound(qi, qf, rhomin=1.0, pitchlims=[-pi / 4, pi / 2])
    F = promote_type(eltype(qi), eltype(qf), typeof(rhomin), eltype(pitchlims))
    qi_sv = SVector{5,F}(qi...)
    qf_sv = SVector{5,F}(qf...)
    return getUpperBound(qi_sv, qf_sv, F(rhomin), SVector{2,F}(pitchlims[1], pitchlims[2]))
end

function getUpperBound(qi::SVector{5,F}, qf::SVector{5,F}, rhomin::F=F(1.0), pitchlims::SVector{2,F}=SVector{2,F}(F(-pi / 4), F(pi / 2))) where {F<:Real}
    maneuver = DubinsManeuver3D(qi, qf, rhomin, pitchlims, Vector{DubinsManeuver2D}(), F(-1.0))

    safeRadius = sqrt(2) * maneuver.rhomin

    p1 = qi[1:2]
    p2 = qf[1:2]
    diff = p2 - p1
    dist = sqrt(diff[1]^2 + diff[2]^2)
    if dist < 4.0 * safeRadius
        maneuver.length = Inf
        return maneuver
    end

    qi2D = @SVector F[maneuver.qi[1], maneuver.qi[2], maneuver.qi[4]]
    qf2D = @SVector F[maneuver.qf[1], maneuver.qf[2], maneuver.qf[4]]
    Dlat = DubinsManeuver2D(qi2D, qf2D; rhomin=safeRadius)

    qi3D = @SVector F[0.0, maneuver.qi[3], maneuver.qi[5]]
    qf3D = @SVector F[Dlat.maneuver.length, maneuver.qf[3], maneuver.qf[5]]

    Dlon = Vertical(qi3D, qf3D, safeRadius, maneuver.pitchlims)

    if Dlon.maneuver.case == "XXX"
        # TODO - update Vertical1D such that it compute the shortest prolongation
        maneuver.length = Inf
        return maneuver
    end

    maneuver.path = [Dlat, Dlon]
    maneuver.length = Dlon.maneuver.length
    return maneuver
end