using StaticArrays

# include("dubinsmaneuver2d.jl")
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
mutable struct DubinsManeuver3D
    qi::SVector{5,Float64}  # Fixed-size array for initial configuration
    qf::SVector{5,Float64}  # Fixed-size array for final configuration

    rhomin::Float64
    pitchlims::SVector{2,Float64}  # Fixed-size array for pitch limits

    path::Vector{DubinsManeuver2D}  # Keeping this dynamic as it depends on the number of paths
    length::Float64
end

function DubinsManeuver3D(qi::Vector, qf::Vector, rhomin, pitchlims::Vector)
    @assert length(qi) >= 5
    @assert length(qf) >= 5
    @assert length(pitchlims) >= 2
    return DubinsManeuver3D(
        SVector{5,Float64}(qi[1:5]),
        SVector{5,Float64}(qf[1:5]),
        rhomin,
        SVector{2,Float64}(pitchlims[1:2])
    )
end

"""
    DubinsManeuver3D(qi, qf, rhomin, pitchlims)

Create 3D Dubins path between two configurations qi, qf
    * qi - initial configuration (x, y, z, heading, pitch)
    * qf - final configuration (x, y, z, heading, pitch)
    * rhomin - minimum turning radius
    * pitchlims - limits of the pitch angle [pitch_min, pitch_max] 
"""
function DubinsManeuver3D(qi::SVector{5,Float64}, qf::SVector{5,Float64},
    rhomin::Float64, pitchlims::SVector{2,Float64})
    maneuver = DubinsManeuver3D(qi, qf, rhomin, pitchlims, Vector{DubinsManeuver2D}(undef, 2), -1.0)

    # Pre-allocate paths for reuse
    dlat_best = DubinsManeuver2D(SVector{3,Float64}(0, 0, 0), SVector{3,Float64}(0, 0, 0), 1.0,
        DubinsStruct(0.0, 0.0, 0.0, Inf, ""))
    dlon_best = DubinsManeuver2D(SVector{3,Float64}(0, 0, 0), SVector{3,Float64}(0, 0, 0), 1.0,
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
    step = 0.1
    while abs(step) > 1e-10
        c = b + step
        c = max(c, 1.0)
        valid = try_to_construct!(result, maneuver, maneuver.rhomin * c)

        if valid && result[2].maneuver.length < best_length
            b = c
            best_length = result[2].maneuver.length
            maneuver.length = best_length
            step *= 2.0
            continue
        end

        step *= -0.1
    end

    # Assign the best paths to maneuver.path
    maneuver.path[1] = result[1]
    maneuver.path[2] = result[2]

    return maneuver
end

function compute_sampling(self::DubinsManeuver3D; numberOfSamples::Integer=1000)
    Dlat, Dlon = self.path
    points = Vector{SVector{5,Float64}}(undef, numberOfSamples)
    lena = Dlon.maneuver.length

    step = lena / (numberOfSamples - 1)
    for i in 1:numberOfSamples
        offsetLon = (i - 1) * step
        qSZ = getCoordinatesAt(Dlon, offsetLon)
        qXY = getCoordinatesAt(Dlat, qSZ[1])
        points[i] = SVector{5,Float64}(qXY[1], qXY[2], qSZ[2], qXY[3], qSZ[3])
    end

    points
end

function try_to_construct!(result::Vector{DubinsManeuver2D}, self::DubinsManeuver3D, horizontal_radius::Float64)
    qi2D = SVector{3,Float64}(self.qi[1], self.qi[2], self.qi[4])
    qf2D = SVector{3,Float64}(self.qf[1], self.qf[2], self.qf[4])

    Dlat = DubinsManeuver2D(qi2D, qf2D; rhomin=horizontal_radius)

    vertical_curvature = sqrt(1.0 / self.rhomin^2 - 1.0 / horizontal_radius^2)
    if vertical_curvature < 1e-5
        return false
    end
    vertical_radius = 1.0 / vertical_curvature

    qi3D = SVector{3,Float64}(0.0, self.qi[3], self.qi[5])
    qf3D = SVector{3,Float64}(Dlat.maneuver.length, self.qf[3], self.qf[5])

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

function getLowerBound(qi, qf, rhomin::Float64=1.0, pitchlims::SVector{2,Float64}=@SVector Float64[-pi/4, pi/2])
    qi = SVector{5,Float64}(qi...)
    qf = SVector{5,Float64}(qf...)
    maneuver = DubinsManeuver3D(qi, qf, rhomin, pitchlims, [], -1.0)

    spiral_radius = rhomin * ((cos(max(-pitchlims[1], pitchlims[2])))^2)

    qi2D = @SVector Float64[maneuver.qi[1], maneuver.qi[2], maneuver.qi[4]]
    qf2D = @SVector Float64[maneuver.qf[1], maneuver.qf[2], maneuver.qf[4]]

    Dlat = DubinsManeuver2D(qi2D, qf2D; rhomin=spiral_radius)

    qi3D = @SVector Float64[0.0, maneuver.qi[3], maneuver.qi[5]]
    qf3D = @SVector Float64[Dlat.maneuver.length, maneuver.qf[3], maneuver.qf[5]]

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

function getUpperBound(qi, qf, rhomin=1, pitchlims=[-pi / 4, pi / 2])
    qi = SVector{5,Float64}(qi...)
    qf = SVector{5,Float64}(qf...)

    maneuver = DubinsManeuver3D(qi, qf, rhomin, pitchlims, [], -1.0)

    safeRadius = sqrt(2) * maneuver.rhomin

    p1 = qi[1:2]
    p2 = qf[1:2]
    diff = p2 - p1
    dist = sqrt(diff[1]^2 + diff[2]^2)
    if dist < 4.0 * safeRadius
        maneuver.length = Inf
        return maneuver
    end

    qi2D = @SVector Float64[maneuver.qi[1], maneuver.qi[2], maneuver.qi[4]]
    qf2D = @SVector Float64[maneuver.qf[1], maneuver.qf[2], maneuver.qf[4]]
    Dlat = DubinsManeuver2D(qi2D, qf2D; rhomin=safeRadius)

    qi3D = @SVector Float64[0.0, maneuver.qi[3], maneuver.qi[5]]
    qf3D = @SVector Float64[Dlat.maneuver.length, maneuver.qf[3], maneuver.qf[5]]

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