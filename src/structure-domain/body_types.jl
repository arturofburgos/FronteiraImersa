abstract type Body{T <: Motion} end

"""
    RigidBody{T} <: Body{T}

This struct represents a rigid body in a simulation, parameterized by the motion type `T`.

# Fields
- `motion::T`: The motion type of the rigid body, which is a subtype of `Motion`.
- `xb::Array{Float64, 2}`: A 2D array containing the (x, y) locations of the body points.
- `x0::Array{Float64, 2}`: A 2D array containing the reference (x, y) locations, used for moving bodies.
- `ub::Array{Float64, 2}`: A 2D array containing the (ẋ, ẏ) velocities of the body points.
- `ds::Array{Float64, 1}`: A 1D array containing the distances between consecutive body points.

# Details
The `RigidBody` struct models a rigid body in a fluid-structure interaction system. It extends the abstract type `Body{T}`, where `T` is a specific type of motion. The struct includes the locations and velocities of the body points, as well as the distances between consecutive points, which are essential for accurately simulating the dynamics and interactions of the rigid body.

"""
struct RigidBody{T} <: Body{T}
    motion::T               # motion type
    xb::Array{Float64, 2}   # (x,y) locations of the body points
    x0::Array{Float64, 2}   # Reference (x,y) location, for moving bodies
    ub::Array{Float64, 2}   # (ẋ, ẏ) velocities of the body points
    ds::Array{Float64, 1}   # distance between consecutive body points
end


# TODO: later on implement FSI bodies
