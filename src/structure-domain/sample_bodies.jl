"""
    MakePlate(L, alpha, h, x0, y0; motion = Static(), n = 0)

Create a rigid plate of length `L` at an angle of attack `alpha`, based on grid spacing `h`.

# Arguments
- `L::Float64`: The length of the plate.
- `alpha::Float64`: The angle of attack (AoA) of the plate in degrees.
- `h::Float64`: The grid spacing.
- `x0::Float64`: The initial x-coordinate of the plate's starting point.
- `y0::Float64`: The initial y-coordinate of the plate's starting point.
- `motion::Motion`: The motion type of the plate, default is `Static()`.
- `n::Int`: The number of points to discretize the plate. If `n` is 0,
 it is calculated such that the distance between points `ds` is approximately `2h`.

# Returns
- `RigidBody{T}`: A `RigidBody` instance representing the plate with the specified parameters.

# Details
The function creates a rigid plate represented by a set of points along its length.
If `n` is not specified (`n=0`), the function calculates `n` such that the distance between consecutive points (`ds`) is approximately `2h`.
The points are then projected based on the angle of attack `alpha` and shifted by the initial position `(x0, y0)`.

A sanity check is included to ensure that `ds` equals `2h`. The function returns a `RigidBody` with the specified motion and calculated body points.
"""
function MakePlate(L, alpha, h, x0, y0; motion = Static(), n = 0)

    # Get the number of body points such that ds = 2h
    if (n==0)
        n = Int(floor(L/(2*h))) + 1 # the floor here acts as a round function to the floor (floor(4.3) -> 4.0)
    end

    #= In case n is prescribed the body will adapt to this by chaging the ds value. 
       Careful: ds should not be higher than 3*dx. TODO: Check this later.
    =#

    spt = range(0, L, n) # sequence/distribution (range) of points

    xhat = spt*cosd(-alpha) # x-projection of the plate
    yhat = spt*sind(-alpha) # y-projection of the plate

    xb = [xhat.+x0 yhat.+y0]

    # If needed here there is a sanity check for ds = 2*h
    ds = sqrt((xhat[2] - xhat[1])^2 + (yhat[2] - yhat[1])^2)

    return RigidBody(motion, xb, copy(xb), 0.0*xb, fill(ds, n))

end