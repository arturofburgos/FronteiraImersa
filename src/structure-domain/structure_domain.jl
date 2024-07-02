"""
    MakeBody(body_idx, h)

Creates a collection of rigid bodies based on provided specifications and grid spacing `h`.

# Arguments
- `body_idx::Vector{Dict{Symbol, Any}}`: A vector of dictionaries where each dictionary contains the parameters for a body.
 Each dictionary should have the keys `:type`, `:lengthscale`, `:alpha`, and optionally `:motion` and `:center`.
- `h::Float64`: The grid spacing to be used in the simulation.

# Returns
- `Vector{RigidBody}`: A vector of `RigidBody` instances created based on the specifications provided in `body_idx`.

# Details
This function processes each body specification provided in `body_idx` to create the corresponding `RigidBody` instances.
If the motion type is not specified, it defaults to `Static`. 
The center of the body is also defaulted to `[0.0, 0.0]` if not provided.
Currently, only the `Static` motion type and `plate` body type are supported.
"""
function MakeBody(body_idx, h)

    # TODO: add functionality for the other types of Motion
    if (:motion in keys(body_idx[1])) == true
        if body_idx[1].motion == :static
            bodies = RigidBody{Static}[] # resource to create multiple bodies on only one single body (either way)
        end
    end
    
    # maybe more efficient than above, TODO: check later
    # if :motion in keys(body_j) && body_j.motion == :static
    #     motion = Static()
    # else
    #     motion = Static()
    # end

    # TODO: add functionality for the other types of Motion
    for body_j in body_idx
        if (:motion in keys(body_j)) == true
            if body_j.motion == :static
                motion = Static()
            end
        else 
            motion = Static()
        end

        if (:center in keys(body_j)) == true
            center = body_j.center
        else
            center = [0.0; 0.0]
        end

        if body_j.type == :plate
            body = MakePlate(body_j.lenghtscale, body_j.alpha, h, center[1], center[2]; motion = motion)
        end

        push!(bodies, body) # TODO: this is a bit different than the version from Jared and Goza. Depending on the output I will have to change

    end

    return bodies
end

# TODO: Check the same function displayed in a more concise and efficient way
# function MakeBody(body_idx, h)
#     # Initialize an empty array to store the bodies
#     bodies = RigidBody{Static}[]

#     # Iterate over each body specification
#     for body_j in body_idx
#         # Determine the motion type, default to Static if not specified
#         motion = get(body_j, :motion, :static) == :static ? Static() : Static() # Placeholder for future motion types

#         # Determine the center of the body, default to [0.0, 0.0] if not specified
#         center = get(body_j, :center, [0.0, 0.0])

#         # Create the body based on the type and other specifications
#         if body_j.type == :plate
#             body = MakePlate(body_j[:lengthscale], body_j[:alpha], h, center[1], center[2]; motion=motion)
#         else
#             throw(ArgumentError("Unsupported body type: $(body_j.type)"))
#         end

#         # Add the created body to the bodies array
#         push!(bodies, body)
#     end

#     return bodies
# end