abstract type Motion end

#= From the type Motion we define 2 subtypes:

    - InertialMotion - for cases where the fluid grid is not moving. We also define more two structs, those are
        ∘ Static
        ∘ MotionFunction (TODO)
        
    - MovingGrid - for cases where the fluid grid is moving (TODO)


    TODO(add the move_body_utils.jl -> for moving body or moving grid cases).
    As for the first implementation of the Immersed Body Method, here I focused only in the Static body.
    In future commits, bring the MotionFunction and MovingGrid functionalities. 
=#

abstract type InertialMotion <: Motion end


"""
    Static <: InertialMotion

This struct represents a static body in a simulation.
"""
struct Static <: InertialMotion
end




