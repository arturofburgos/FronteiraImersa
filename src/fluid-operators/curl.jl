"""
    Curl!(q, s, grid::MultiGrid)

Computes the discrete curl operator of a scalar field `s`, in this case streamfunction and stores the result in the vector field `q`.

Mathematical expression:
    q = Curl(s)

This is wrapped into a LinearMap, so call with:

    q = C * s

where C is in the IBMatrices struct. Note that this doesn't actualy get used in the MultiGrid formulation, but it is the adjoint to `CurlT`, which does get called.  See below for the MultiGrid version of curl.
# Arguments
- `q`: Vector field (flux) array to store the computed curl. It will be modified in place.
- `s`: Streamfunction Scalar field array representing the potential function. It will be reshaped within the function.
- `grid`: The grid on which the curl operation is performed.

# Functionality
1. Extracts the grid dimensions `nx` and `ny`.
2. Reshapes the scalar field `s` to match the grid dimensions.
3. Splits the vector field `q` into its x-component `qx` and y-component `qy`.
4. Computes the curl in the inner domain and on the boundaries for both x and y fluxes.
5. Reshapes the scalar field `s` back to its original shape.

# Returns
- Nothing. The function modifies `q` and `s` in place.

"""
function Curl!(q, s, grid)
    nx = grid.nx
    ny = grid.ny
    s = reshape(s, nx-1, ny-1)
    qx, qy = grid.SplitFlux(q)


    # x-fluxes

    # Inner domain
    i = 2:nx
    j = 2:ny-1
    @views broadcast!(-, qx[i, j], s[i.-1, j], s[i.-1, j.-1])

    # Bottom boundary
    j = 1
    @views @. qx[i, j] = s[i-1, j]
    # Top boundary
    j = ny
    @views @. qx[i, j] = -s[i-1, j-1]


    # y-fluxes

    # Inner domain
    i = 2:nx-1
    j = 2:ny
    @views broadcast!(-, qy[i, j], s[i.-1, j.-1], s[i, j.-1])

    # Left boundary
    i = 1
    @views @. qy[i, j] = -s[i, j-1]
    # Right boundary
    i = nx
    @views @. qy[i, j] = s[i-1, j-1]

    s = reshape(s, grid.nΓ, 1)
    
    return nothing
end


#= 
As mentioned in the definition of the previous Curl function, below there is defined the full Curl that takes into account the boundary conitions for each domain.
This is the function called in IBMatrices.
=#
"""
    Curl!(q, s, sbc, grid::MultiGrid)

Computes the curl of the streamfunction `s` to obtain the grid flux `q` in a fluid dynamics simulation. This function accounts for boundary conditions and updates the flux components in place.

# Arguments
- `q`: grid flux array to be updated. It will be modified in place.
- `s`: Streamfunction array. It is reshaped within the function to a matrix of size `(nx-1, ny-1)`.
- `sbc`: Boundary condition values for the streamfunction on the domain boundaries.
- `grid::MultiGrid`: The fluid grid object containing grid dimensions and boundary indices.

# Functionality
1. Extracts grid dimensions and boundary indices from the `grid` object.
2. Reshapes the streamfunction array `s` to match the inner grid dimensions `(nx-1, ny-1)`.
3. Splits the flux array `q` into its x-component `qx` and y-component `qy`.
4. Computes the x-fluxes:
   - Updates inner domain x-fluxes.
   - Handles bottom and top boundary x-fluxes using boundary conditions.
   - Handles left and right boundary x-fluxes using boundary conditions.
5. Computes the y-fluxes:
   - Updates inner domain y-fluxes.
   - Handles left and right boundary y-fluxes using boundary conditions.
   - Handles bottom and top boundary y-fluxes using boundary conditions.
6. Reshapes the streamfunction array `s` back to its original shape.
7. Returns `nothing`.

# Returns
- Nothing. The function modifies `q` and `s` in place.

"""
function Curl!(q, s, sbc, grid::MultiGrid)
    nx = grid.nx
    ny = grid.ny
    T = grid.TOP
    B = grid.BOTTOM
    L = grid.LEFT
    R = grid.RIGHT
    s = reshape(s, nx-1, ny-1)
    qx, qy = grid.SplitFlux(q)

    
    # x-fluxes
    
    # Inner domain
    i = 2:nx
    j = 2:ny-1
    @views broadcast!(-, qx[i, j], s[i.-1, j], s[i.-1, j.-1])
    
    # Bottom boundary
    j = 1
    @views broadcast!(-, qx[i, j], s[i.-1, j], sbc[B.+i]) # TODO: here is different than the old formulation (sbc[(B-1).+i]). I believe mine is correct, since I think I have fixxed the boundary conditions in the grid (fluid_dommain.jl). Delete this once the first teste in done and you checked the result.
    # Top boundary
    j = ny
    @views broadcast!(-, qx[i, j], sbc[T.+i], s[i.-1, j.-1]) 

    # Left boundary
    j = 1:ny
    i = 1
    @views broadcast!(-, qx[i, j], sbc[(L+1).+j], sbc[(L).+j]) # remember that the streamfunction quantity is stored in the vertices and there are ny+1. Therefore L+1 reaches ny+1 at the maximum value.
    # Right boundary
    i = nx+1
    @views broadcast!(-, qx[i, j], sbc[(R+1).+j], sbc[(R).+j])


    # y-fluxes
    
    # Inner domain
    i = 2:nx-1
    j = 2:ny
    @views broadcast!(-, qy[i, j], s[i.-1, j.-1], s[i, j.-1])
    
    # Left boundary
    i = 1
    @views broadcast!(-, qy[i, j], sbc[L.+j], s[i, j.-1]) 
    # Right boundary
    i = nx
    @views broadcast!(-, qy[i, j], s[i.-1, j.-1], sbc[R.+j])

    # Bottom boundary
    i = 1:nx
    j = 1
    @views broadcast!(-, qy[i, j], sbc[B.+i], sbc[(B+1).+i])
    # Top boundary
    j = ny+1
    @views broadcast!(-, qy[i, j], sbc[T.+i], sbc[(T+1).+i])


    s = reshape(s, grid.nΓ, 1)
    
    return nothing
end

#=
Now, there is the definition of the CurlT

    Γ = ∇×q --> Mathematical expression

    Γ = C^T * q

Hence the circulation (or "vorticity" in 2D).
=#

"""
    CurlT!(Γ, q, grid, Γwork)

Computes the circulation (or vorticity in 2D) from the grid forces in a fluid-structure interaction simulation.

# Arguments
- `Γ`: Circulation array to be updated. It will be reshaped within the function to match the grid dimensions.
- `q`: Grid flux array. It is split within the function into its x-component and y-component.
- `grid`: Grid object containing the fluid grid dimensions and properties.
- `Γwork`: Working array for intermediate calculations.

# Functionality
1. Reshapes the circulation array `Γ` from a vector to a matrix with dimensions `(nx-1, ny-1)`.
2. Splits the grid flux array `q` into its x-component `qx` and y-component `qy`.
3. Iterates over the grid interior to update the circulation based on the grid flux.
4. Uses the working array `Γwork` to store intermediate results and combines them with `Γ`.
5. Reshapes the circulation array back to a single vector.

# Internal Details
- The grid indices around each point are calculated considering the grid dimensions.
- The function uses `@views` to avoid unnecessary array allocations when updating `Γ` and `Γwork`.
- The `@.` macro is used to broadcast operations element-wise across arrays.

# Returns
- Nothing. The function modifies `Γ` in place.
"""
function CurlT!(Γ, q, grid, Γwork)
    nx = grid.nx
    ny = grid.ny

    Γ = reshape(Γ, nx-1, ny-1)

    qx, qy = grid.SplitFlux(q)

    i = 2:nx
    j = 2:ny

    @views broadcast!(-, Γ[i.-1, j.-1], qx[i, j.-1], qx[i, j])
    @views broadcast!(-, Γwork[i.-1, j.-1], qy[i, j], qy[i.-1, j])
    #=
        The below @. macro in Julia is a shorthand for broadcasting all the operations in an expression. 
        Broadcasting allows operations to be applied element-wise across arrays or collections, automatically handling different shapes and sizes.
        When you use @., it will apply the dot (.) to each function call and operator in the expression.
    =#
    @views @. Γ[i-1, j-1] = Γ[i-1, j-1] + Γwork[i-1, j-1]
    Γ = reshape(Γ, grid.nΓ, 1)

    return nothing
end


# TODO: build Vort2Flux! function. Maybe move from here to it's own julia script. 
function Vort2Flux!(s, q, Γ, model::IBModel{MultiGrid, <:Body}, ngrids::Int64)

end