"""
    Curl!(q, s, grid::MultiGrid)

Computes the discrete curl operator of a scalar field `s`, in this case streamfunction and stores the result in the vector field `q`.

Mathematical expression:
    q = Curl(s)

This is wrapped into a LinearMap, so call with:

    q = C * ψ

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
    # inner domain
    i = 2:nx
    j = 2:ny-1
    @views broadcast!(-, qx[i, j], s[i.-1, j], s[i.-1, j.-1])
    # bottom boundary
    j = 1
    @views @. qx[i, j] = s[i-1, j]
    # top boundary
    j = ny
    @views @. qx[i, j] = -s[i-1, j-1]


    # y-fluxes
    # inner domain
    i = 2:nx-1
    j = 2:ny
    @views broadcast!(-, qy[i, j], s[i.-1, j.-1], s[i, j.-1])
    # left boundary
    i = 1
    @views @. qy[i, j] = -s[i, j-1]
    # right boundary
    i = nx
    @views @. qy[i, j] = s[i-1, j-1]

    s = reshape(s, grid.nΓ, 1)
    
    return nothing
end


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
    # inner domain
    i = 2:nx
    j = 2:ny-1
    @views broadcast!(-, qx[i, j], s[i.-1, j], s[i.-1, j.-1])
    # bottom boundary
    j = 1
    @views broadcast!(-, qx[i, j], s[i.-1, j], sbc[B.+i])
    # top boundary
    j = ny
    @views broadcast!(-, qx[i, j], sbc[T.+i], s[i.-1, j.-1]) 

    # TODO: Continue Curl in MultiGrid

    # y-fluxes
    # inner domain
    i = 2:nx-1
    j = 2:ny
    @views broadcast!(-, qy[i, j], s[i.-1, j.-1], s[i, j.-1])
    # left boundary
    i = 1
    @views @. qy[i, j] = -s[i, j-1]
    # right boundary
    i = nx
    @views @. qy[i, j] = s[i-1, j-1]

    s = reshape(s, grid.nΓ, 1)
    
    return nothing
end