"""
Definition of fluid grid domain
"""

abstract type Grid end

struct MultiGrid <: Grid
    nx::Int64
    ny::Int64
    nΓ::Int64
    nu::Int64
    nv::Int64
    nq::Int64
    mg::Int64
    offx::Float64
    offy::Float64
    len::Float64
    h::Float64
    SplitFlux::Any
    LEFT::Int64
    RIGHT::Int64
    BOTTOM::Int64
    TOP::Int64
end

"""
    MakeGrid(h::Float64, boundary::NTuple{4,Float64}; mg=5::Int64)

This function creates a grid based on the provided parameters and returns a `MultiGrid` struct.

# Arguments
- `h::Float64`: Grid spacing.
- `boundary::NTuple{4,Float64}`: Tuple defining the grid boundaries as (xmin, xmax, ymin, ymax).
- `mg::Int64=5`: Number of multigrid levels (default is 5).

# Returns
- `MultiGrid`: A `MultiGrid` object with the grid parameters and utility function for splitting fluxes.

# Details
1. Computes grid offsets (`offx`, `offy`) based on the boundary.
2. Determines grid dimensions (`xlen`, `ylen`) and the number of grid points in x and y directions (`nx`, `ny`).
3. Calculates the number of fluxes in the x (`nu`) and y (`nv`) directions, and the total number of fluxes (`nq`).
4. Computes the number of circulation points (`nΓ`), corresponding to cell vertices.
5. Defines a utility function `SplitFlux` to split the flux vector into 2D arrays for u-velocity (`qu`) and v-velocity (`qv`).
6. Predefines constant offsets for indexing boundary conditions (`left`, `right`, `bottom`, `top`).
7. Returns a `MultiGrid` object initialized with the computed parameters and the `SplitFlux` function.

"""
function MakeGrid(h::Float64, boundary::NTuple{4,Float64}; mg=5::Int64)
    
    offx = -boundary[1]
    offy = -boundary[3]
    xlen = boundary[2] - boundary[1]
    ylen = boundary[4] - boundary[3]

    nx = Int64(round(xlen/h)) # Defined based in the centered points
    ny = Int64(round(ylen/h)) # Defined based in the centered points
    
    nu = ny*(nx+1) # "u-velocity" fluxes quantity number
    nv = nx*(ny+1) # "v-velocity" fluxes quantity number
    
    nq = nu+nv # total number of fluxes


    nΓ = (nx-1)*(ny-1) # total number of circulation("vorticity" in 3D) point - corresponding to the cell vertices
    # TODO - check if it is not nΓ = (nx-2)*(ny-2). I believe the two TODO's in this script are related to one another


    """
    SplitFlux(q; lev = 1)
    
    Return views to 2D arrays of fluxes. Initially we start with a big matrix where the rows represent the fluxes
    in the x, y order [qu (x), qv (y)] and the columns represents the level where this flux is located

    . . . . .                .               . . . 
    . . . . .                .     qu        . . . 
    . . . . .                .      
    . . . . .                .
    . . . . . --> 1st column . --> 
    . . . . .                .
    . . . . .                .
    . . . . .                .               . . . 
    . . . . .                .     qv        . . .
    . . . . .                .   
    . . . . .                .
    . . . . .                .
    """
    SplitFlux(q; lev = 1) = reshape(@view(q[1:nu, lev]), nx+1, ny), reshape(@view(q[nu+1:end, lev]), nx, ny+1) # qu, qv

    # Predefine constant offsets for indexing boundary conditions
    # TODO - Still need to understand here, why those values for them - Most likely to understand when dealing with the CURL operator
    left = 0
    right =  ny+1
    bottom = 2*(ny+1)
    top = 2*(ny+1) + nx+1

    return MultiGrid(nx, ny, nΓ, nu, nv, nq, mg, offx, offy, xlen, h, SplitFlux, left, right, bottom, top)
end

