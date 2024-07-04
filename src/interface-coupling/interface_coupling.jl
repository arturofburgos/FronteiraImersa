# Here we  define the delta function, responsible for the coupling between the fluid and the structure.
# This is the core of the IBM framework. TODO: check with Professor Andres how this delta function was constructed

"""
    δh(rf, rb, dr)

Discrete delta function used to relate flow to structure quantities.

# Arguments
- `rf`: Position in the flow domain.
- `rb`: Position of the immersed boundary (IB) point.
- `dr`: Grid spacing in the r direction.

# Returns
- `del_h`: Value of the discrete delta function at the given positions.

# Details
This function computes a discrete delta function using the Yang smooth delta function (Yang et al, JCP, 2009 - https://doi.org/10.1016/j.jcp.2009.07.023). 
The support of this delta function is 6*h (3*h on each side). This is crucial to make sure that there is no fluid leaking inside the immersed body.

The computation is based on the distance between `rf` and `rb`. Different polynomial expressions are used depending 
on the distance range:

1. If the distance `r1 = abs(rf - rb) / dr` is less than or equal to 1, a specific polynomial expression is used.
2. If `r1` is between 1 and 2, another polynomial expression is used.
3. If `r1` is between 2 and 3, yet another polynomial expression is used.
4. If `r1` is greater than 3, the delta function value is zero (maximum support for this delta function in the stagerred grid).

Disregard the following information if you are an user. 
Note: There are slight differences in the floating point arithmetic which can lead to results differing by around 1e-4 
from those computed using Fortran.
"""
function δh(rf, rb, dr)
    # Calculate the absolute distance between rf and rb
    r = abs(rf - rb)
    r1 = r / dr
    r2 = r1 * r1
    r3 = r2 * r1
    r4 = r3 * r1

    if r1 <= 1.0
        # Coefficients for r1 <= 1.0
        a5 = asin(0.5 * sqrt(3.0) * (2.0 * r1 - 1.0))
        a8 = sqrt(1.0 - 12.0 * r2 + 12.0 * r1)

        del_h = 4.166666667e-2 * r4 + (-0.1388888889 + 3.472222222e-2 * a8) * r3 +
                (-7.121664902e-2 - 5.208333333e-2 * a8 + 0.2405626122 * a5) * r2 +
                (-0.2405626122 * a5 - 0.3792313933 + 0.1012731481 * a8) * r1 +
                8.0187537413e-2 * a5 - 4.195601852e-2 * a8 + 0.6485698427

    elseif r1 <= 2.0
        # Coefficients for 1.0 < r1 <= 2.0
        a6 = asin(0.5 * sqrt(3.0) * (-3.0 + 2.0 * r1))
        a9 = sqrt(-23.0 + 36.0 * r1 - 12.0 * r2)

        del_h = -6.250000000e-2 * r4 + (0.4861111111 - 1.736111111e-2 * a9) * r3 +
                (-1.143175026 + 7.812500000e-2 * a9 - 0.1202813061 * a6) * r2 +
                (0.8751991178 + 0.3608439183 * a6 - 0.1548032407 * a9) * r1 -
                0.2806563809 * a6 + 8.22848104e-3 + 0.1150173611 * a9

    elseif r1 <= 3.0
        # Coefficients for 2.0 < r1 <= 3.0
        a1 = asin(0.5 * sqrt(3.0) * (2.0 * r1 - 5.0))
        a7 = sqrt(-71.0 - 12.0 * r2 + 60.0 * r1)

        del_h = 2.083333333e-2 * r4 + (3.472222222e-3 * a7 - 0.2638888889) * r3 +
                (1.214391675 - 2.604166667e-2 * a7 + 2.405626122e-2 * a1) * r2 +
                (-0.1202813061 * a1 - 2.449273192 + 7.262731481e-2 * a7) * r1 +
                0.1523563211 * a1 + 1.843201677 - 7.306134259e-2 * a7

    else
        del_h = 0.0
    end

    return del_h
end

# The delta function is used in the formulation of the discrete interpolation and regularization operators
# Note that one is the transpose of the other, hence H = -E^t or E = -H^t. 

#=
TODO: Check with Professor Andres is the variable supp stands for support (in this case 6 as it seems to be for the δ function support).
=#
function SetupReg(grid::T, bodies::Array{<:Body, 1}; supp = 6) where T <:Grid
    
    # Extract the sizes of the fluid staggered grid
    nx = grid.nx
    ny = grid.ny
    h = grid.h

    # Stack all the body points, don't need to distinguish them here
    xb = vcat([body.xb for body in bodies]...)
    nb = size(xb,1) # total number of body points (considering all of the bodies since we stack them in the previous code line)
    nf = 2*nb # total number of forces in the body points. Remember that each body will have a x and y force component
 
    supp_idx = -supp:supp # TODO: ask Nick if shouldn't be -supp


