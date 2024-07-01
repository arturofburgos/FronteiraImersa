include("../src/IBPM.jl")
using .IBPM


# User requiered variables

Re = 100.0


# User optional variables (uncomment below if you want to define them)

# dx = 0.02
# dt = 0.004
# T  = 1.0

dx = missing
dt = missing
T = missing

Ux = t -> t^0.0
freestream = (Ux = Ux,)

dx, freestream, dt, T = read_user_vars((Ux = t-> t^0.0,), missing, missing, Re, missing)