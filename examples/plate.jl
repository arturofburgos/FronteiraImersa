include("../src/FronteiraImersa.jl")
using .FronteiraImersa


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

dx, freestream, dt, T = ReadUserVars((Ux = t-> t^0.0,), missing, missing, Re, missing)

dx = 0.1
boundary = (0.0, 1.0, 0.0, 1.0)
grid = MakeGrid(dx, boundary)