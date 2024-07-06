include("../src/FronteiraImersa.jl")
using .FronteiraImersa
using Plots

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

# dx = 0.1
# boundary = (0.0, 1.0, 0.0, 1.0)

boundary = (-1.0, 2.0, -1.0, 1.0)
grid = MakeGrid(dx, boundary)
L = 1.0
alpha = 10
plate = MakePlate(L, alpha, dx, 0.0, 0.0)

type = :plate

lenghtscale = 1.0
motion = :static

center = [-0.5; 0.1]

body = [(type = type, lenghtscale = lenghtscale, alpha = alpha, motion = motion, center = center)]

body = MakeBody(body, dx)

xb = body[1].xb

E = SetupReg(grid, body)

plot(xb[:,1], xb[:,2])