module FronteiraImersa

using LinearMaps
using Plots


include("fluid-domain/include_fluid_domain.jl")
include("structure-domain/include_structure_domain.jl")
include("interface-coupling/include_interface_coupling.jl")
include("pre-processing/include_pre_processing.jl")
include("fluid-operators/include_fluid_operators.jl")

export ReadUserVars
export MakeGrid
export MakePlate, MakeBody
export SetupReg

# dx, freestream, dt, T = read_user_vars(freestream, dx, dt, Re, T)



end