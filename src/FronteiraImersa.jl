module FronteiraImersa

include("pre-processing/include_pre_processing.jl")
include("fluid-domain/include_fluid_domain.jl")
include("structure-domain/include_structure_domain.jl")
include("interface-coupling/include_interface_coupling.jl")

export ReadUserVars
export MakeGrid
export MakePlate, MakeBody


# dx, freestream, dt, T = read_user_vars(freestream, dx, dt, Re, T)



end