module CharacterVarieties


using Chevie

# Export structs
export GType, gType


# Export functions
export group_types, algebra_types
export group_type_table, algebra_type_table
export EX, EY

export group_type_data, algebra_type_data
export Mtau, Stau, qdtau, Htau

export isisolated, islevi, ispalindromic, isnonnegative
export myorbit, orderpol, mobius, pi0, nu
export plorbit_reps, plorbits, plevis
export iplorbit_reps, iplorbits, iplevis
export lorbit_reps, lorbits, levis

export palindrome_X, euler_X, nonnegative_Y


include("algebratypes.jl")
include("checks.jl")
include("epolys.jl")
include("grouptypes.jl")
include("gtypes.jl")
include("helpers.jl")
include("plevis.jl")
include("tables.jl")
include("testing.jl")



end # End of module CharacterVarieties
