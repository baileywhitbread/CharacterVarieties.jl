module CharacterVarieties

using Reexport
@reexport using Chevie





# Export structs
export GType, gType

# Export functions
## These are for users
export group_type_table, algebra_type_table
export EX, EY

## These support the above
export plorbit_reps, plorbits, plevis
export iplorbit_reps, iplorbits, iplevis
export lorbit_reps, lorbits, levis

export group_types, algebra_types
export group_type_data, algebra_type_data

export isisolated, islevi, myorbit, orderpol, mobius, pi0, nu

## These are for testing
export ispalindromic, isnonnegative
export palindrome_X, euler_X, nonnegative_Y
export bigint_EX, bigint_EY


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
