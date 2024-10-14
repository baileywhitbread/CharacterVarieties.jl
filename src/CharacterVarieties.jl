module CharacterVarieties

using Reexport
@reexport using Chevie
@reexport using Logging
@reexport using JLD2
@reexport using Random

# Export structs
export GType, gType

# Export functions
## These are for users
export group_type_table, algebra_type_table
export group_types, algebra_types
export EX, EY
export ESp4
export plorbit_reps
export plorbits
export plevis
export iplorbit_reps
export iplorbits
export iplevis
export lorbit_reps
export lorbits
export levis

## These are for testing
export palindrome_X, euler_X
export nonnegative_Y, nonnegative_X
export log_nonnegative_X, log_nonnegative_Y
export fast_log_nonnegative_Y
export check_dim_X, check_dim_Y

## Eventually stop exporting these
export group_type_data
export algebra_type_data
export fast_EX
export fast_EY
export Stau_EX
export Mtau_EX
export orderpol

# Needed for main files
include("checks.jl")
include("helpers.jl")
include("plevis.jl")

# Main files
include("grouptypes.jl")
include("algebratypes.jl")
include("epolys.jl")

# Extras
include("cambo.jl")
include("tables.jl")
include("testing.jl")
include("logging.jl")

# file = dirname(dirname(pathof(CharacterVarieties)))*"\\data\\E6_algebra_type_data.jld2"

end # End of module CharacterVarieties