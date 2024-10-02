module CharacterVarieties

using Reexport
@reexport using Chevie, Logging, Random

# Export structs
export GType, gType

# Export user functions
export group_type_table, algebra_type_table
export group_types, algebra_types
export EX, EY, ESp4
export plorbit_reps, plorbits, plevis
export iplorbit_reps, iplorbits, iplevis
export lorbit_reps, lorbits, levis

# Export functions for testing
export palindrome_X, euler_X, nonnegative_Y, nonnegative_X
export log_nonnegative_Y
export check_dim_X, check_dim_Y

# These may eventually stop being exported
export dimension_XY
export group_type_data
export algebra_type_data

# Include additional files
include("checks.jl")
include("helpers.jl")
include("plevis.jl")
include("grouptypes.jl")
include("algebratypes.jl")
include("epolys.jl")
include("cambo.jl")
include("tables.jl")
include("testing.jl")

end # End of module CharacterVarieties