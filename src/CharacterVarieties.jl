module CharacterVarieties

using Reexport
@reexport using Chevie





# Export structs
export GType, gType

# Export functions
## These are for users
export group_type_table, algebra_type_table
export group_types, algebra_types
export EX, EY
export bigint_EX, bigint_EY

## These are for testing
export palindrome_X, euler_X, nonnegative_Y, nonnegative_X



include("checks.jl")
include("helpers.jl")
include("plevis.jl")

include("grouptypes.jl")
include("algebratypes.jl")

include("epolys.jl")
include("tables.jl")

include("testing.jl")



end # End of module CharacterVarieties
