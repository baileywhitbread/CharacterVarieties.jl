module CharacterVarieties

using Reexport
@reexport using Chevie





# Export structs
export GType, gType

# Export functions
## These are for users
export group_type_table, algebra_type_table
export EX, EY

## These are for testing
export palindrome_X, euler_X, nonnegative_Y, nonnegative_X
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
