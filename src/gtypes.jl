struct GType
	# A G-type is a pair [L,ρ] where 
	# L is an endoscopy group of G containing T
	# ρ is a principal unipotent character of L(Fq)
	
	# To record ρ, we record a string representation and its degree
	# eg. the Steinberg character of GL3 is recorded as 
	# "3" (for the partition (3,0,...) of 3) and Pol(:q)^3
	
	endoscopy::FiniteCoxeterGroup
	character::String
	degree::Pol{Rational{Int64}}
	
end # End of struct GType

# Make GTypes display nicely on the REPL
Base.show(io::IO, tau::GType) = print(io,
"[",tau.endoscopy,",",tau.character,"]"
)

struct gType
	# A g-type is a pair [L,N] where 
	# L is a Levi subgroup of G containing T
	# N is an L(Fq)-orbit of a nilpotent element of L(Fq)
	
	# To record N, we record a string representation and its size
	# eg. the regular nilpotent orbit of GL3 is recorded as 
	# "3" (for the partition (3,0,...) of 3) and q⁶-q⁴-q³+q
	
	levi::FiniteCoxeterGroup
	orbit::String
	size::Pol{Rational{Int64}}
	green::Pol{Rational{Int64}}
	
end # End of struct GType

# Make gTypes display nicely on the REPL
Base.show(io::IO, tau::gType) = print(io,
"[",tau.levi,",",tau.orbit,"]"
)