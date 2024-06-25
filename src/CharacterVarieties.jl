# module CharacterVarieties

using Chevie

# Export structs
export GType

# Export functions
export dual 
export plorbit_reps, plorbits, plevis
export iplorbit_reps, iplorbits, iplevis

###############################################################################
###############################################################################

# Define structs
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

###############################################################################
###############################################################################

# Define functions

## Dualising
function dual(L::FiniteCoxeterGroup) 
	# Not sure what I want this function to actually do...
	# Returns the Langlands dual group of L
	# The output will have a parent group iff the input had one
	try
		L_parent = L.parent
		L_parent_dual = rootdatum(simplecoroots(L_parent),simpleroots(L_parent))
		return reflection_subgroup(L_parent_dual,inclusiongens(L))
	catch err
		return rootdatum(simplecoroots(L),simpleroots(L))
	end
End

###############################################################################

## Checks
function isisolated(L::FiniteCoxeterGroup)
	# Returns true if L is isolated in its parent
	# If no parent, returns true
	try
		return length(gens(L)) == length(gens(L.parent))
	catch err
		return true
	end
end

###############################################################################

## Calculating pseudo-Levis
function plorbit_reps(G::FiniteCoxeterGroup)
	# Returns pseudo-Levi orbit representatives as a vector of FiniteCoxeterSubGroup's
	return reflection_subgroup.(Ref(G),sscentralizer_reps(G))
end

function iplorbit_reps(G::FiniteCoxeterGroup)
	# Returns isolated pseudo-Levi orbit representatives as a vector of FiniteCoxeterSubGroup's
	return filter(isisolated,plorbit_reps(G))
end

function plorbits(G::FiniteCoxeterGroup)
	# Returns pseudo-Levi orbits as a vector of vectors of FiniteCoxeterSubGroup's
	# orbits(G,plorbit_reps(G)) is the obvious solution but it creates duplicates for some reason
	# Duplicates are killed by converting the vector to a set then back to a vector
	map(orbits(G,plorbit_reps(G))) do x
		return reflection_subgroup.(Ref(G),collect(Set(sort.(inclusion.(x)))))
	end
end

function iplorbits(G::FiniteCoxeterGroup)
	# Returns isolated pseudo-Levi orbits as a vector of vectors of FiniteCoxeterSubGroup's
	map(orbits(G,iplorbit_reps(G))) do x
		return reflection_subgroup.(Ref(G),collect(Set(sort.(inclusion.(x)))))
	end
end

function plevis(G::FiniteCoxeterGroup)
	# Returns all pseudo-Levis as a vector of FiniteCoxeterSubGroup's
	return reduce(vcat,plorbits(G))
end

function iplevis(G::FiniteCoxeterGroup)
	# Returns all pseudo-Levis as a vector of FiniteCoxeterSubGroup's
	return reduce(vcat,iplorbits(G))
end	

###############################################################################

## Calculations
function orderpol(L::FiniteCoxeterGroup)
	# Returns ||L||(q) = |L(Fq)|
	return PermRoot.generic_order(L,Pol(:q))
end

function pi0(L::FiniteCoxeterGroup)
	# Returns |pi_0(Z(L))|
	return length(algebraic_center(L).AZ)
end

function mobius(A::FiniteCoxeterGroup,B::FiniteCoxeterGroup,P::Vector)
	# Returns the Mobius function of the poset P evaluated at (A,B) in PxP
	# Throws ArgumentError if the interval [A,B] is empty
	# Only intended to be used with P being a poset of pseudo-Levi subgroups
	if isequal(A,B)
		return 1
	elseif issubset(A,B)
		mobius_value = 0
		for element in P
			if issubset(A,element) && issubset(element,B) && !isequal(element,B)
				mobius_value += mobius(A,element,P)
			end
		end
		return (-1)*mobius_value
	else
		ArgumentError("First argument is not contained in the second argument")
	end
end

function nu(L::FiniteCoxeterGroup)
	nu_value = 0
	try
		L_parent = L.parent
		G_plevis = plevis(L_parent)
		G_iplevis = iplevis(L_parent)
		for iplevi in G_iplevis
			if issubset(L,iplevi)
				nu_value += mobius(L,iplevi,G_plevis)*pi0(iplevi)
			end
		end
		return nu_value
	catch err
		G_plevis = plevis(L_parent)
		G_iplevis = iplevis(L_parent)
		for iplevi in G_iplevis
			if issubset(L,iplevi)
				nu_value += mobius(L,iplevi,G_plevis)*pi0(iplevi)
			end
		end
		return nu_value
	end
end

###############################################################################

## G-type functions
function group_types(G::FiniteCoxeterGroup)
	# Returns a vector of GTypes, ie. the G-types of G
	types = [];
	for plevi in plorbit_reps(G)
		plevi_uc = UnipotentCharacters(plevi)
		plevi_uc_names = charnames(plevi_uc,limit=true)
		plevi_uc_degs = degrees(plevi_uc)
		for i in 1:length(plevi_uc)
			# Check if unipotent character is principal
			if Int(plevi_uc_degs[i](1))!=0
				append!(types,[GType(plevi,plevi_uc_names[i],plevi_uc_degs[i])])
			end
		end
	end
	return types
end

function type_data(G::FiniteCoxeterGroup)
	d = Array{Any}(nothing,0,8)
	for gtype in group_types(G)
		type_row = Array{Any}(nothing,1,0)
		type_row = hcat(type_row,[gtype])
		type_row = hcat(type_row,[Int64(length(roots(gtype.endoscopy))/2)])
		type_row = hcat(type_row,[orderpol(gtype.endoscopy)])
		type_row = hcat(type_row,[gtype.degree])
		type_row = hcat(type_row,[Int64(gtype.degree(1))])
		type_row = hcat(type_row,[length(gtype.endoscopy)])
		type_row = hcat(type_row,["Orbit size"])
		type_row = hcat(type_row,[nu(gtype.endoscopy)])		
		d = vcat(d,type_row)
	end
	return sortslices(d,dims=1,by = x -> x[2],rev=true)
end
# end # End of module CharacterVarieties
