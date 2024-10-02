## Helper functions
function myorbit(L::FiniteCoxeterGroup)
	# This is a temporary fix because PermGroups.jl orbit() function creates duplicates
	try
		return reflection_subgroup.(Ref(L.parent),collect(Set(sort.(inclusion.(orbit(L.parent,L))))))
	catch err 
		# I am assuming L has no parent iff L=G
		return [G]
	end
end

function orderpol(L::FiniteCoxeterGroup)
	# Returns ||L||(q) = |L(Fq)|
	return PermRoot.generic_order(L,Pol(:q))
end

function mobius(A::FiniteCoxeterGroup,B::FiniteCoxeterGroup,P::Vector)
	# Returns the Mobius function of the poset P evaluated at (A,B) in PxP
	# Note the relation on P must be A <= B iff A is a subset of B
	# Throws ArgumentError if the interval [A,B] is empty
	# Only intended to be used with P a set of subgroups of G
	if isequal(sort(inclusion(A)),sort(inclusion(B)))
		return 1
	elseif issubset(sort(inclusion(A)),sort(inclusion(B)))
		mobius_value = 0
		for element in P
			if issubset(sort(inclusion(A)),sort(inclusion(element))) && issubset(sort(inclusion(element)),sort(inclusion(B))) && !isequal(sort(inclusion(element)),sort(inclusion(B)))
				mobius_value += mobius(A,element,P)
			end
		end
		return (-1)*mobius_value
	else
		ArgumentError("First argument is not contained in the second argument")
	end
end

function pi0(L::FiniteCoxeterGroup)
	# Returns |pi_0(Z(L))|
	return length(algebraic_center(L).AZ)
end

function nu(L::FiniteCoxeterGroup,iso_plevis::Vector,all_plevis::Vector)
	nu_value = 0
	for iplevi in iso_plevis
		if issubset(L,iplevi)
			nu_value += mobius(L,iplevi,all_plevis)*pi0(iplevi)
		end
	end
	return nu_value
end

function dimension_XY(G::FiniteCoxeterGroup,genus::Number,puncture::Number)
	return (2*genus-2+puncture)*dimension(G) + 2*(rank(G)-semisimplerank(G)) - puncture*rank(G)
end