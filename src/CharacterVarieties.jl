module CharacterVarieties

using Chevie

# Export structs
export GType

# Export functions
export dual 
export plorbit_reps, plorbits, plevis
export iplorbit_reps, iplorbits, iplevis





# Define structs
struct GType
	endoscopy::FiniteCoxeterGroup
	character::String
	degree::Pol{Rational{Int64}}
end # End of struct GType








# Define functions
function dual(G::FiniteCoxeterGroup)
	return rootdatum(simplecoroots(G),simpleroots(G))
end

function plorbit_reps(G::FiniteCoxeterGroup)
	return reflection_subgroup.(Ref(G),sscentralizer_reps(G))
end

function iplorbit_reps(G::FiniteCoxeterGroup)
	iplorbit_reps = [];
	for plorbit_rep in reflection_subgroup.(Ref(G),sscentralizer_reps(G))
		if length(gens(plorbit_rep)) == length(gens(G))
			append!(iplorbit_reps,[plorbit_rep])
		end
	end
	return iplorbit_reps
end

function plorbits(G::FiniteCoxeterGroup)
	# orbits(G,plorbit_reps(G)) is the obvious solution but it creates duplicates
	# Duplicates are killed by converting the vector to a set then back to a vector
	return map(plorbit -> reflection_subgroup.(Ref(G),collect(Set(sort.(inclusion.(plorbit))))), orbits(G,plorbit_reps(G)))
end

function iplorbits(G::FiniteCoxeterGroup)
	return map(plorbit -> reflection_subgroup.(Ref(G),collect(Set(sort.(inclusion.(plorbit))))), orbits(G,iplorbit_reps(G)))
end


function plevis(G::FiniteCoxeterGroup)
	plevis = [];
	iplevis = [];
	for plorbit in plorbits(G)
		for plevi in plorbit
			append!(plevis,[plevi])
		end
	end
	return plevis
end

function iplevis(G::FiniteCoxeterGroup)
	iplevis = [];
	for iplorbit in iplorbits(G)
		for iplevi in iplorbit
			append!(iplevis,[iplevi])
		end
	end
	return iplevis
end

function orderpol(L::FiniteCoxeterGroup)
	return PermRoot.generic_order(L,Pol(:q))
end

function pi0(L::FiniteCoxeterGroup)
	return length(algebraic_center(L).AZ)
end

function subset(L::FiniteCoxeterGroup,M::FiniteCoxeterGroup)
	return issubset(inclusion(L),inclusion(M))
end

function equal(L::FiniteCoxeterGroup,M::FiniteCoxeterGroup)
	return subset(L,M) && subset(M,L)
end

function mobius(A::FiniteCoxeterGroup,B::FiniteCoxeterGroup,poset::Vector)
	if equal(A,B)
		return 1
	elseif subset(A,B)
		mobius_value = 0
		for element in poset
			if subset(A,element) && subset(element,B) && !equal(element,B)
				mobius_value += mobius(A,element,poset)
			end
		end
		return (-1)*mobius_value
	else
		error("First argument is not contained in the second argument")
	end
end

function nu(L::FiniteCoxeterGroup)
	nu_value = 0
	G_plevis = plevis(L.parent)
	G_iplevis = iplevis(L.parent)
	for iplevi in G_iplevis
		if subset(L,iplevi)
			nu_value += mobius(L,iplevi,G_plevis)*pi0(iplevi)
		end
	end
	return nu_value
end







end # End of module CharacterVarieties



