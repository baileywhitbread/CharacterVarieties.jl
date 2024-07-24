## Calculating Levi orbit representatives
function plorbit_reps(G::FiniteCoxeterGroup)
	# Returns pseudo-Levi orbit representatives as a vector of FiniteCoxeterSubGroup's
	return reflection_subgroup.(Ref(G),sscentralizer_reps(G))
end

function iplorbit_reps(G::FiniteCoxeterGroup)
	# Returns isolated pseudo-Levi orbit representatives as a vector of FiniteCoxeterSubGroup's
	return filter(isisolated,plorbit_reps(G))
end

function lorbit_reps(G::FiniteCoxeterGroup)
	# Returns Levi orbit representatives as a vector of FiniteCoxeterSubGroup's
	return filter(islevi,plorbit_reps(G))
end





## Calculating Levi orbits
function plorbits(G::FiniteCoxeterGroup)
	# Returns pseudo-Levi orbits as a vector of vectors of FiniteCoxeterSubGroup's
	# orbits(G,plorbit_reps(G)) is obvious solution but creates duplicates for some reason
	# Kill dupes by converting vector to set and back to vector
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

function lorbits(G::FiniteCoxeterGroup)
	# Returns Levi orbits as a vector of vectors of FiniteCoxeterSubGroup's
	map(orbits(G,lorbit_reps(G))) do x
		return reflection_subgroup.(Ref(G),collect(Set(sort.(inclusion.(x)))))
	end
end





## Calculating all Levis
function plevis(G::FiniteCoxeterGroup)
	# Returns all pseudo-Levis as a vector of FiniteCoxeterSubGroup's
	return reduce(vcat,plorbits(G))
end

function iplevis(G::FiniteCoxeterGroup)
	# Returns all isolated pseudo-Levis as a vector of FiniteCoxeterSubGroup's
	return reduce(vcat,iplorbits(G))
end	

function levis(G::FiniteCoxeterGroup)
	# Returns all Levis as a vector of FiniteCoxeterSubGroup's
	return reduce(vcat,lorbits(G))
end	