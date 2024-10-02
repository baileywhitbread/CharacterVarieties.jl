struct gType
	# A g-type is a pair [L,N] where 
	# L is a Levi subgroup of G containing T
	# N is an L(Fq)-orbit of a nilpotent element of L(Fq)
	
	# To record N, we record a string representation and its size
	# eg. the regular nilpotent orbit of GL3 is recorded as 
	# "3" (for the partition (3,0,...) of 3) and q⁶-q⁴-q³+q
	
	levi::FiniteCoxeterGroup
	orbit::String
	size::Pol{Rational{BigInt}}
	green::Pol{Rational{BigInt}}
	
end # End of struct GType

# Make gTypes display nicely on the REPL
Base.show(io::IO, tau::gType) = print(io,
"[",tau.levi,",",tau.orbit,"]"
)


## g-type functions
function algebra_types(G::FiniteCoxeterGroup)
	# Returns a vector of gTypes, ie. the g-types of G
	types = gType[];
	for levi in lorbit_reps(G)
		levi_uc = UnipotentClasses(levi); # This is a dictionary object, not the actual classes
		levi_xt = XTable(levi_uc;classes=true) # Another dictionary object containing cardinality of centralisers and conjugacy classes
		levi_gt = GreenTable(levi_uc;classes=true) # Another dictionary containing Greens functions
		levi_uc_classes = levi_uc.classes; # Grab unipotent classes over algebraic closure
		levi_rational_orbit_labels = levi_xt.classes # Grab rational unipotent classes using XTable
		# These keep track of rational orbits: they are pairs of integers [n,m] where 
		# n is the nth unipotent class over algebraic closure and
		# m is the mth component after taking Fq points
		# eg. If G=G2 then G2(a1) is the 4th class over algebraic closure which splits into 3 after taking Fq points
		# so the rational orbits associated to G2(a1) are represented by [4,1], [4,2] and [4,3]
		levi_rational_orbit_TeX_names = map(label->name(TeX(rio();class=label[2]),levi_uc_classes[label[1]]),levi_rational_orbit_labels)
		levi_rational_orbit_names = fromTeX.(Ref(rio()),levi_rational_orbit_TeX_names)
		levi_class_sizes = Pol{Rational{BigInt}}.(Pol.(levi_xt.cardClass)) # Contains "Mvp"s, we convert to single-variable polys
		levi_greens_functions = Pol{Rational{BigInt}}.(Pol.(levi_gt.scalar[1,:])) # Contains "Mvp"s, we convert to single-variable polys
		for i in 1:length(levi_rational_orbit_labels)
			append!(types,[gType(levi,levi_rational_orbit_names[i],levi_class_sizes[i],levi_greens_functions[i])])
		end
	end
	return types
end

function algebra_type_data(G::FiniteCoxeterGroup)
	d = Array{Any}(nothing,0,6)
	for type in algebra_types(G)
		type_row = Array{Any}(nothing,1,0)
		type_row = hcat(type_row,[type]) # type_row[1] = type
		type_row = hcat(type_row,[degree(orderpol(type.levi))-degree(type.size)]) # type_row[2] = d(tau)
		type_row = hcat(type_row,[Pol{Rational{BigInt}}(type.size)]) # type_row[3] = N size
		type_row = hcat(type_row,[type.green]) # type_row[4] = green
		type_row = hcat(type_row,[length(myorbit(type.levi))]) # type_row[5] = |[L]|
		type_row = hcat(type_row,[mobius(type.levi,type.levi.parent,levis(G))])	# type_row[6] = µ(L,G)
		d = vcat(d,type_row)
	end
	return d
end
