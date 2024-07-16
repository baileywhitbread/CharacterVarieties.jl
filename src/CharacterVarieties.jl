# module CharacterVarieties

using Chevie

# Export structs


# Export functions


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

###############################################################################
###############################################################################

# Define functions

## Dualising
function old_dual(L::FiniteCoxeterGroup) 
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
end

function dual(L::FiniteCoxeterGroup) 
	# Not sure what I want this function to actually do...
	# Returns the Langlands dual group of L
	# The output will have a parent group iff the input had one
	try
		L_parent = L.parent
		L_parent_dual = rootdatum(simplecoroots(L_parent),simpleroots(L_parent))
		L_parent_dual_plevis = plevis(L_parent_dual)
	catch err
		return rootdatum(simplecoroots(L),simpleroots(L))
	end
end

## Orbits
function myorbit(L::FiniteCoxeterGroup)
	# This is a temporary fix because PermGroups.jl orbit() function creates duplicates
	try
		return reflection_subgroup.(Ref(L.parent),collect(Set(sort.(inclusion.(orbit(L.parent,L))))))
	catch err
		return [G]
	end
end


###############################################################################
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

function islevi(L::FiniteCoxeterGroup)
	# Returns true if L is a Levi subgroup of its parent
	# If no parent, returns true
	try
		return issubset(inclusiongens(L),inclusiongens(L.parent))
	catch err
		return true
	end
end

function ispalindromic(f::Union{Pol{BigInt},Pol{Int64}})
	# Returns true iff f palindromic
	return ((f(0)!= 0) && (f.c == f.c[end:-1:1])) || (f(exp(1))==0)
end

function isnonnegative(f::Union{Pol{BigInt},Pol{Int64}})
	# Returns true iff f has non-negative coefficients
	return length(filter(x -> x<0, f.c)) == 0
end

###############################################################################
###############################################################################

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

###############################################################################
###############################################################################

## Helper functions
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
	# Note the relation on P must be A <= B iff A is a subset of B
	# Throws ArgumentError if the interval [A,B] is empty
	# Only intended to be used with P a set of subgroups of G
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

function nu(L::FiniteCoxeterGroup,G::FiniteCoxeterGroup)
	nu_value = 0
	G_dual = dual(G)
	G_dual_plevis = plevis(G_dual)
	for iplevi in iplevis(G_dual)
		if issubset(dual(L),iplevi)
			nu_value += mobius(dual(L),iplevi,G_dual_plevis)*pi0(iplevi)
		end
	end
	return nu_value
end

###############################################################################
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

function group_type_data(G::FiniteCoxeterGroup)
	d = Array{Any}(nothing,0,8)
	for type in group_types(G)
		type_row = Array{Any}(nothing,1,0)
		type_row = hcat(type_row,[type])
		type_row = hcat(type_row,[Int64(length(roots(type.endoscopy))/2)])
		type_row = hcat(type_row,[orderpol(type.endoscopy)])
		type_row = hcat(type_row,[type.degree])
		type_row = hcat(type_row,[Int64(type.degree(1))])
		type_row = hcat(type_row,[length(type.endoscopy)])
		type_row = hcat(type_row,[length(myorbit(type.endoscopy))])
		type_row = hcat(type_row,[nu(type.endoscopy,G)])		
		d = vcat(d,type_row)
	end
	return sortslices(d,dims=1,by = x -> x[2],rev=true)
end

###############################################################################
###############################################################################

## g-type functions
function algebra_types(G::FiniteCoxeterGroup)
	# Returns a vector of gTypes, ie. the g-types of G
	types = [];
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
		levi_class_sizes = Pol{Rational{Int64}}.(Pol.(levi_xt.cardClass)) # Contains "Mvp"s, we convert to single-variable polys
		levi_greens_functions = Pol{Rational{Int64}}.(Pol.(levi_gt.scalar[1,:])) # Contains "Mvp"s, we convert to single-variable polys
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
		type_row = hcat(type_row,[Pol{Rational{Int64}}(Pol{Rational{Int64}}((type.size))*(orderpol(G)//orderpol(type.levi)))]) # type_row[3] = N size
		type_row = hcat(type_row,[type.green]) # type_row[4] = green
		type_row = hcat(type_row,[length(myorbit(type.levi))]) # type_row[5] = |[L]|
		type_row = hcat(type_row,[mobius(type.levi,type.levi.parent,levis(G))])	# type_row[6] = µ(L,G)
		d = vcat(d,type_row)
	end
	return d
end


function fast_algebra_type_data(G::FiniteCoxeterGroup,type_data)
	d = Array{Any}(nothing,0,6)
	for type in type_data
		type_row = Array{Any}(nothing,1,0)
		type_row = hcat(type_row,[type]) # type_row[1] = type
		type_row = hcat(type_row,[degree(orderpol(type.levi))-degree(type.size)]) # type_row[2] = d(tau)
		type_row = hcat(type_row,[Pol{Rational{Int64}}(Pol{Rational{Int64}}((type.size))*(orderpol(G)//orderpol(type.levi)))]) # type_row[3] = N size
		type_row = hcat(type_row,[type.green]) # type_row[4] = green
		type_row = hcat(type_row,[length(myorbit(type.levi))]) # type_row[5] = |[L]|
		type_row = hcat(type_row,[mobius(type.levi,type.levi.parent,levis(G))])	# type_row[6] = µ(L,G)
		d = vcat(d,type_row)
	end
	return d
end

###############################################################################
###############################################################################

## Build E-polynomial of X
function Mtau(G::FiniteCoxeterGroup,i::Int64)
	# Returns Mτ(q) = q^(|Φ(G)+|-|Φ(L)+|) |L(Fq)|/ρ(1) 
	# where τ = [L,ρ] is the ith GType
	row_data = group_type_data(G)[i,:]
	L_pos_root_size = row_data[2]
	L_size = row_data[3]
	rho_deg = row_data[4]
	return Pol{Rational{Int64}}(Pol(:q)^(Int64(length(roots(G))/2)-L_pos_root_size)*L_size//rho_deg)
end

function Stau(G::FiniteCoxeterGroup,n::Int64,i::Int64)
	# Returns Sτ(q) = |Z(Fq)| * χᵨ(1)^n * |[L]| * (|W|/|W(L)|)^(n-1) * ν(L)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	row_data = group_type_data(G)[i,:]
	chi_rho_deg = row_data[5]
	G_weyl_size = length(G)
	L_weyl_size = row_data[6]
	orbit_size = row_data[7]
	nu_L = row_data[8]
	return Pol{Rational{Int64}}(Z_size * chi_rho_deg^n * orbit_size * (Int64(G_weyl_size//L_weyl_size))^(n-1) * nu_L)
end

function fast_Mtau(G::FiniteCoxeterGroup,i::Int64,type_data)
	# Returns Mτ(q) = q^(|Φ(G)+|-|Φ(L)+|) |L(Fq)|/ρ(1) 
	# where τ = [L,ρ] is the ith GType
	L_pos_root_size = type_data[i,:][2]
	L_size = type_data[i,:][3]
	rho_deg = type_data[i,:][4]
	return Pol{Rational{Int64}}(Pol(:q)^(Int64(length(roots(G))/2)-L_pos_root_size)*L_size//rho_deg)
end

function fast_Stau(G::FiniteCoxeterGroup,n::Int64,i::Int64,type_data)
	# Returns Sτ(q) = |Z(Fq)| * χᵨ(1)^n * |[L]| * (|W|/|W(L)|)^(n-1) * ν(L)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	chi_rho_deg = type_data[i,:][5]
	G_weyl_size = length(G)
	L_weyl_size = type_data[i,:][6]
	orbit_size = type_data[i,:][7]
	nu_L = type_data[i,:][8]
	return Pol{Rational{Int64}}(Z_size * chi_rho_deg^n * orbit_size * (Int64(G_weyl_size//L_weyl_size))^(n-1) * nu_L)
end

function EX(G::FiniteCoxeterGroup,g::Int64,n::Int64)
	# Returns the E-polynomial E(X;q) associated to the group G and a genus g surface with n punctures
	d = group_type_data(G)
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	T_size = orderpol(torus(rank(G)))
	type_sum = 0
	for i in 1:size(d)[1]
		type_sum += fast_Mtau(G,i,d)^(2g-2+n)*fast_Stau(G,n,i,d)
	end
	return (Z_size//T_size^n)*type_sum
end

function fast_EX(G::FiniteCoxeterGroup,g::Int64,n::Int64,type_data)
	# Returns the E-polynomial E(X;q) associated to the group G and a genus g surface with n punctures
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	T_size = orderpol(torus(rank(G)))
	type_sum = 0
	for i in 1:size(type_data)[1]
		type_sum += fast_Mtau(G,i,type_data)^(2g-2+n)*fast_Stau(G,n,i,type_data)
	end
	return Pol{Int64}((Z_size//T_size^n)*type_sum)
end

###############################################################################
###############################################################################

## Build E-polynomial of Y
function qdtau(G::FiniteCoxeterGroup,i::Int64)
	# Returns q^(d(τ)) where τ = [L,ρ] is the ith GType
	row_data = algebra_type_data(G)[i,:]
	return Pol(:q)^(row_data[2])
end

function Htau(G::FiniteCoxeterGroup,n::Int64,i::Int64)
	# Returns Hτ(q) = q^(n|Φ(G)+| + dim(Z)) * (|G(Fq)|/|L(Fq)|) * |N| * Q_L^T(N)^n * |[L]| * (|W|/|W(L)|)^(n-1) * µ(L,G)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	G_weyl_size = BigInt(length(G))
	G_pos_root_size = BigInt(length(G)/2)
	Z_dim = BigInt(rank(G)-semisimplerank(G))
	row_data = algebra_type_data(G)[i,:]
	L_size = orderpol(row_data[1].levi)
	N_size = row_data[1].size
	L_green = row_data[1].green
	orbit_size = row_data[5]
	L_weyl_size = BigInt(length(row_data[1].levi))
	mu_L = row_data[6]
	return Pol{Rational{BigInt}}(Pol(:q)^(n*G_pos_root_size + Z_dim) * Pol{BigInt}(orderpol(G)//L_size) * N_size * L_green^n * BigInt(orbit_size * BigInt((G_weyl_size/L_weyl_size)^(n-1)) * mu_L))
end

function fast_qdtau(G::FiniteCoxeterGroup,i::Int64,type_data)
	# Returns q^(d(τ)) where τ = [L,ρ] is the ith GType
	return Pol{BigInt}(Pol(:q)^(type_data[i,:][2]))
end

function fast_Htau(G::FiniteCoxeterGroup,n::Int64,i::Int64,type_data)
	# Returns Hτ(q) = q^(n|Φ(G)+| + dim(Z)) * (|G(Fq)|/|L(Fq)|) * |N| * Q_L^T(N)^n * |[L]| * (|W|/|W(L)|)^(n-1) * µ(L,G)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	G_weyl_size = BigInt(length(G))
	G_pos_root_size = BigInt(Int64(length(G)/2))
	Z_dim = BigInt(rank(G)-semisimplerank(G))
	L_size = Pol{BigInt}(orderpol(type_data[i,:][1].levi))
	N_size = Pol{Rational{BigInt}}(type_data[i,:][1].size)
	L_green = Pol{Rational{BigInt}}(type_data[i,:][1].green)
	orbit_size = BigInt(type_data[i,:][5])
	L_weyl_size = BigInt(length(type_data[i,:][1].levi))
	mu_L = BigInt(type_data[i,:][6])
	return Pol{Rational{BigInt}}(Pol(:q)^(n*G_pos_root_size + Z_dim) * Pol{BigInt}(orderpol(G)//L_size) * N_size * L_green^n * BigInt(orbit_size * BigInt((G_weyl_size/L_weyl_size)^(n-1)) * mu_L))
end

function EY(G::FiniteCoxeterGroup,g::Int64,n::Int64)
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	d = algebra_type_data(G)
	Z_size = Pol{BigInt}(orderpol(torus(rank(G)-semisimplerank(G))))
	G_size = Pol{BigInt}(orderpol(G))
	g_size = Pol{Int64}(Pol(:q)^(degree(G_size)))

	type_sum = Pol{BigInt}(0)
	for i in 1:size(d)[1]
		type_sum += fast_qdtau(G,i,d)^(g)*fast_Htau(G,n,i,d)
	end
	return Pol{BigInt}((Z_size//G_size) * g_size^(g-1) * type_sum ) 
end

function fast_EY(G::FiniteCoxeterGroup,g::Int64,n::Int64,type_data)
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	Z_size = Pol{BigInt}(orderpol(torus(rank(G)-semisimplerank(G))))
	G_size = Pol{BigInt}(orderpol(G))
	g_size = Pol{Int64}(Pol(:q)^(degree(G_size)))

	type_sum = Pol{BigInt}(0)
	for i in 1:size(type_data)[1]
		type_sum += fast_qdtau(G,i,type_data)^(g)*fast_Htau(G,n,i,type_data)
	end
	return Pol{BigInt}((Z_size//G_size) * g_size^(g-1) * type_sum ) 
end

###############################################################################
###############################################################################

## Display human-readable tables
function group_type_table(G::FiniteCoxeterGroup;summands=false,n=1)
	if summands == false
		d = group_type_data(G)
		num_of_types = size(d)[1]
		# The next lines make the |L(Fq)| and ρ(1) columns readable
		# By converting the entries from Pol to CycPol
		for i in [2,3]
			d[i*num_of_types+1:(i+1)*num_of_types] = CycPol.(d[i*num_of_types+1:(i+1)*num_of_types])
		end
		clabels = ["|Φ(L)+|","|L(Fq)|","ρ(1)","χᵨ(1)","|W(L)|","|[L]|","ν(L)"];
		rlabels = xrepr.(Ref(rio()),d[:,1]); 
		# xrepr(rio(), __ ) is a string of __ when printed on the REPR
		repr_d = xrepr.(Ref(rio()),d[:,2:size(d)[2]]);
		println("A G-type is a W-orbit [L,ρ] where ")
		println("L is an endoscopy group of G containing T")
		println("ρ is a principal unipotent character of L(Fq)")
		println("Φ(L)+ is the set of positive roots of L")
		println("|L(Fq)| is the size of L(Fq)")
		println("ρ(1) is the degree of the unipotent character ρ")
		println("χᵨ(1) is the degree of the Weyl group character associated to ρ")
		println("W(L) is the Weyl group of L")
		println("[L] is the orbit of L under the W-action")
		println("ν(L) is an integer only depending on L")
		println("")
		return showtable(repr_d;col_labels=clabels,rows_label="Types [L,ρ]",row_labels=rlabels)
	elseif summands == true
		d = group_type_data(G)
		num_of_types = size(d)[1]
		# Add Mτ and Sτ to this data
		Mtau_column = Array{Any}(nothing,num_of_types,1)
		Stau_column = Array{Any}(nothing,num_of_types,1)
		for i in 1:num_of_types
			Mtau_column[i] = fast_Mtau(G,i,d)
			Stau_column[i] = fast_Stau(G,n,i,d)
		end
		d = hcat(d,Mtau_column)
		d = hcat(d,Stau_column)
		# The next lines make the |L(Fq)|, ρ(1), Mτ and Sτ columns readable
		# By converting the entries from Pol to CycPol
		for i in [2,3,8,9]
			d[i*num_of_types+1:(i+1)*num_of_types] = CycPol.(d[i*num_of_types+1:(i+1)*num_of_types])
		end
		clabels = ["|Φ(L)+|","|L(Fq)|","ρ(1)","χᵨ(1)","|W(L)|","|[L]|","ν(L)","Mτ","Sτ"];
		rlabels = xrepr.(Ref(rio()),d[:,1]); 
		# xrepr(rio(), __ ) is a string of __ when printed on the REPR
		repr_d = xrepr.(Ref(rio()),d[:,2:size(d)[2]]);
		println("A G-type is a W-orbit [L,ρ] where ")
		println("L is an endoscopy group of G containing T")
		println("ρ is a principal unipotent character of L(Fq)")
		println("Φ(L)+ is the set of positive roots of L")
		println("|L(Fq)| is the size of L(Fq)")
		println("ρ(1) is the degree of the unipotent character ρ")
		println("χᵨ(1) is the degree of the Weyl group character associated to ρ")
		println("W(L) is the Weyl group of L")
		println("[L] is the orbit of L under the W-action")
		println("ν(L) is an integer only depending on L")
		println("Mτ(q) is the q-mass of τ")
		println("Sτ(q) is the character sum of τ")
		println("")
		return showtable(repr_d;col_labels=clabels,rows_label="Types [L,ρ]",row_labels=rlabels)
	else 
		ArgumentError("The optional summands argument must be true or false; default is false")
	end
end


function algebra_type_table(G::FiniteCoxeterGroup;summands=false,g=1,n=1)
	if summands == false
		d = algebra_type_data(G)
		num_of_types = size(d)[1]
		# The next lines make the |N| and Q_T^L(N) columns readable
		# By converting the entries from Pol to CycPol
		for i in [2,3]
			d[i*num_of_types+1:(i+1)*num_of_types] = CycPol.(d[i*num_of_types+1:(i+1)*num_of_types])
		end
		clabels = ["d(τ)","|N|","Q_T^L(N)","|[L]|","µ(L,G)"];
		rlabels = xrepr.(Ref(rio()),d[:,1]); 
		# xrepr(rio(), __ ) is a string of __ when printed on the REPR
		repr_d = xrepr.(Ref(rio()),d[:,2:size(d)[2]]);
		println("A g-type is a W-orbit [L,N] where ")
		println("L is a Levi subgroup of G containing T")
		println("N is an L(Fq)-orbit of a nilpotent element of Lie(L)(Fq)")
		println("d(τ) is dim(L)-degree(|N|)")
		println("Q_T^L(N) is the Greens function evaluated at the unipotent L(Fq)-conjugacy class associated to N")
		println("[L] is the orbit of L under the W-action")
		println("µ(L,G) is the Mobius function of the poset of Levi subgroups of G containing T")
		println("")
		return showtable(repr_d;col_labels=clabels,rows_label="Types [L,N]",row_labels=rlabels)
	elseif summands == true
		d = algebra_type_data(G)
		num_of_types = size(d)[1]
		# Add qdτ and Hτ to this data
		qdtau_column = Array{Any}(nothing,num_of_types,1)
		Htau_column = Array{Any}(nothing,num_of_types,1)
		for i in 1:num_of_types
			qdtau_column[i] = fast_qdtau(G,i,d)
			Htau_column[i] = fast_Htau(G,n,i,d)
		end
		d = hcat(d,qdtau_column)
		d = hcat(d,Htau_column)
		# The next lines make the |N| and Q_T^L(N) columns readable
		# By converting the entries from Pol to CycPol
		for i in [2,3]
			d[i*num_of_types+1:(i+1)*num_of_types] = CycPol.(d[i*num_of_types+1:(i+1)*num_of_types])
		end
		clabels = ["d(τ)","|N|","Q_T^L(N)","|[L]|","µ(L,G)","qdτ","Hτ"];
		rlabels = xrepr.(Ref(rio()),d[:,1]); 
		# xrepr(rio(), __ ) is a string of __ when printed on the REPR
		repr_d = xrepr.(Ref(rio()),d[:,2:size(d)[2]]);
		println("A g-type is a W-orbit [L,N] where ")
		println("L is a Levi subgroup of G containing T")
		println("N is an L(Fq)-orbit of a nilpotent element of Lie(L)(Fq)")
		println("d(τ) is dim(L)-degree(|N|)")
		println("Q_T^L(N) is the Greens function evaluated at the unipotent L(Fq)-conjugacy class associated to N")
		println("[L] is the orbit of L under the W-action")
		println("µ(L,G) is the Mobius function of the poset of Levi subgroups of G containing T")
		println("qdτ is q^d(τ)")
		println("Hτ(q) is the character sum of τ")
		println("")
		return showtable(repr_d;col_labels=clabels,rows_label="Types [L,N]",row_labels=rlabels)
	else 
		ArgumentError("The optional summands argument must be true or false; default is false")
	end
end

###############################################################################
###############################################################################

## Testing E-polynomials
function palindrome_X(G,genus,puncture_min,puncture_max)
	# Checks palindromicity of X with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking E(X;q) palindromic when g=",genus," and n=",n,": ")
			if ispalindromic(fast_EX(G,genus,n,d))
				println("Yes")
			else
				println("No")
			end
		catch err
			if isa(err,OverflowError)
				println("Overflow error")
				break
			elseif isa(err,ErrorException)
				println("Error exception")
				break
			else
				println(err)
				break
			end
		end
	end
end

function euler_X(G,genus,puncture_min,puncture_max)
	# Checks Euler characteristic of X with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Computing χ(X) when g=",genus," and n=",n,": ")
			if fast_EX(G,genus,n,d)(1)==0
				println("χ(X)=0")
			else
				println("χ(X)=",fast_EX(G,genus,n,d)(1))
			end
		catch err
			if isa(err,OverflowError)
				println("Overflow error")
				break
			elseif isa(err,ErrorException)
				println("Error exception")
				break
			else
				println(err)
				break
			end
		end
	end
end

function nonnegative_Y(G,genus,puncture_min,puncture_max)
	# Checks negativity of coefficients of E(Y;q) with g=genus and n=puncture_min,...,puncture_max
	d=algebra_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking coefficients of E(Y;q) when g=",genus," and n=",n,": ")
			if isnonnegative(fast_EY(G,genus,n,d))
				println("All non-negative")
			else
				println("Negative coefficients found")
			end
		catch err
			if isa(err,OverflowError)
				println("Overflow error")
				break
			elseif isa(err,ErrorException)
				println("Error exception")
				break
			else
				println(err)
				break
			end
		end
	end
end



# end # End of module CharacterVarieties
