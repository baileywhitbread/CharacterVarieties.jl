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
	G_plevis = plevis(G)
	for iplevi in iplevis(G)
		if issubset(L,iplevi)
			nu_value += mobius(L,iplevi,G_plevis)*pi0(iplevi)
		end
	end
	return nu_value
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

function group_type_data(G::FiniteCoxeterGroup)
	d = Array{Any}(nothing,0,8)
	for gtype in group_types(G)
		type_row = Array{Any}(nothing,1,0)
		type_row = hcat(type_row,[gtype])
		type_row = hcat(type_row,[Int64(length(roots(gtype.endoscopy))/2)])
		type_row = hcat(type_row,[orderpol(gtype.endoscopy)])
		type_row = hcat(type_row,[gtype.degree])
		type_row = hcat(type_row,[Int64(gtype.degree(1))])
		type_row = hcat(type_row,[length(gtype.endoscopy)])
		type_row = hcat(type_row,[length(myorbit(gtype.endoscopy))])
		type_row = hcat(type_row,[nu(gtype.endoscopy)])		
		d = vcat(d,type_row)
	end
	return sortslices(d,dims=1,by = x -> x[2],rev=true)
end

###############################################################################

## Build E-polynomial
function Mtau(G::FiniteCoxeterGroup,i::Int64)
	# Returns Mτ(q) = q^(|Φ(G)+|-|Φ(L)+|) |L(Fq)|/ρ(1) 
	# where τ = [L,ρ] is the ith GType
	row_data = group_type_data(G)[i,:]
	L_pos_root_size = row_data[2]
	L_size = row_data[3]
	rho_deg = row_data[4]
	return Pol(:q)^(Int64(length(roots(G))/2)-L_pos_root_size)*L_size//rho_deg
end

function Stau(G::FiniteCoxeterGroup,n::Int64,i::Int64)
	# Returns S_τ(q) = |Z(Fq)| * χᵨ(1)^n * |[L]| * (|W|/|W(L)|)^(n-1) * ν(L)
	# where τ = [L,ρ] is the ith GType
	# and n is the number of punctures
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	row_data = group_type_data(G)[i,:]
	chi_rho_deg = row_data[5]
	G_weyl_size = length(G)
	L_weyl_size = row_data[6]
	orbit_size = row_data[7]
	nu_L = row_data[8]
	return Z_size * chi_rho_deg^n * orbit_size * (Int64(G_weyl_size//L_weyl_size))^(n-1) * nu_L
end

function EX(G::FiniteCoxeterGroup,g::Int64,n::Int64)
	# Returns the E-polynomial E(X;q) associated to the group G and a genus g surface with n punctures
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	T_size = orderpol(torus(rank(G)))
	type_sum = 0
	for i in 1:length(group_types(G))
		type_sum += Mtau(G,i)^(2g-2+n)*Stau(G,n,i)
	end
	return Pol{Int64}((Z_size//T_size^n)*type_sum)
end

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
			Mtau_column[i] = Mtau(G,i)
			Stau_column[i] = Stau(G,n,i)
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

###############################################################################


## Testing counting functions
function is_palindromic(f)
	return ((f(0)!= 0) && (f.c == f.c[end:-1:1])) || (f(exp(1))==0)
end

function palindrome_X(G,genus_max,puncture_max)
	# Checks palindromicity of X with g=0,1,2,..,genus_max and n=1,2,...,puncture_max
	for g in 0:genus_max
		for n in 1:puncture_max
			try 
				if is_palindromic(EX(G,g,n))
					println("EX palindromic when g=",g," and n=",n)
				else
					println("EX not palindromic when g=",g," and n=",n)
				end
			catch err
				if isa(err,OverflowError)
					println("Overflow error when g=",g," and n=",n)
				else
					println(err," when g=",g," and n=",n)
				end
			end
		end
	end
end

function euler_X(G,genus_max,puncture_max)
	# Checks Euler characteristic of X with g=0,1,2,..,genus_max and n=1,2,...,puncture_max
	for g in 0:genus_max
		for n in 1:puncture_max
			try 
				if EX(G,g,n)(1)==0
					println("χ(X)=0 when g=",g," and n=",n)
				elseif EX(G,g,n)(1)!=0
					println("χ(X) non-zero when g=",g," and n=",n)
				end
			catch err
				if isa(err,OverflowError)
					println("Overflow error when g=",g," and n=",n)
				else
					println(err," when g=",g," and n=",n)
				end
			end
		end
	end
end



# end # End of module CharacterVarieties
