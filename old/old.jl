# The purpose of this script is to compute E-polynomials associated to multiplicative and additive character varieties

# Load Jean Michel's package
using Chevie

#################################################################################
#################################################################################
#################################################################################

# Helper functions
function orderpol(G)
	return PermRoot.generic_order(G,Pol(:q))
end

function pi0(L)
	return length(algebraic_center(L).AZ)
end

function subset(L,M)
	return issubset(inclusion(L),inclusion(M))
end

function equal(L,M)
	return subset(L,M) && subset(M,L)
end

function mob(A,B,poset)
	if equal(A,B)
		return 1
	elseif subset(A,B)
		mob_value = 0
		for element in poset
			if subset(A,element) && subset(element,B) && !equal(element,B)
				mob_value += mob(A,element,poset)
			end
		end
		return (-1)*mob_value
	else
		error("First argument must be a subset of the second argument")
	end
end

function nu(L,iplevis,plevis)
	nu_value = 0
	for iplevi in iplevis
		if subset(L,iplevi)
			nu_value += mob(L,iplevi,plevis)*pi0(iplevi)
		end
	end
	return nu_value
end

#################################################################################
#################################################################################
#################################################################################

# Functions for computing E-polynomials
# The two important functions are EX(G,g,n) and EY(G,g,n)



# First let's create EX(G,g,n)
function gptypes(G)
	# This creates a matrix of data (called gptypes) required to compute E-polynomial
	# gptypes[i,:][1] = ith G-type, ie. τ=[L,ρ]
	# gptypes[i,:][2] = |Φ(L)+| of ith type
	# gptypes[i,:][3] = ρ(1) of ith type (unipotent character degree)
	# gptypes[i,:][4] = |L(Fq)| of ith type
	# gptypes[i,:][5] = χᵨ(1) of ith type (Weyl gp character degree)
	# gptypes[i,:][6] = |W(L)| of ith type
	# gptypes[i,:][7] = |[L]| of ith type
	# gptypes[i,:][8] = ν(L) of ith type
	
	# Gather required info about G
	G_dual = rootdatum(simplecoroots(G),simpleroots(G));
	G_positive_root_size = Int(length(roots(G))/2);
	
	# Compute pseudo Levi orbit representatives and pseudo Levi orbits
	plorbit_reps = reflection_subgroup.(Ref(G_dual),sscentralizer_reps(G_dual)); 
	plorbits = orbits(G_dual,plorbit_reps); 
	# A PROBLEM: |[G]| > 1 IFF G=SO5
	# THEREFORE WE KILL THE DUPLICATE THAT OCCURS
	# Probably a better way of fixing this
	if G==rootdatum(:so,5)
		plorbits[4] = [plorbits[4][1]];
	end
	
	# Compute all pseudo Levis and all isolated pseudo Levis	
	plevis = [];
	iplevis = [];
	for plorbit in plorbits
		for plevi in plorbit
			append!(plevis,[plevi])
			# Check is pseudo Levi is isolated
			if length(gens(plevi)) == length(gens(G))
				append!(iplevis,[plevi])
			end
		end
	end
	
	# Create gptypes
	gptypes = Array{Any}(nothing,0,8);
	for plevi in plorbit_reps
		plevi_order = orderpol(plevi);
		plevi_positive_root_size = Int(length(roots(plevi))/2);
		plevi_orbit_size = length(orbit(G_dual,plevi)); 
		# Above line's calculation has already been performed, room for optimisation here
		plevi_weyl_size = length(plevi);
		plevi_nu = nu(plevi,iplevis,plevis);
		plevi_uc = UnipotentCharacters(plevi);
		plevi_uc_names = charnames(plevi_uc,limit=true);
		plevi_uc_degs = degrees(plevi_uc);
		# Pick unipotent character
		for i in 1:length(plevi_uc)
			# Check if unipotent character is principal
			if Int(plevi_uc_degs[i](1))!=0
				gptypes = vcat(gptypes,
				[(plevi,plevi_uc_names[i]) plevi_positive_root_size plevi_uc_degs[i] plevi_order Int(plevi_uc_degs[i](1)) plevi_weyl_size plevi_orbit_size plevi_nu]
				);
			end
		end
	end
	return gptypes
end


function gp_row_term(G,i,genus_num,puncture_num)
	# Returns (||Z||(q)/||T||(q)^n) * ||τ||(q)^(2g-2+n) * S_τ(q)
	# where τ is the ith type, g = genus_num and n = puncture_num
	
	# Grab necessary data
	G_rank = rank(G);
	G_ssrank = semisimplerank(G);
	T = torus(G_rank);
	Z = torus(G_rank-G_ssrank);
	# The above four lines could be quicker
	G_positive_root_size = Int(length(roots(G))/2);
	G_weyl_group_size = length(G);
	gptype_data = gptypes(G)
	# For readability:
	# gptypes[i,:][1] = ith G-type, ie. τ=[L,ρ]
	# gptypes[i,:][2] = |Φ(L)+| of ith type
	# gptypes[i,:][3] = ρ(1) of ith type (unipotent character degree)
	# gptypes[i,:][4] = |L(Fq)| of ith type
	# gptypes[i,:][5] = χᵨ(1) of ith type (Weyl gp character degree)
	# gptypes[i,:][6] = |W(L)| of ith type
	# gptypes[i,:][7] = |[L]| of ith type
	# gptypes[i,:][8] = ν(L) of ith type
	
	coeff = orderpol(Z)//(orderpol(T)^puncture_num);
	tau_q = ((Pol(:q)^(G_positive_root_size-gptype_data[i,:][2]))*gptype_data[i,:][4])//gptype_data[i,:][3];
	S_tau_q = orderpol(Z)*(gptype_data[i,:][5]^puncture_num)*gptype_data[i,:][7]*((G_weyl_group_size//gptype_data[i,:][6])^(puncture_num-1))*gptype_data[i,:][8]
	return coeff*(tau_q)^(2*genus_num-2+puncture_num)*S_tau_q
end

function EX(G,genus_num,puncture_num)
	# Returns the E-polynomial of the character variety 
	# associated to (G,g,n) where g = genus_num and n = puncture_num
	gptype_data = gptypes(G);
	epol = Pol([0]);
	for row_number in 1:size(gptype_data)[1]
		epol += gp_row_term(G,row_number,genus_num,puncture_num)
	end
	return Pol{Int64}(epol)
end

### Integer overflows occur if G,g,n are big
### Overflow occurs when G=SO7 and g=n=2 but not when G=SO7, g=1 and n=2
### Overflow occurs when G=G2 and g=n=2 but not when G=G2, g=1 and n=2
### Overflow does not occur when G=GL4 and g=n=5 but does when g>5 and n>5
### Solution: Replace Pol{Int64} with Pol{Int128} or something larger, or rewrite computations


# And now let's create EY(G,g,n)
function algtypes(G)
	# This creates a matrix of data (called gptypes) required to compute E-polynomial
	# algtypes[i,:][1] = ith g-type, ie. τ=[L,Oᴸ(Fq)]
	# algtypes[i,:][2] = dim(τ) of ith type
	# algtypes[i,:][3] = |Oᴳ(Fq)| of ith type
	# algtypes[i,:][4] = Qᴸ(Oᴸ(Fq)) of ith type
	# algtypes[i,:][5] = |[L]| of ith type
	# algtypes[i,:][6] = µ(L,G) of ith type
	
	# Gather required info about G
	G_orderpol = orderpol(G)
	G_dual = rootdatum(simplecoroots(G),simpleroots(G));
	
	# Compute Levis in G_dual
	lorbit_reps = map(L -> L.W, split_levis(G_dual));
	# A PROBLEM: TWO LORBITS REPS OF G IFF G = GLN
	# THEREFORE WE KILL THE DUPLICATE THAT OCCURS
	# Definitely a better way of fixing this
	if G.TeXcallname[1:2] == "gl"
		lorbit_reps = lorbit_reps[2:end]
	end
	lorbits = orbits(G_dual,lorbit_reps)
	levis = [];
	for lorbit in lorbits
		for levi in lorbit
			append!(levis,[levi])
		end
	end
	
	# Create algtypes
	algtypes = Array{Any}(nothing,0,6);
	for levi in lorbit_reps		
		# First we compute g-types
		# Create required dictionary objects containing the data we need
		levi_uc = UnipotentClasses(levi); # This is a dictionary object, not the actual classes
		levi_xt = XTable(levi_uc;classes=true) # Another dictionary object containing cardinality of centralisers and conjugacy classes
		levi_gt = GreenTable(levi_uc;classes=true); # Another dictionary containing Greens functions
		
		# Grab unipotent classes over algebraic closure
		levi_uc_classes = levi_uc.classes;
		
		# Grab rational unipotent classes using XTable
		levi_rational_orbit_labels = levi_xt.classes 
		# These keep track of rational orbits: they are pairs of integers [n,m] where 
		# n is the nth unipotent class over algebraic closure and
		# m is the mth component after taking Fq points
		# eg. If G=G2 then G2(a1) is the 4th class over algebraic closure which splits into 3 after taking Fq points
		# so the rational orbits associated to G2(a1) are represented by [4,1], [4,2] and [4,3]
		levi_rational_orbit_TeX_names = map(label->name(TeX(rio();class=label[2]),levi_uc_classes[label[1]]),levi_rational_orbit_labels)
		levi_rational_orbit_names = fromTeX.(Ref(rio()),levi_rational_orbit_TeX_names)
		
		# Now we have computed g-types, we can compute the required data
		levi_cent_sizes = Pol{Rational{Int64}}.(Pol.(levi_xt.cardCent))
		levi_class_sizes = Pol{Rational{Int64}}.(Pol.(levi_xt.cardClass))
		levi_greens_functions = Pol{Rational{Int64}}.(Pol.(levi_gt.scalar[1,:]))
		# .cardCent, .cardClass and .scalar returns multi-variable polys (specifically the "Mvp" data type) so we converted them to single-variable polys
		
		levi_orbit_size = length(orbit(G_dual,levi));
		levi_mobius_value = mob(levi,G_dual,levis);
		levi_orderpol = orderpol(levi)
		
		for i = 1:length(levi_rational_orbit_labels)
				algtypes = vcat(algtypes,
				[(levi,levi_rational_orbit_names[i]) degree(orderpol(levi))-degree(levi_class_sizes[i]) Pol{Int64}((levi_class_sizes[i])*(G_orderpol//levi_orderpol)) levi_greens_functions[i] levi_orbit_size levi_mobius_value]
				);
		end
	end
	return algtypes
end


function alg_row_term(G,i,genus_num,puncture_num)
	# Returns (||Z||(q)/||G||(q)) * q^\xi(G) * q^(g*dim(τ)) * k_τ(q) 
	
	# Grab necessary data
	algtype_data = algtypes(G)
	G_rank = rank(G);
	G_ssrank = semisimplerank(G);
	T = torus(G_rank);
	Z = torus(G_rank-G_ssrank);
	G_root_size = length(roots(G));
	Z_orderpol = orderpol(Z)
	G_orderpol = orderpol(G)
	Z_dim = degree(Z_orderpol)
	G_dim = degree(G_orderpol)
	
	# For readability:
	# algtype_data[i,:][1] = ith type, ie. τ=[L,Oᴸ(Fq)]
	# algtype_data[i,:][2] = dim(τ) of ith type
	# algtype_data[i,:][3] = |Oᴳ(Fq)| of ith type
	# algtype_data[i,:][4] = Qᴸ(Oᴸ(Fq)) of ith type
	# algtype_data[i,:][5] = |[L]| of ith type
	# algtype_data[i,:][6] = µ(L,G) of ith type
	
	coeff = (Z_orderpol*Pol(:q)^(G_root_size+Z_dim+(genus_num)*G_dim))//G_orderpol
	q_dim_tau = Pol(:q)^(algtype_data[i,:][2])
	k_tau_q = algtype_data[i,:][3]*algtype_data[i,:][4]*algtype_data[i,:][5]*algtype_data[i,:][6]
	
	return coeff*q_dim_tau^genus_num*k_tau_q
end



function EY(G,genus_num,puncture_num)
	# Returns the E-polynomial of the additive character variety 
	# associated to (G,g,n) where g = genus_num and n = puncture_num
	algtype_data = algtypes(G);
	epol = Pol([0]);
	for row_number in 1:size(algtype_data)[1]
		epol += alg_row_term(G,row_number,genus_num,puncture_num)
	end
	return Pol{Int64}(epol)
end

#################################################################################
#################################################################################
#################################################################################


# Display human-readable table
function gptable(G)
	gptype_data = gptypes(G)
	clabels = ["|Φ(L)+|","|L(Fq)|","ρ(1)","χᵨ(1)","|W(L)|","|[L]|","ν(L)"];
	rlabels = xrepr.(Ref(rio()),gptype_data[:,1]); 
	# xrepr(rio(), __ ) is a string of __ when printed on the REPR
	repr_gptype_data = xrepr.(Ref(rio()),gptype_data[:,2:size(gptype_data)[2]]);
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
	return showtable(repr_gptype_data;col_labels=clabels,rows_label="Types [L,ρ]",row_labels=rlabels)
end

function algtable(G)
	algtype_data = algtypes(G)
	clabels = ["dim(τ)","|Oᴳ(Fq)|","Qᴸ(Oᴸ(Fq))","|[L]|","µ(L,G)"];
	rlabels = xrepr.(Ref(rio()),algtype_data[:,1]); 
	# xrepr(rio(), __ ) is a string of __ when printed on the REPR
	repr_algtype_data = xrepr.(Ref(rio()),algtype_data[:,2:size(algtype_data)[2]]);
	println("A g-type is a W-orbit τ=[L,Oᴸ(Fq)] where ")
	println("L is a Levi subgroup of G containing T")
	println("Oᴸ(Fq) is the rational unipotent orbit of the L(Fq)-action on Lie(L)(Fq)")
	println("dim(τ) is dim(Lie(L)) - dim(Oᴸ(Fq)) = degree |L(Fq)| - degree |Oᴸ(Fq)|")
	println("Oᴳ(Fq) is the rational unipotent orbit of the G(Fq)-action on Lie(G)(Fq) associated to τ")
	println("Qᴸ(Oᴸ(Fq)) is the Greens function associated to (L,T) evaluated at Oᴸ(Fq)")
	println("[L] is the orbit of L under the W-action")
	println("µ(L,G) is the Mobius function of the poset of Levi subgroups of G containing T evaluated at L and G")
	println("")
	return showtable(repr_algtype_data;col_labels=clabels,rows_label="Types [L,Oᴸ(Fq)]",row_labels=rlabels)
end


#################################################################################
#################################################################################
#################################################################################

# Testing counting functions
function is_palindromic(f)
	return (f(0)!= 0) && (f.c == f.c[end:-1:1])
end

function euler_zero_X(G,genus_max,puncture_max)
	# Checks if χ(X) is zero or non-zero from g=1,2,..,genus_max and n=1,2,...,puncture_max
	for g in 1:genus_max
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

