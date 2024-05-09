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
	# gptypes[i,:][1] = a human-readable representation of the ith type, ie. a pair [L,ρ]
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
	# A problem: |[G]| > 1 iff G=SO5
	# Therefore we kill the duplicate that occurs
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
		# pick unipotent character
		for i in 1:length(plevi_uc)
			# check if unipotent character is principal
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
	# Returns ||tau||(q)^(2g-2+n)*S_tau(q)
	# where tau is the ith type, g = genus_num and n = puncture_num
	
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
	# gptypes[i,:][1] = a human-readable representation of the ith type, ie. a pair [L,ρ]
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
	# algtypes[i,:][1] = ith type, ie. [L,O]
	# algtypes[i,:][2] = dim(ith type) of ith type
	# algtypes[i,:][3] = |O(Fq)| of ith type
	# algtypes[i,:][4] = Q_T^L(O) of ith type
	# algtypes[i,:][5] = |[L]| of ith type
	# algtypes[i,:][6] = mu(L,G) of ith type
	
	# Gather required info about G
	G_dual = rootdatum(simplecoroots(G),simpleroots(G));
	
	# Compute Levis in G_dual
	lorbit_reps = map(L -> L.W, split_levis(G_dual));
	lorbits = orbits(G_dual,lorbit_reps)
	levis = [];
	for lorbit in lorbits
		for levi in lorbit
			append!(levis,[levi])
		end
	end
	
	# Create algtypes
	# I am grabbing the unipotent orbits over the algebraic closure... How do I grab the rational orbits?
	algtypes = Array{Any}(nothing,0,6);
	for levi in lorbit_reps
		levi_orbit_size = length(orbit(G_dual,levi));
		levi_mobius_value = mob(levi,G_dual,levis);
		levi_uc = UnipotentClasses(levi);
		levi_uc_classes = UnipotentClasses(levi).classes;
		levi_xt = XTable(levi_uc;classes=true);
		for class in levi_uc_classes
			algtypes = vcat(algtypes,
			[(levi,class) "?" "?" "?" levi_orbit_size levi_mobius_value]
			);
		end
	end
	return algtypes
end


function alg_row_term(G,i,genus_num,puncture_num)
	# Returns (formula from paper)
	# where tau is the ith type, g = genus_num and n = puncture_num
	
	# Grab necessary data
	algtype_data = algtypes(G)
	# For readability:
	# algtypes[i,:][1] = ith type, ie. [L,O]
	# algtypes[i,:][2] = dim(ith type) of ith type
	# algtypes[i,:][3] = |O(Fq)| of ith type
	# algtypes[i,:][4] = Q_T^L(O) of ith type
	# algtypes[i,:][5] = |[L]| of ith type
	# algtypes[i,:][6] = mu(L,G) of ith type
	
	return 0
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


# Display human-readable table (work in progress, many columns display poorly)
function gptable(G)
	gptype_data = gptypes(G)
	clabels = ["|Φ(L)+|","|L(Fq)|","ρ(1)","χᵨ(1)","|W(L)|","|[L]|","ν(L)"];
	rlabels = gptype_data[:,1];
	alignment = "llllllll";
	return showtable(gptype_data[:,2:size(gptype_data)[2]];col_labels=clabels,row_labels=rlabels,align=alignment)
end

function algtable(G)
	algtype_data = algtypes(G)
	clabels = [];
	rlabels = algtype_data[:,1];
	alignment = "";
	return showtable(algtype_data[:,2:size(algtype_data)[2]];col_labels=clabels,row_labels=rlabels,align=alignment)
end


#################################################################################
#################################################################################
#################################################################################

# Testing counting functions
function is_palindromic(f)
	return (f(0)!= 0) && (f.c == f.c[end:-1:1])
end

function euler_zero(G,genus_max,puncture_max)
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





















