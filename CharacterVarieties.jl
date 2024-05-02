# Load Jean Michel's package
#import Pkg
#Pkg.add("Chevie")
using Chevie

# Choose group G
G = rootdatum(:gl,2);
G_dual = rootdatum(:gl,2);
# Used to be
# G_dual = rootdatum(simplecoroots(G),simpleroots(G));

# Gather info about G
T = torus(rank(G));
Z = torus(rank(G)-semisimplerank(G));
G_positive_root_size = Int(length(roots(G))/2);
G_weyl_group_size = length(G);


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

function nu(L)
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

# Display readable table
function gptable(types)
	clabels = ["|Φ(L)+|","|L(Fq)|","ρ(1)","χᵨ(1)","|W(L)|","|[L]|","ν(L)"];
	rlabels = types[:,1];
	alignment = "llllllll";
	return showtable(types[:,2:size(types)[2]];col_labels=clabels,row_labels=rlabels,align=alignment)
end

function algtable(types)
	clabels = [];
	rlabels = types[:,1];
	alignment = "";
	return showtable(types[:,2:size(types)[2]];col_labels=clabels,row_labels=rlabels,align=alignment)
end

#################################################################################
#################################################################################
#################################################################################

# EX() turns gptypes matrix into an E-polynomial
# gptypes[i,:][2] = |Φ(L)+| of ith type
# gptypes[i,:][3] = ρ(1) (unipotent character degree) of ith type
# gptypes[i,:][4] = |L(Fq)| of ith type
# gptypes[i,:][5] = χᵨ(1) (Weyl gp character degree) of ith type
# gptypes[i,:][6] = |W(L)| of ith type
# gptypes[i,:][7] = |[L]| of ith type
# gptypes[i,:][8] = ν(L) of ith type

function gp_row_term(data,type_num,genus_num,puncture_num)
	coeff = orderpol(Z)//(orderpol(T)^puncture_num);
	tau_q = ((Pol(:q)^(G_positive_root_size-data[type_num,:][2]))*data[type_num,:][4])//data[type_num,:][3];
	S_tau_q = orderpol(Z)*(data[type_num,:][5]^puncture_num)*data[type_num,:][7]*((G_weyl_group_size//data[type_num,:][6])^(puncture_num-1))*data[type_num,:][8]
	return coeff*(tau_q)^(2*genus_num-2+puncture_num)*S_tau_q
end

function EX(data,genus_num,puncture_num)
	epol = Pol([0]);
	for row_number in 1:size(data)[1]
		epol += gp_row_term(data,row_number,genus_num,puncture_num)
	end
	return Pol{Int128}(epol) 
end

### Integer overflows occur if G,g,n are big
### overflow occurs when G=SO7 and g=n=2 but not when G=SO7, g=1 and n=2
### overflow occurs when G=G2 and g=n=2 but not when G=G2, g=1 and n=2
### Overflow does not occur when G=GL4 and g=n=5 but does when g>5 and n>5

### Solution: Replace Pol{Int64} with Pol{Int128}

#################################################################################
#################################################################################
#################################################################################

# Testing counting functions
function is_palindromic(f)
	return (f(0)!= 0) && (f.c == f.c[end:-1:1])
end

#################################################################################
#################################################################################
#################################################################################


# Compute pseudo Levi orbit representatives, pseudo Levi orbits, pseudo Levis and isolated pseudo Levis
plorbit_reps = reflection_subgroup.(Ref(G_dual),sscentralizer_reps(G_dual)); 
plorbits = orbits(G_dual,plorbit_reps); # duplicates. eg. G=SO5
plevis = [];
iplevis = [];
for plorbit in plorbits
	for plevi in plorbit
		append!(plevis,[plevi])
		if length(gens(plevi)) == length(gens(G))
			append!(iplevis,[plevi])
		end
	end
end

# Compute group G-types
# gptypes[i,:][1] = ith type, ie. [L,ρ]
# gptypes[i,:][2] = |Φ(L)+| of ith type
# gptypes[i,:][3] = ρ(1) (unipotent character degree) of ith type
# gptypes[i,:][4] = |L(Fq)| of ith type
# gptypes[i,:][5] = χᵨ(1) (Weyl gp character degree) of ith type
# gptypes[i,:][6] = |W(L)| of ith type
# gptypes[i,:][7] = |[L]| of ith type
# gptypes[i,:][8] = ν(L) of ith type

# I am building this row-by-row... Should change to column-by-column
gptypes = Array{Any}(nothing,0,8);
for plevi in plorbit_reps
	plevi_order = orderpol(plevi);
	plevi_positive_root_size = Int(length(roots(plevi))/2);
	plevi_orbit_size = length(orbit(G_dual,plevi));
	plevi_weyl_size = length(plevi);
	plevi_nu = nu(plevi);
	plevi_uc = UnipotentCharacters(plevi);
	plevi_uc_names = charnames(plevi_uc,limit=true);
	plevi_uc_degs = degrees(plevi_uc);
	# pick unipotent character
	for i in 1:length(plevi_uc)
		# check if unipotent character is principal
		if Int(plevi_uc_degs[i](1))!=0
			global gptypes = vcat(gptypes,
			[(plevi,plevi_uc_names[i]) plevi_positive_root_size plevi_uc_degs[i] plevi_order Int(plevi_uc_degs[i](1)) plevi_weyl_size plevi_orbit_size plevi_nu]
			);
		end
	end
end
sortslices(gptypes,dims=1,by=x->x[2]); # Sort according to number of positive roots of L


#################################################################################
#################################################################################
#################################################################################

# Compute Levis in G_dual
lorbit_reps = map(L -> L.W, split_levis(G_dual));
lorbits = orbits(G_dual,lorbit_reps)
levis = [];
for lorbit in lorbits
	for levi in lorbit
		append!(levis,[levi])
	end
end

# Compute algebra G-types
# algtypes[i,:][1] = ith type, ie. [L,O]
# algtypes[i,:][2] = dim(ith type) of ith type
# algtypes[i,:][3] = |O(Fq)| of ith type
# algtypes[i,:][4] = Q_T^L(O) of ith type
# algtypes[i,:][5] = |[L]| of ith type
# algtypes[i,:][6] = mu(L,G) of ith type


# I am grabbing the unipotent orbits over the algebraic closure... How do I grab the rational orbits?
algtypes = Array{Any}(nothing,0,3);
for levi in lorbit_reps
	levi_orbit_size = length(orbit(G_dual,levi));
	levi_mobius_value = mob(levi,G_dual,levis);
	levi_uc = UnipotentClasses(levi).classes;
	for class in levi_uc
		global algtypes = vcat(algtypes,
		[(levi,class) levi_orbit_size levi_mobius_value]
		);
	end
end

