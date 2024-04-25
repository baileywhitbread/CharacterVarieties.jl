# Load Jean Michel's package
#import Pkg
#Pkg.add("Chevie")
using Chevie

# Helper functions
function orderpol(G)
	return CycPol(PermRoot.generic_order(G,Pol(:q)))
end

function pi0(L)
	return length(algebraic_center(L).AZ)
end

function root_diff(L)
	Int((length(roots(G)) - length(roots(L)))/2)
end

function subset(L,M)
	return issubset(inclusion(L),inclusion(M))
end

function mob(A,B,poset)
	if A == B
		return 1
	elseif subset(A,B)
		mob_value = 0
		for element in poset
			if subset(A,element) && subset(element,B) && element != B
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

function orbit_size(L)
	return length(orbit(G_dual,L))
end

# Choose group G
G = rootdatum(:G2);
G_dual = rootdatum(:G2);
# Used to be
# G_dual = rootdatum(simplecoroots(G),simpleroots(G));

# Gather info about G
T = torus(rank(G));
Z = torus(rank(G)-semisimplerank(G));

# Compute Levis in G_dual
levis = map(L -> L.W, split_levis(G_dual)); 

# Compute pseudo Levi orbit representatives, pseudo Levi orbits and all pseudo Levis
plorbit_reps = reflection_subgroup.(Ref(G_dual),sscentralizer_reps(G_dual)); 
plorbits = orbits(G_dual,plorbit_reps);
plevis = [];
for plorbit in plorbits
	for plevi in plorbit
		append!(plevis,[plevi])
	end
end

# Compute isolated pseudo Levi orbit representatives, isolated pseudo Levi orbits and all isolated pseudo Levis
iplorbit_reps = unique(map(L -> L.group, centralizer.(Ref(G_dual),quasi_isolated_reps(G_dual))));
iplorbits = orbits(G_dual, iplorbit_reps);
iplevis = [];
for plevi in plevis
	for iplorbit in iplorbits
		for iplevi in iplorbit
			if sort(inclusion(plevi)) == sort(inclusion(iplevi))
				append!(iplevis,[plevi])
			end
		end
	end
end

# Compute group G-types
gptypes = Array{Any}(nothing,0,7);
for plevi in plorbit_reps
	plevi_order = orderpol(plevi);
	plevi_root_diff = root_diff(plevi);
	plevi_orbit_size = orbit_size(plevi);
	plevi_nu = nu(plevi);
	plevi_uc = UnipotentCharacters(plevi);
	plevi_uc_names = charnames(plevi_uc,limit=true);
	plevi_uc_degs = CycPoldegrees(plevi_uc);
	# pick unipotent character
	for i in 1:length(plevi_uc)
		# check if unipotent character is principal
		if Int(plevi_uc_degs[i](1))!=0
			global gptypes = vcat(gptypes,
			[(plevi,plevi_uc_names[i]) plevi_root_diff plevi_uc_degs[i] plevi_order Int(plevi_uc_degs[i](1)) plevi_orbit_size plevi_nu]
			);
		end
	end
end

# Sort according to ssrank(G)-ssrank(L)
sortslices(gptypes,dims=1,by=x->x[2]);

#labels = ["([L],ρ)","r(G)-r(L)","|L(Fq)|","χᵨ(1)","ρ(1)","|[L]|","nu(L)"];
#showtable(gptypes;col_labels=labels,align="|r|r|r|r|r|r|r|r|")

