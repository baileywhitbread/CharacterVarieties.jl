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

function ssrank_diff(L)
	Int((length(roots(G)) - length(roots(L)))/2)
end

function subset(L,M)
	return issubset(inclusion(L),inclusion(M))
end

function mob(A,B,poset)
	if A == B
		return 1
	elseif subset(A,B)
		mobvalue = 0
		for element in poset
			if subset(A,element) && subset(element,B) && element != B
				mobvalue += mob(A,element,poset)
			end
		end
		return (-1)*mobvalue
	else
		error("First argument must be a subset of the second argument")
	end
end

function orbit(L)
	return 1
end

function nu(L)
	nuvalue = 0
	for iplevi in iplevis
		if subset(L,iplevi)
			nuvalue += mob(L,iplevi,plevis)*pi0(iplevi)
		end
	end
	return nuvalue
end



# Choose group G
G = rootdatum(:G2);

# Gather info about G
G_dual = rootdatum(simplecoroots(G),simpleroots(G));
G_rank = size(simpleroots(G),2);
G_ssrank = size(simpleroots(G),1);
W = rootdatum(cartan(G));
T_order = orderpol(torus(G_rank));
Z_order = orderpol(torus(G_rank-G_ssrank));

# Compute pseudo Levis and isolated pseudo Levis in G_dual
plevis = reflection_subgroup.(Ref(G_dual),sscentralizer_reps(G_dual)); 
iplevis = [];
for plevi in plevis
	for iplevi in unique(map(L -> L.group, centralizer.(Ref(G_dual),quasi_isolated_reps(G_dual))))
		if sort(inclusion(plevi)) == sort(inclusion(iplevi))
			append!(iplevis,[plevi])
		end
	end
end

# Compute Levis in G_dual
levis = map(L -> L.W, split_levis(G_dual)); 

# Compute number of G-types
number_of_gptypes = 0;
for plevi in plevis
	for unipotent_degree in CycPoldegrees(UnipotentCharacters(plevi))
		if Int(unipotent_degree(1)) != 0
			global number_of_gptypes += 1;
		end
	end
end

# Compute group G-types
gptypes = Array{Any}(nothing,0,7);
for plevi in plevis
	order_plevi = orderpol(plevi);
	ssrank_diff_plevi = ssrank_diff(plevi);
	orbit_plevi = "???";
	nu_plevi = nu(plevi);
	uc_plevi = UnipotentCharacters(plevi);
	name = charnames(uc_plevi,limit=true);
	ucdegs = CycPoldegrees(uc_plevi);
	for i in 1:length(uc_plevi)
		if Int(ucdegs[i](1))!=0
			global gptypes = vcat(gptypes,
			[(plevi,name[i]) ssrank_diff_plevi ucdegs[i] order_plevi Int(ucdegs[i](1)) orbit_plevi nu_plevi]
			);
		end
	end
end

# Sort according to ssrank(G)-ssrank(L)
sortslices(gptypes,dims=1,by=x->x[2]);

#labels = ["([L],ρ)","r(G)-r(L)","|L(Fq)|","χᵨ(1)","ρ(1)","|[L]|","nu(L)"];
#showtable(gptypes;col_labels=labels,align="|r|r|r|r|r|r|r|r|")

