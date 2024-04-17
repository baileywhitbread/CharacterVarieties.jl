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
Gdual = rootdatum(simplecoroots(G),simpleroots(G));
rank = size(simpleroots(G),2);
ssrank = size(simpleroots(G),1);
W = rootdatum(cartan(G));
orderT = orderpol(torus(rank));
orderZ = orderpol(torus(rank-ssrank));

# Compute pseudo Levis and isolated pseudo Levis in Gdual
plevis = reflection_subgroup.(Ref(Gdual),sscentralizer_reps(Gdual)); 
iplevis = [];
for plevi in plevis
	for iplevi in unique(map(L -> L.group, centralizer.(Ref(Gdual),quasi_isolated_reps(Gdual))))
		if sort(inclusion(plevi)) == sort(inclusion(iplevi))
			append!(iplevis,[plevi])
		end
	end
end



# Compute Levis in Gdual
levis = map(L -> L.W, split_levis(Gdual)); 

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
	uc_plevi = UnipotentCharacters(plevi);
	name = charnames(uc_plevi,limit=true);
	ucdegs = CycPoldegrees(uc_plevi);
	for i in 1:length(uc_plevi)
		if Int(ucdegs[i](1))!=0
			global gptypes = vcat(gptypes,[(plevi,name[i]) Int((length(roots(G)) - length(roots(plevi)))/2) ucdegs[i] orderpol(plevi) Int(ucdegs[i](1)) 0 nu(plevi)]);
		end
	end
end

# Sort according to semisimple rank
sortslices(gptypes,dims=1,by=x->x[2]);
