# Load Jean Michel's package
using Chevie

# Helper functions
function orderpol(G)
	return CycPol(PermRoot.generic_order(G,Pol(:q)))
end



# Choose group G. Some options are: 
# gl, pgl, so, psp, G2, F4, E6, E7, E8 
G=rootdatum(:so,5)

# Gather info about G
Gdual = rootdatum(simplecoroots(G),simpleroots(G))
rank = size(simpleroots(G),2)
ssrank = size(simpleroots(G),1)
W = rootdatum(cartan(G))
QIR = quasi_isolated_reps(G)
UC = UnipotentCharacters(G)
Z = algebraic_center(G)[1]
T = torus(G,position_class(G,G(1,2))).W 

# Compute order polynomials
orderG = orderpol(G)
orderT = orderpol(T)
orderZ = orderpol(torus(rank-ssrank))

# Compute pseudo levis in Gdual
plevis = reflection_subgroup.(Ref(Gdual),sscentralizer_reps(Gdual)); 
# For a one-argument function f, we have f.([a b c]) = [f(a) f(b) f(c)]
# For a two-argument function g, we have g.(Ref(x),[a b c]) = [g(x,a) g(x,b) g(x,c)]

# Compute levis in Gdual
levis = map(L -> L.W, split_levis(Gdual)); 
# split_levis(Gdual) is a list of spets objects, L -> L.W makes them groups

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
gptypes = Array{Any}(nothing,0,6);
for plevi in plevis
	uc_plevi = UnipotentCharacters(plevi);
	name = charnames(uc_plevi,limit=true);
	ucdegs = CycPoldegrees(uc_plevi);
	for i in 1:length(uc_plevi)
		if Int(ucdegs[i](1))!=0
			global gptypes = vcat(gptypes,[plevi Int((length(roots(G)) - length(roots(plevi)))/2) orderpol(plevi) name[i] ucdegs[i] Int(ucdegs[i](1))]);
		end
	end
end

# Compute alg G-types
algtypes = Array{Any}(nothing,0,6);
for levi in levis
	uc_levi = UnipotentCharacters(levi);
	name = charnames(uc_levi,limit=true);
	ucdegs = CycPoldegrees(uc_levi);
	levi_size = CycPol(PermRoot.generic_order(levi,Pol(:q)))
	for i in 1:length(uc_levi)
		if Int(ucdegs[i](1))!=0
			global algtypes = vcat(algtypes,[levi Int((length(roots(G)) - length(roots(levi)))/2) levi_size name[i] ucdegs[i] Int(ucdegs[i](1))]);
		end
	end
end







