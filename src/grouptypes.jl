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




## G-type functions
function group_types(G::FiniteCoxeterGroup)
	# Returns a vector of GTypes, ie. the G-types of G
	G_dual = rootdatum(simplecoroots(G),simpleroots(G))
	types = []
	for plevi in plorbit_reps(G_dual)
		# I am grabbing pseudo-Levis of G rather than endoscopies of G...
		# So far this has not caused a problem because I only need data
		# preserved by Langlands duality, eg. unipotent character degrees
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
	G_dual = rootdatum(simplecoroots(G),simpleroots(G))
	for type in group_types(G)
		type_row = Array{Any}(nothing,1,0)
		type_row = hcat(type_row,[type])
		type_row = hcat(type_row,[Int64(length(roots(type.endoscopy))/2)])
		type_row = hcat(type_row,[orderpol(type.endoscopy)])
		type_row = hcat(type_row,[type.degree])
		type_row = hcat(type_row,[Int64(type.degree(1))])
		type_row = hcat(type_row,[length(type.endoscopy)])
		type_row = hcat(type_row,[length(myorbit(type.endoscopy))])
		type_row = hcat(type_row,[nu(type.endoscopy,iplevis(G_dual),plevis(G_dual))])		
		d = vcat(d,type_row)
	end
	return sortslices(d,dims=1,by = x -> x[2],rev=true)
end