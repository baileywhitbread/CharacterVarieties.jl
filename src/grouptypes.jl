struct GType
	# A G-type is a pair [L,ρ] where 
	# L is an endoscopy group of G containing T
	# ρ is a principal unipotent character of L(Fq)
	
	# To record ρ, we record a string representation and its degree
	# eg. the Steinberg character of GL3 is recorded as 
	# "3" (for the partition (3,0,...) of 3) and Pol(:q)^3

	# I also force the computation of |[L]| and nu(L) 
	# at the creation of a type so I don't have to 
	# compute them many times
	
	endoscopy::FiniteCoxeterGroup
	character::String
	degree::Pol{Rational{BigInt}}
	orbit_size::Union{BigInt,String}
	nu::Union{BigInt,String}

	
end # End of struct GType

# Make GTypes display nicely on the REPL
Base.show(io::IO, tau::GType) = print(io,
"[",tau.endoscopy,",",tau.character,"]"
)




## G-type functions
function group_types(G::FiniteCoxeterGroup)
	# Returns a vector of GTypes, ie. the G-types of G
	types = GType[]
	G_dual = rootdatum(simplecoroots(G),simpleroots(G))
	G_dual_iplevis = iplevis(G_dual)
	G_dual_plevis = plevis(G_dual)
	for plevi in plorbit_reps(G_dual)
		# I am grabbing pseudo-Levis of G rather than endoscopies of G...
		# So far this has not caused a problem because I only need data
		# preserved by Langlands duality, eg. unipotent character degrees
		plevi_uc = UnipotentCharacters(plevi)
		plevi_uc_names = charnames(plevi_uc,limit=true)
		plevi_uc_degs = degrees(plevi_uc)
		plevi_orbit_size = length(myorbit(plevi))
		plevi_nu = nu(plevi,G_dual_iplevis,G_dual_plevis)
		for i in 1:length(plevi_uc)
			# Check if unipotent character is principal
			if BigInt(plevi_uc_degs[i](1))!=0
				append!(types,[GType(plevi,plevi_uc_names[i],plevi_uc_degs[i],plevi_orbit_size,plevi_nu)])
			end
		end

	end
	return types
end

function group_types_no_data(G::FiniteCoxeterGroup)
	# This is only intended for quick inspection of GTypes
	G_dual = rootdatum(simplecoroots(G),simpleroots(G))
	types = GType[]
	for plevi in plorbit_reps(G_dual)
		# I am grabbing pseudo-Levis of G rather than endoscopies of G...
		# So far this has not caused a problem because I only need data
		# preserved by Langlands duality, eg. unipotent character degrees
		plevi_uc = UnipotentCharacters(plevi)
		plevi_uc_names = charnames(plevi_uc,limit=true)
		plevi_uc_degs = degrees(plevi_uc)
		for i in 1:length(plevi_uc)
			# Check if unipotent character is principal
			if BigInt(plevi_uc_degs[i](1))!=0
				append!(types,[GType(plevi,plevi_uc_names[i],plevi_uc_degs[i],"???","???")])
			end
		end
	end
	return types
end

function group_type_data(G::FiniteCoxeterGroup)
	# Packing GTypes and their data into an array for epolys.jl and tables.jl
	d = Array{Any}(nothing,0,8)
	for type in group_types(G)
		type_row = Array{Any}(nothing,1,0)
		type_row = hcat(type_row,[type])									# Type
		type_row = hcat(type_row,[BigInt(length(roots(type.endoscopy))/2)])	# |Phi(L)+|
		type_row = hcat(type_row,[orderpol(type.endoscopy)])				# |L(Fq)|
		type_row = hcat(type_row,[type.degree])								# rho(1)
		type_row = hcat(type_row,[BigInt(type.degree(1))])					# phi(1)
		type_row = hcat(type_row,[length(type.endoscopy)])					# |W(L)|
		type_row = hcat(type_row,[type.orbit_size])							# |[L]|
		type_row = hcat(type_row,[type.nu])									# nu(L)
		d = vcat(d,type_row)
	end
	return sortslices(d,dims=1,by = x -> x[2],rev=true)
end