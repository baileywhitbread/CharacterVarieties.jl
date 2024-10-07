### "Fast" = no group_type_data(G) or algebra_type_data(G) calls
function fast_Mtau(G::FiniteCoxeterGroup,i::Union{BigInt,Integer},type_data::Any)
	row = type_data[i,:] # type_data[i,:] is the ith row
	term = Pol{BigInt}(1)
	term *= Pol(:q)^(BigInt((length(roots(G))/2)-(length(roots(row[1].endoscopy))/2))) # q^(|Phi(G)+|-|Phi(L)+|)
	term *= Pol{BigInt}(orderpol(row[1].endoscopy)//row[1].degree) # |L(Fq)|/rho(q)
	return Pol{BigInt}(term)
end

function fast_Stau(G::FiniteCoxeterGroup,n::Union{BigInt,Integer},i::Union{BigInt,Integer},type_data::Any)
	row = type_data[i,:] # type_data[i,:] is the ith row
	term = Pol{BigInt}(1)
	term *= (Pol(:q)-1)^BigInt(rank(G)-semisimplerank(G)) # |Z(Fq)|
	term *= BigInt(row[5])^BigInt(n) # dim(rho)^n
	term *= BigInt(length(G)//row[6])^BigInt(n-1) # (|W|/|W(L)|)^(n-1)
	term *= row[7] # |[L]|
	term *= row[8] # nu(L)
	return Pol{BigInt}(term)
end

function fast_EX(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer},type_data::Any)
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	sum = Pol{BigInt}(0)
	for i in 1:size(type_data)[1]
		term = Pol{BigInt}(1)
		term *= fast_Mtau(G,i,type_data)^BigInt(2*g-2+n)
		term *= fast_Stau(G,n,i,type_data)
		sum += term
	end
	sum *= ((Pol(:q)-1)^BigInt(rank(G)-semisimplerank(G))//(((Pol(:q)-1)^(BigInt(rank(G))))^BigInt(n))) # |Z(Fq)|/|T(Fq)|^n
	return Pol{BigInt}(sum)
end

function EX(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer})
	# Returns the E-polynomial E(X;q) associated to the group G and a genus g surface with n punctures
	type_data = group_type_data(G)
	return fast_EX(G,g,n,type_data)
end

#############################################################################################################################

function fast_qdtau(G::FiniteCoxeterGroup,i::Union{BigInt,Integer},type_data::Any)
	# Returns q^(d(τ)) where τ = [L,ρ] is the ith GType
	row = type_data[i,:] # type_data[i,:] is the ith row
	return Pol{BigInt}(Pol(:q)^(BigInt(row[2])))
end

function fast_Htau(G::FiniteCoxeterGroup,n::Union{BigInt,Integer},i::Union{BigInt,Integer},type_data::Any)
	row = type_data[i,:] # type_data[i,:] is the ith row
	term = Pol{BigInt}(1)
	term *= Pol(:q)^(BigInt(n*(length(roots(G))/2)+(rank(G)-semisimplerank(G))))
	term *= Frac{Pol{BigInt}}(orderpol(G)//orderpol(row[1].levi))
	term *= Pol{Rational{BigInt}}(row[1].size)
	term *= (Pol{BigInt}(row[1].green))^BigInt(n)
	term *= (BigInt(BigInt(length(G))//BigInt(length(row[1].levi))))^BigInt(n-1)
	term *= BigInt(length(myorbit(row[1].levi)))
	term *= BigInt(mobius(row[1].levi,row[1].levi.parent,levis(G)))
	return Frac{Pol{Rational{BigInt}}}(term)
end

function fast_EY(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer},type_data::Any)
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	sum = Pol{BigInt}(0)
	for i in 1:size(type_data)[1]
		term = Pol{BigInt}(1)
		term *= fast_qdtau(G,i,type_data)^BigInt(g)
		term *= fast_Htau(G,n,i,type_data)
		sum += term
	end
	sum *= (((Pol(:q)-1)^(rank(G)-semisimplerank(G)))//(orderpol(G)))
	sum *= ((Pol(:q)^(BigInt(g)*BigInt(degree(orderpol(G)))))//(Pol(:q)^(BigInt(1)*BigInt(degree(orderpol(G))))))
	return Pol{BigInt}(sum)
end

function EY(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer})
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	type_data = algebra_type_data(G)
	return fast_EY(G,g,n,type_data)
end