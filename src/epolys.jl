### "Fast" = no group_type_data(G) or algebra_type_data(G) calls
function fast_Mtau(G::FiniteCoxeterGroup,i::Union{BigInt,Integer},type_data::Any)
	# Returns Mτ(q) = q^(|Φ(G)+|-|Φ(L)+|) |L(Fq)|/ρ(1) 
	# where τ = [L,ρ] is the ith GType
	L_pos_root_size = BigInt(type_data[i,:][2])
	L_size = Pol{Rational{BigInt}}(type_data[i,:][3])
	rho_deg = Pol{Rational{BigInt}}(type_data[i,:][4])
	return Pol{Rational{BigInt}}(Pol(:q)^(BigInt(length(roots(G))/2)-L_pos_root_size)*Pol{Rational{BigInt}}(L_size//rho_deg))
end

function fast_Stau(G::FiniteCoxeterGroup,n::Union{BigInt,Integer},i::Union{BigInt,Integer},type_data::Any)
	# Returns Sτ(q) = |Z(Fq)| * χᵨ(1)^n * |[L]| * (|W|/|W(L)|)^(n-1) * ν(L)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	Z_size = Pol{Rational{BigInt}}(orderpol(torus(rank(G)-semisimplerank(G))))
	chi_rho_deg = BigInt(type_data[i,:][5])
	G_weyl_size = BigInt(length(G))
	L_weyl_size = BigInt(type_data[i,:][6])
	orbit_size = BigInt(type_data[i,:][7])
	nu_L = BigInt(type_data[i,:][8])
	return Pol{Rational{BigInt}}(Z_size * chi_rho_deg^n * orbit_size * (BigInt(G_weyl_size//L_weyl_size))^(n-1) * BigInt(nu_L))
end

function EX(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer})
	# Returns the E-polynomial E(X;q) associated to the group G and a genus g surface with n punctures
	type_data = group_type_data(G)
	Z_size = Pol{BigInt}(orderpol(torus(rank(G)-semisimplerank(G))))
	T_size = Pol{BigInt}(orderpol(torus(rank(G))))
	type_sum = BigInt(0)
	for i in 1:size(type_data)[1]
		type_sum += fast_Mtau(G,i,type_data)^(2g-2+n)*fast_Stau(G,n,i,type_data)
	end
	return Pol{BigInt}(((Z_size//T_size^n))*type_sum)
end

function fast_EX(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer},type_data::Any)
	# Returns the E-polynomial E(X;q) associated to the group G and a genus g surface with n punctures
	Z_size = Pol{BigInt}(orderpol(torus(rank(G)-semisimplerank(G))))
	T_size = Pol{BigInt}(orderpol(torus(rank(G))))
	type_sum = BigInt(0)
	for i in 1:size(type_data)[1]
		type_sum += fast_Mtau(G,i,type_data)^(2g-2+n)*fast_Stau(G,n,i,type_data)
	end
	return Pol{BigInt}((Z_size//T_size^n)*type_sum)
end


function new_EX(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer})
	type_data = group_type_data(G)
	sum = Pol{BigInt}(0)
	for i in 1:size(type_data)[1]
		row = type_data[i,:] # type_data[i,:] is the ith row
		term = Pol{BigInt}(1)
		term *= Pol(:q)^(BigInt(length(roots(G))/2-length(roots(row[1].endoscopy))/2)) # q^(|Phi(G)+|-|Phi(L)+|)
		term *= Pol{Rational{BigInt}}(orderpol(row[1].endoscopy)//row[1].degree) # |L(Fq)|/rho(q)
		term *= (Pol(:q)-1)^BigInt(rank(G)-semisimplerank(G)) # |Z(Fq)|
		term *= BigInt(row[5])^BigInt(n) # dim(rho)^n
		term *= BigInt(length(G)//row[6])^BigInt(n-1) # (|W|/|W(L)|)^(n-1)
		term *= row[7] # |[L]|
		term *= row[8] # nu(L)
		sum += term
	end
	sum *= (Pol(:q)-1)^BigInt(rank(G)-semisimplerank(G))//((sizeofT)^BigInt(n)) # |Z(Fq)|/|T(Fq)|^n
	return Pol{BigInt}(sum)
end

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

function EY(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer})
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	type_data = algebra_type_data(G)
	sum = Pol{BigInt}(0)
	for i in 1:size(type_data)[1]
		row = type_data[i,:] # type_data[i,:] is the ith row
		term = Pol{BigInt}(1)
		term *= Pol(:q)^(BigInt(g*row[2]))
		term *= Pol(:q)^(BigInt(n*(length(roots(G))/2)+(rank(G)-semisimplerank(G))))
		term *= Frac{Pol{BigInt}}(orderpol(G)//orderpol(row[1].levi))
		term *= Pol{Rational{BigInt}}(row[1].size)
		term *= (Pol{BigInt}(row[1].green))^BigInt(n)
		term *= (BigInt(BigInt(length(G))//BigInt(length(row[1].levi))))^BigInt(n-1)
		term *= BigInt(length(myorbit(row[1].levi)))
		term *= BigInt(mobius(row[1].levi,row[1].levi.parent,levis(G)))
		sum += term
	end
	sum *= (((Pol(:q)-1)^BigInt(rank(G)-semisimplerank(G)))//(orderpol(G)))
	sum *= ((Pol(:q)^(BigInt(g)*BigInt(degree(orderpol(G)))))//(Pol(:q)^(BigInt(1)*BigInt(degree(orderpol(G))))))
	return Pol{BigInt}(sum)
end

function fast_EY(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},n::Union{BigInt,Integer},type_data::Any)
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	sum = Pol{BigInt}(0)
	for i in 1:size(type_data)[1]
		row = type_data[i,:] # type_data[i,:] is the ith row
		term = Pol{BigInt}(1)
		term *= Pol(:q)^(BigInt(g*row[2]))
		term *= Pol(:q)^(BigInt(n*(length(roots(G))/2)+(rank(G)-semisimplerank(G))))
		term *= Frac{Pol{BigInt}}(orderpol(G)//orderpol(row[1].levi))
		term *= Pol{Rational{BigInt}}(row[1].size)
		term *= (Pol{BigInt}(row[1].green))^BigInt(n)
		term *= (BigInt(BigInt(length(G))//BigInt(length(row[1].levi))))^BigInt(n-1)
		term *= BigInt(length(myorbit(row[1].levi)))
		term *= BigInt(mobius(row[1].levi,row[1].levi.parent,levis(G)))
		sum += term
	end
	sum *= (((Pol(:q)-1)^(rank(G)-semisimplerank(G)))//(orderpol(G)))
	sum *= ((Pol(:q)^(BigInt(g)*BigInt(degree(orderpol(G)))))//(Pol(:q)^(BigInt(1)*BigInt(degree(orderpol(G))))))
	return Pol{BigInt}(sum)
end































# I don't think these are called anywhere...
function Mtau(G::FiniteCoxeterGroup,i::Int64)
	# Returns Mτ(q) = q^(|Φ(G)+|-|Φ(L)+|) |L(Fq)|/ρ(1) 
	# where τ = [L,ρ] is the ith GType
	type_data = group_type_data(G)
	L_pos_root_size = BigInt(type_data[i,:][2])
	L_size = Pol{Rational{BigInt}}(type_data[i,:][3])
	rho_deg = Pol{Rational{BigInt}}(type_data[i,:][4])
	return Pol{Rational{BigInt}}(Pol(:q)^(BigInt(length(roots(G))/2)-L_pos_root_size)*Pol{Rational{BigInt}}(L_size//rho_deg))
end

function Stau(G::FiniteCoxeterGroup,n::Int64,i::Int64)
	# Returns Sτ(q) = |Z(Fq)| * χᵨ(1)^n * |[L]| * (|W|/|W(L)|)^(n-1) * ν(L)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	type_data = group_type_data(G)
	Z_size = Pol{Rational{BigInt}}(orderpol(torus(rank(G)-semisimplerank(G))))
	chi_rho_deg = BigInt(type_data[i,:][5])
	G_weyl_size = BigInt(length(G))
	L_weyl_size = BigInt(type_data[i,:][6])
	orbit_size = BigInt(type_data[i,:][7])
	nu_L = BigInt(type_data[i,:][8])
	return Pol{Rational{BigInt}}(Z_size * chi_rho_deg^n * orbit_size * (BigInt(G_weyl_size//L_weyl_size))^(n-1) * BigInt(nu_L))
end

function qdtau(G::FiniteCoxeterGroup,i::Int64)
	# Returns q^(d(τ)) where τ = [L,ρ] is the ith GType
	type_data = group_type_data(G)
	return Pol{BigInt}(Pol(:q)^(type_data[i,:][2]))
end

function Htau(G::FiniteCoxeterGroup,n::Int64,i::Int64)
	# Returns Hτ(q) = q^(n|Φ(G)+| + dim(Z)) * (|G(Fq)|/|L(Fq)|) * |N| * Q_L^T(N)^n * |[L]| * (|W|/|W(L)|)^(n-1) * µ(L,G)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	type_data = group_type_data(G)
	G_weyl_size = BigInt(length(G))
	G_pos_root_size = BigInt(length(G)/2)
	Z_dim = BigInt(rank(G)-semisimplerank(G))
	L_size = Pol{Rational{BigInt}}(orderpol(type_data[i,:][1].levi))
	N_size = Pol{Rational{BigInt}}(type_data[i,:][1].size)
	L_green = Pol{Rational{BigInt}}(type_data[i,:][1].green)
	orbit_size = BigInt(type_data[i,:][5])
	L_weyl_size = BigInt(length(type_data[i,:][1].levi))
	mu_L = BigInt(type_data[i,:][6])
	return Pol(:q)^(BigInt(n*G_pos_root_size + Z_dim)) * (orderpol(G)//L_size) * N_size * L_green^n * BigInt(orbit_size) * BigInt((BigInt(G_weyl_size)/BigInt(L_weyl_size))^(n-1)) * BigInt(mu_L)
end