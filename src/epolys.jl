## E-polynomial of X
function fast_Mtau(G::FiniteCoxeterGroup,i::Int64,type_data)
	# Returns Mτ(q) = q^(|Φ(G)+|-|Φ(L)+|) |L(Fq)|/ρ(1) 
	# where τ = [L,ρ] is the ith GType
	L_pos_root_size = type_data[i,:][2]
	L_size = type_data[i,:][3]
	rho_deg = type_data[i,:][4]
	return Pol{Rational{Int64}}(Pol(:q)^(Int64(length(roots(G))/2)-L_pos_root_size)*L_size//rho_deg)
end

function fast_Stau(G::FiniteCoxeterGroup,n::Int64,i::Int64,type_data)
	# Returns Sτ(q) = |Z(Fq)| * χᵨ(1)^n * |[L]| * (|W|/|W(L)|)^(n-1) * ν(L)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	chi_rho_deg = type_data[i,:][5]
	G_weyl_size = length(G)
	L_weyl_size = type_data[i,:][6]
	orbit_size = type_data[i,:][7]
	nu_L = type_data[i,:][8]
	return Pol{Rational{Int64}}(Z_size * chi_rho_deg^n * orbit_size * (Int64(G_weyl_size//L_weyl_size))^(n-1) * nu_L)
end

function EX(G::FiniteCoxeterGroup,g::Int64,n::Int64)
	# Returns the E-polynomial E(X;q) associated to the group G and a genus g surface with n punctures
	d = group_type_data(G)
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	T_size = orderpol(torus(rank(G)))
	type_sum = 0
	for i in 1:size(d)[1]
		type_sum += fast_Mtau(G,i,d)^(2g-2+n)*fast_Stau(G,n,i,d)
	end
	return Pol{Int64}((Z_size//T_size^n)*type_sum)
end

## E-polynomial of Y
function fast_qdtau(G::FiniteCoxeterGroup,i::Int64,type_data)
	# Returns q^(d(τ)) where τ = [L,ρ] is the ith GType
	return Pol(:q)^(type_data[i,:][2])
end

function fast_Htau(G::FiniteCoxeterGroup,n::Int64,i::Int64,type_data)
	# Returns Hτ(q) = q^(n|Φ(G)+| + dim(Z)) * (|G(Fq)|/|L(Fq)|) * |N| * Q_L^T(N)^n * |[L]| * (|W|/|W(L)|)^(n-1) * µ(L,G)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
	G_weyl_size = length(G)
	G_pos_root_size = Int64(length(G)/2)
	Z_dim = rank(G)-semisimplerank(G)
	L_size = orderpol(type_data[i,:][1].levi)
	N_size = type_data[i,:][1].size
	L_green = type_data[i,:][1].green
	orbit_size = type_data[i,:][5]
	L_weyl_size = length(type_data[i,:][1].levi)
	mu_L = type_data[i,:][6]
	return Pol(:q)^(n*G_pos_root_size + Z_dim) * (orderpol(G)//L_size) * N_size * L_green^n * orbit_size * (G_weyl_size/L_weyl_size)^(n-1) * mu_L
end

function EY(G::FiniteCoxeterGroup,g::Int64,n::Int64)
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	d = algebra_type_data(G)
	Z_size = orderpol(torus(rank(G)-semisimplerank(G)))
	G_size = orderpol(G)
	g_size = Pol(:q)^(degree(G_size))
	type_sum = Pol{Int64}(0)
	for i in 1:size(d)[1]
		type_sum += fast_qdtau(G,i,d)^(g)*fast_Htau(G,n,i,d)
	end
	return Pol{Int64}((Z_size//G_size) * (g_size^(g)//g_size) * type_sum ) 
end








### BigInt version for testing.jl
function bigint_Mtau(G::FiniteCoxeterGroup,i::Int64,type_data)
	# Returns Mτ(q) = q^(|Φ(G)+|-|Φ(L)+|) |L(Fq)|/ρ(1) 
	# where τ = [L,ρ] is the ith GType
	L_pos_root_size = BigInt(type_data[i,:][2])
	L_size = Pol{Rational{BigInt}}(type_data[i,:][3])
	rho_deg = Pol{Rational{BigInt}}(type_data[i,:][4])
	return Pol{Rational{BigInt}}(Pol(:q)^(BigInt(length(roots(G))/2)-L_pos_root_size)*Pol{Rational{BigInt}}(L_size//rho_deg))
end

function bigint_Stau(G::FiniteCoxeterGroup,n::Int64,i::Int64,type_data)
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

function bigint_EX(G::FiniteCoxeterGroup,g::Int64,n::Int64,type_data)
	# Returns the E-polynomial E(X;q) associated to the group G and a genus g surface with n punctures
	Z_size = Pol{BigInt}(orderpol(torus(rank(G)-semisimplerank(G))))
	T_size = Pol{BigInt}(orderpol(torus(rank(G))))
	type_sum = BigInt(0)
	for i in 1:size(type_data)[1]
		type_sum += bigint_Mtau(G,i,type_data)^(2g-2+n)*bigint_Stau(G,n,i,type_data)
	end
	return Pol{BigInt}(Pol{Rational{BigInt}}((Z_size//T_size^n))*type_sum)
end

function bigint_qdtau(G::FiniteCoxeterGroup,i::Int64,type_data)
	# Returns q^(d(τ)) where τ = [L,ρ] is the ith GType
	return Pol{BigInt}(Pol(:q)^(type_data[i,:][2]))
end

function bigint_Htau(G::FiniteCoxeterGroup,n::Int64,i::Int64,type_data)
	# Returns Hτ(q) = q^(n|Φ(G)+| + dim(Z)) * (|G(Fq)|/|L(Fq)|) * |N| * Q_L^T(N)^n * |[L]| * (|W|/|W(L)|)^(n-1) * µ(L,G)
	# where τ = [L,ρ] is the ith GType and n is the number of punctures
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



function bigint_EY(G::FiniteCoxeterGroup,g::Int64,n::Int64,type_data)
	# Returns the E-polynomial E(Y;q) associated to the group G and a genus g surface with n punctures
	Z_size = Pol{Rational{BigInt}}(orderpol(torus(rank(G)-semisimplerank(G))))
	G_size = Pol{Rational{BigInt}}(orderpol(G))
	g_size = Pol{Rational{BigInt}}(Pol(:q)^(BigInt(degree(G_size))))

	type_sum = Pol{BigInt}(0)
	for i in 1:size(type_data)[1]
		type_sum += bigint_qdtau(G,i,type_data)^(g)*bigint_Htau(G,n,i,type_data)
	end
	return Pol{BigInt}((Z_size//G_size) * (g_size^(g)//g_size) * type_sum)
end