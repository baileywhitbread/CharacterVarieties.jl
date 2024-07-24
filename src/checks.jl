function isisolated(L::FiniteCoxeterGroup)
	# Returns true if L is isolated in its parent
	# If no parent, returns true
	try
		return length(gens(L)) == length(gens(L.parent))
	catch err
		return true
	end
end

function islevi(L::FiniteCoxeterGroup)
	# Returns true if L is a Levi subgroup of its parent
	# If no parent, returns true
	try
		return issubset(inclusiongens(L),inclusiongens(L.parent))
	catch err
		return true
	end
end

function ispalindromic(f::Union{Pol{BigInt},Pol{Int64}})
	# Returns true <=> f palindromic
	# The zero poly is palindromic and f=0 <=> f(e)=0
	# A non-zero poly is palindromic <=> it has non-zero constant term and coeffs are symmetric
	return (f(exp(1))==0) || ((f(0)!= 0) && (f.c == f.c[end:-1:1]))
end

function isnonnegative(f::Union{Pol{BigInt},Pol{Int64}})
	# Returns true <=> f has non-negative coefficients
	return length(filter(x -> x<0, f.c)) == 0
end