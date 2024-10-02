"""
    isisolated(L::FiniteCoxeterGroup) -> Bool

Check if a Coxeter group is isolated in its parent group.
A group is considered isolated if it has the same number of generators as its parent.
Returns true if L has no parent.

# Examples
```julia
G = CoxeterGroup(:A, 4)
H = reflection_subgroup(G, [1, 2, 3])
isisolated(H)  # returns false
```
"""
function isisolated(L::FiniteCoxeterGroup)
    parent_gens_length = try
        length(gens(L.parent))
    catch
        return true  # No parent
    end
    
    return length(gens(L)) == parent_gens_length
end

"""
    islevi(L::FiniteCoxeterGroup) -> Bool

Check if a Coxeter group is a Levi subgroup of its parent.
A group is a Levi subgroup if its inclusion generators are a subset of its parent's.
Returns true if L has no parent.

# Examples
```julia
G = CoxeterGroup(:A, 4)
H = reflection_subgroup(G, [1, 2])
islevi(H)  # returns true
```
"""
function islevi(L::FiniteCoxeterGroup)
    parent_inclusion_gens = try
        inclusiongens(L.parent)
    catch
        return true  # No parent
    end
    
    return issubset(inclusiongens(L), parent_inclusion_gens)
end

"""
    ispalindromic(f::Pol{BigInt}) -> Bool

Check if a polynomial is palindromic.
A polynomial is palindromic if:
1. It is zero (f(e) = 0), or
2. It has a non-zero constant term and its coefficients are symmetric

"""
function ispalindromic(f::Pol{BigInt})
    iszero = f(exp(1)) == 0
    has_nonzero_constant = f(0) != 0
    coeffs_symmetric = f.c == f.c[end:-1:1]
    
    return iszero || (has_nonzero_constant && coeffs_symmetric)
end

"""
    isnonnegative(f::Pol{BigInt}) -> Bool

Check if a polynomial has all non-negative coefficients.

"""
function isnonnegative(f::Pol{BigInt})
    return all(x -> x >= 0, f.c)
end