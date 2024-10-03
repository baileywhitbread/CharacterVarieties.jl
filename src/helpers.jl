"""
    myorbit(L::FiniteCoxeterGroup) -> Vector{FiniteCoxeterGroup}

Compute the orbit of a Coxeter group L under its parent group's action.
Handles duplicate elimination that may occur in PermGroups.jl orbit() function.

# Returns
- If L has a parent: Vector of unique subgroups in the orbit
- If L has no parent: Single-element vector containing L
"""
function myorbit(L::FiniteCoxeterGroup)
    parent_group = try
        L.parent
    catch
        return [L]  # L has no parent, return singleton orbit
    end
    
    orbit_elements = orbit(parent_group, L)
    unique_inclusions = collect(Set(sort.(inclusion.(orbit_elements))))
    return reflection_subgroup.(Ref(parent_group), unique_inclusions)
end

"""
    orderpol(L::FiniteCoxeterGroup) -> Pol{BigInt}

Compute the order polynomial ||L||(q) = |L(Fq)|.
"""
function orderpol(L::FiniteCoxeterGroup)
    return PermRoot.generic_order(L, Pol(:q))
end

"""
    mobius(A::FiniteCoxeterGroup, B::FiniteCoxeterGroup, P::Vector{FiniteCoxeterGroup}) -> Int

Compute the Möbius function value for the interval [A,B] in the poset P.
P must be ordered by subgroup inclusion: A ≤ B iff A is a subgroup of B.

# Arguments
- A, B: Coxeter groups in P
- P: Vector of Coxeter groups forming a poset

# Returns
The Möbius function value μ(A,B)

# Throws
ArgumentError if A is not a subgroup of B
"""
function mobius(A::FiniteCoxeterGroup, B::FiniteCoxeterGroup, P::Vector{FiniteCoxeterGroup})
    inc_A = sort(inclusion(A))
    inc_B = sort(inclusion(B))
    
    if !issubset(inc_A, inc_B)
        throw(ArgumentError("First argument must be a subgroup of the second argument"))
    end
    
    if inc_A == inc_B
        return 1
    end
    
    mobius_value = 0
    for element in P
        inc_elem = sort(inclusion(element))
        if issubset(inc_A, inc_elem) && issubset(inc_elem, inc_B) && inc_elem != inc_B
            mobius_value += mobius(A, element, P)
        end
    end
    
    return -mobius_value
end

"""
    pi0(L::FiniteCoxeterGroup) -> Int

Compute |π₀(Z(L))|, the number of connected components of the algebraic center of L.
"""
function pi0(L::FiniteCoxeterGroup)
    return length(algebraic_center(L).AZ)
end

"""
    nu(L::FiniteCoxeterGroup, iso_plevis::Vector{FiniteCoxeterGroup}, 
       all_plevis::Vector{FiniteCoxeterGroup}) -> Int

Compute the ν invariant for a Coxeter group L.

# Arguments
- L: The Coxeter group to compute ν for
- iso_plevis: Vector of isolated pseudo-Levi subgroups
- all_plevis: Vector of all pseudo-Levi subgroups
"""
function nu(L::FiniteCoxeterGroup, 
            iso_plevis::Vector{FiniteCoxeterGroup}, 
            all_plevis::Vector{FiniteCoxeterGroup})
    return sum(
        mobius(L, iplevi, all_plevis) * pi0(iplevi)
        for iplevi in iso_plevis
        if issubset(L, iplevi)
    )
end

"""
    dimension_XY(G::FiniteCoxeterGroup, genus::Int, puncture::Int) -> Int

Compute the dimension of the moduli space XY for a given Coxeter group G,
genus, and number of punctures.
"""
function dimension_XY(G::FiniteCoxeterGroup, genus::Int, puncture::Int)
    return (2*genus - 2 + puncture) * dimension(G) + 
           2 * (rank(G) - semisimplerank(G)) - 
           puncture * rank(G)
end