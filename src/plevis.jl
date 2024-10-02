# Types of Levi subgroups
@enum LeviType begin
    Standard
    Pseudo
    IsolatedPseudo
end

"""
    orbit_representatives(G::FiniteCoxeterGroup, levi_type::LeviType=Pseudo)

Returns orbit representatives for the specified Levi type.
"""
function orbit_representatives(G::FiniteCoxeterGroup, levi_type::LeviType=Pseudo)
    all_reps = reflection_subgroup.(Ref(G), sscentralizer_reps(G))
    
    return filter(all_reps) do rep
        if levi_type == Standard
            return islevi(rep)
        elseif levi_type == IsolatedPseudo
            return isisolated(rep)
        else
            return true
        end
    end
end

"""
    calculate_orbits(G::FiniteCoxeterGroup, levi_type::LeviType=Pseudo)

Calculates orbits for the specified Levi type.
"""
function calculate_orbits(G::FiniteCoxeterGroup, levi_type::LeviType=Pseudo)
    reps = orbit_representatives(G, levi_type)
    raw_orbits = orbits(G, reps)
    
    return map(raw_orbits) do orbit
        unique_inclusions = collect(Set(sort.(inclusion.(orbit))))
        return reflection_subgroup.(Ref(G), unique_inclusions)
    end
end

"""
    get_all_levis(G::FiniteCoxeterGroup, levi_type::LeviType=Pseudo)

Returns all Levi subgroups of the specified type.
"""
function get_all_levis(G::FiniteCoxeterGroup, levi_type::LeviType=Pseudo)
    return reduce(vcat, calculate_orbits(G, levi_type))
end

# Functions to maintain compatibility
plorbit_reps(G::FiniteCoxeterGroup) = orbit_representatives(G, Pseudo)
iplorbit_reps(G::FiniteCoxeterGroup) = orbit_representatives(G, IsolatedPseudo)
lorbit_reps(G::FiniteCoxeterGroup) = orbit_representatives(G, Standard)

plorbits(G::FiniteCoxeterGroup) = calculate_orbits(G, Pseudo)
iplorbits(G::FiniteCoxeterGroup) = calculate_orbits(G, IsolatedPseudo)
lorbits(G::FiniteCoxeterGroup) = calculate_orbits(G, Standard)

plevis(G::FiniteCoxeterGroup) = get_all_levis(G, Pseudo)
iplevis(G::FiniteCoxeterGroup) = get_all_levis(G, IsolatedPseudo)
levis(G::FiniteCoxeterGroup) = get_all_levis(G, Standard)