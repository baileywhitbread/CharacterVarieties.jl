"""
    GType

Represents a G-type, which is a pair [L,ρ] where:
- L is an endoscopy group of G containing T
- ρ is a principal unipotent character of L(Fq)

The character ρ is recorded as a string representation and its degree.
For example, the Steinberg character of GL3 is recorded as "3" 
(for the partition (3,0,...) of 3) and q^3.
"""
struct GType
    endoscopy::FiniteCoxeterGroup
    character::String
    degree::Pol{Rational{BigInt}}
    orbit_size::BigInt
    nu::BigInt
    
    """
        GType(endoscopy, character, degree, orbit_size, nu)
    
    Construct a new G-type.
    """
    function GType(endoscopy::FiniteCoxeterGroup, 
                  character::String, 
                  degree::Pol{Rational{BigInt}}, 
                  orbit_size::Union{BigInt,String}, 
                  nu::Union{BigInt,String})
        # Convert String "???" to BigInt for compatibility with group_types_no_data
        actual_orbit_size = orbit_size isa String ? BigInt(0) : orbit_size
        actual_nu = nu isa String ? BigInt(0) : nu
        
        new(endoscopy, character, degree, actual_orbit_size, actual_nu)
    end
end

# Make GTypes display nicely on the REPL
Base.show(io::IO, tau::GType) = print(io, "[$(tau.endoscopy),$(tau.character)]")

"""
    group_types(G::FiniteCoxeterGroup) -> Vector{GType}

Compute all G-types of the given Coxeter group G.
"""
function group_types(G::FiniteCoxeterGroup)
    G_dual = rootdatum(simplecoroots(G), simpleroots(G))
    G_dual_iplevis = iplevis(G_dual)
    G_dual_plevis = plevis(G_dual)
    
    return mapreduce(vcat, plorbit_reps(G_dual)) do plevi
        get_plevi_types(plevi, G_dual_iplevis, G_dual_plevis)
    end
end

"""
    group_types_no_data(G::FiniteCoxeterGroup) -> Vector{GType}

Compute G-types without computing orbit sizes and nu values.
Intended for quick inspection of GTypes.
"""
function group_types_no_data(G::FiniteCoxeterGroup)
    G_dual = rootdatum(simplecoroots(G), simpleroots(G))
    
    return mapreduce(vcat, plorbit_reps(G_dual)) do plevi
        get_plevi_types(plevi, nothing, nothing)
    end
end

"""
    get_plevi_types(plevi::FiniteCoxeterGroup, 
                    iplevis::Union{Vector{FiniteCoxeterGroup}, Nothing}, 
                    plevis::Union{Vector{FiniteCoxeterGroup}, Nothing}) -> Vector{GType}

Helper function to compute G-types for a given pseudo-Levi subgroup.
"""
function get_plevi_types(plevi::FiniteCoxeterGroup, 
                        iplevis::Union{Vector{FiniteCoxeterGroup}, Nothing}, 
                        plevis::Union{Vector{FiniteCoxeterGroup}, Nothing})
    uc = UnipotentCharacters(plevi)
    names = charnames(uc, limit=true)
    degs = degrees(uc)
    
    orbit_size = isnothing(iplevis) ? "???" : BigInt(length(myorbit(plevi)))
    nu_value = isnothing(plevis) ? "???" : BigInt(nu(plevi, iplevis, plevis))
    
    return [GType(plevi, names[i], degs[i], orbit_size, nu_value)
            for i in 1:length(uc)
            if BigInt(degs[i](1)) != 0]
end

"""
    group_type_data(G::FiniteCoxeterGroup) -> Matrix{Any}

Compute a matrix containing detailed data about all G-types of G.
Used by epolys.jl and tables.jl.

# Returns
A matrix where each row contains:
1. GType
2. |Φ(L)⁺| (number of positive roots)
3. |L(Fq)| (order polynomial)
4. ρ(1) (character degree)
5. φ(1) (evaluation at 1)
6. |W(L)| (order of Weyl group)
7. |[L]| (orbit size)
8. ν(L) (nu invariant)
"""
function group_type_data(G::FiniteCoxeterGroup)
    return mapreduce(vcat, group_types(G)) do type
        [type                                      # Type
         BigInt(length(roots(type.endoscopy))/2)  # |Φ(L)⁺|
         orderpol(type.endoscopy)                 # |L(Fq)|
         type.degree                              # ρ(1)
         BigInt(type.degree(1))                   # φ(1)
         length(type.endoscopy)                   # |W(L)|
         type.orbit_size                          # |[L]|
         type.nu]                                 # ν(L)
    end |> data -> sortslices(data, dims=1, by=x->x[2], rev=true)
end