using CharacterVarieties
using Test



@testset "CharacterVarieties.jl" begin
    @test EX(rootdatum(:gl,2),0,4)==Pol([1,4,1])
end
