using CharacterVarieties
using Test



@testset "CharacterVarieties.jl" begin
    @test length(plorbit_reps(rootdatum(:gl,2))) == 2
end
