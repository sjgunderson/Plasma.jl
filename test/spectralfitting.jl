using Plasma, SpectralFitting, Test

include("fixtures.jl")

@testset "SpectralFitting APEC model" begin
    mktempdir() do dir
        db = AtomDB(write_test_apec(joinpath(dir, "apec")))
        edges = [1.0, 4.0, 8.0, 13.0]
        K = FitParam(0.5)
        model = Apec(db; K, kT = 2.0, abundance = FitParam(0.4))
        expected = apec(edges, db; norm = 0.5, kT = 2.0, abundance = 0.4)
        @test model.K !== K
        @test typeof(Apec(db; K = FitParam(0.5f0)).K) === typeof(K)
        @test Apec(db).abundance.frozen && Apec(db).redshift.frozen
        @test SpectralFitting.parameter_names(model) == (:K, :kT, :abundance, :redshift)
        @test isapprox(invokemodel(edges, model), expected)
        @test isapprox(invokemodel(edges, model + model), 2 .* expected)
    end
end
