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

@testset "SpectralFitting APEC variants" begin
    mktempdir() do dir
        db = AtomDB(write_test_apec(joinpath(dir, "variants")))
        edges = collect(range(1.0, 13.0; length = 401))
        for (constructor, evaluator, count, kwargs) in (
                (Apec, apec, 4, (; abundance = 0.4)),
                (Bapec, bapec, 5, (; abundance = 0.4, velocity = 250.0)),
                (Vapec, vapec, 16, (; Fe = 0.4, trace_abundance = :Fe)),
                (Vvapec, vvapec, 33, (; Fe = 0.4, Na = 0.5)),
                (Bvapec, bvapec, 17, (; Fe = 0.4, velocity = 250.0)),
                (Bvvapec, bvvapec, 34, (; Fe = 0.4, Na = 0.5, velocity = 250.0)))
            model = constructor(db; K = 0.5, kT = 2.0, abundance_table = :wilm, kwargs...)
            expected = evaluator(edges, db; norm = 0.5, kT = 2.0, abundance_table = :wilm, kwargs...)
            names = SpectralFitting.parameter_names(model)
            @test length(names) == count
            @test names[1:2] == (:K, :kT)
            copied = copy(model)
            @test copied.db === model.db && copied.K !== model.K
            @test isapprox(invokemodel(edges, model), expected)
            @test isapprox(invokemodel(edges, model + model), 2 .* expected)
        end
        model = Bvvapec(db; K = 0.5, K_abundance = 0.6, Na = FitParam(0.5), velocity = FitParam(200.0))
        @test !model.Na.frozen && !model.velocity.frozen
        @test model.K.value == 0.5 && model.K_abundance.value == 0.6
        @test SpectralFitting.parameter_names(model)[3:32] ==
            map(name -> name === :K ? :K_abundance : name, keys(solar_abundances()))
        function redshifted_bin(z)
            return bvapec(edges, db; kT = 2.0, velocity = 500.0, redshift = z, Fe = 0.4)[70]
        end
        @test isfinite(SpectralFitting.ForwardDiff.derivative(redshifted_bin, 0.002))
    end
end
