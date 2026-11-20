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

@testset "SpectralFitting DEM families" begin
    mktempdir() do dir
        db = AtomDB(write_test_apec(joinpath(dir, "dem")))
        edges = collect(range(1.0, 13.0; length = 201))
        for (constructor, evaluator, names, kwargs) in (
                (Wdem, wdem, (:Tmax, :beta), (; abundance = 0.6)),
                (Vwdem, vwdem, (:Tmax, :beta), (; Fe = 0.4, Na = 0.3)),
                (Vvwdem, vvwdem, (:Tmax, :beta), (; H = 0.7, Zn = 0.0)),
                (Bwdem, bwdem, (:Tmax, :beta), (; velocity = 200.0)),
                (Bvwdem, bvwdem, (:Tmax, :beta), (; velocity = 200.0, Na = 0.3)),
                (Bvvwdem, bvvwdem, (:Tmax, :beta), (; velocity = 200.0, Na = 0.3)),
                (Cempow, cempow, (:alpha, :Tmax), (; alpha = 0.4, abundance = 0.6)),
                (Vcempow, vcempow, (:alpha, :Tmax), (; Fe = 0.4, Na = 0.3)),
                (Bcempow, bcempow, (:alpha, :Tmax), (; velocity = 200.0)),
                (Bvcempow, bvcempow, (:alpha, :Tmax), (; velocity = 200.0, Na = 0.3)))
            model = constructor(db; Tmax = 2.5, K = 0.4, abundance_table = :wilm, kwargs...)
            expected = evaluator(edges, db; Tmax = 2.5, norm = 0.4, abundance_table = :wilm, kwargs...)
            @test SpectralFitting.parameter_names(model)[1:3] == (:K, names...)
            @test isapprox(invokemodel(edges, model), expected)
            @test isapprox(invokemodel(edges, model + model), 2 .* expected)
        end
        @test Bvvwdem(db; K = 0.4, K_abundance = 0.3).K_abundance.value == 0.3
        @test isfinite(SpectralFitting.ForwardDiff.derivative(0.37) do beta
            sum(wdem(edges, db; beta, Tmax = 2.5, inv_slope = -1.0))
        end)
        @test isfinite(SpectralFitting.ForwardDiff.derivative(0.7) do alpha
            sum(cempow(edges, db; alpha, Tmax = 2.5))
        end)
    end
end

@testset "SpectralFitting cooling-flow families" begin
    mktempdir() do dir
        db = AtomDB(write_test_apec(joinpath(dir, "cooling")))
        edges = collect(range(1.0, 13.0; length = 201))
        distance_mpc, redshift, K = 11.0, 0.02, 0.4
        for (constructor, evaluator, names, kwargs) in (
                (Coolflow, coolflow, (:lowT, :highT), (; abundance = 0.6)),
                (Vcoolflow, vcoolflow, (:lowT, :highT), (; Fe = 0.4, Na = 0.3)),
                (Bcoolflow, bcoolflow, (:lowT, :highT), (; velocity = 200.0)),
                (Bvcoolflow, bvcoolflow, (:lowT, :highT), (; velocity = 200.0, Na = 0.3)),
                (Cph, cph, (:peakT,), (; peakT = 1.5)),
                (Vcph, vcph, (:peakT,), (; Fe = 0.4, Na = 0.3)),
                (Bcph, bcph, (:peakT,), (; velocity = 200.0)),
                (Bvcph, bvcph, (:peakT,), (; velocity = 200.0, Na = 0.3)))
            model = constructor(db; K, distance_mpc, redshift, abundance_table = :wilm, kwargs...)
            expected = evaluator(edges, db; norm = K, distance_mpc, redshift, abundance_table = :wilm, kwargs...)
            @test SpectralFitting.parameter_names(model)[1:(1 + length(names))] == (:K, names...)
            copied = copy(model)
            @test copied.db === model.db && copied.K !== model.K
            @test isapprox(invokemodel(edges, model), expected)
            @test isapprox(invokemodel(edges, model + model), 2 .* expected)
        end
        @test Mkcflow === Coolflow && Vmcflow === Vcoolflow
        @test isfinite(SpectralFitting.ForwardDiff.derivative(0.08) do z
            sum(coolflow(edges, db; redshift = z))
        end)
        @test isfinite(SpectralFitting.ForwardDiff.derivative(0.7) do lowT
            sum(vcoolflow(edges, db; lowT))
        end)
    end
end
