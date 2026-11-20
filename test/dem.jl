@testset "XSPEC DEM temperature sampling" begin
    t, w = Plasma.wdem_sampling(3.1, 0.1, 0.25)
    @test length(t) == length(w) == 51
    @test isapprox(t[[1, 2, 3, 50, 51]], [0.3379, 0.3937, 0.4495, 3.0721, 3.1279])
    @test isapprox(w[[1, 2, 3, 50, 51]],
        [1.1778751e-5, 2.163042e-5, 3.667029e-5, 0.079409753, 0.085337783]; rtol = 1e-7)
    @test isapprox(sum(w), 1)
    t, w = Plasma.wdem_sampling(3.1, 0.37, -1.0)
    @test isapprox(t[[1, 51]], [1.16653, 3.11953])
    @test isapprox(w[[1, 51]], [0.033261751, 0.012437038]; rtol = 1e-7)
    @test Plasma.wdem_sampling(3.1, 1.0, 0.25) == ((3.1,), (1.0,))
    @test all(isfinite, last(Plasma.wdem_sampling(3.1, 0.1, 1e-12)))
    @test all(isfinite, last(Plasma.wdem_sampling(3.1, 0.1, -1e-12)))
    @test all(isfinite, last(Plasma.wdem_sampling(3.1, 0.1, 0.0)))
    @test_throws DomainError Plasma.wdem_sampling(3.1, 0.0, 1.0)

    t, w = Plasma.cempow_sampling(0.4, 3.1)
    @test length(t) == length(w) == 21
    @test isapprox(t[[1, 2, 21]], [3.1, 2.4624175, 0.031])
    @test isapprox(sum(w), 0.9722288707647349)
    @test isapprox(sum(last(Plasma.cempow_sampling(2.5, 3.1))), 0.22848727424897336)
    @test Plasma.cempow_sampling(1.0, 0.02725) == ([0.02725], [0.1])
    @test_throws DomainError Plasma.cempow_sampling(0.0, 1.0)
end

@testset "Multi-temperature AtomDB evaluation" begin
    mktempdir() do dir
        db = AtomDB(write_test_apec(joinpath(dir, "dem"); variable_continuum = true))
        edges = collect(range(1.0, 13.0; length = 241))
        t, w = Plasma.wdem_sampling(2.7, 0.3, 0.5)
        expected = sum(measure .* bapec(edges, db; kT) for (kT, measure) in zip(t, w))
        @test isapprox(bwdem(edges, db; Tmax = 2.7, beta = 0.3, inv_slope = 0.5), expected)
        @test isapprox(bvwdem(edges, db; Tmax = 2.7, beta = 0.3, inv_slope = 0.5), expected)
        @test isapprox(bvvwdem(edges, db; Tmax = 2.7, beta = 0.3, inv_slope = 0.5), expected)

        @test isapprox(wdem(edges, db; Tmax = 2.0, beta = 1.0), apec(edges, db; kT = 2.0))
        @test isapprox(wdem(edges, db; Tmax = 2.0, abundance = 0.3),
            vwdem(edges, db; Tmax = 2.0, Fe = 0.3, Na = 0.3))
        @test isapprox(vwdem(edges, db; Tmax = 2.0, Fe = 0.3, Na = 0.5, abundance_table = :wilm),
            vvwdem(edges, db; Tmax = 2.0, Fe = 0.3, Na = 0.5, abundance_table = :wilm))
        @test !isapprox(vvwdem(edges, db; Tmax = 2.0, H = 0.0),
            vvwdem(edges, db; Tmax = 2.0, H = 1.0))

        t, w = Plasma.cempow_sampling(0.6, 2.7)
        expected = sum(measure .* apec(edges, db; kT) for (kT, measure) in zip(t, w))
        @test isapprox(cempow(edges, db; alpha = 0.6, Tmax = 2.7), expected)
        @test isapprox(vcempow(edges, db; alpha = 0.6, Tmax = 2.7), expected)
        @test isapprox(bcempow(edges, db; alpha = 0.6, Tmax = 2.7, thermal_broadening = false), expected)
        @test isapprox(bcempow(edges, db; alpha = 0.6, Tmax = 2.7, velocity = 300.0),
            bvcempow(edges, db; alpha = 0.6, Tmax = 2.7, velocity = 300.0))
        @test isapprox(cempow(edges, db; nH = 11.0), cempow(edges, db))
        @test_throws ArgumentError cempow(edges, db; switch = 3)
        @test_throws ArgumentError cempow(edges, db; Fe = 0.3)
        @test_throws ArgumentError vwdem(edges, db; H = 0.3)
    end
end
