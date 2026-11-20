using Plasma
using Test

include("fixtures.jl")

@testset "Plain AtomDB CIE APEC" begin
    mktempdir() do dir
        prefix = write_test_apec(joinpath(dir, "apec"))
        db = AtomDB(prefix)
        @test db.temperatures == [1.0, 3.0]
        @test all(isnothing, db.slices)

        edges = [1.0, 4.0, 8.0, 13.0]
        fixed = [11.0, 12.0, 20.0]
        metals = [3.0, 7.0, 5.0]
        solar = 1e14 .* (fixed .+ metals)
        @test isapprox(apec(edges, db; kT = 1.0), solar)
        @test !isnothing(db.slices[1]) && isnothing(db.slices[2])
        first_slice = db.slices[1]
        @test isapprox(apec(edges, db; kT = 2.0, abundance = 0.25), 2e14 .* (fixed .+ 0.25 .* metals))
        @test db.slices[1] === first_slice
        @test isapprox(apec(edges, db; kT = 3.0), 3 .* solar)
        @test isapprox(apec(edges, db; kT = 0.5), solar)
        @test isapprox(apec(edges, db; kT = 5.0), 3 .* solar)
        @test isapprox(apec(edges, db; abundance = 0), 1e14 .* fixed)
        @test isapprox(apec(edges ./ 2, db; redshift = 1), solar ./ 2)
        @test apec(edges, db; norm = 0.0) == zeros(3)

        out = fill(NaN, length(edges) - 1)
        @test apec!(out, edges, db; norm = 2) === out
        @test isapprox(out, 2 .* solar)

        # A line coincident with an internal edge belongs to the lower XSPEC bin.
        line_edges = [2.0, Plasma.KEV_ANGSTROM / 4, 4.0]
        @test isapprox(apec(line_edges, db; abundance = 0), 1e14 .* (3 .* diff(line_edges) .+ [2, 0]))
        @test apec([20.0, 30.0], db) == [0.0]

        @test_throws DimensionMismatch apec!(zeros(2), edges, db)
        @test_throws ArgumentError apec([2.0, 2.0], db)
        @test_throws ArgumentError apec([2.0], db)
        @test_throws DomainError apec(edges, db; redshift = -1)
        @test_throws DomainError apec(edges, db; abundance = -1)
        @test_throws DomainError apec(edges, db; kT = NaN)
        @test_throws ArgumentError AtomDB(joinpath(dir, "missing"))
        mismatched = write_test_apec(joinpath(dir, "mismatched"); continuum_temperatures = Float32[1, 4])
        @test_throws ArgumentError AtomDB(mismatched)
    end
end

@testset "Continuum sampled at bin midpoints" begin
    # The first midpoint equals the first node. The second bin straddles an
    # internal node, where midpoint sampling differs from an exact integral.
    continuum = Plasma.Continuum(Float32[2, 3, 4], Float32[1, 3, 1])
    flux = zeros(3)
    Plasma.add_continuum!(flux, [2.0, 3.5, 5.0], [2.0, 2.0, 1.0], continuum, 1.0)
    @test isapprox(flux, [2.0, 4.0, 0.0])
end
