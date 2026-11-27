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
        line_edges = [2.0, db.constants.kev_angstrom / 4, 4.0]
        @test isapprox(apec(line_edges, db; abundance = 0), 1e14 .* (3 .* diff(line_edges) .+ [2, 0]))
        @test apec([20.0, 30.0], db) == [0.0]

        @test_throws DimensionMismatch apec!(zeros(2), edges, db)
        @test_throws ArgumentError apec([2.0, 2.0], db)
        @test_throws ArgumentError apec([2.0], db)
        @test_throws DomainError apec(edges, db; redshift = -1)
        @test_throws DomainError apec(edges, db; abundance = -1)
        @test_throws DomainError apec(edges, db; kT = NaN)
        missing_prefix = joinpath(dir, "missing")
        @test_throws ArgumentError AtomDB(missing_prefix)
        @test_throws ArgumentError Plasma.load(AtomDB, missing_prefix)
        mismatched = write_test_apec(joinpath(dir, "mismatched"); continuum_temperatures = Float32[1, 4])
        @test_throws ArgumentError AtomDB(mismatched)
        @test_throws ArgumentError Plasma.load(AtomDB, mismatched * "_line.fits", mismatched * "_coco.fits")
    end
end

@testset "Variable-length AtomDB continuum rows" begin
    mktempdir() do dir
        fixed = AtomDB(write_test_apec(joinpath(dir, "fixed")))
        variable = AtomDB(write_test_apec(joinpath(dir, "variable"); variable_continuum = true))
        edges = [1.0, 4.0, 8.0, 13.0]
        for kT in (1.0, 2.0, 3.0), abundance in (0.0, 1.0)
            @test isapprox(apec(edges, variable; kT, abundance), apec(edges, fixed; kT, abundance))
        end
    end
end

@testset "Parametric AtomDB storage" begin
    mktempdir() do dir
        edges = [1.0, 4.0, 8.0, 13.0]
        for variable_continuum in (false, true)
            prefix = write_test_apec(joinpath(dir, string(variable_continuum)); variable_continuum)
            reference = apec(edges, AtomDB(prefix); kT = 2.0)
            for (db, C, E, F) in (
                    (AtomDB(prefix), Float32, Float64, Float32),
                    (AtomDB{Float32,Float32,Float32}(prefix), Float32, Float32, Float32),
                    (AtomDB{Float64,Float64,Float64}(prefix * "_line.fits", prefix * "_coco.fits"),
                        Float64, Float64, Float64))
                @test all(isnothing, db.slices) && all(isnothing, db.line_lists)
                lines = Plasma.temperature_lines(db, 1)
                @test lines isa Plasma.Lines{E,F}
                @test eltype(db.line_lists) === Union{Nothing,Plasma.Lines{E,F}}
                @test isnothing(db.slices[1])
                slice = Plasma.temperature_slice(db, 1)
                @test slice.lines === lines
                @test eltype(slice.continua) === Tuple{Int,Plasma.Continuum{C}}
                @test eltype(db.slices) === Union{Nothing,typeof(slice)}
                @test apec(edges, db; kT = 2.0) ≈ reference
            end
        end
    end
end

@testset "Load AtomDB by type" begin
    mktempdir() do dir
        edges = [1.0, 4.0, 8.0, 13.0]
        for variable_continuum in (false, true)
            prefix = write_test_apec(joinpath(dir, string(variable_continuum)); variable_continuum)
            for db_type in (AtomDB, AtomDB{Float64,Float64,Float64})
                reference = db_type(prefix)
                expected = apec(edges, reference; kT = 2.0)
                from_prefix = Plasma.load(db_type, prefix)
                from_paths = Plasma.load(db_type, reference.linefile, reference.cocofile)
                for db in (from_prefix, from_paths)
                    @test typeof(db) === typeof(reference)
                    @test db.temperatures == reference.temperatures
                    @test all(isnothing, db.slices) && all(isnothing, db.line_lists)
                    @test apec(edges, db; kT = 2.0) == expected
                end
            end
        end
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

@testset "APEC variants and elemental abundances" begin
    mktempdir() do dir
        db = AtomDB(write_test_apec(joinpath(dir, "variants"); variable_continuum = true))
        edges = [1.0, 4.0, 8.0, 13.0]
        H, Fe, Na = [8.0, 8.0, 10.0], [3.0, 7.0, 5.0], [3.0, 4.0, 10.0]
        @test isapprox(vapec(edges, db; kT = 1.0, Fe = 0.0), 1e14 .* (H .+ Na))
        @test isapprox(vvapec(edges, db; kT = 1.0, H = 0.0, Na = 0.0), 1e14 .* Fe)
        @test isapprox(vapec(edges, db; kT = 1.0, Fe = 0.3, trace_abundance = :Fe),
            1e14 .* (H .+ 0.3 .* (Fe .+ Na)))
        @test isapprox(apec(edges, db; kT = 1.0, abundance = 0.3, trace_abundance = :Fe),
            1e14 .* (H .+ 0.3 .* (Fe .+ Na)))
        @test_throws ArgumentError vapec(edges, db; Na = 0.2)
        @test_throws ArgumentError vapec(edges, db; trace_abundance = :Na)
        @test_throws DomainError bapec(edges, db; velocity = -1.0)

        angr, wilm = solar_abundances(:angr), solar_abundances(:wilm)
        expected = 1e14 .* (H .+ Fe .* (wilm.Fe / angr.Fe) .+ Na .* (wilm.Na / angr.Na))
        @test isapprox(apec(edges, db; kT = 1.0, abundance_table = :wilm), expected)
        @test isapprox(vvapec(edges, db; kT = 1.0, abundance_table = :wilm), expected)
        path = joinpath(dir, "abundances.txt")
        write(path, join(values(wilm), '\n') * "\n")
        @test solar_abundances(path) == wilm
        @test isapprox(apec(edges, db; kT = 1.0, abundance_table = path), expected)

        fine = collect(range(1.0, 13.0; length = 12_001))
        plain = apec(fine, db; kT = 2.0)
        @test bapec(fine, db; kT = 2.0, thermal_broadening = false) == plain
        thermal = bapec(fine, db; kT = 2.0)
        broadened = bapec(fine, db; kT = 2.0, velocity = 300.0)
        @test !isapprox(thermal, plain)
        @test !isapprox(broadened, thermal)
        @test isapprox(sum(broadened), sum(plain))
        @test isapprox(bvapec(fine, db; kT = 2.0, velocity = 300.0), broadened)
        @test isapprox(bvvapec(fine, db; kT = 2.0, velocity = 300.0), broadened)
    end
end

@testset "XSPEC abundance tables and Gaussian line integral" begin
    for table in (:angr, :aspl, :feld, :aneb, :grsa, :wilm, :lodd, :lpgp, :lpgs, :felc)
        row = solar_abundances(table)
        @test length(row) == 30 && row.H == 1 && all(>=(0), row)
    end
    @test solar_abundances(:angr).Fe == 4.68e-5
    @test solar_abundances(:wilm).Zn == 0
    @test solar_abundances(:felc).Fe == 1.26e-4
    @test_throws ArgumentError solar_abundances(:unknown)

    constants = Plasma.LEGACY_CONSTANTS
    thermal = Plasma.line_width_fractions(1.0, true, 0.0, constants)
    mixed = Plasma.line_width_fractions(1.0, true, 300.0, constants)
    @test isapprox(thermal[1] / thermal[26], sqrt(Plasma.ATOMIC_MASSES[26] / Plasma.ATOMIC_MASSES[1]))
    @test isapprox(mixed[26], hypot(thermal[26], 300.0 / constants.light_speed_km_s))
    lines = Plasma.Lines([5.0], Float32[2.0], UInt8[26])
    out = zeros(4)
    Plasma.add_broadened_lines!(out, [0.0, 4.5, 5.0, 5.5, 10.0], lines,
        ntuple(_ -> 1.0, 30), 1.0, ntuple(_ -> 0.1, 30))
    @test isapprox(out, [0.3173105078629141, 0.6826894921370859, 0.6826894921370859, 0.3173105078629141])
end

include("constants.jl")
include("dem.jl")
include("cooling.jl")
include("line_search.jl")
