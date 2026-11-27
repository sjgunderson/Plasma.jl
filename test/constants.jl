@testset "AtomDB physical constants" begin
    mktempdir() do dir
        prefix = write_test_apec(joinpath(dir, "constants"))
        legacy = AtomDB(prefix; constants = :legacy)
        db_type = AtomDB{Float32,Float64,Float32}
        codata = Plasma.load(db_type, legacy.linefile, legacy.cocofile; constants = :codata2022)
        @test AtomDB(prefix).constants === codata.constants
        @test Plasma.load(AtomDB, prefix).constants === codata.constants
        @test Plasma.load(db_type, legacy.linefile, legacy.cocofile).constants === codata.constants
        @test Plasma.load(AtomDB, prefix; constants = :legacy).constants === legacy.constants
        @test_throws ArgumentError AtomDB(prefix; constants = :unknown)

        @test legacy.constants.kev_angstrom == 12.39841974
        @test isapprox(codata.constants.kev_angstrom, 12.398419843320026; rtol = 1e-14)
        @test isapprox(codata.constants.kev_per_kelvin, 8.61733326214518e-8; rtol = 1e-14)
        @test isapprox(codata.constants.kev_per_atomic_mass, 931494.103718569; rtol = 1e-14)
        @test isapprox(codata.constants.light_speed_km_s, 299792.458; rtol = 1e-14)
        thermal_fe = Plasma.line_width_fractions(1.0, true, 0.0, codata.constants)[26]
        @test isapprox(thermal_fe, sqrt(1 / (55.845 * 931494.103718569)); rtol = 1e-14)

        matches = map((legacy, codata)) do db
            row = only(find_lines(db, 1.99, 2.01; element = :Fe, kT = 1.0))
            @test row.wavelength ≈ 2.0
            @test only(find_lines(db, row.energy, row.energy + 0.01; unit = :keV, kT = 1.0)).id == row.id
            @test line_emissivity(db, row, 1.0 / 8.61733326214518e-8; teunit = :K) ≈ 3.0
            row
        end
        old_line, new_line = matches
        @test old_line.id == new_line.id
        @test old_line.energy < new_line.energy
        edges = [5.5, (old_line.energy + new_line.energy) / 2, 7.0]
        old_flux, new_flux = map(db -> apec(edges, db; kT = 1.0), (legacy, codata))
        @test old_flux .- new_flux ≈ [3e14, -3e14]
        @test sum(old_flux) ≈ sum(new_flux)

        edges = [1.0, 4.0, 8.0, 13.0]
        flow = (; lowT = 1.0, highT = 3.0, ntemps = 2)
        distance_mpc = 10.0
        norm = inv(1e-14 * (1e6 * distance_mpc)^2)
        spectra = map((legacy, codata)) do db
            flux = coolflow(edges, db; flow..., distance_mpc, redshift = 0.0)
            @test localcoolflow(edges, db; flow..., norm) ≈ flux
            flux
        end
        mass_ratio = 1.66053906660 / 1.66053906892
        @test isapprox(spectra[2], mass_ratio .* spectra[1]; rtol = 1e-12)
    end
end
