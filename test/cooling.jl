@testset "AtomDB cooling-flow weights" begin
    mktempdir() do dir
        db = AtomDB(write_test_apec(joinpath(dir, "cooling"); variable_continuum = true))
        edges = collect(range(1.0, 14.0; length = 301))
        distance_mpc = 10.0
        t1, t2 = 1.0, 3.0
        area = (14^2 - 1^2) / 2
        e1 = Plasma.KEV_ANGSTROM / 4
        lum = 4 * area + (2 + 2 * 3 + 4 * 5) * e1
        @test all(isnothing, db.bolometric) && all(isnothing, db.direct_bolometric)
        totals = Plasma.bolometric_elements(db, 1, :direct)
        @test isapprox(sum(totals), lum)
        @test isapprox(totals[11], area + 5 * 4 * e1)
        @test db.direct_bolometric[1] === totals && isnothing(db.direct_bolometric[2])
        binned_flux = apec(Plasma.BOLOMETRIC_EDGES, db; kT = t1)
        binned = Plasma.bolometric_elements(db, 1)
        @test sum(binned) ≈ sum(Plasma.BOLOMETRIC_CENTERS .* binned_flux) / 1e14
        @test db.bolometric[1] === binned && db.direct_bolometric[1] === totals

        normalization = 1e-14 * Plasma.FLOW_PARTICLES_PER_SECOND / (4π * (distance_mpc * Plasma.MPC_CM)^2)
        expected = apec(edges, db; kT = t1) .* (4 * normalization * log(3) / (2 * lum))
        args = (; lowT = t1, highT = t2, redshift = 0.0, distance_mpc, ntemps = 2)
        plain = coolflow(edges, db; args...)
        @test coolflow(edges, db; args..., bolometric = :xspec) == plain
        @test isapprox(coolflow(edges, db; args..., bolometric = :direct), expected; atol = 0.0)
        @test isapprox(coolflow(edges, db; args..., norm = 4.5), 4.5 .* plain; atol = 0.0)
        @test coolflow(edges, db; args..., highT = t1) == zeros(length(edges) - 1)
        @test coolflow(edges, db; args..., norm = 0.0) == zeros(length(edges) - 1)
        @test mkcflow === coolflow && vmcflow === vcoolflow
        @test isapprox(vcoolflow(edges, db; args...), plain; atol = 0.0)
        @test isapprox(bcoolflow(edges, db; args..., thermal_broadening = false), plain; atol = 0.0)
        @test isapprox(bcoolflow(edges, db; args..., velocity = 300.0),
            bvcoolflow(edges, db; args..., velocity = 300.0); atol = 0.0)
        @test !isapprox(vcoolflow(edges, db; args..., Na = 0.0), plain; atol = 0.0)
        @test_throws ArgumentError coolflow(edges, db; switch = 1)
        @test_throws ArgumentError coolflow(edges, db; ntemps = 1)
        @test_throws DomainError coolflow(edges, db; redshift = 0.0)
        @test_throws DomainError coolflow(edges, db; lowT = 5.0, highT = 4.0)
        @test_throws DomainError coolflow(edges, db; distance_mpc = -1.0)
        @test_throws ArgumentError coolflow(edges, db; bolometric = :unknown)
        @test_throws DomainError coolflow(edges, db; normalization_redshift = 0.0)

        z = 0.08
        d = Plasma.cooling_luminosity_distance(z, CoolingCosmology())
        @test isapprox(d, 363.91869824; rtol = 1e-7)
        @test isapprox(Plasma.cooling_luminosity_distance(z, CoolingCosmology(lambda0 = 0.0)),
            Plasma.LIGHT_SPEED_KM_S / 70 * z * (1 + z / 2))
        @test isapprox(coolflow(edges, db; redshift = z, distance_mpc = d), coolflow(edges, db; redshift = z))
        @test isapprox(coolflow(edges, db; redshift = z, cosmology = CoolingCosmology(H0 = 140.0)),
            4 .* coolflow(edges, db; redshift = z))

        shift = 0.25
        shifted_edges = edges ./ (1 + shift)
        @test coolflow(shifted_edges, db; args..., redshift = shift, normalization_redshift = 0.0) ≈ plain
        @test coolflow(shifted_edges, db; args..., redshift = shift) ≈ plain ./ (1 + shift)
        @test coolflow(shifted_edges, db; redshift = shift, normalization_redshift = z) ≈
            coolflow(edges ./ (1 + z), db; redshift = z)

        peakT = first(Plasma.cooling_temperatures(Plasma.CPH_LIMITS..., 10))[6]
        heating = cph(edges, db; peakT, distance_mpc, redshift = 0.0)
        @test all(isfinite, heating) && sum(heating) > 0
        @test isapprox(vcph(edges, db; peakT, distance_mpc, redshift = 0.0), heating)
        @test isapprox(bcph(edges, db; peakT, distance_mpc, redshift = 0.0, velocity = 300.0),
            bvcph(edges, db; peakT, distance_mpc, redshift = 0.0, velocity = 300.0))
        @test_throws DomainError cph(edges, db; peakT = 0.0)
    end
end
