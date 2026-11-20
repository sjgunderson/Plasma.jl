@testset "Search loaded AtomDB lines" begin
    mktempdir() do dir
        db = AtomDB(write_test_apec(joinpath(dir, "lines")))
        fe = LineID(:Fe, 17, 5, 1)
        rows = find_lines(db, 1.0, 4.0; kT = 1.0)
        @test [row.wavelength for row in rows] ≈ [1.0, 2.0, 4.0]
        @test [row.emissivity for row in rows] == [5.0, 3.0, 2.0]
        @test db.line_lists[1] !== nothing && db.line_lists[2] === nothing
        @test all(isnothing, db.slices)

        match = only(find_lines(db, 0.1, 0.4; unit = :nm, element = "fe", ion = 17, kT = 2.0))
        @test match.id == fe
        @test match.emissivity == 6.0
        @test only(find_lines(db, match.energy, match.energy + 0.01; unit = :keV, kT = 2.0)).id == fe
        @test isempty(find_lines(db, 1.0, 4.0; kT = 2.0, element = :Fe, ion = 16))
        @test only(find_lines(db, 1.0, 4.0; kT = 2.0, brightest = 1)).id.element == 11
        @test length(find_lines(db, 1.0, 4.0; kT = 2.0, min_emissivity = 6.5)) == 1

        @test line_emissivity(db, fe, [1.0, 2.0, 3.0]) == [3.0, 6.0, 9.0]
        @test line_emissivity(db, fe, 2.0; abundances = (Fe = 0.25,)) == 1.5
        h = LineID("H", 1, 2, 1)
        @test ratio_emissivity(db, fe, h, [1.0, 2.0]) == [1.5, 1.5]
        @test isnan(ratio_emissivity(db, fe, LineID(:O, 7, 2, 1), 2.0))
        @test_throws ArgumentError find_lines(db, 1.0, 4.0)
    end
end
