using FITSIO: FITS

function write_test_apec(prefix; continuum_temperatures = Float32[1, 3], variable_continuum = false)
    temperatures = Float32[1, 3]
    FITS(prefix * "_line.fits", "w") do fits
        write(fits, ["kT", "Nline"], Any[temperatures, Int32[3, 3]]; name = "PARAMETERS")
        for (index, scale) in enumerate(Float32[1, 3])
            write(fits, ["Lambda", "Epsilon", "Element"],
                Any[Float32[4, 2, 1], scale .* Float32[2, 3, 5], Int32[1, 26, 11]];
                name = "EMISSIVITY", ver = index)
        end
    end
    FITS(prefix * "_coco.fits", "w") do fits
        write(fits, ["kT", "NCont"], Any[continuum_temperatures, Int32[3, 3]]; name = "PARAMETERS")
        for (index, scale) in enumerate(Float32[1, 3])
            # The CIE rows are H, Fe, and Na.
            if variable_continuum
                continuum_energy = [Float32[1, 14], Float32[1, 7, 14], Float32[]]
                continuum = [fill(2 * scale, 2), fill(scale, 3), Float32[]]
                pseudo_energy = [Float32[], Float32[1, 14], Float32[1, 4, 9, 14]]
                pseudo = [Float32[], zeros(Float32, 2), fill(scale, 4)]
                continuum_counts = Int32.(length.(continuum_energy))
                pseudo_counts = Int32.(length.(pseudo_energy))
                varcols = ["E_Cont", "Continuum", "E_Pseudo", "Pseudo"]
            else
                # The third entry of each FITS row is padding.
                continuum_energy = pseudo_energy = repeat(Float32[1, 14, 0], 1, 3)
                continuum_counts = pseudo_counts = Int32[2, 2, 2]
                continuum = scale .* Float32[2 1 0; 2 1 0; -1 -1 -1]
                pseudo = scale .* Float32[0 0 1; 0 0 1; -1 -1 -1]
                varcols = nothing
            end
            write(fits, ["Z", "rmJ", "N_Cont", "E_Cont", "Continuum", "N_Pseudo", "E_Pseudo", "Pseudo"],
                Any[Int32[1, 26, 11], zeros(Int32, 3), continuum_counts, continuum_energy, continuum,
                    pseudo_counts, pseudo_energy, pseudo]; name = "EMISSIVITY", ver = index, varcols)
        end
    end
    return prefix
end
