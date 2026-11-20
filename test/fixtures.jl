using FITSIO: FITS

function write_test_apec(prefix; continuum_temperatures = Float32[1, 3])
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
            # The third entry of each FITS row is padding. The CIE rows are H, Fe, and Na.
            energy = repeat(Float32[1, 14, 0], 1, 3)
            counts = Int32[2, 2, 2]
            continuum = scale .* Float32[2 1 0; 2 1 0; -1 -1 -1]
            pseudo = scale .* Float32[0 0 1; 0 0 1; -1 -1 -1]
            write(fits, ["Z", "rmJ", "N_Cont", "E_Cont", "Continuum", "N_Pseudo", "E_Pseudo", "Pseudo"],
                Any[Int32[1, 26, 11], zeros(Int32, 3), counts, energy, continuum, counts, energy, pseudo];
                name = "EMISSIVITY", ver = index)
        end
    end
    return prefix
end
