# XSPEC's wavelength-to-energy conversion in keV Å.
const KEV_ANGSTROM = 12.39841974
const APEC_METALS = (6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28)

struct Continuum
    energies::Vector{Float32}
    density::Vector{Float32}
end

struct Lines
    energies::Vector{Float64}
    emissivities::Vector{Float32}
end

struct TemperatureSlice
    fixed_continua::Vector{Continuum}
    metal_continua::Vector{Continuum}
    fixed_lines::Lines
    metal_lines::Lines
end

"""
    AtomDB(linefile, cocofile)
    AtomDB(prefix)

Open an AtomDB CIE emissivity table pair on the Anders & Grevesse (`angr`)
abundance scale. A prefix is expanded to `prefix_line.fits` and
`prefix_coco.fits`. The constructor reads the temperature grid in keV.
Each temperature's emissivities are read on demand and cached in this object.
"""
struct AtomDB
    linefile::String
    cocofile::String
    temperatures::Vector{Float64}
    slices::Vector{Union{Nothing,TemperatureSlice}}
    cache_lock::ReentrantLock

    function AtomDB(linefile::AbstractString, cocofile::AbstractString)
        linepath, cocopath = abspath(linefile), abspath(cocofile)
        line_temperatures = read_temperatures(linepath, "Nline")
        continuum_temperatures = read_temperatures(cocopath, "NCont")
        line_temperatures == continuum_temperatures ||
            throw(ArgumentError("AtomDB line and continuum files have different temperature grids"))
        slices = Vector{Union{Nothing,TemperatureSlice}}(nothing, length(line_temperatures))
        return new(linepath, cocopath, line_temperatures, slices, ReentrantLock())
    end
end

AtomDB(prefix::AbstractString) = AtomDB(prefix * "_line.fits", prefix * "_coco.fits")

function read_temperatures(path, count_column)
    isfile(path) || throw(ArgumentError("AtomDB file does not exist: $path"))
    return FITS(path, "r") do fits
        parameters = fits["PARAMETERS"]
        temperatures = Float64.(read(parameters, "kT"; case_sensitive = false))
        counts = read(parameters, count_column; case_sensitive = false)
        if isempty(temperatures) || !all(t -> isfinite(t) && t > 0, temperatures) ||
                !all(i -> temperatures[i] > temperatures[i - 1], 2:length(temperatures)) ||
                length(fits) != length(temperatures) + 2 || any(<(0), counts)
            throw(ArgumentError("Invalid AtomDB temperature table in $path"))
        end
        return temperatures
    end
end

function read_lines(path, index)
    return FITS(path, "r") do fits
        hdu = fits["EMISSIVITY", index]
        wavelength = read(hdu, "Lambda"; case_sensitive = false)
        emissivities = read(hdu, "Epsilon"; case_sensitive = false)
        elements = read(hdu, "Element"; case_sensitive = false)
        if !all(x -> isfinite(x) && x > 0, wavelength) ||
                !all(f -> isfinite(f) && f >= 0, emissivities) ||
                !all(z -> 1 <= z <= 30, elements)
            throw(ArgumentError("Invalid CIE line record at temperature index $index in $path"))
        end

        energies = KEV_ANGSTROM ./ wavelength
        order = sortperm(energies)
        fixed, metals = Lines(Float64[], Float32[]), Lines(Float64[], Float32[])
        for i in order
            lines = elements[i] in APEC_METALS ? metals : fixed
            push!(lines.energies, energies[i])
            push!(lines.emissivities, emissivities[i])
        end
        return fixed, metals
    end
end

# FITSIO returns a vector of row vectors for variable-length columns, and
# a numeric array for fixed-width columns (including scalar-width columns).
continuum_rows(column::AbstractVector{<:AbstractVector}, _nrows) = column
continuum_rows(column::AbstractArray{<:Real}, nrows) = eachcol(reshape(column, :, nrows))

function read_continua(path, index)
    return FITS(path, "r") do fits
        hdu = fits["EMISSIVITY", index]
        elements = read(hdu, "Z"; case_sensitive = false)
        ions = read(hdu, "rmJ"; case_sensitive = false)
        if !all(iszero, ions) || !all(z -> 1 <= z <= 30, elements) || !allunique(elements)
            throw(ArgumentError("Expected one CIE continuum row per element at temperature index $index in $path"))
        end
        fixed, metals = Continuum[], Continuum[]
        isempty(elements) && return fixed, metals
        for (count_column, energy_column, density_column) in (
                ("N_Cont", "E_Cont", "Continuum"), ("N_Pseudo", "E_Pseudo", "Pseudo"))
            counts = read(hdu, count_column; case_sensitive = false)
            energies = continuum_rows(read(hdu, energy_column; case_sensitive = false), length(elements))
            densities = continuum_rows(read(hdu, density_column; case_sensitive = false), length(elements))
            for i in eachindex(elements)
                n = counts[i]
                energy_row, density_row = energies[i], densities[i]
                if !(0 <= n <= length(energy_row) && n <= length(density_row))
                    throw(ArgumentError("Invalid $count_column for Z=$(elements[i]) at temperature index $index in $path: " *
                        "declared $n points, found $(length(energy_row)) energies and $(length(density_row)) emissivities"))
                end
                iszero(n) && continue
                e = Vector{Float32}(view(energy_row, 1:n))
                f = Vector{Float32}(view(density_row, 1:n))
                if !all(x -> isfinite(x) && x > 0, e) ||
                        !all(j -> e[j] > e[j - 1], 2:length(e)) ||
                        !all(x -> isfinite(x) && x >= 0, f)
                    throw(ArgumentError("Invalid $density_column data at temperature index $index in $path"))
                end
                push!(elements[i] in APEC_METALS ? metals : fixed, Continuum(e, f))
            end
        end
        return fixed, metals
    end
end

function temperature_slice(db::AtomDB, index::Int)
    return lock(db.cache_lock) do
        slice = db.slices[index]
        if isnothing(slice)
            fixed_lines, metal_lines = read_lines(db.linefile, index)
            fixed_continua, metal_continua = read_continua(db.cocofile, index)
            slice = TemperatureSlice(fixed_continua, metal_continua, fixed_lines, metal_lines)
            db.slices[index] = slice
        end
        return slice
    end
end
