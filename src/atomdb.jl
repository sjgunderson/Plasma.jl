# XSPEC's wavelength-to-energy conversion in keV Å.
const KEV_ANGSTROM = 12.39841974

struct Continuum
    energies::Vector{Float32}
    density::Vector{Float32}
end

struct Lines
    energies::Vector{Float64}
    emissivities::Vector{Float32}
    elements::Vector{UInt8}
end

"""
    LineID(element, ion, upper_level, lower_level)

Identify an AtomDB transition by atomic number or element symbol, spectroscopic
ion stage, and one-based level indices. For example, O VII has ion stage 7.
"""
struct LineID
    element::Int16
    ion::Int16
    upper_level::Int32
    lower_level::Int32

    function LineID(element, ion::Integer, upper_level::Integer, lower_level::Integer)
        z = line_atomic_number(element)
        1 <= ion <= z + 1 || throw(ArgumentError("Invalid spectroscopic ion stage $ion for element $z"))
        1 <= upper_level <= typemax(Int32) && 1 <= lower_level <= typemax(Int32) ||
            throw(ArgumentError("AtomDB level indices must be positive 32-bit integers"))
        return new(z, ion, upper_level, lower_level)
    end
end

struct TemperatureSlice
    continua::Vector{Tuple{Int,Continuum}}
    lines::Lines
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
    line_lists::Vector{Union{Nothing,Lines}}
    line_ids::Vector{Union{Nothing,Vector{LineID}}}
    bolometric::Vector{Union{Nothing,NTuple{30,Float64}}}
    direct_bolometric::Vector{Union{Nothing,NTuple{30,Float64}}}
    cache_lock::ReentrantLock

    function AtomDB(linefile::AbstractString, cocofile::AbstractString)
        linepath, cocopath = abspath(linefile), abspath(cocofile)
        line_temperatures = read_temperatures(linepath, "Nline")
        continuum_temperatures = read_temperatures(cocopath, "NCont")
        line_temperatures == continuum_temperatures ||
            throw(ArgumentError("AtomDB line and continuum files have different temperature grids"))
        slices = Vector{Union{Nothing,TemperatureSlice}}(nothing, length(line_temperatures))
        line_lists = Vector{Union{Nothing,Lines}}(nothing, length(line_temperatures))
        line_ids = Vector{Union{Nothing,Vector{LineID}}}(nothing, length(line_temperatures))
        bolometric = Vector{Union{Nothing,NTuple{30,Float64}}}(nothing, length(line_temperatures))
        return new(linepath, cocopath, line_temperatures, slices, line_lists, line_ids,
            bolometric, copy(bolometric), ReentrantLock())
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
        read_lines(fits, path, index)
    end
end

function read_lines(fits, path, index)
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
    return Lines(energies[order], Float32.(emissivities[order]), UInt8.(elements[order]))
end

# FITSIO returns a vector of row vectors for variable-length columns, and
# a numeric array for fixed-width columns (including scalar-width columns).
continuum_rows(column::AbstractVector{<:AbstractVector}, _nrows) = column
continuum_rows(column::AbstractArray{<:Real}, nrows) = eachcol(reshape(column, :, nrows))

function read_continua(fits, path, index)
    hdu = fits["EMISSIVITY", index]
    elements = read(hdu, "Z"; case_sensitive = false)
    ions = read(hdu, "rmJ"; case_sensitive = false)
    if !all(iszero, ions) || !all(z -> 1 <= z <= 30, elements) || !allunique(elements)
        throw(ArgumentError("Expected one CIE continuum row per element at temperature index $index in $path"))
    end
    continua = Tuple{Int,Continuum}[]
    isempty(elements) && return continua
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
            push!(continua, (Int(elements[i]), Continuum(e, f)))
        end
    end
    return continua
end

function temperature_lines(db::AtomDB, index::Int)
    return lock(db.cache_lock) do
        lines = db.line_lists[index]
        if isnothing(lines)
            lines = read_lines(db.linefile, index)
            db.line_lists[index] = lines
        end
        return lines
    end
end

function temperature_slice(db::AtomDB, index::Int)
    return lock(db.cache_lock) do
        slice = db.slices[index]
        if isnothing(slice)
            load_temperature_slices!(db, (index,))
            slice = db.slices[index]
        end
        return slice::TemperatureSlice
    end
end

function load_temperature_slices!(db, indices)
    return lock(db.cache_lock) do
        any(i -> isnothing(db.slices[i]), indices) || return nothing
        FITS(db.linefile, "r") do linefile
            FITS(db.cocofile, "r") do cocofile
                for i in indices
                    isnothing(db.slices[i]) || continue
                    lines = db.line_lists[i]
                    if isnothing(lines)
                        lines = read_lines(linefile, db.linefile, i)
                        db.line_lists[i] = lines
                    end
                    continua = read_continua(cocofile, db.cocofile, i)
                    db.slices[i] = TemperatureSlice(continua, lines)
                end
            end
        end
        return nothing
    end
end
