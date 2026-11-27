struct Continuum{T<:AbstractFloat}
    energies::Vector{T}
    density::Vector{T}
end

struct Lines{E<:AbstractFloat,F<:AbstractFloat}
    energies::Vector{E}
    emissivities::Vector{F}
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

struct TemperatureSlice{C<:AbstractFloat,E<:AbstractFloat,F<:AbstractFloat}
    continua::Vector{Tuple{Int,Continuum{C}}}
    lines::Lines{E,F}
end

"""
    AtomDB(linefile, cocofile; constants=:codata2022)
    AtomDB(prefix; constants=:codata2022)
    AtomDB{C,E,F}(linefile, cocofile; constants=:codata2022)
    AtomDB{C,E,F}(prefix; constants=:codata2022)

Open an AtomDB CIE emissivity table pair on the Anders & Grevesse (`angr`)
abundance scale. A prefix is expanded to `prefix_line.fits` and
`prefix_coco.fits`. The constructor reads the temperature grid in keV.
Each temperature's emissivities are read on demand and cached in this object.

The default storage types are `Float32` for continuum energies and densities,
`Float64` for line energies, and `Float32` for line emissivities. Use the
typed constructors to select other concrete floating-point types for these
three quantities, respectively. Storage precision does not change the
temperature grid or the precision of model parameters and outputs.

CODATA 2022 is the default for line conversion, broadening, and cooling-flow
normalization. Pass `constants=:legacy` to reproduce the original values.
"""
struct AtomDB{C<:AbstractFloat,E<:AbstractFloat,F<:AbstractFloat}
    linefile::String
    cocofile::String
    temperatures::Vector{Float64}
    slices::Vector{Union{Nothing,TemperatureSlice{C,E,F}}}
    line_lists::Vector{Union{Nothing,Lines{E,F}}}
    line_ids::Vector{Union{Nothing,Vector{LineID}}}
    bolometric::Vector{Union{Nothing,NTuple{30,Float64}}}
    direct_bolometric::Vector{Union{Nothing,NTuple{30,Float64}}}
    constants::PlasmaConstants
    cache_lock::ReentrantLock

    function AtomDB{C,E,F}(linefile::AbstractString, cocofile::AbstractString;
            constants::Symbol = DEFAULT_CONSTANTS_PRESET) where {C,E,F}
        all(isconcretetype, (C, E, F)) ||
            throw(ArgumentError("AtomDB storage types must be concrete floating-point types"))
        physical = plasma_constants(constants)
        linepath, cocopath = abspath(linefile), abspath(cocofile)
        line_temperatures = read_temperatures(linepath, "Nline")
        continuum_temperatures = read_temperatures(cocopath, "NCont")
        line_temperatures == continuum_temperatures ||
            throw(ArgumentError("AtomDB line and continuum files have different temperature grids"))
        slices = Vector{Union{Nothing,TemperatureSlice{C,E,F}}}(nothing, length(line_temperatures))
        line_lists = Vector{Union{Nothing,Lines{E,F}}}(nothing, length(line_temperatures))
        line_ids = Vector{Union{Nothing,Vector{LineID}}}(nothing, length(line_temperatures))
        bolometric = Vector{Union{Nothing,NTuple{30,Float64}}}(nothing, length(line_temperatures))
        return new{C,E,F}(linepath, cocopath, line_temperatures, slices, line_lists, line_ids,
            bolometric, copy(bolometric), physical, ReentrantLock())
    end
end

AtomDB(linefile::AbstractString, cocofile::AbstractString; constants::Symbol = DEFAULT_CONSTANTS_PRESET) =
    AtomDB{Float32,Float64,Float32}(linefile, cocofile; constants)

AtomDB(prefix::AbstractString; constants::Symbol = DEFAULT_CONSTANTS_PRESET) =
    AtomDB{Float32,Float64,Float32}(prefix; constants)

AtomDB{C,E,F}(prefix::AbstractString; constants::Symbol = DEFAULT_CONSTANTS_PRESET) where {C,E,F} =
    AtomDB{C,E,F}(prefix * "_line.fits", prefix * "_coco.fits"; constants)

"""
    Plasma.load(AtomDB, prefix; constants=:codata2022)
    Plasma.load(AtomDB, linefile, cocofile; constants=:codata2022)

Open an AtomDB CIE database from a common prefix or two explicit FITS paths.
Pass `AtomDB{C,E,F}` to select storage precision or `constants=:legacy`
to reproduce the original constants, as described in [`AtomDB`](@ref).
Temperature slices are loaded on demand.
"""
load(db_type::Type{<:AtomDB}, prefix::AbstractString; constants::Symbol = DEFAULT_CONSTANTS_PRESET) =
    db_type(prefix; constants)

load(db_type::Type{<:AtomDB}, linefile::AbstractString, cocofile::AbstractString;
        constants::Symbol = DEFAULT_CONSTANTS_PRESET) =
    db_type(linefile, cocofile; constants)

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

function read_lines(path, index, line_type::Type{<:Lines}, constants::PlasmaConstants)
    return FITS(path, "r") do fits
        read_lines(fits, path, index, line_type, constants)
    end
end

function ordered_line_energies(wavelengths, constants::PlasmaConstants)
    energies = constants.kev_angstrom ./ wavelengths
    return energies, sortperm(energies)
end

function read_lines(fits, path, index, ::Type{Lines{E,F}}, constants::PlasmaConstants) where {E,F}
    hdu = fits["EMISSIVITY", index]
    wavelength = read(hdu, "Lambda"; case_sensitive = false)
    emissivities = read(hdu, "Epsilon"; case_sensitive = false)
    elements = read(hdu, "Element"; case_sensitive = false)
    if !all(x -> isfinite(x) && x > 0, wavelength) ||
            !all(f -> isfinite(f) && f >= 0, emissivities) ||
            !all(z -> 1 <= z <= 30, elements)
        throw(ArgumentError("Invalid CIE line record at temperature index $index in $path"))
    end

    # Sort before rounding so the line records and transition IDs stay aligned.
    energies, order = ordered_line_energies(wavelength, constants)
    return Lines{E,F}(E.(view(energies, order)), F.(emissivities[order]), UInt8.(elements[order]))
end

# FITSIO returns a vector of row vectors for variable-length columns, and
# a numeric array for fixed-width columns (including scalar-width columns).
continuum_rows(column::AbstractVector{<:AbstractVector}, _nrows) = column
continuum_rows(column::AbstractArray{<:Real}, nrows) = eachcol(reshape(column, :, nrows))

function read_continua(fits, path, index, ::Type{T}) where {T<:AbstractFloat}
    hdu = fits["EMISSIVITY", index]
    elements = read(hdu, "Z"; case_sensitive = false)
    ions = read(hdu, "rmJ"; case_sensitive = false)
    if !all(iszero, ions) || !all(z -> 1 <= z <= 30, elements) || !allunique(elements)
        throw(ArgumentError("Expected one CIE continuum row per element at temperature index $index in $path"))
    end
    continua = Tuple{Int,Continuum{T}}[]
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
            e = Vector{T}(view(energy_row, 1:n))
            f = Vector{T}(view(density_row, 1:n))
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

function temperature_lines(db::AtomDB{C,E,F}, index::Int) where {C,E,F}
    return lock(db.cache_lock) do
        lines = db.line_lists[index]
        if isnothing(lines)
            lines = read_lines(db.linefile, index, Lines{E,F}, db.constants)
            db.line_lists[index] = lines
        end
        return lines::Lines{E,F}
    end
end

function temperature_slice(db::AtomDB{C,E,F}, index::Int) where {C,E,F}
    return lock(db.cache_lock) do
        slice = db.slices[index]
        if isnothing(slice)
            load_temperature_slices!(db, (index,))
            slice = db.slices[index]
        end
        return slice::TemperatureSlice{C,E,F}
    end
end

function load_temperature_slices!(db::AtomDB{C,E,F}, indices) where {C,E,F}
    return lock(db.cache_lock) do
        any(i -> isnothing(db.slices[i]), indices) || return nothing
        FITS(db.linefile, "r") do linefile
            FITS(db.cocofile, "r") do cocofile
                for i in indices
                    isnothing(db.slices[i]) || continue
                    lines = db.line_lists[i]
                    if isnothing(lines)
                        lines = read_lines(linefile, db.linefile, i, Lines{E,F}, db.constants)
                        db.line_lists[i] = lines
                    end
                    continua = read_continua(cocofile, db.cocofile, i, C)
                    db.slices[i] = TemperatureSlice(continua, lines)
                end
            end
        end
        return nothing
    end
end
