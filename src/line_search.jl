"""
    LineMatch

A CIE line search result. `id` identifies the transition, `wavelength` is in
Å, `energy` is in keV, and `emissivity` is in photons cm³ s⁻¹.
"""
struct LineMatch
    id::LineID
    wavelength::Float64
    energy::Float64
    emissivity::Float64
end

function line_atomic_number(element)
    if element isa Integer && 1 <= element <= length(ELEMENT_NAMES)
        return Int(element)
    end
    if element isa Symbol || element isa AbstractString
        symbol = lowercase(String(element))
        z = findfirst(name -> lowercase(String(name)) == symbol, ELEMENT_NAMES)
        !isnothing(z) && return z
    end
    throw(ArgumentError("Unknown AtomDB element $element"))
end

function line_ion_stage(value)
    value isa Integer && 1 <= value <= 31 || throw(ArgumentError("Ion must be a spectroscopic stage from 1 to 31"))
    return Int(value)
end

function line_selection(value, parse)
    isnothing(value) && return nothing
    selection = value isa Union{Tuple,AbstractVector} ? value : (value,)
    return Set(map(parse, selection))
end

function temperature_line_ids(db::AtomDB, index::Int, lines::Lines)
    return lock(db.cache_lock) do
        ids = db.line_ids[index]
        if isnothing(ids)
            ids = FITS(db.linefile, "r") do fits
                hdu = fits["EMISSIVITY", index]
                wavelength = read(hdu, "Lambda"; case_sensitive = false)
                ions = read(hdu, "Ion"; case_sensitive = false)
                upper = read(hdu, "UpperLev"; case_sensitive = false)
                lower = read(hdu, "LowerLev"; case_sensitive = false)
                _, order = ordered_line_energies(wavelength)
                return [LineID(lines.elements[i], ions[j], upper[j], lower[j]) for (i, j) in enumerate(order)]
            end
            db.line_ids[index] = ids
        end
        return ids
    end
end

function line_search_unit(unit)
    name = lowercase(String(unit))
    name in ("angstrom", "å", "a") && return (KEV_ANGSTROM, true)
    name == "nm" && return (KEV_ANGSTROM / 10, true)
    name == "kev" && return (1.0, false)
    name == "ev" && return (1e3, false)
    throw(ArgumentError("Line search unit must be Angstrom, nm, keV, or eV"))
end

function line_search_temperature(kT::Real, teunit)
    name = lowercase(String(teunit))
    name in ("kev", "k", "kelvin") || throw(ArgumentError("Temperature unit must be keV or K"))
    isfinite(kT) && kT > 0 || throw(DomainError(kT, "Temperature must be positive and finite"))
    # Boltzmann's constant in keV/K.
    return Float64(kT) * (name == "kev" ? 1.0 : 1.380649e-23 / 1.602176634e-16)
end

function line_abundance_scales(abundance_table, abundances)
    scales = collect(apec_options(something(abundance_table, :angr), 1.0, false).ratios)
    if !isnothing(abundances)
        for (element, value) in pairs(abundances)
            value isa Real && isfinite(value) && value >= 0 ||
                throw(DomainError(value, "Abundance multipliers must be nonnegative and finite"))
            scales[line_atomic_number(element)] *= value
        end
    end
    return scales
end

line_model_state(_db, _model) = throw(ArgumentError("Line ranking requires a native Plasma.jl SpectralFitting model"))

"""
    find_lines(db, lo, hi; kT, unit=:angstrom, teunit=:keV, element=nothing,
               ion=nothing, brightest=nothing, min_emissivity=0,
               abundance_table=:angr, abundances=nothing)
    find_lines(model, lo, hi; unit=:angstrom, element=nothing, ion=nothing,
               brightest=nothing, min_emissivity=0)

Search the loaded AtomDB CIE emissivity tables over an inclusive energy or
wavelength interval. `unit` accepts Angstrom, nm, keV, or eV. `element` accepts
a symbol, name such as `"Fe"`, atomic number, or a collection of these.
`ion` accepts a spectroscopic stage or collection of stages.

At a single `kT`, use Plasma.jl's temperature interpolation. A native
SpectralFitting component can be passed as the first argument or as the
`model` keyword with its database. Model searches report intrinsic rest-frame
emissivities using the component's abundances and relative temperature weights.

Return a vector of [`LineMatch`](@ref) values, ordered by wavelength unless
`brightest` requests the top lines by emissivity.
"""
function find_lines(db::AtomDB, lo::Real, hi::Real; kT = nothing, model = nothing,
        unit = :angstrom, teunit = :keV, element = nothing, ion = nothing,
        brightest = nothing, min_emissivity::Real = 0.0, abundance_table = nothing, abundances = nothing)
    isfinite(lo) && isfinite(hi) && 0 < lo < hi ||
        throw(ArgumentError("Line search bounds must be finite and satisfy 0 < lo < hi"))
    isfinite(min_emissivity) && min_emissivity >= 0 ||
        throw(ArgumentError("min_emissivity must be nonnegative and finite"))
    isnothing(brightest) || (brightest isa Integer && brightest > 0) ||
        throw(ArgumentError("brightest must be a positive integer"))

    scale, reciprocal = line_search_unit(unit)
    lower, upper = reciprocal ? (scale / hi, scale / lo) : (lo / scale, hi / scale)
    elements = line_selection(element, line_atomic_number)
    ions = line_selection(ion, line_ion_stage)
    weights, scales = if isnothing(model)
        isnothing(kT) && throw(ArgumentError("Specify kT or a native Plasma.jl model to search the CIE tables"))
        temperature = line_search_temperature(kT, teunit)
        dem_node_weights(db, ((temperature,), (1.0,))), line_abundance_scales(abundance_table, abundances)
    else
        all(isnothing, (kT, abundance_table, abundances)) ||
            throw(ArgumentError("Model searches use the model's temperatures and abundances"))
        line_model_state(db, model)
    end

    matches = Dict{LineID,LineMatch}()
    for index in eachindex(weights)
        weight = weights[index]
        iszero(weight) && continue
        lines = temperature_lines(db, index)
        # Recheck bounds in the requested units after converting the search window.
        start = searchsortedfirst(lines.energies, prevfloat(Float64(lower)))
        stop = searchsortedlast(lines.energies, nextfloat(Float64(upper)))
        start > stop && continue
        ids = temperature_line_ids(db, index, lines)
        for j in start:stop
            id = ids[j]
            !isnothing(elements) && !(id.element in elements) && continue
            !isnothing(ions) && !(id.ion in ions) && continue
            energy = lines.energies[j]
            coordinate = reciprocal ? scale / energy : energy * scale
            lo <= coordinate <= hi || continue
            emissivity = Float64(weight * scales[id.element] * lines.emissivities[j])
            prior = get(matches, id, nothing)
            matches[id] = isnothing(prior) ? LineMatch(id, KEV_ANGSTROM / energy, energy, emissivity) :
                LineMatch(id, prior.wavelength, prior.energy, prior.emissivity + emissivity)
        end
    end

    rows = [row for row in values(matches) if row.emissivity >= min_emissivity]
    sort!(rows; by = row -> (row.wavelength, row.id.element, row.id.ion, row.id.upper_level, row.id.lower_level))
    if !isnothing(brightest)
        sort!(rows; by = row -> -row.emissivity, alg = Base.Sort.MergeSort)
        resize!(rows, min(length(rows), brightest))
    end
    return rows
end

function line_emissivity_at_node(db, id, index)
    lines = temperature_lines(db, index)
    ids = temperature_line_ids(db, index, lines)
    value = 0.0
    for i in eachindex(ids)
        ids[i] == id && (value += lines.emissivities[i])
    end
    return value
end

"""
    line_emissivity(db, line, kT; teunit=:keV, abundance_table=:angr, abundances=nothing)

Return the CIE emissivity of a [`LineID`](@ref) or [`LineMatch`](@ref) in
photons cm³ s⁻¹. A match supplies the ID for a new evaluation at the requested
scalar or vector `kT`. Missing transitions contribute zero.
"""
function line_emissivity(db::AtomDB, id::LineID, kT::Real; teunit = :keV, abundance_table = :angr, abundances = nothing)
    temperature = line_search_temperature(kT, teunit)
    scales = line_abundance_scales(abundance_table, abundances)
    lower, upper, fraction = temperature_bracket(db.temperatures, temperature)
    value = line_emissivity_at_node(db, id, lower)
    if lower != upper
        value = (1 - fraction) * value + fraction * line_emissivity_at_node(db, id, upper)
    end
    return scales[id.element] * value
end

line_emissivity(db::AtomDB, id::LineID, temperatures::AbstractVector; options...) =
    map(kT -> line_emissivity(db, id, kT; options...), temperatures)

line_emissivity(db::AtomDB, match::LineMatch, kT::Union{Real,AbstractVector}; options...) =
    line_emissivity(db, match.id, kT; options...)

"""
    ratio_emissivity(db, numerator, denominator, kT; kwargs...)

Return the emissivity ratio for two transitions at a scalar or vector `kT`.
Each transition can be a [`LineID`](@ref) or [`LineMatch`](@ref). A zero
denominator yields `NaN`.
"""
function ratio_emissivity(db::AtomDB, numerator::Union{LineID,LineMatch},
        denominator::Union{LineID,LineMatch}, kT::Real; options...)
    value = line_emissivity(db, denominator, kT; options...)
    return iszero(value) ? NaN : line_emissivity(db, numerator, kT; options...) / value
end

function ratio_emissivity(db::AtomDB, numerator::Union{LineID,LineMatch},
        denominator::Union{LineID,LineMatch}, temperatures::AbstractVector; options...)
    return map(kT -> ratio_emissivity(db, numerator, denominator, kT; options...), temperatures)
end
