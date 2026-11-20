"""
    apec(edges, db; kT = 1.0, abundance = 1.0, redshift = 0.0, norm = 1.0)

Evaluate the unbroadened single-temperature XSPEC APEC spectrum in Julia.
`edges` is a strictly increasing vector of observed-frame bin edges in keV.
The output is integrated photon flux per bin in photons cm⁻² s⁻¹.
`kT` is in keV, `abundance` scales the 12 fitted metals relative to `angr`,
and `norm` is the XSPEC APEC normalization. H, He, and the trace elements
retain their reference abundances by default. Choose a different XSPEC
table with `abundance_table`, or a trace-element multiplier or tie with
`trace_abundance`. Temperatures outside the table use the nearest tabulated
temperature.
"""
function apec(edges::AbstractVector, db::AtomDB; abundance::Real = 1.0, abundance_table = :angr,
        trace_abundance = 1.0, thermal_broadening::Bool = false, kwargs...)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    return evaluate_apec(edges, db, tied_abundances(abundance, options.trace), options; kwargs...)
end

"""
    apec!(output, edges, db; kwargs...)

Evaluate [`apec`](@ref) into a preallocated output vector of length
`length(edges) - 1`, overwriting the previous spectrum.
"""
function apec!(output::AbstractVector, edges::AbstractVector, db::AtomDB; abundance::Real = 1.0,
        abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = false, kwargs...)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    return evaluate_apec!(output, edges, db, tied_abundances(abundance, options.trace), options; kwargs...)
end

"""Evaluate APEC with thermal broadening and a Gaussian `velocity` sigma in km/s."""
function bapec(edges, db::AtomDB; thermal_broadening::Bool = true, kwargs...)
    return apec(edges, db; thermal_broadening, kwargs...)
end

"""Evaluate APEC with independent He and 12 metal abundances, e.g. `vapec(edges, db; Fe=0.4)`."""
function vapec(edges, db::AtomDB; kT::Real = 6.5, redshift::Real = 0.0, norm::Real = 1.0,
        velocity::Real = 0.0, abundance_table = :angr, trace_abundance = 1.0,
        thermal_broadening::Bool = false, elements...)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    abundances = partial_abundances(element_parameters(elements, VAPEC_NAMES), options.trace)
    return evaluate_apec(edges, db, abundances, options; kT, redshift, norm, velocity)
end

"""Evaluate APEC with independent abundances from H through Zn, all defaulting to 1."""
function vvapec(edges, db::AtomDB; kT::Real = 6.5, redshift::Real = 0.0, norm::Real = 1.0,
        velocity::Real = 0.0, abundance_table = :angr, thermal_broadening::Bool = false, elements...)
    options = apec_options(abundance_table, 1.0, thermal_broadening)
    abundances = element_parameters(elements, ELEMENT_NAMES)
    return evaluate_apec(edges, db, abundances, options; kT, redshift, norm, velocity)
end

"""Evaluate Vapec with thermal broadening and a Gaussian `velocity` sigma in km/s."""
function bvapec(edges, db::AtomDB; thermal_broadening::Bool = true, kwargs...)
    return vapec(edges, db; thermal_broadening, kwargs...)
end

"""Evaluate Vvapec with thermal broadening and a Gaussian `velocity` sigma in km/s."""
function bvvapec(edges, db::AtomDB; thermal_broadening::Bool = true, kwargs...)
    return vvapec(edges, db; thermal_broadening, kwargs...)
end

function evaluate_apec(edges::AbstractVector, db, abundances, options;
        kT::Real = 1.0, redshift::Real = 0.0, norm::Real = 1.0, velocity::Real = 0.0)
    return evaluate_cie(edges, db, abundances, options, kT; redshift, norm, velocity)
end

function evaluate_apec!(output::AbstractVector, edges::AbstractVector, db, abundances, options;
        kT::Real = 1.0, redshift::Real = 0.0, norm::Real = 1.0, velocity::Real = 0.0)
    return evaluate_cie!(output, edges, db, abundances, options, kT; redshift, norm, velocity)
end

emission_type(temperature::Real) = typeof(temperature)
emission_type(node_weights::AbstractVector) = eltype(node_weights)

function evaluate_cie(edges, db, abundances, options, emission; redshift, norm, velocity)
    T = promote_type(Float64, eltype(edges), emission_type(emission), typeof(redshift), typeof(norm),
        typeof(velocity), map(typeof, abundances)...)
    output = Vector{T}(undef, max(0, length(edges) - 1))
    return evaluate_cie!(output, edges, db, abundances, options, emission; redshift, norm, velocity)
end

function evaluate_cie!(output::AbstractVector, edges::AbstractVector, db, abundances, options, emission;
        redshift, norm = 1.0, velocity)
    Base.require_one_based_indexing(output, edges)
    length(edges) >= 2 || throw(ArgumentError("At least two energy edges are required"))
    length(output) == length(edges) - 1 ||
        throw(DimensionMismatch("The output must have one entry per energy bin"))
    if !all(value -> value isa Real && 0 <= value && isfinite(value), abundances)
        throw(DomainError(abundances, "Element abundances must be nonnegative and finite"))
    end
    if !(-1 < redshift && isfinite(redshift))
        throw(DomainError(redshift, "Redshift must be finite and greater than -1"))
    end
    if !(0 <= norm && isfinite(norm))
        throw(DomainError(norm, "Normalization must be nonnegative and finite"))
    end
    if !(0 <= velocity && isfinite(velocity))
        throw(DomainError(velocity, "Velocity broadening must be nonnegative and finite"))
    end

    z1 = 1.0 + redshift
    source_edges = edges .* z1
    if !all(x -> isfinite(x) && x >= 0, source_edges) ||
            !all(i -> source_edges[i] > source_edges[i - 1], 2:length(source_edges))
        throw(ArgumentError("Energy edges must be finite, nonnegative, and strictly increasing in the source frame"))
    end
    widths = diff(source_edges)
    centers = [source_edges[i] + widths[i] / 2 for i in eachindex(output)]
    fill!(output, zero(eltype(output)))
    scales = map(*, promote(abundances...), options.ratios)

    add_emission!(output, source_edges, centers, widths, db, emission, scales, options.thermal, velocity)
    output .*= 1.0e14 * norm / z1
    return output
end

function temperature_bracket(temperatures, kT)
    lower = searchsortedlast(temperatures, kT)
    if lower == 0 || lower == length(temperatures) || kT == temperatures[lower]
        index = max(1, lower)
        return index, index, zero(kT)
    end
    upper = lower + 1
    fraction = (kT - temperatures[lower]) / (temperatures[upper] - temperatures[lower])
    return lower, upper, fraction
end

function add_emission!(output, edges, centers, widths, db, kT::Real, scales, thermal, velocity)
    if !(0 < kT && isfinite(kT))
        throw(DomainError(kT, "Plasma temperature must be positive and finite"))
    end
    lower, upper, fraction = temperature_bracket(db.temperatures, kT)
    add_slice!(output, edges, centers, widths, db, lower, scales, 1 - fraction, thermal, velocity)
    if lower != upper
        add_slice!(output, edges, centers, widths, db, upper, scales, fraction, thermal, velocity)
    end
    return nothing
end

function add_emission!(output, edges, centers, widths, db, weights::AbstractVector, scales, thermal, velocity)
    for index in eachindex(weights)
        weight = weights[index]
        iszero(weight) && continue
        add_slice!(output, edges, centers, widths, db, index, scales, weight, thermal, velocity)
    end
    return nothing
end

function add_slice!(output, edges, centers, widths, db, index, scales, weight, thermal, velocity)
    slice = temperature_slice(db, index)
    for (z, continuum) in slice.continua
        add_continuum!(output, centers, widths, continuum, weight * scales[z])
    end
    if thermal || !iszero(velocity)
        fractions = line_width_fractions(db.temperatures[index], thermal, velocity)
        add_broadened_lines!(output, edges, slice.lines, scales, weight, fractions)
    else
        add_lines!(output, edges, slice.lines, scales, weight)
    end
    return nothing
end

function add_continuum!(output, centers, widths, continuum, weight)
    energies, density = continuum.energies, continuum.density
    firstbin = searchsortedfirst(centers, first(energies))
    lastbin = searchsortedlast(centers, last(energies))
    firstbin > lastbin && return nothing
    if length(energies) == 1
        # A single tabulation has no interpolation interval.
        for i in firstbin:lastbin
            output[i] += weight * density[1] * widths[i]
        end
        return nothing
    end
    j = clamp(searchsortedfirst(energies, centers[firstbin]), 2, length(energies))
    @inbounds for i in firstbin:lastbin
        energy = centers[i]
        while j < length(energies) && energy > energies[j]
            j += 1
        end
        lo, hi = Float64(energies[j - 1]), Float64(energies[j])
        fraction = (energy - lo) / (hi - lo)
        flux_density = (1 - fraction) * density[j - 1] + fraction * density[j]
        output[i] += weight * flux_density * widths[i]
    end
    return nothing
end

function add_lines!(output, edges, lines, scales, weight)
    energies, emissivities = lines.energies, lines.emissivities
    firstline = searchsortedfirst(energies, first(edges))
    lastline = searchsortedlast(energies, last(edges))
    firstline > lastline && return nothing
    i = max(1, searchsortedfirst(edges, energies[firstline]) - 1)
    @inbounds for j in firstline:lastline
        # XSPEC's line histogram assigns a line on an internal edge to the lower bin.
        while i < length(output) && energies[j] > edges[i + 1]
            i += 1
        end
        output[i] += weight * scales[lines.elements[j]] * emissivities[j]
    end
    return nothing
end
