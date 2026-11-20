"""
    apec(edges, db; kT = 1.0, abundance = 1.0, redshift = 0.0, norm = 1.0)

Evaluate the unbroadened single-temperature XSPEC APEC spectrum in Julia.
`edges` is a strictly increasing vector of observed-frame bin edges in keV.
The output is integrated photon flux per bin in photons cm⁻² s⁻¹.
`kT` is in keV, `abundance` scales the 12 fitted metals relative to `angr`,
and `norm` is the XSPEC APEC normalization. H, He, and the trace elements
retain their `angr` abundances. Temperatures outside the table use the
nearest tabulated temperature.
"""
function apec(edges::AbstractVector, db::AtomDB;
        kT::Real = 1.0, abundance::Real = 1.0, redshift::Real = 0.0, norm::Real = 1.0)
    T = promote_type(Float64, eltype(edges), typeof(kT), typeof(abundance), typeof(redshift), typeof(norm))
    output = Vector{T}(undef, max(0, length(edges) - 1))
    return apec!(output, edges, db; kT, abundance, redshift, norm)
end

"""
    apec!(output, edges, db; kT = 1.0, abundance = 1.0, redshift = 0.0, norm = 1.0)

Evaluate [`apec`](@ref) into a preallocated output vector of length
`length(edges) - 1`, overwriting the previous spectrum.
"""
function apec!(output::AbstractVector, edges::AbstractVector, db::AtomDB;
        kT::Real = 1.0, abundance::Real = 1.0, redshift::Real = 0.0, norm::Real = 1.0)
    Base.require_one_based_indexing(output, edges)
    length(edges) >= 2 || throw(ArgumentError("At least two energy edges are required"))
    length(output) == length(edges) - 1 ||
        throw(DimensionMismatch("The output must have one entry per energy bin"))
    if !(0 < kT && isfinite(kT))
        throw(DomainError(kT, "Plasma temperature must be positive and finite"))
    end
    if !(0 <= abundance && isfinite(abundance))
        throw(DomainError(abundance, "Metal abundance must be nonnegative and finite"))
    end
    if !(-1 < redshift && isfinite(redshift))
        throw(DomainError(redshift, "Redshift must be finite and greater than -1"))
    end
    if !(0 <= norm && isfinite(norm))
        throw(DomainError(norm, "Normalization must be nonnegative and finite"))
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

    temperatures = db.temperatures
    lower = searchsortedlast(temperatures, kT)
    if lower == 0 || lower == length(temperatures) || kT == temperatures[lower]
        index = max(1, lower)
        add_slice!(output, source_edges, centers, widths, temperature_slice(db, index), abundance, 1.0)
    else
        upper = lower + 1
        fraction = (kT - temperatures[lower]) / (temperatures[upper] - temperatures[lower])
        add_slice!(output, source_edges, centers, widths, temperature_slice(db, lower), abundance, 1 - fraction)
        add_slice!(output, source_edges, centers, widths, temperature_slice(db, upper), abundance, fraction)
    end
    output .*= 1.0e14 * norm / z1
    return output
end

function add_slice!(output, edges, centers, widths, slice, abundance, weight)
    for continuum in slice.fixed_continua
        add_continuum!(output, centers, widths, continuum, weight)
    end
    for continuum in slice.metal_continua
        add_continuum!(output, centers, widths, continuum, weight * abundance)
    end
    add_lines!(output, edges, slice.fixed_lines, weight)
    add_lines!(output, edges, slice.metal_lines, weight * abundance)
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

function add_lines!(output, edges, lines, weight)
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
        output[i] += weight * emissivities[j]
    end
    return nothing
end
