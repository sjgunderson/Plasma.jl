# Keep one working grid per fitting component. Trial temperatures are used
# only to select and weight the tabulated spectra on that grid.
const SPECTRUM_CACHE_BYTES = 64 * 1024^2

mutable struct SpectrumBasis
    variable_groups::UInt32
    columns::NTuple{30,Int}
    scales::NTuple{30,Float64}
    representatives::Vector{Int}
    spectra::Vector{Union{Nothing,Matrix{Float64}}}
    capacity::Int
    count::Int
end

mutable struct SpectrumGrid
    observed_edges::Vector{Float64}
    source_edges::Vector{Float64}
    centers::Vector{Float64}
    widths::Vector{Float64}
    redshift::Float64
    velocity::Float64
    options::ApecOptions
    reference_abundances::NTuple{30,Float64}
    basis::Union{Nothing,SpectrumBasis}
end

mutable struct SpectrumCache
    # Zero identifies abundances that are constant by construction.
    abundance_groups::NTuple{30,Int}
    lock::ReentrantLock
    grid::Union{Nothing,SpectrumGrid}
end

function spectrum_cache(enabled, kind, options, indices)
    enabled || return nothing
    groups = spectrum_abundance_groups(kind, options.trace, indices)
    return SpectrumCache(groups, ReentrantLock(), nothing)
end

spectrum_abundance_groups(::Val{:full}, _trace, _indices) = ntuple(identity, 30)
spectrum_abundance_groups(::Val{:tied}, trace, indices) =
    map(abundance_group, tied_abundances(-1, trace, indices))

function spectrum_abundance_groups(::Val{:partial}, trace, indices)
    markers = ntuple(-, maximum(indices))
    return map(abundance_group, partial_abundances(markers, trace, indices))
end

abundance_group(marker) = marker < 0 ? Int(-marker) : 0

# The fitting extension also recognizes dual numbers with constant geometry.
cache_value(x::Float64) = x
cache_value(_x) = nothing
cache_constant(x) = cache_value(x)

function add_cached_emission!(_output, _edges, _db, _abundances, _options, _emission, _z, _v, ::Nothing)
    return false
end

function add_cached_emission!(output, edges, db, abundances, options, emission, redshift, velocity, cache::SpectrumCache)
    eltype(edges) === Float64 || return false
    cache_value(zero(eltype(output))) === 0.0 || return false
    z = cache_constant(redshift)
    v = cache_constant(velocity)
    (isnothing(z) || isnothing(v)) && return false
    values = map(cache_value, abundances)
    any(isnothing, values) && return false

    return lock(cache.lock) do
        grid = cache.grid
        if isnothing(grid) || z != grid.redshift || v != grid.velocity || options !== grid.options ||
                edges != grid.observed_edges
            # Build templates only after two evaluations share a geometry.
            # A fit varying redshift or velocity can often use direct synthesis.
            observed = collect(edges)
            source, centers, widths = source_grid(observed, z)
            cache.grid = SpectrumGrid(observed, source, centers, widths, z, v, options, values, nothing)
            fill!(output, zero(eltype(output)))
            add_emission!(output, source, centers, widths, db, emission,
                abundance_scales(abundances, options), options.thermal, v)
            return true
        end

        basis = grid.basis
        variable = isnothing(basis) ? UInt32(0) : basis.variable_groups
        for (element, group) in enumerate(cache.abundance_groups)
            if group != 0 && cache_constant(abundances[element]) !== grid.reference_abundances[element]
                variable |= UInt32(1) << (group - 1)
            end
        end
        if isnothing(basis) || variable != basis.variable_groups
            basis = spectrum_basis(grid, cache.abundance_groups, variable, length(db.temperatures))
            grid.basis = basis
        end
        fill!(output, zero(eltype(output)))
        add_binned_emission!(output, grid, basis, db, emission, abundances)
        return true
    end
end

function spectrum_basis(grid, groups, variable, ntemperatures)
    representatives = Int[]
    group_columns = zeros(Int, 30)
    columns = ntuple(30) do z
        group = groups[z]
        if group == 0 || iszero(variable & (UInt32(1) << (group - 1)))
            return 1
        end
        if group_columns[group] == 0
            push!(representatives, z)
            group_columns[group] = length(representatives) + 1
        end
        return group_columns[group]
    end
    scales = ntuple(30) do z
        return grid.options.ratios[z] * (columns[z] == 1 ? grid.reference_abundances[z] : 1.0)
    end
    bytes_per_node = length(grid.widths) * (1 + length(representatives)) * sizeof(Float64)
    capacity = min(ntemperatures, div(SPECTRUM_CACHE_BYTES, bytes_per_node))
    slots = Vector{Union{Nothing,Matrix{Float64}}}(nothing, ntemperatures)
    return SpectrumBasis(variable, columns, scales, representatives, slots, capacity, 0)
end

function add_binned_emission!(output, grid, basis, db, kT::Real, abundances)
    lower, upper, fraction = emission_bracket(db, kT)
    if isnothing(basis.spectra[lower]) || isnothing(basis.spectra[upper])
        load_temperature_slices!(db, (lower, upper))
    end
    add_binned_slice!(output, grid, basis, db, lower, 1 - fraction, abundances)
    if lower != upper
        add_binned_slice!(output, grid, basis, db, upper, fraction, abundances)
    end
    return nothing
end

function add_binned_emission!(output, grid, basis, db, weights::AbstractVector, abundances)
    missing = (i for i in eachindex(weights) if !iszero(weights[i]) && isnothing(basis.spectra[i]))
    isempty(missing) || load_temperature_slices!(db, missing)
    for i in eachindex(weights)
        weight = weights[i]
        iszero(weight) && continue
        add_binned_slice!(output, grid, basis, db, i, weight, abundances)
    end
    return nothing
end

function add_binned_slice!(output, grid, basis, db, index, weight, abundances)
    spectrum = basis.spectra[index]
    if isnothing(spectrum)
        if basis.count == basis.capacity
            add_slice!(output, grid.source_edges, grid.centers, grid.widths, db, index,
                abundance_scales(abundances, grid.options), weight, grid.options.thermal, grid.velocity)
            return nothing
        end
        spectrum = zeros(length(grid.widths), 1 + length(basis.representatives))
        add_slice!(spectrum, grid.source_edges, grid.centers, grid.widths, db, index,
            basis.scales, 1.0, grid.options.thermal, grid.velocity, basis.columns)
        basis.spectra[index] = spectrum
        basis.count += 1
    end
    accumulate_spectrum!(output, spectrum, 1, weight)
    for (i, z) in enumerate(basis.representatives)
        accumulate_spectrum!(output, spectrum, i + 1, weight * abundances[z])
    end
    return nothing
end

function accumulate_spectrum!(output, spectrum, column, coefficient)
    @inbounds @simd for i in eachindex(output)
        output[i] += coefficient * spectrum[i, column]
    end
    return nothing
end
