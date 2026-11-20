const WDEM_VARIANTS = (
    (:wdem, :Wdem, :WdemModel, :tied, false),
    (:vwdem, :Vwdem, :VwdemModel, :partial, false),
    (:vvwdem, :Vvwdem, :VvwdemModel, :full, false),
    (:bwdem, :Bwdem, :BwdemModel, :tied, true),
    (:bvwdem, :Bvwdem, :BvwdemModel, :partial, true),
    (:bvvwdem, :Bvvwdem, :BvvwdemModel, :full, true),
)

const CEMPOW_VARIANTS = (
    (:cempow, :Cempow, :CempowModel, :tied, false),
    (:vcempow, :Vcempow, :VcempowModel, :partial, false),
    (:bcempow, :Bcempow, :BcempowModel, :tied, true),
    (:bvcempow, :Bvcempow, :BvcempowModel, :partial, true),
)

const WDEM_INTERVALS = 50
const CEMPOW_MINIMUM = 0.02725 # keV, corresponding to 10^5.5 K
const CEMPOW_LOG_STEP = 0.1

function wdem_sampling(Tmax, beta, inv_slope)
    isfinite(Tmax) && Tmax > 0 || throw(DomainError(Tmax, "Tmax must be positive and finite"))
    isfinite(beta) && 0 < beta <= 1 || throw(DomainError(beta, "beta must be in (0, 1]"))
    isfinite(inv_slope) || throw(DomainError(inv_slope, "inv_slope must be finite"))
    beta == 1 && return ((Tmax,), (one(Tmax),))

    # XSPEC steps by (Tmax - Tmin)/50 and keeps 51 intervals, including the
    # one immediately above Tmax. Integrate T^(1/inv_slope) over each bin.
    step = (1 - beta) / WDEM_INTERVALS
    samples = WDEM_INTERVALS + 1
    exponent = 1.0 / inv_slope + 1
    if !isfinite(exponent)
        # A vertical power law collapses to the cold or hot endpoint.
        i = inv_slope < 0 ? 1 : samples
        return ((Tmax * (beta + (i - 0.5) * step),), (one(Tmax),))
    end
    upper = beta + samples * step
    weights = map(1:samples) do i
        lo = beta + (i - 1) * step
        logwidth = log1p(step / lo)
        if abs(exponent) < 1e-5
            # Smooth logarithmic limit near a T^-1 distribution.
            x = exponent * logwidth
            logwidth * exp(exponent * log(lo)) * (1 + x * (0.5 + x * (1 / 6 + x / 24)))
        elseif exponent > 0
            hi = beta + i * step
            exp(exponent * log(hi / upper)) * -expm1(-exponent * logwidth)
        else
            exp(exponent * (log(lo) - log(beta))) * -expm1(exponent * logwidth)
        end
    end
    temperatures = [Tmax * (beta + (i - 0.5) * step) for i in 1:samples]
    return temperatures, weights ./ sum(weights)
end

function cempow_sampling(alpha, Tmax)
    isfinite(alpha) && alpha > 0 || throw(DomainError(alpha, "alpha must be positive and finite"))
    isfinite(Tmax) && Tmax >= CEMPOW_MINIMUM ||
        throw(DomainError(Tmax, "Cempow Tmax must be at least $CEMPOW_MINIMUM keV"))
    last = 0
    while Tmax * 10.0^(-CEMPOW_LOG_STEP * (last + 1)) >= CEMPOW_MINIMUM
        last += 1
    end
    ratios = [10.0^(-CEMPOW_LOG_STEP * i) for i in 0:last]
    # Preserve the absolute XSPEC weights. Their sum changes with alpha.
    return Tmax .* ratios, CEMPOW_LOG_STEP .* (ratios .^ alpha)
end

function dem_options(table, trace, thermal, nH, switch)
    switch == 2 || throw(ArgumentError("Plasma.jl implements the AtomDB CIE backend (switch=2)"))
    isfinite(nH) && nH > 0 || throw(DomainError(nH, "nH must be positive and finite"))
    return apec_options(table, trace, thermal, DEM_INDEX)
end

dem_abundance_names(::Val{:tied}) = (:abundance,)
dem_abundance_names(::Val{:partial}) = DEM_NAMES
dem_abundance_names(::Val{:full}) = ELEMENT_NAMES

dem_abundances(::Val{:tied}, values, options) = tied_abundances(only(values), options.trace, DEM_INDEX)
dem_abundances(::Val{:partial}, values, options) = partial_abundances(values, options.trace, DEM_INDEX)
dem_abundances(::Val{:full}, values, _options) = values

# Collapse the DEM sampling onto the AtomDB temperature nodes. Temperature
# interpolation, binning, and broadening are all linear in these weights.
function dem_node_weights(db::AtomDB, sampling)
    sample_temperatures, emission_measures = sampling
    T = promote_type(Float64, eltype(sample_temperatures), eltype(emission_measures))
    temperatures = db.temperatures
    weights = zeros(T, length(temperatures))
    for (kT, measure) in zip(sample_temperatures, emission_measures)
        lower, upper, fraction = temperature_bracket(temperatures, kT)
        weights[lower] += measure * (1 - fraction)
        if lower != upper
            weights[upper] += measure * fraction
        end
    end
    return weights
end

for (evaluator, _, _, mode, broadened) in WDEM_VARIANTS
    @eval function $evaluator(edges::AbstractVector, db::AtomDB;
            Tmax::Real = 1.0, beta::Real = 0.1, inv_slope::Real = 0.25, redshift::Real = 0.0,
            norm::Real = 1.0, velocity::Real = 0.0, nH::Real = 1.0, switch::Integer = 2,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, nH, switch)
        kind = Val($(QuoteNode(mode)))
        values = element_parameters(elements, dem_abundance_names(kind))
        abundances = dem_abundances(kind, values, options)
        weights = dem_node_weights(db, wdem_sampling(Tmax, beta, inv_slope))
        return evaluate_cie(edges, db, abundances, options, weights; redshift, norm, velocity)
    end
end

for (evaluator, _, _, mode, broadened) in CEMPOW_VARIANTS
    @eval function $evaluator(edges::AbstractVector, db::AtomDB;
            alpha::Real = 1.0, Tmax::Real = 1.0, redshift::Real = 0.0,
            norm::Real = 1.0, velocity::Real = 0.0, nH::Real = 1.0, switch::Integer = 2,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, nH, switch)
        kind = Val($(QuoteNode(mode)))
        values = element_parameters(elements, dem_abundance_names(kind))
        abundances = dem_abundances(kind, values, options)
        weights = dem_node_weights(db, cempow_sampling(alpha, Tmax))
        return evaluate_cie(edges, db, abundances, options, weights; redshift, norm, velocity)
    end
end
