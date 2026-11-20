const COOLFLOW_VARIANTS = (
    (:coolflow, :Coolflow, :CoolflowModel, :tied, false),
    (:vcoolflow, :Vcoolflow, :VcoolflowModel, :partial, false),
    (:bcoolflow, :Bcoolflow, :BcoolflowModel, :tied, true),
    (:bvcoolflow, :Bvcoolflow, :BvcoolflowModel, :partial, true),
)

const CPH_VARIANTS = (
    (:cph, :Cph, :CphModel, :tied, false),
    (:vcph, :Vcph, :VcphModel, :partial, false),
    (:bcph, :Bcph, :BcphModel, :tied, true),
    (:bvcph, :Bvcph, :BvcphModel, :partial, true),
)

function local_cooling_variants(variants)
    return map(variants) do variant
        evaluator, constructor, model, mode, broadened = variant
        (Symbol("local", evaluator), Symbol("Local", constructor), Symbol("Local", model), mode, broadened)
    end
end

const LOCAL_COOLFLOW_VARIANTS = local_cooling_variants(COOLFLOW_VARIANTS)
const LOCAL_CPH_VARIANTS = local_cooling_variants(CPH_VARIANTS)

const COOLING_NTEMPS = 10
const CPH_LIMITS = (0.01, 50.0)
const MPC_CM = 3.0856775814913673e24
const PC_CM = MPC_CM / 1e6
const LOCAL_COOLING_NORM_UNIT = 1e-14 # solar masses per year per square parsec
const FLOW_PARTICLES_PER_SECOND = (5 / 2) * (1.98847e33 / 31_557_600) / (0.6 * 1.66053906660e-24)

# XSPEC 12.15's AtomDB grid from plasmaFileInfo, in keV.
const BOLOMETRIC_EDGES = let lo = 0.0105, hi = 1000.0
    edges = exp.(range(log(lo), log(hi); length = 1000))
    edges[1], edges[end] = lo, hi
    edges
end
const BOLOMETRIC_WIDTHS = diff(BOLOMETRIC_EDGES)
const BOLOMETRIC_CENTERS = (BOLOMETRIC_EDGES[1:end-1] .+ BOLOMETRIC_EDGES[2:end]) ./ 2

"""
    CoolingCosmology(; H0=70.0, q0=0.0, lambda0=0.73)

Distance parameters for the cooling-flow models. `H0` is in km/s/Mpc.
A nonzero `lambda0` uses a flat matter–Λ cosmology. With `lambda0=0`
the distance follows the matter-only Mattig relation with deceleration `q0`.
"""
struct CoolingCosmology
    H0::Float64
    q0::Float64
    lambda0::Float64

    function CoolingCosmology(; H0::Real = 70.0, q0::Real = 0.0, lambda0::Real = 0.73)
        isfinite(H0) && H0 > 0 || throw(DomainError(H0, "H0 must be positive and finite"))
        isfinite(q0) || throw(DomainError(q0, "q0 must be finite"))
        isfinite(lambda0) && lambda0 < 1 ||
            throw(DomainError(lambda0, "lambda0 must be finite and less than 1"))
        return new(H0, q0, lambda0)
    end
end

struct CoolingOptions
    cosmology::CoolingCosmology
    distance_mpc::Union{Nothing,Float64}
    normalization_redshift::Union{Nothing,Float64}
    ntemps::Int
    bolometric::Symbol

    function CoolingOptions(cosmology; distance_mpc, normalization_redshift, ntemps::Integer, bolometric::Symbol)
        validate_cooling_quadrature(ntemps, bolometric)
        if !isnothing(distance_mpc)
            isfinite(distance_mpc) && distance_mpc > 0 ||
                throw(DomainError(distance_mpc, "Luminosity distance must be positive and finite"))
        end
        if !isnothing(normalization_redshift)
            isfinite(normalization_redshift) && normalization_redshift > -1 ||
                throw(DomainError(normalization_redshift, "normalization_redshift must be finite and greater than -1"))
            if isnothing(distance_mpc) && normalization_redshift <= 0
                throw(DomainError(normalization_redshift, "A nonpositive normalization_redshift requires distance_mpc"))
            end
        end
        return new(cosmology, distance_mpc, normalization_redshift, ntemps, bolometric)
    end
end

struct LocalCoolingOptions
    ntemps::Int
    bolometric::Symbol

    function LocalCoolingOptions(ntemps::Integer, bolometric::Symbol)
        validate_cooling_quadrature(ntemps, bolometric)
        return new(ntemps, bolometric)
    end
end

function validate_cooling_quadrature(ntemps, bolometric)
    ntemps >= 2 || throw(ArgumentError("ntemps must be at least 2"))
    bolometric in (:xspec, :direct) || throw(ArgumentError("bolometric must be :xspec or :direct"))
    return nothing
end

cooling_norm_redshift(settings::CoolingOptions, redshift) = something(settings.normalization_redshift, redshift)
cooling_norm_redshift(::LocalCoolingOptions, _redshift) = 0.0

function cooling_luminosity_distance(z, cosmology)
    isfinite(z) && z > 0 ||
        throw(DomainError(z, "Cooling flows require positive redshift or an explicit distance_mpc"))
    (; H0, q0, lambda0) = cosmology
    hubble_distance = LIGHT_SPEED_KM_S / H0
    if iszero(lambda0)
        term = 1 + 2 * q0 * z
        term >= 0 || throw(DomainError(z, "The Mattig distance is undefined at this redshift"))
        r = sqrt(term) + 1
        return hubble_distance * (2 * z / r) * (1 + z / r)
    end

    # Simpson's rule for the flat cosmology in sqrt(scale factor), whose
    # integrand remains smooth up to the high redshifts allowed by the models.
    n = 128
    start = inv(sqrt(1 + z))
    step = z / ((1 + sqrt(1 + z)) * sqrt(1 + z) * n)
    integrand(a) = inv(sqrt(1 + lambda0 * (a^6 - 1)))
    integral = integrand(start) + integrand(1.0)
    for i in 1:(n - 1)
        integral += (isodd(i) ? 4 : 2) * integrand(start + i * step)
    end
    return (1 + z) * hubble_distance * (2 / 3) * step * integral
end

function cooling_normalization(z_norm, settings::CoolingOptions)
    isfinite(z_norm) && z_norm > -1 ||
        throw(DomainError(z_norm, "Normalization redshift must be finite and greater than -1"))
    distance = settings.distance_mpc
    if isnothing(distance)
        distance = cooling_luminosity_distance(z_norm, settings.cosmology)
    end
    distance_cm = distance * MPC_CM
    return 1e-14 * FLOW_PARTICLES_PER_SECOND / (4π * distance_cm^2)
end

# Convert the scaled local K or norm to a mass accretion rate per square parsec.
cooling_normalization(_z_norm, ::LocalCoolingOptions) =
    1e-14 * LOCAL_COOLING_NORM_UNIT * FLOW_PARTICLES_PER_SECOND / (4π * PC_CM^2)

function continuum_luminosity(continuum)
    energy, density = continuum.energies, continuum.density
    luminosity = 0.0
    for i in 2:length(energy)
        lo, hi = Float64(energy[i - 1]), Float64(energy[i])
        luminosity += (hi - lo) / 6 *
            ((2 * lo + hi) * density[i - 1] + (lo + 2 * hi) * density[i])
    end
    return luminosity
end

function direct_bolometric_elements(slice)
    totals = zeros(30)
    for (element, continuum) in slice.continua
        totals[element] += continuum_luminosity(continuum)
    end
    lines = slice.lines
    for i in eachindex(lines.energies)
        totals[lines.elements[i]] += lines.energies[i] * lines.emissivities[i]
    end
    return Tuple(totals)
end

function binned_bolometric_elements(slice, edges, centers, widths)
    totals = zeros(30)
    photons = zeros(length(centers))
    for (element, continuum) in slice.continua
        fill!(photons, 0.0)
        add_continuum!(photons, centers, widths, continuum, 1.0)
        totals[element] += sum(i -> centers[i] * photons[i], eachindex(photons))
    end

    lines = slice.lines
    firstline = searchsortedfirst(lines.energies, first(edges))
    lastline = searchsortedlast(lines.energies, last(edges))
    for i in firstline:lastline
        # As in add_lines!, a line on an internal edge goes in the lower bin.
        bin = max(1, searchsortedfirst(edges, lines.energies[i]) - 1)
        totals[lines.elements[i]] += centers[bin] * lines.emissivities[i]
    end
    return Tuple(totals)
end

# Store each element's bolometric emissivity in keV cm³/s. A fit can change
# the elemental scales without re-integrating the underlying table.
function bolometric_elements(db, index, method = :xspec)
    cache = method === :xspec ? db.bolometric : db.direct_bolometric
    return lock(db.cache_lock) do
        cached = cache[index]
        if isnothing(cached)
            slice = temperature_slice(db, index)
            cached = method === :xspec ?
                binned_bolometric_elements(slice, BOLOMETRIC_EDGES, BOLOMETRIC_CENTERS, BOLOMETRIC_WIDTHS) :
                direct_bolometric_elements(slice)
            cache[index] = cached
        end
        return cached
    end
end

function bolometric_emissivity(db, temperature, scales, method)
    lower, upper, fraction = temperature_bracket(db.temperatures, temperature)
    lo = sum(map(*, bolometric_elements(db, lower, method), scales))
    if lower == upper
        return lo
    end
    hi = sum(map(*, bolometric_elements(db, upper, method), scales))
    return (1 - fraction) * lo + fraction * hi
end

function cooling_temperatures(lowT, highT, ntemps)
    isfinite(lowT) && isfinite(highT) && 0 < lowT <= highT ||
        throw(DomainError((lowT, highT), "Cooling flow requires 0 < lowT ≤ highT in keV"))
    log_low = log(lowT)
    logstep = (log(highT) - log_low) / (ntemps - 1)
    temperatures = [exp(log_low + i * logstep) for i in 0:(ntemps - 1)]
    return temperatures, logstep
end

function cooling_weights(db, abundances, options, settings, z_norm, lowT, highT; peakT = nothing)
    if !all(value -> value isa Real && isfinite(value) && value >= 0, abundances)
        throw(DomainError(abundances, "Element abundances must be nonnegative and finite"))
    end
    norm = cooling_normalization(z_norm, settings)
    temperatures, logstep = cooling_temperatures(lowT, highT, settings.ntemps)
    if iszero(logstep)
        return zeros(typeof(lowT + highT + norm + first(abundances)), length(db.temperatures))
    end

    scales = map(*, promote(abundances...), options.ratios)
    if !isnothing(peakT)
        isfinite(peakT) && peakT > 0 || throw(DomainError(peakT, "peakT must be positive and finite"))
    end
    # XSPEC 12.15 evaluates its heating reference at 1 keV for peakT >= 1.
    # Below 1 keV its heating reference is empty.
    heat_reference = !isnothing(peakT) && peakT >= 1 ?
        bolometric_emissivity(db, 1.0, scales, settings.bolometric) : nothing
    emission_measures = map(eachindex(temperatures)) do i
        temperature = temperatures[i]
        emissivity = bolometric_emissivity(db, temperature, scales, settings.bolometric)
        cooling = emissivity
        if !isnothing(heat_reference)
            heating = heat_reference * (temperature / peakT)^2.5
            cooling = abs(emissivity - heating)
            # At exact balance, approximate the denominator over a finite
            # logarithmic cell instead of dividing by floating-point zero.
            if cooling <= eps(Float64) * max(emissivity, heating)
                cooling = emissivity * min(logstep, 1.0)
            end
        end
        iszero(cooling) && return zero(cooling * norm)
        endpoint = i == 1 || i == length(temperatures)
        return norm * temperature * (endpoint ? logstep / 2 : logstep) / cooling
    end
    return dem_node_weights(db, (temperatures, emission_measures))
end

function evaluate_cooling(edges, db, kind, options, settings, elements, lowT, highT;
        redshift, norm, velocity, peakT = nothing)
    values = element_parameters(elements, dem_abundance_names(kind))
    abundances = dem_abundances(kind, values, options)
    z_norm = cooling_norm_redshift(settings, redshift)
    weights = cooling_weights(db, abundances, options, settings, z_norm, lowT, highT; peakT)
    return evaluate_cie(edges, db, abundances, options, weights;
        redshift, normalization_redshift = z_norm, norm, velocity)
end

for (evaluator, _, _, mode, broadened) in COOLFLOW_VARIANTS
    @eval function $evaluator(edges::AbstractVector, db::AtomDB;
            lowT::Real = 0.1, highT::Real = 4.0, redshift::Real = 0.1, norm::Real = 1.0,
            velocity::Real = 0.0, switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            cosmology::CoolingCosmology = CoolingCosmology(), distance_mpc::Union{Nothing,Real} = nothing,
            normalization_redshift::Union{Nothing,Real} = nothing, bolometric::Symbol = :xspec,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        settings = CoolingOptions(cosmology; distance_mpc, normalization_redshift, ntemps, bolometric)
        return evaluate_cooling(edges, db, Val($(QuoteNode(mode))), options, settings, elements, lowT, highT;
            redshift, norm, velocity)
    end
end

for (evaluator, _, _, mode, broadened) in CPH_VARIANTS
    @eval function $evaluator(edges::AbstractVector, db::AtomDB;
            peakT::Real = 2.2, redshift::Real = 0.1, norm::Real = 1.0, velocity::Real = 0.0,
            switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            cosmology::CoolingCosmology = CoolingCosmology(), distance_mpc::Union{Nothing,Real} = nothing,
            normalization_redshift::Union{Nothing,Real} = nothing, bolometric::Symbol = :xspec,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        settings = CoolingOptions(cosmology; distance_mpc, normalization_redshift, ntemps, bolometric)
        return evaluate_cooling(edges, db, Val($(QuoteNode(mode))), options, settings, elements, CPH_LIMITS...;
            redshift, norm, velocity, peakT)
    end
end

for (evaluator, _, _, mode, broadened) in LOCAL_COOLFLOW_VARIANTS
    @eval function $evaluator(edges::AbstractVector, db::AtomDB;
            lowT::Real = 0.1, highT::Real = 4.0, redshift::Real = 0.0, norm::Real = 1.0,
            velocity::Real = 0.0, switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            bolometric::Symbol = :xspec, abundance_table = :angr, trace_abundance = 1.0,
            thermal_broadening::Bool = $broadened, elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        settings = LocalCoolingOptions(ntemps, bolometric)
        return evaluate_cooling(edges, db, Val($(QuoteNode(mode))), options, settings, elements, lowT, highT;
            redshift, norm, velocity)
    end
end

for (evaluator, _, _, mode, broadened) in LOCAL_CPH_VARIANTS
    @eval function $evaluator(edges::AbstractVector, db::AtomDB;
            peakT::Real = 2.2, redshift::Real = 0.0, norm::Real = 1.0, velocity::Real = 0.0,
            switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS, bolometric::Symbol = :xspec,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        settings = LocalCoolingOptions(ntemps, bolometric)
        return evaluate_cooling(edges, db, Val($(QuoteNode(mode))), options, settings, elements, CPH_LIMITS...;
            redshift, norm, velocity, peakT)
    end
end

@doc """
    localcoolflow(edges, db::AtomDB; redshift=0.0, norm=1.0, ...)

Evaluate cooling flow with `norm` equal to the mass accretion rate divided by
the squared source distance, in units of 10⁻¹⁴ solar masses per year per
square parsec.
The signed spectral redshift is independent of the normalization.
""" localcoolflow

const mkcflow = coolflow
const vmcflow = vcoolflow
