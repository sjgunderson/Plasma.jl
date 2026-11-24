using Plasma: COOLFLOW_VARIANTS, CPH_VARIANTS, LOCAL_COOLFLOW_VARIANTS, LOCAL_CPH_VARIANTS,
    COOLING_NTEMPS, CPH_LIMITS, CoolingCosmology, CoolingOptions, LocalCoolingOptions,
    cooling_norm_redshift, cooling_weights
import Plasma: Coolflow, Vcoolflow, Bcoolflow, Bvcoolflow, Cph, Vcph, Bcph, Bvcph
import Plasma: LocalCoolflow, LocalVcoolflow, LocalBcoolflow, LocalBvcoolflow,
    LocalCph, LocalVcph, LocalBcph, LocalBvcph

abstract type AbstractCoolingModel{T} <: AbstractCIEModel{T} end
abstract type AbstractCoolflowModel{T} <: AbstractCoolingModel{T} end
abstract type AbstractCphModel{T} <: AbstractCoolingModel{T} end

for (variants, supertype, temperature_parameters, cooling_options) in (
        (COOLFLOW_VARIANTS, :AbstractCoolflowModel, (:lowT, :highT), CoolingOptions),
        (LOCAL_COOLFLOW_VARIANTS, :AbstractCoolflowModel, (:lowT, :highT), LocalCoolingOptions),
        (CPH_VARIANTS, :AbstractCphModel, (:peakT,), CoolingOptions),
        (LOCAL_CPH_VARIANTS, :AbstractCphModel, (:peakT,), LocalCoolingOptions))
    for (_, _, model, mode, broadened) in variants
        kind = Val(mode)
        abundance_parameters = dem_parameter_names(kind)
        names = (:K, temperature_parameters..., abundance_parameters..., :redshift)
        broadened && (names = (names..., :velocity))
        supertype === :AbstractCoolflowModel && (names = (names..., :slope))
        fields = [:($parameter::T) for parameter in names]
        @eval struct $model{T,D<:AtomDB} <: $supertype{T}
            db::D
            options::ApecOptions
            cache::Union{Nothing,SpectrumCache}
            cooling::$cooling_options
            $(fields...)
        end
        @eval function model_abundances(model::$model)
            values = map(name -> getfield(model, name), $abundance_parameters)
            return dem_abundances($kind, values, model.options)
        end
        if broadened
            @eval model_velocity(model::$model) = model.velocity
        end
    end
end

cooling_redshift_lower_limit(settings::CoolingOptions) =
    isnothing(settings.distance_mpc) && isnothing(settings.normalization_redshift) ? 0.0 : -0.999
cooling_redshift_lower_limit(::LocalCoolingOptions) = -0.999

function fit_cooling_tail(kind, elements, redshift, velocity, broadened, settings)
    if !broadened && !isnothing(velocity)
        throw(ArgumentError("Use a broadened cooling-flow model to fit velocity"))
    end
    values = element_parameters(elements, dem_parameter_names(kind))
    abundances = map(value -> fit_parameter(value, 0.0, 1000.0, kind isa Val{:partial}), values)
    z = fit_parameter(redshift, cooling_redshift_lower_limit(settings), 50.0, true)
    return broadened ? (abundances..., z, fit_parameter(something(velocity, 0.0), 0.0, 1e6, true)) :
        (abundances..., z)
end

for (_, constructor, model, mode, broadened) in COOLFLOW_VARIANTS
    @eval function $constructor(db::AtomDB; K = 1.0, lowT = 0.1, highT = 4.0, redshift = 0.1,
            slope = 0.0, velocity = nothing, switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            cosmology::CoolingCosmology = CoolingCosmology(), distance_mpc::Union{Nothing,Real} = nothing,
            normalization_redshift::Union{Nothing,Real} = nothing, bolometric::Symbol = :xspec,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            cache_spectra::Bool = true, elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        cache = spectrum_cache(cache_spectra, Val($(QuoteNode(mode))), options, DEM_INDEX)
        settings = CoolingOptions(cosmology; distance_mpc, normalization_redshift, ntemps, bolometric)
        tail = fit_cooling_tail(Val($(QuoteNode(mode))), elements, redshift, velocity, $broadened, settings)
        return $model(db, options, cache, settings, fit_parameter(K, 0.0, 1e24, false),
            fit_parameter(lowT, 0.0808, 79.9, false), fit_parameter(highT, 0.0808, 79.9, false), tail...,
            fit_parameter(slope, -5.0, 5.0, true))
    end
end

for (_, constructor, model, mode, broadened) in CPH_VARIANTS
    @eval function $constructor(db::AtomDB; K = 1.0, peakT = 2.2, redshift = 0.1,
            velocity = nothing, switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            cosmology::CoolingCosmology = CoolingCosmology(), distance_mpc::Union{Nothing,Real} = nothing,
            normalization_redshift::Union{Nothing,Real} = nothing, bolometric::Symbol = :xspec,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            cache_spectra::Bool = true, elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        cache = spectrum_cache(cache_spectra, Val($(QuoteNode(mode))), options, DEM_INDEX)
        settings = CoolingOptions(cosmology; distance_mpc, normalization_redshift, ntemps, bolometric)
        tail = fit_cooling_tail(Val($(QuoteNode(mode))), elements, redshift, velocity, $broadened, settings)
        return $model(db, options, cache, settings, fit_parameter(K, 0.0, 1e24, false),
            fit_parameter(peakT, 0.1, 100.0, false), tail...)
    end
end

for (_, constructor, model, mode, broadened) in LOCAL_COOLFLOW_VARIANTS
    @eval function $constructor(db::AtomDB; K = 1.0, lowT = 0.1, highT = 4.0, slope = 0.0,
            redshift = 0.0, velocity = nothing, switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            bolometric::Symbol = :xspec, abundance_table = :angr, trace_abundance = 1.0,
            thermal_broadening::Bool = $broadened, cache_spectra::Bool = true, elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        cache = spectrum_cache(cache_spectra, Val($(QuoteNode(mode))), options, DEM_INDEX)
        settings = LocalCoolingOptions(ntemps, bolometric)
        tail = fit_cooling_tail(Val($(QuoteNode(mode))), elements, redshift, velocity, $broadened, settings)
        return $model(db, options, cache, settings, fit_parameter(K, 0.0, 1e24, false),
            fit_parameter(lowT, 0.0808, 79.9, false), fit_parameter(highT, 0.0808, 79.9, false), tail...,
            fit_parameter(slope, -5.0, 5.0, true))
    end
end

for (_, constructor, model, mode, broadened) in LOCAL_CPH_VARIANTS
    @eval function $constructor(db::AtomDB; K = 1.0, peakT = 2.2, redshift = 0.0,
            velocity = nothing, switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            bolometric::Symbol = :xspec, abundance_table = :angr, trace_abundance = 1.0,
            thermal_broadening::Bool = $broadened, cache_spectra::Bool = true, elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        cache = spectrum_cache(cache_spectra, Val($(QuoteNode(mode))), options, DEM_INDEX)
        settings = LocalCoolingOptions(ntemps, bolometric)
        tail = fit_cooling_tail(Val($(QuoteNode(mode))), elements, redshift, velocity, $broadened, settings)
        return $model(db, options, cache, settings, fit_parameter(K, 0.0, 1e24, false),
            fit_parameter(peakT, 0.1, 100.0, false), tail...)
    end
end

function model_cooling_weights(model::AbstractCoolflowModel, abundances, z_norm)
    return cooling_weights(model.db, abundances, model.options, model.cooling, z_norm, model.lowT, model.highT;
        slope = model.slope)
end

function model_cooling_weights(model::AbstractCphModel, abundances, z_norm)
    return cooling_weights(model.db, abundances, model.options, model.cooling, z_norm, CPH_LIMITS...;
        peakT = model.peakT)
end

function SpectralFitting.invoke!(output, edges, model::AbstractCoolingModel)
    abundances = model_abundances(model)
    z_norm = cooling_norm_redshift(model.cooling, model.redshift)
    weights = model_cooling_weights(model, abundances, z_norm)
    return evaluate_cie!(output, edges, model.db, abundances, model.options, weights;
        redshift = model.redshift, normalization_redshift = z_norm, velocity = model_velocity(model), cache = model.cache)
end
