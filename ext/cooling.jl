using Plasma: COOLFLOW_VARIANTS, CPH_VARIANTS, COOLING_NTEMPS, CPH_LIMITS, CoolingCosmology,
    CoolingOptions, cooling_weights
import Plasma: Coolflow, Vcoolflow, Bcoolflow, Bvcoolflow, Cph, Vcph, Bcph, Bvcph

abstract type AbstractCoolingModel{T} <: AbstractCIEModel{T} end
abstract type AbstractCoolflowModel{T} <: AbstractCoolingModel{T} end
abstract type AbstractCphModel{T} <: AbstractCoolingModel{T} end

for (variants, supertype, temperature_parameters) in (
        (COOLFLOW_VARIANTS, :AbstractCoolflowModel, (:lowT, :highT)),
        (CPH_VARIANTS, :AbstractCphModel, (:peakT,)))
    for (_, _, model, mode, broadened) in variants
        kind = Val(mode)
        abundance_parameters = dem_parameter_names(kind)
        names = (:K, temperature_parameters..., abundance_parameters..., :redshift)
        broadened && (names = (names..., :velocity))
        fields = [:($parameter::T) for parameter in names]
        @eval struct $model{T} <: $supertype{T}
            db::AtomDB
            options::ApecOptions
            cooling::CoolingOptions
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

function fit_cooling_tail(kind, elements, redshift, velocity, broadened, distance_mpc)
    if !broadened && !isnothing(velocity)
        throw(ArgumentError("Use a broadened cooling-flow model to fit velocity"))
    end
    values = element_parameters(elements, dem_parameter_names(kind))
    abundances = map(value -> fit_parameter(value, 0.0, 1000.0, kind isa Val{:partial}), values)
    z = fit_parameter(redshift, isnothing(distance_mpc) ? 0.0 : -0.999, 50.0, true)
    return broadened ? (abundances..., z, fit_parameter(something(velocity, 0.0), 0.0, 1e6, true)) :
        (abundances..., z)
end

for (_, constructor, model, mode, broadened) in COOLFLOW_VARIANTS
    @eval function $constructor(db::AtomDB; K = 1.0, lowT = 0.1, highT = 4.0, redshift = 0.1,
            velocity = nothing, switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            cosmology::CoolingCosmology = CoolingCosmology(), distance_mpc::Union{Nothing,Real} = nothing,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        settings = CoolingOptions(cosmology, distance_mpc, ntemps)
        tail = fit_cooling_tail(Val($(QuoteNode(mode))), elements, redshift, velocity, $broadened, distance_mpc)
        return $model(db, options, settings, fit_parameter(K, 0.0, 1e24, false),
            fit_parameter(lowT, 0.0808, 79.9, false), fit_parameter(highT, 0.0808, 79.9, false), tail...)
    end
end

for (_, constructor, model, mode, broadened) in CPH_VARIANTS
    @eval function $constructor(db::AtomDB; K = 1.0, peakT = 2.2, redshift = 0.1,
            velocity = nothing, switch::Integer = 2, ntemps::Integer = COOLING_NTEMPS,
            cosmology::CoolingCosmology = CoolingCosmology(), distance_mpc::Union{Nothing,Real} = nothing,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, 1.0, switch)
        settings = CoolingOptions(cosmology, distance_mpc, ntemps)
        tail = fit_cooling_tail(Val($(QuoteNode(mode))), elements, redshift, velocity, $broadened, distance_mpc)
        return $model(db, options, settings, fit_parameter(K, 0.0, 1e24, false),
            fit_parameter(peakT, 0.1, 100.0, false), tail...)
    end
end

function model_cooling_weights(model::AbstractCoolflowModel, abundances)
    return cooling_weights(model.db, abundances, model.options, model.cooling, model.redshift, model.lowT, model.highT)
end

function model_cooling_weights(model::AbstractCphModel, abundances)
    return cooling_weights(model.db, abundances, model.options, model.cooling, model.redshift, CPH_LIMITS...;
        peakT = model.peakT)
end

function SpectralFitting.invoke!(output, edges, model::AbstractCoolingModel)
    abundances = model_abundances(model)
    weights = model_cooling_weights(model, abundances)
    return evaluate_cie!(output, edges, model.db, abundances, model.options, weights;
        redshift = model.redshift, velocity = model_velocity(model))
end
