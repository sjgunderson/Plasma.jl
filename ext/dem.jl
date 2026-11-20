using Plasma: WDEM_VARIANTS, CEMPOW_VARIANTS, dem_options, dem_abundance_names, dem_abundances,
    dem_node_weights, wdem_sampling, cempow_sampling, evaluate_cie!
import Plasma: Wdem, Vwdem, Vvwdem, Bwdem, Bvwdem, Bvvwdem, Cempow, Vcempow, Bcempow, Bvcempow

abstract type AbstractDemModel{T} <: AbstractCIEModel{T} end
abstract type AbstractWdemModel{T} <: AbstractDemModel{T} end
abstract type AbstractCempowModel{T} <: AbstractDemModel{T} end

dem_parameter_names(::Val{:tied}) = (:abundance,)
dem_parameter_names(::Val{:partial}) = dem_abundance_names(Val(:partial))
dem_parameter_names(::Val{:full}) = VVAPEC_PARAMETER_NAMES

for (variants, supertype, temperature_parameters) in (
        (WDEM_VARIANTS, :AbstractWdemModel, (:Tmax, :beta, :inv_slope)),
        (CEMPOW_VARIANTS, :AbstractCempowModel, (:alpha, :Tmax)))
    for (_, _, model, mode, broadened) in variants
        kind = Val(mode)
        abundance_parameters = dem_parameter_names(kind)
        names = (:K, temperature_parameters..., abundance_parameters..., :redshift)
        broadened && (names = (names..., :velocity))
        fields = [:($parameter::T) for parameter in names]
        @eval struct $model{T} <: $supertype{T}
            db::AtomDB
            options::ApecOptions
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

function fit_dem_tail(kind, elements, redshift, velocity)
    limit = kind isa Val{:full} ? 1000.0 : 10.0
    values = element_parameters(elements, dem_parameter_names(kind))
    abundances = map(value -> fit_parameter(value, 0.0, limit, true), values)
    z = fit_parameter(redshift, -0.999, 10.0, true)
    return isnothing(velocity) ? (abundances..., z) : (abundances..., z, fit_parameter(velocity, 0.0, 1e6, true))
end

for (_, constructor, model, mode, broadened) in WDEM_VARIANTS
    @eval function $constructor(db::AtomDB; K = 1.0, Tmax = 1.0, beta = 0.1, inv_slope = 0.25,
            redshift = 0.0, velocity = nothing, nH::Real = 1.0, switch::Integer = 2,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        if !$broadened && !isnothing(velocity)
            throw(ArgumentError("Use a broadened Wdem model to fit velocity"))
        end
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, nH, switch)
        tail = fit_dem_tail(Val($(QuoteNode(mode))), elements, redshift, $broadened ? something(velocity, 0.0) : nothing)
        return $model(db, options, common_parameters(K, Tmax, (0.01, 20.0))...,
            fit_parameter(beta, 0.01, 1.0, false), fit_parameter(inv_slope, -1.0, 10.0, false), tail...)
    end
end

for (_, constructor, model, mode, broadened) in CEMPOW_VARIANTS
    @eval function $constructor(db::AtomDB; K = 1.0, alpha = 1.0, Tmax = 1.0,
            redshift = 0.0, velocity = nothing, nH::Real = 1.0, switch::Integer = 2,
            abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = $broadened,
            elements...)
        if !$broadened && !isnothing(velocity)
            throw(ArgumentError("Use a broadened Cempow model to fit velocity"))
        end
        options = dem_options(abundance_table, trace_abundance, thermal_broadening, nH, switch)
        tail = fit_dem_tail(Val($(QuoteNode(mode))), elements, redshift, $broadened ? something(velocity, 0.0) : nothing)
        return $model(db, options, fit_parameter(K, 0.0, 1e24, false), fit_parameter(alpha, 0.01, 20.0, true),
            fit_parameter(Tmax, 0.02725, 100.0, false), tail...)
    end
end

dem_sampling(model::AbstractWdemModel) = wdem_sampling(model.Tmax, model.beta, model.inv_slope)
dem_sampling(model::AbstractCempowModel) = cempow_sampling(model.alpha, model.Tmax)

function SpectralFitting.invoke!(output, edges, model::AbstractDemModel)
    weights = dem_node_weights(model.db, dem_sampling(model))
    return evaluate_cie!(output, edges, model.db, model_abundances(model), model.options, weights;
        redshift = model.redshift, velocity = model_velocity(model))
end
