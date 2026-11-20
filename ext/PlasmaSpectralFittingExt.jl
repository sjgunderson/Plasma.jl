module PlasmaSpectralFittingExt

using Plasma: AtomDB, ApecOptions, ELEMENT_NAMES, VAPEC_NAMES, apec_options,
    element_parameters, evaluate_apec!, tied_abundances, partial_abundances
import Plasma: Apec, Bapec, Vapec, Vvapec, Bvapec, Bvvapec
using SpectralFitting: AbstractSpectralModel, Additive, FitParam
import SpectralFitting

abstract type AbstractApecModel{T} <: AbstractSpectralModel{T,Additive} end

# SpectralFitting reserves K for additive normalization. Keep potassium separate.
const VVAPEC_PARAMETER_NAMES = map(name -> name === :K ? :K_abundance : name, ELEMENT_NAMES)

# SpectralFitting discovers FitParams by field name and rebuilds models with
# numerical parameters during evaluation. Keep the database and options first.
for (name, abundance_names, broadened) in (
        (:ApecModel, (:abundance,), false), (:BapecModel, (:abundance,), true),
        (:VapecModel, VAPEC_NAMES, false), (:BvapecModel, VAPEC_NAMES, true),
        (:VvapecModel, VVAPEC_PARAMETER_NAMES, false), (:BvvapecModel, VVAPEC_PARAMETER_NAMES, true))
    names = (:K, :kT, abundance_names..., :redshift)
    broadened && (names = (names..., :velocity))
    fields = [:($parameter::T) for parameter in names]
    @eval struct $name{T} <: AbstractApecModel{T}
        db::AtomDB
        options::ApecOptions
        $(fields...)
    end
end

function Apec(db::AtomDB; K = 1.0, kT = 1.0, abundance = 1.0, redshift = 0.0,
        abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = false)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    return ApecModel(db, options, common_parameters(K, kT, (0.008, 64.0))...,
        fit_parameter(abundance, 0.0, 5.0, true), fit_parameter(redshift, -0.999, 10.0, true))
end

function Bapec(db::AtomDB; K = 1.0, kT = 1.0, abundance = 1.0, redshift = 0.0, velocity = 0.0,
        abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = true)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    return BapecModel(db, options, common_parameters(K, kT, (0.008, 64.0))...,
        fit_parameter(abundance, 0.0, 5.0, true), fit_parameter(redshift, -0.999, 10.0, true),
        fit_parameter(velocity, 0.0, 1e6, true))
end

function Vapec(db::AtomDB; K = 1.0, kT = 6.5, redshift = 0.0, abundance_table = :angr,
        trace_abundance = 1.0, thermal_broadening::Bool = false, elements...)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    return VapecModel(db, options, variable_parameters(K, kT, redshift, elements, VAPEC_NAMES)...)
end

function Vvapec(db::AtomDB; K = 1.0, kT = 6.5, redshift = 0.0, abundance_table = :angr,
        thermal_broadening::Bool = false, elements...)
    options = apec_options(abundance_table, 1.0, thermal_broadening)
    return VvapecModel(db, options, variable_parameters(K, kT, redshift, elements, VVAPEC_PARAMETER_NAMES)...)
end

function Bvapec(db::AtomDB; K = 1.0, kT = 6.5, redshift = 0.0, velocity = 0.0,
        abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = true, elements...)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    return BvapecModel(db, options, variable_parameters(K, kT, redshift, elements, VAPEC_NAMES)...,
        fit_parameter(velocity, 0.0, 1e6, true))
end

function Bvvapec(db::AtomDB; K = 1.0, kT = 6.5, redshift = 0.0, velocity = 0.0,
        abundance_table = :angr, thermal_broadening::Bool = true, elements...)
    options = apec_options(abundance_table, 1.0, thermal_broadening)
    return BvvapecModel(db, options, variable_parameters(K, kT, redshift, elements, VVAPEC_PARAMETER_NAMES)...,
        fit_parameter(velocity, 0.0, 1e6, true))
end

function common_parameters(K, kT, limits)
    return fit_parameter(K, 0.0, 1e24, false), fit_parameter(kT, limits..., false)
end

function variable_parameters(K, kT, redshift, elements, names)
    abundances = map(value -> fit_parameter(value, 0.0, 1000.0, true), element_parameters(elements, names))
    return (common_parameters(K, kT, (0.0808, 68.447))..., abundances...,
        fit_parameter(redshift, -0.999, 10.0, true))
end

function fit_parameter(value::Real, lower_limit, upper_limit, frozen)
    return FitParam(Float64(value); lower_limit, upper_limit, frozen)
end

function fit_parameter(param::FitParam, _lower_limit, _upper_limit, _frozen)
    # Each component owns its fit parameters. Use one numeric type for all parameters.
    return FitParam(Float64(param.value); error = param.error,
        lower_limit = param.lower_limit, upper_limit = param.upper_limit, frozen = param.frozen)
end

model_abundances(model::Union{ApecModel,BapecModel}) = tied_abundances(model.abundance, model.options.trace)

function model_abundances(model::Union{VapecModel,BvapecModel})
    values = map(name -> getfield(model, name), VAPEC_NAMES)
    return partial_abundances(values, model.options.trace)
end

model_abundances(model::Union{VvapecModel,BvvapecModel}) =
    map(name -> getfield(model, name), VVAPEC_PARAMETER_NAMES)
model_velocity(::AbstractApecModel) = 0.0
model_velocity(model::Union{BapecModel,BvapecModel,BvvapecModel}) = model.velocity

function Base.copy(model::AbstractApecModel{<:FitParam})
    parameters = copy.(SpectralFitting.parameter_vector(model))
    # Reuse the AtomDB cache and immutable options across copied components.
    return SpectralFitting.remake_with_parameters(model, parameters)
end

Base.copy(model::AbstractApecModel{<:Number}) = model

function SpectralFitting.invoke!(output, edges, model::AbstractApecModel)
    # SpectralFitting applies K when invoking an additive model.
    return evaluate_apec!(output, edges, model.db, model_abundances(model), model.options;
        kT = model.kT, redshift = model.redshift, velocity = model_velocity(model))
end

end
