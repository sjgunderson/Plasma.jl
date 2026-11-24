module PlasmaSpectralFittingExt

using Plasma: AtomDB, ApecOptions, ELEMENT_NAMES, VAPEC_NAMES, apec_options,
    element_parameters, evaluate_cie!, tied_abundances, partial_abundances, line_search_temperature,
    SpectrumCache, spectrum_cache, VAPEC_INDEX, DEM_INDEX
import Plasma: Apec, Bapec, Vapec, Vvapec, Bvapec, Bvvapec
import Plasma: find_lines, line_model_state, cache_value, cache_constant
using SpectralFitting: AbstractSpectralModel, Additive, FitParam
import SpectralFitting

abstract type AbstractCIEModel{T} <: AbstractSpectralModel{T,Additive} end
abstract type AbstractApecModel{T} <: AbstractCIEModel{T} end

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
    @eval struct $name{T,D<:AtomDB} <: AbstractApecModel{T}
        db::D
        options::ApecOptions
        cache::Union{Nothing,SpectrumCache}
        $(fields...)
    end
end

function Apec(db::AtomDB; K = 1.0, kT = 1.0, abundance = 1.0, redshift = 0.0,
        abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = false,
        cache_spectra::Bool = true)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    cache = spectrum_cache(cache_spectra, Val(:tied), options, VAPEC_INDEX)
    return ApecModel(db, options, cache, common_parameters(K, kT, (0.008, 64.0))...,
        fit_parameter(abundance, 0.0, 5.0, true), fit_parameter(redshift, -0.999, 10.0, true))
end

function Bapec(db::AtomDB; K = 1.0, kT = 1.0, abundance = 1.0, redshift = 0.0, velocity = 0.0,
        abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = true,
        cache_spectra::Bool = true)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    cache = spectrum_cache(cache_spectra, Val(:tied), options, VAPEC_INDEX)
    return BapecModel(db, options, cache, common_parameters(K, kT, (0.008, 64.0))...,
        fit_parameter(abundance, 0.0, 5.0, true), fit_parameter(redshift, -0.999, 10.0, true),
        fit_parameter(velocity, 0.0, 1e6, true))
end

function Vapec(db::AtomDB; K = 1.0, kT = 6.5, redshift = 0.0, abundance_table = :angr,
        trace_abundance = 1.0, thermal_broadening::Bool = false,
        cache_spectra::Bool = true, elements...)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    cache = spectrum_cache(cache_spectra, Val(:partial), options, VAPEC_INDEX)
    return VapecModel(db, options, cache, variable_parameters(K, kT, redshift, elements, VAPEC_NAMES)...)
end

function Vvapec(db::AtomDB; K = 1.0, kT = 6.5, redshift = 0.0, abundance_table = :angr,
        thermal_broadening::Bool = false,
        cache_spectra::Bool = true, elements...)
    options = apec_options(abundance_table, 1.0, thermal_broadening)
    cache = spectrum_cache(cache_spectra, Val(:full), options, VAPEC_INDEX)
    return VvapecModel(db, options, cache, variable_parameters(K, kT, redshift, elements, VVAPEC_PARAMETER_NAMES)...)
end

function Bvapec(db::AtomDB; K = 1.0, kT = 6.5, redshift = 0.0, velocity = 0.0,
        abundance_table = :angr, trace_abundance = 1.0, thermal_broadening::Bool = true,
        cache_spectra::Bool = true, elements...)
    options = apec_options(abundance_table, trace_abundance, thermal_broadening)
    cache = spectrum_cache(cache_spectra, Val(:partial), options, VAPEC_INDEX)
    return BvapecModel(db, options, cache, variable_parameters(K, kT, redshift, elements, VAPEC_NAMES)...,
        fit_parameter(velocity, 0.0, 1e6, true))
end

function Bvvapec(db::AtomDB; K = 1.0, kT = 6.5, redshift = 0.0, velocity = 0.0,
        abundance_table = :angr, thermal_broadening::Bool = true,
        cache_spectra::Bool = true, elements...)
    options = apec_options(abundance_table, 1.0, thermal_broadening)
    cache = spectrum_cache(cache_spectra, Val(:full), options, VAPEC_INDEX)
    return BvvapecModel(db, options, cache, variable_parameters(K, kT, redshift, elements, VVAPEC_PARAMETER_NAMES)...,
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
model_velocity(::AbstractCIEModel) = 0.0
model_velocity(model::Union{BapecModel,BvapecModel,BvvapecModel}) = model.velocity

function Base.copy(model::AbstractCIEModel{<:FitParam})
    parameters = copy.(SpectralFitting.parameter_vector(model))
    # Reuse the AtomDB cache and immutable options across copied components.
    return SpectralFitting.remake_with_parameters(model, parameters)
end

Base.copy(model::AbstractCIEModel{<:Number}) = model

const FD = SpectralFitting.ForwardDiff
cache_value(x::FD.Dual) = cache_value(FD.value(x))
function cache_constant(x::FD.Dual)
    all(partial -> cache_constant(partial) == 0.0, FD.partials(x)) || return nothing
    return cache_constant(FD.value(x))
end

function SpectralFitting.invoke!(output, edges, model::AbstractApecModel)
    # SpectralFitting applies K when invoking an additive model.
    return evaluate_cie!(output, edges, model.db, model_abundances(model), model.options, model.kT;
        redshift = model.redshift, velocity = model_velocity(model), cache = model.cache)
end

include("dem.jl")
include("cooling.jl")

find_lines(model::AbstractCIEModel, lo::Real, hi::Real; options...) =
    find_lines(model.db, lo, hi; options..., model)

function line_model_state(db::AtomDB, model::AbstractCIEModel{<:FitParam})
    values = map(parameter -> parameter.value, SpectralFitting.parameter_vector(model))
    return line_model_state(db, SpectralFitting.remake_with_parameters(model, values))
end

function line_model_state(db::AtomDB, model::AbstractCIEModel{<:Number})
    model.db === db || throw(ArgumentError("Search with the model's AtomDB object"))
    abundances = model_abundances(model)
    all(value -> isfinite(value) && value >= 0, abundances) ||
        throw(DomainError(abundances, "Element abundances must be nonnegative and finite"))
    scales = map(*, promote(abundances...), model.options.ratios)
    weights = line_model_weights(model, abundances)
    total = sum(weights)
    return iszero(total) ? (weights, scales) : (weights ./ total, scales)
end

line_model_weights(model::AbstractApecModel, _abundances) =
    dem_node_weights(model.db, ((line_search_temperature(model.kT, :keV),), (1.0,)))
line_model_weights(model::AbstractDemModel, _abundances) =
    dem_node_weights(model.db, dem_sampling(model))
line_model_weights(model::AbstractCoolingModel, abundances) =
    model_cooling_weights(model, abundances, cooling_norm_redshift(model.cooling, model.redshift))

end
