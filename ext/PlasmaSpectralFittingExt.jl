module PlasmaSpectralFittingExt

using Plasma: AtomDB, apec!
import Plasma: Apec
using SpectralFitting: AbstractSpectralModel, Additive, FitParam
import SpectralFitting

struct ApecModel{T} <: AbstractSpectralModel{T,Additive}
    db::AtomDB
    K::T
    kT::T
    abundance::T
    redshift::T
end

function Apec(db::AtomDB; K = 1.0, kT = 1.0, abundance = 1.0, redshift = 0.0)
    return ApecModel(db,
        fit_parameter(K, 0.0, 1e24, false),
        fit_parameter(kT, 0.008, 64.0, false),
        fit_parameter(abundance, 0.0, 5.0, true),
        fit_parameter(redshift, -0.999, 10.0, true))
end

function fit_parameter(value::Real, lower_limit, upper_limit, frozen)
    return FitParam(Float64(value); lower_limit, upper_limit, frozen)
end

function fit_parameter(param::FitParam, _lower_limit, _upper_limit, _frozen)
    # Each component owns its fit parameters. Use one numeric type for all four.
    return FitParam(Float64(param.value); error = param.error,
        lower_limit = param.lower_limit, upper_limit = param.upper_limit, frozen = param.frozen)
end

function SpectralFitting.invoke!(output, edges, model::ApecModel)
    # SpectralFitting applies K when invoking an additive model.
    return apec!(output, edges, model.db; kT = model.kT, abundance = model.abundance,
        redshift = model.redshift)
end

end
