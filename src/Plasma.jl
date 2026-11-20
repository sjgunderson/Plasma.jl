module Plasma

using FITSIO: FITS
using SpecialFunctions: erf

export AtomDB, solar_abundances
export Apec, Bapec, Vapec, Vvapec, Bvapec, Bvvapec
export apec, apec!, bapec, vapec, vvapec, bvapec, bvvapec
export Wdem, Vwdem, Vvwdem, Bwdem, Bvwdem, Bvvwdem
export Cempow, Vcempow, Bcempow, Bvcempow
export wdem, vwdem, vvwdem, bwdem, bvwdem, bvvwdem
export cempow, vcempow, bcempow, bvcempow
export Coolflow, Vcoolflow, Bcoolflow, Bvcoolflow, Mkcflow, Vmcflow
export Cph, Vcph, Bcph, Bvcph, CoolingCosmology
export coolflow, vcoolflow, bcoolflow, bvcoolflow, mkcflow, vmcflow
export cph, vcph, bcph, bvcph

include("abundances.jl")
include("atomdb.jl")
include("broadening.jl")
include("apec.jl")
include("dem.jl")
include("cooling.jl")

"""
    Apec(db::AtomDB; K, kT, abundance, redshift)

Construct a fit-ready additive APEC model after loading SpectralFitting.jl.
`K` is the XSPEC normalization, `kT` is in keV, and `abundance` is the
common abundance of the 12 metals fitted by plain APEC. The initial values
are `1.0`, `1.0`, `1.0`, and `0.0`, respectively. Pass a number to set a
starting value or a `SpectralFitting.FitParam` to customize the fit bounds.
"""
function Apec end

"""Fit-ready APEC model with thermal and velocity broadening. `velocity` is in km/s."""
function Bapec end

"""Fit-ready APEC model with independent He and 12 metal abundances, e.g. `Vapec(db; Fe=0.4)`."""
function Vapec end

"""Fit-ready APEC model with independent abundances for H through Zn."""
function Vvapec end

"""Fit-ready Vapec model with thermal and velocity broadening. `velocity` is in km/s."""
function Bvapec end

"""Fit-ready Vvapec model with thermal and velocity broadening. `velocity` is in km/s."""
function Bvvapec end

# The SpectralFitting extension provides the constructors for these model types.
for (_, constructor, _, _, _) in (WDEM_VARIANTS..., CEMPOW_VARIANTS..., COOLFLOW_VARIANTS..., CPH_VARIANTS...)
    @eval function $constructor end
end

const Mkcflow = Coolflow
const Vmcflow = Vcoolflow

end
