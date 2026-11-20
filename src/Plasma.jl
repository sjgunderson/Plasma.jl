module Plasma

using FITSIO: FITS
using SpecialFunctions: erf

export AtomDB, solar_abundances
export Apec, Bapec, Vapec, Vvapec, Bvapec, Bvvapec
export apec, apec!, bapec, vapec, vvapec, bvapec, bvvapec

include("abundances.jl")
include("atomdb.jl")
include("broadening.jl")
include("apec.jl")

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

end
