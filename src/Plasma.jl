module Plasma

using FITSIO: FITS
using SpecialFunctions: erf

export AtomDB, solar_abundances
export LineID, LineMatch, find_lines, line_emissivity, ratio_emissivity
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
export LocalCoolflow, LocalVcoolflow, LocalBcoolflow, LocalBvcoolflow
export LocalCph, LocalVcph, LocalBcph, LocalBvcph
export localcoolflow, localvcoolflow, localbcoolflow, localbvcoolflow
export localcph, localvcph, localbcph, localbvcph

include("abundances.jl")
include("atomdb.jl")
include("broadening.jl")
include("apec.jl")
include("dem.jl")
include("cooling.jl")
include("line_search.jl")

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
for variants in (WDEM_VARIANTS, CEMPOW_VARIANTS, COOLFLOW_VARIANTS, CPH_VARIANTS,
        LOCAL_COOLFLOW_VARIANTS, LOCAL_CPH_VARIANTS)
    for (_, constructor, _, _, _) in variants
        @eval function $constructor end
    end
end

@doc """
    Coolflow(db::AtomDB; lowT=0.1, highT=4.0, slope=0.0, ...)

Fit-ready AtomDB cooling flow. The dimensionless `slope` multiplies the
emission measure at each temperature by `(T / highT)^slope`.
Pass a `FitParam` to fit the slope.
""" Coolflow

@doc """
    LocalCoolflow(db::AtomDB; slope=0.0, redshift=0.0, K=1.0, ...)

Fit-ready cooling flow with `K` equal to the mass accretion rate divided by
the squared source distance, in units of 10⁻¹⁴ solar masses per year per
square parsec.
The signed spectral `redshift` can be fitted independently.
""" LocalCoolflow

const Mkcflow = Coolflow
const Vmcflow = Vcoolflow

end
