module Plasma

using FITSIO: FITS

export AtomDB, Apec, apec, apec!

include("atomdb.jl")
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

end
