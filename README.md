# Plasma.jl

Plasma.jl evaluates the plain single-temperature, unbroadened XSPEC APEC
model in Julia using AtomDB's precomputed CIE emissivities. Supply the matching
`*_line.fits` and `*_coco.fits` files from AtomDB or an XSPEC installation.
The initial target is the AtomDB 3.0.9 pair with 201 temperatures.

## Evaluate a spectrum

```julia
using Plasma

# Reads /path/to/apec_v3.0.9_line.fits and /path/to/apec_v3.0.9_coco.fits.
db = AtomDB("/path/to/apec_v3.0.9")
edges = range(0.3, 8.0; length = 2001)  # observed energy edges in keV
flux = apec(edges, db; kT = 0.8, abundance = 0.5, redshift = 0.01, norm = 1e-3)
```

`flux` contains photons cm⁻² s⁻¹ integrated over each bin. `apec!` accepts a
preallocated output vector. Reuse `db` across calls to cache loaded temperature
slices. The evaluator uses the `angr` abundance table. The abundance parameter
scales C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, and Ni. H, He, and the trace
elements keep their `angr` values. Temperatures outside the tabulated range use
the nearest tabulated spectrum.

## Fit with SpectralFitting.jl

```julia
using Plasma, SpectralFitting

model = Apec(db; kT = 0.8, K = 1e-3, abundance = FitParam(0.5), redshift = 0.01)
flux = invokemodel(edges, model)
```

The model works as an additive component. `K` is the XSPEC normalization. Use
`Apec(db)` for the plain XSPEC-style parameter defaults, with abundance and
redshift frozen. Use `FitParam` for custom bounds or to thaw a parameter.
The adapter supports SpectralFitting 0.7.3 and 0.7.4. Loading
and evaluating `apec` requires FITSIO.jl. The calculation has no XSPEC runtime
dependency.

## Validation

Run `julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'` to check the
FITS loader, binning, interpolation, redshift, and abundance scaling against
synthetic tables. In an environment with Plasma and SpectralFitting installed,
run `julia test/spectralfitting.jl` to check the additive model interface.

For a direct comparison, run `julia test/xspec_parity.jl /path/to/apec_v3.0.9`
in an environment with Plasma and `LibXSPEC_jll` installed. The comparison pins
both evaluators to the supplied FITS pair, reports checksums and the XSPEC
library package version, and checks the bin-by-bin agreement on broad and
fine grating grids. The first release covers the standard CIE path with
unbroadened lines, default temperature interpolation and trace abundances,
and standard XSPEC normalization.
