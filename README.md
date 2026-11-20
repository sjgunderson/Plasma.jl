# Plasma.jl

Plasma.jl evaluates the single-temperature XSPEC APEC family in Julia using
AtomDB's precomputed CIE emissivities. Supply the matching
`*_line.fits` and `*_coco.fits` files from AtomDB or an XSPEC installation.
The loader accepts fixed-width and variable-length continuum columns.

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
slices. The default reference table is `angr`. In plain APEC, the abundance parameter
scales C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, and Ni. H, He, and the trace
elements keep their reference-table values by default. Temperatures outside the tabulated range use
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

## APEC variants and abundance tables

| XSPEC | Julia spectrum | SpectralFitting model | Abundance parameters |
| --- | --- | --- | --- |
| `apec` | `apec` | `Apec` | One multiplier for 12 metals |
| `bapec` | `bapec` | `Bapec` | One multiplier for 12 metals |
| `vapec` | `vapec` | `Vapec` | He plus 12 metals |
| `vvapec` | `vvapec` | `Vvapec` | All 30 elements, H through Zn |
| `bvapec` | `bvapec` | `Bvapec` | He plus 12 metals |
| `bvvapec` | `bvvapec` | `Bvvapec` | All 30 elements, H through Zn |

The `b` models broaden resolved lines by default. `velocity` is the additional
Gaussian sigma in km/s. It adds in quadrature to the thermal width computed
from each temperature slice and element mass. Use `thermal_broadening=false`
for velocity-only broadening. The continuum and pseudo-continuum retain their
standard interpolation. The Julia implementation integrates the Gaussian over
each bin and includes wings down to eight sigma.

Use `abundance_table=:wilm` to set the reference per evaluation or model.
The built-in XSPEC tables are `angr`, `aspl`, `feld`, `aneb`, `grsa`, `wilm`,
`lodd`, `lpgp`, `lpgs`, and `felc`. `solar_abundances(:wilm).Fe` returns the Fe/H
number ratio. For a custom table, pass a path with 30 whitespace-separated
ratios in atomic-number order. Model abundance parameters are multipliers of
the selected table. For the `apec` and `vapec` families, set `trace_abundance`
to a multiplier or tie it to an element such as `:Fe`.

```julia
model = Bvapec(db; kT = 0.8, K = 1e-3, Fe = FitParam(0.5),
               velocity = FitParam(200.0), abundance_table = :wilm)
flux = invokemodel(edges, model)
```

Element and velocity parameters start frozen when supplied as numbers. Pass a
`FitParam` to thaw them or set their bounds. `Vvapec` and `Bvvapec` use
`K_abundance` for potassium because `K` is the additive normalization.
The direct `vvapec` and `bvvapec` functions accept `K` for potassium and
`norm` for normalization.

## Validation

Run `julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'` to check the
FITS loader, binning, interpolation, redshift, and abundance scaling against
synthetic tables. In an environment with Plasma and SpectralFitting installed,
run `julia test/spectralfitting.jl` to check the additive model interface.

For a direct comparison, run `julia test/xspec_parity.jl /path/to/apec_v3.0.9`
in an environment with Plasma and `LibXSPEC_jll` installed. The comparison pins
both evaluators to the supplied FITS pair, reports checksums and the XSPEC
library package version, and checks the bin-by-bin agreement on broad and
fine grating grids for plain APEC. The package tests also cover element
scaling, reference tables, and the broadened variants with synthetic tables.
