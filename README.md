# Plasma.jl

Plasma.jl evaluates single-temperature APEC, multi-temperature emission-measure
distributions, and cooling-flow models in Julia using AtomDB's precomputed
collisional ionization equilibrium (CIE) emissivities. Supply the matching
`*_line.fits` and `*_coco.fits` files from AtomDB or an XSPEC installation.
The loader accepts fixed-width and variable-length continuum columns.

## Available models

| Base model | What it describes | Temperature parameters | Normalization (`K` or `norm`) |
| --- | --- | --- | --- |
| [`Apec`](#apec-variants-and-abundance-tables) | Plasma at one temperature | `kT` | XSPEC APEC emission-measure normalization |
| [`Wdem`](#multi-temperature-models) | Power-law differential emission measure with a low-temperature cutoff | `Tmax`, `beta`, `inv_slope` | APEC scale applied to weights that sum to one |
| [`Cempow`](#multi-temperature-models) | Power-law emission measure sampled at logarithmic temperature intervals | `alpha`, `Tmax` | APEC scale applied to the absolute weights |
| [`Coolflow`](#cooling-flow-models) | Isobaric cooling with optional power-law emissivity weighting | `lowT`, `highT`, `slope` | Solar masses per year |
| [`Cph`](#cooling-flow-models) | Cooling and heating in a multi-temperature plasma | `peakT` | Solar masses per year |
| [`LocalCoolflow`](#cooling-flow-models) | Cooling flow with a signed spectral shift for a nearby source | `lowT`, `highT`, `slope` | 10⁻¹⁴ solar masses per year per square parsec |
| [`LocalCph`](#cooling-flow-models) | Cooling and heating with a signed spectral shift for a nearby source | `peakT` | 10⁻¹⁴ solar masses per year per square parsec |

Temperatures are in keV. The lowercase functions evaluate spectra directly
and take `norm`. The uppercase constructors create additive models for
SpectralFitting.jl and take `K`. The family tables below list the available
abundance and line-broadening variants.

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
slices. Temperatures outside the tabulated range use the nearest tabulated
spectrum.

## Fit with SpectralFitting.jl

```julia
using Plasma, SpectralFitting

model = Apec(db; kT = 0.8, K = 1e-3, abundance = FitParam(0.5), redshift = 0.01)
flux = invokemodel(edges, model)
```

For `Apec(db)`, abundance and redshift start frozen. Use `FitParam` for custom
bounds or to thaw a parameter. The adapter supports SpectralFitting 0.7.3 and
0.7.4. The models use FITSIO.jl to read the tables and run without an XSPEC
runtime dependency.

Fit-ready models cache AtomDB spectra on the current energy grid once the grid,
redshift, and velocity repeat. Temperature weights and normalization are still
evaluated at each trial. Free abundances use separate contributions for the
affected elements or tied abundance groups. The cache retains up to 64 MiB of
spectral templates per component, and trials with active redshift or velocity
derivatives use direct synthesis. With multiple Julia threads, direct synthesis
can split broadened models across temperature nodes.

Use `cache_spectra=false` on a fit-ready model to compare with direct synthesis.
The [fit-grid cache notebook](benchmarks/fit_grid_cache.ipynb) checks spectrum
and gradient agreement and measures warm objective and gradient times.

## APEC variants and abundance tables

| XSPEC | Julia spectrum | SpectralFitting model | Abundance parameters |
| --- | --- | --- | --- |
| `apec` | `apec` | `Apec` | One multiplier for 12 metals |
| `bapec` | `bapec` | `Bapec` | One multiplier for 12 metals |
| `vapec` | `vapec` | `Vapec` | He plus 12 metals |
| `vvapec` | `vvapec` | `Vvapec` | All 30 elements, H through Zn |
| `bvapec` | `bvapec` | `Bvapec` | He plus 12 metals |
| `bvvapec` | `bvvapec` | `Bvvapec` | All 30 elements, H through Zn |

In plain APEC, `abundance` scales C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe,
and Ni. H, He, and the trace elements keep their reference-table values by
default.

The `b` models broaden resolved lines by default. `velocity` is the additional
Gaussian sigma in km/s. It adds in quadrature to the thermal width computed
from each temperature slice and element mass. Use `thermal_broadening=false`
for velocity-only broadening. The continuum and pseudo-continuum retain their
standard interpolation. The Julia implementation integrates the Gaussian over
each bin and includes wings down to eight sigma.

The default abundance table is `angr`. Use `abundance_table=:wilm` to set the
reference per evaluation or model.
The built-in XSPEC tables are `angr`, `aspl`, `feld`, `aneb`, `grsa`, `wilm`,
`lodd`, `lpgp`, `lpgs`, and `felc`. `solar_abundances(:wilm).Fe` returns the Fe/H
number ratio. For a custom table, pass a path with 30 whitespace-separated
ratios in atomic-number order. Model abundance parameters are multipliers of
the selected table. In variants with tied or selected-element abundances, set
`trace_abundance` to a multiplier or tie it to an element such as `:Fe`.

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

## Multi-temperature models

| Family | SpectralFitting constructors | Standalone evaluators |
| --- | --- | --- |
| Wdem | `Wdem`, `Vwdem`, `Vvwdem`, `Bwdem`, `Bvwdem`, `Bvvwdem` | `wdem`, `vwdem`, `vvwdem`, `bwdem`, `bvwdem`, `bvvwdem` |
| Cempow | `Cempow`, `Vcempow`, `Bcempow`, `Bvcempow` | `cempow`, `vcempow`, `bcempow`, `bvcempow` |

`Wdem` takes `Tmax`, `beta = Tmin/Tmax`, and `inv_slope`, which default to
1.0 keV, 0.1, and 0.25, respectively. The XSPEC
sampling integrates `T^(1/inv_slope)` in 51 intervals of width
`(Tmax - Tmin)/50`, including the final interval above `Tmax`. Its weights sum
to 1. When `beta=1`, Plasma uses the single-temperature limit at `Tmax`.

`Cempow` takes the power-law exponent `alpha` and `Tmax`, which default to
1.0 and 1.0 keV. It steps down from `Tmax` by 0.1 dex and keeps samples
at or above 0.02725 keV. It retains the absolute `0.1*(T/Tmax)^alpha`
weights, whose sum varies with `alpha`.

The standalone `norm` and SpectralFitting `K` scale the resulting spectrum.
The `V` models fit He and 13 metals, including Na. The `Vv` models fit H
through Zn, using `K_abundance` for potassium in SpectralFitting. The `B`
models enable thermal broadening by default and accept a `velocity` sigma
in km/s. The abundance and trace settings follow the APEC interface.
These models use the AtomDB CIE backend (`switch=2`). `nH` is accepted as a
compatibility keyword and does not enter the AtomDB calculation.

```julia
source = Bvwdem(db; Tmax=1.5, beta=FitParam(0.15), inv_slope=0.25,
                Fe=FitParam(0.5), velocity=FitParam(200.0), K=1e-3)
flux = invokemodel(edges, source)

other_flux = cempow(edges, db; Tmax=1.5, alpha=0.8, norm=1e-3)
```

## Cooling-flow models

| Family | SpectralFitting constructors | Standalone evaluators |
| --- | --- | --- |
| Isobaric cooling flow | `Coolflow`, `Vcoolflow`, `Bcoolflow`, `Bvcoolflow` | `coolflow`, `vcoolflow`, `bcoolflow`, `bvcoolflow` |
| Legacy aliases | `Mkcflow`, `Vmcflow` | `mkcflow`, `vmcflow` |
| Cooling plus heating | `Cph`, `Vcph`, `Bcph`, `Bvcph` | `cph`, `vcph`, `bcph`, `bvcph` |
| Local cooling flow | `LocalCoolflow`, `LocalVcoolflow`, `LocalBcoolflow`, `LocalBvcoolflow` | `localcoolflow`, `localvcoolflow`, `localbcoolflow`, `localbvcoolflow` |
| Local cooling plus heating | `LocalCph`, `LocalVcph`, `LocalBcph`, `LocalBvcph` | `localcph`, `localvcph`, `localbcph`, `localbvcph` |

The cooling-flow models use the AtomDB CIE backend (`switch=2`).
`Mkcflow` and `Vmcflow` are aliases for `Coolflow` and `Vcoolflow`.
The isobaric models integrate from `lowT` to `highT` (default 0.1 to 4.0 keV),
weighting the emission measure by the inverse bolometric emissivity. The
cooling-plus-heating models take `peakT` (default 2.2 keV) and integrate from
0.01 to 50 keV. The temperature integral uses 10 logarithmic samples by default.
Set `ntemps` to change the sample count, analogous to XSPEC's `CFLOW_NTEMPS`.

All `Coolflow` and `LocalCoolflow` variants accept a dimensionless `slope`.
It multiplies each temperature's emission measure by `(T / highT)^slope`.
The factor is one at the upper temperature. Positive slopes suppress the
cooler contributions. The default `slope=0.0` gives ordinary isobaric cooling flow.
Numeric slopes are frozen in SpectralFitting. Use `slope=FitParam(0.5)` to
fit the slope.

Plasma.jl uses AtomDB for both the spectra and the bolometric emissivity.
XSPEC's legacy [`cflow` model](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSmodelCflow.html)
uses a fixed 10⁷ K pivot and a Raymond–Smith table.

The cooling-plus-heating models follow the effective XSPEC 12.15 heating
prescription. At `peakT ≥ 1` keV the heating term uses a 1 keV emissivity
reference and scales as `(T/peakT)^2.5`. For `peakT < 1` keV the model uses
cooling alone. By default, `bolometric=:xspec` calculates the bolometric
emissivity in Julia using XSPEC's 1,000 logarithmic AtomDB energy edges
from 0.0105 to 1,000 keV and bin-midpoint weighting. Use
`bolometric=:direct` to integrate the tabulated continuum and use the
original energies of resolved lines. The choice can have a large effect on
`cph` where heating and cooling nearly cancel.

### Cosmological and local normalization

For the cosmological families, `K` and `norm` are mass accretion rates in
solar masses per year. The default `CoolingCosmology()` uses `H0=70.0`,
`q0=0.0`, and `lambda0=0.73` to determine distance from `redshift`, which
defaults to 0.1. Set `distance_mpc` to supply an independently measured distance.
By default, `redshift` also supplies the normalization's photon arrival-rate
factor. Set `normalization_redshift` to fix both redshift-dependent parts of
the normalization while varying the spectral shift with `redshift`.
A positive `normalization_redshift` can supply the cosmological distance.
The `V` and `B` variants follow the abundance and broadening conventions of
`Cempow`.

For the local families, `K` and `norm` are the mass accretion rate divided by
the square of the source distance, in units of 10⁻¹⁴ solar masses per year per
square parsec. Use `K = (mdot / distance_pc^2) / 1e-14` when converting a
physical accretion rate `mdot` in solar masses per year. The model includes
the geometric factor of `1/(4π)`. The signed `redshift` defaults to zero and
only shifts the spectrum. A negative value blueshifts the lines.

```julia
model = Bvcoolflow(db; lowT = 0.3, highT = 4.0, redshift = 0.08,
    K = 25.0, Fe = FitParam(0.5), velocity = 200.0)
flux = invokemodel(edges, model)

local_model = LocalCoolflow(db; lowT = 0.2, highT = 1.2,
    redshift = FitParam(-0.0002; lower_limit = -0.01, upper_limit = 0.01),
    K = 0.025)  # 1e-9 solar masses per year at 2,000 pc
local_flux = invokemodel(edges, local_model)
```

## Search AtomDB emission lines

Use the same `AtomDB` object that supplies the spectra. Give an energy or
wavelength interval and either a temperature or a native Plasma.jl model.
The search reads the CIE line tables and caches the transition IDs for the
temperatures it visits. A temperature-only search can run without loading
continuum slices.

```julia
# Rank O VII lines near the He-like triplet at a specified temperature.
lines = find_lines(db, 21.5, 22.2; element = :O, ion = 7, kT = 0.2, brightest = 3)

# Track one transition or compare two as a function of temperature.
temperatures = range(0.1, 0.5; length = 50)
emissivity = line_emissivity(db, lines[1].id, temperatures)
ratio = ratio_emissivity(db, lines[1].id, lines[2].id, temperatures)

# Rank by the current temperature distribution and abundances of a native model.
model = LocalBvcoolflow(db; lowT = 0.1, highT = 2.0, Fe = 0.5)
strongest = find_lines(db, 0.5, 2.0; unit = :keV, model, brightest = 10)
```

Bounds are inclusive. Use `unit=:angstrom` (default), `:nm`, `:keV`, or
`:eV`. Temperatures are in keV by default and `teunit=:K` accepts kelvin.
Filter by element symbol or atomic number and by spectroscopic ion stage.
The result exposes `wavelength` in Å, `energy` in keV, and `emissivity` in
photons cm³ s⁻¹. The `id` contains the element, ion, upper level, and lower
level. At a specified temperature, use `abundance_table` and
`abundances=(O=0.5, Fe=1.2)` to adjust the reference composition. The search
uses the same linear temperature interpolation and endpoint clamping as the
spectrum models.

Model ranking uses rest-frame emissivities averaged over the model's
abundance-scaled temperature mixture. Select an individual native APEC,
Wdem, Cempow, or cooling-flow component when the fitted source is composite.
Transition IDs are specific to the loaded AtomDB files. Detailed level
properties and Einstein A values require separate AtomDB atomic data files.

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
