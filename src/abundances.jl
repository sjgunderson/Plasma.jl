# XSPEC abund tables, in linear number ratios relative to H (atomic numbers 1–30).
# https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSabund.html
const ELEMENT_NAMES = (
    :H, :He, :Li, :Be, :B, :C, :N, :O, :F, :Ne,
    :Na, :Mg, :Al, :Si, :P, :S, :Cl, :Ar, :K, :Ca,
    :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn,
)

const XSPEC_ABUNDANCES = (
    angr = (
        1.00e+00, 9.77e-02, 1.45e-11, 1.41e-11, 3.98e-10, 3.63e-04, 1.12e-04, 8.51e-04, 3.63e-08, 1.23e-04,
        2.14e-06, 3.80e-05, 2.95e-06, 3.55e-05, 2.82e-07, 1.62e-05, 3.16e-07, 3.63e-06, 1.32e-07, 2.29e-06,
        1.26e-09, 9.77e-08, 1.00e-08, 4.68e-07, 2.45e-07, 4.68e-05, 8.32e-08, 1.78e-06, 1.62e-08, 3.98e-08,
    ),
    aspl = (
        1.00e+00, 8.51e-02, 1.12e-11, 2.40e-11, 5.01e-10, 2.69e-04, 6.76e-05, 4.90e-04, 3.63e-08, 8.51e-05,
        1.74e-06, 3.98e-05, 2.82e-06, 3.24e-05, 2.57e-07, 1.32e-05, 3.16e-07, 2.51e-06, 1.07e-07, 2.19e-06,
        1.41e-09, 8.91e-08, 8.51e-09, 4.37e-07, 2.69e-07, 3.16e-05, 9.77e-08, 1.66e-06, 1.55e-08, 3.63e-08,
    ),
    feld = (
        1.00e+00, 9.77e-02, 1.26e-11, 2.51e-11, 3.55e-10, 3.98e-04, 1.00e-04, 8.51e-04, 3.63e-08, 1.29e-04,
        2.14e-06, 3.80e-05, 2.95e-06, 3.55e-05, 2.82e-07, 1.62e-05, 3.16e-07, 4.47e-06, 1.32e-07, 2.29e-06,
        1.48e-09, 1.05e-07, 1.00e-08, 4.68e-07, 2.45e-07, 3.24e-05, 8.32e-08, 1.78e-06, 1.62e-08, 3.98e-08,
    ),
    aneb = (
        1.00e+00, 8.01e-02, 2.19e-09, 2.87e-11, 8.82e-10, 4.45e-04, 9.12e-05, 7.39e-04, 3.10e-08, 1.38e-04,
        2.10e-06, 3.95e-05, 3.12e-06, 3.68e-05, 3.82e-07, 1.89e-05, 1.93e-07, 3.82e-06, 1.39e-07, 2.25e-06,
        1.24e-09, 8.82e-08, 1.08e-08, 4.93e-07, 3.50e-07, 3.31e-05, 8.27e-08, 1.81e-06, 1.89e-08, 4.63e-08,
    ),
    grsa = (
        1.00e+00, 8.51e-02, 1.26e-11, 2.51e-11, 3.55e-10, 3.31e-04, 8.32e-05, 6.76e-04, 3.63e-08, 1.20e-04,
        2.14e-06, 3.80e-05, 2.95e-06, 3.55e-05, 2.82e-07, 2.14e-05, 3.16e-07, 2.51e-06, 1.32e-07, 2.29e-06,
        1.48e-09, 1.05e-07, 1.00e-08, 4.68e-07, 2.45e-07, 3.16e-05, 8.32e-08, 1.78e-06, 1.62e-08, 3.98e-08,
    ),
    wilm = (
        1.00e+00, 9.77e-02, 0.00e+00, 0.00e+00, 0.00e+00, 2.40e-04, 7.59e-05, 4.90e-04, 0.00e+00, 8.71e-05,
        1.45e-06, 2.51e-05, 2.14e-06, 1.86e-05, 2.63e-07, 1.23e-05, 1.32e-07, 2.57e-06, 0.00e+00, 1.58e-06,
        0.00e+00, 6.46e-08, 0.00e+00, 3.24e-07, 2.19e-07, 2.69e-05, 8.32e-08, 1.12e-06, 0.00e+00, 0.00e+00,
    ),
    lodd = (
        1.00e+00, 7.92e-02, 1.90e-09, 2.57e-11, 6.03e-10, 2.45e-04, 6.76e-05, 4.90e-04, 2.88e-08, 7.41e-05,
        1.99e-06, 3.55e-05, 2.88e-06, 3.47e-05, 2.88e-07, 1.55e-05, 1.82e-07, 3.55e-06, 1.29e-07, 2.19e-06,
        1.17e-09, 8.32e-08, 1.00e-08, 4.47e-07, 3.16e-07, 2.95e-05, 8.13e-08, 1.66e-06, 1.82e-08, 4.27e-08,
    ),
    lpgp = (
        1.00e+00, 8.41e-02, 1.26e-11, 2.40e-11, 5.01e-10, 2.45e-04, 7.24e-05, 5.37e-04, 3.63e-08, 1.12e-04,
        2.00e-06, 3.47e-05, 2.95e-06, 3.31e-05, 2.88e-07, 1.38e-05, 3.16e-07, 3.16e-06, 1.32e-07, 2.14e-06,
        1.26e-09, 7.94e-08, 1.00e-08, 4.37e-07, 2.34e-07, 2.82e-05, 8.32e-08, 1.70e-06, 1.62e-08, 4.17e-08,
    ),
    lpgs = (
        1.00e+00, 9.69e-02, 2.15e-09, 2.36e-11, 7.26e-10, 2.78e-04, 8.19e-05, 6.06e-04, 3.10e-08, 1.27e-04,
        2.23e-06, 3.98e-05, 3.27e-06, 3.86e-05, 3.20e-07, 1.63e-05, 2.00e-07, 3.58e-06, 1.45e-07, 2.33e-06,
        1.33e-09, 9.54e-08, 1.11e-08, 5.06e-07, 3.56e-07, 3.27e-05, 9.07e-08, 1.89e-06, 2.09e-08, 5.02e-08,
    ),
    felc = (
        1.00e+00, 7.94e-02, 0.00e+00, 0.00e+00, 0.00e+00, 3.89e-04, 1.00e-04, 7.76e-04, 0.00e+00, 1.20e-04,
        8.51e-06, 1.41e-04, 1.10e-05, 1.26e-04, 0.00e+00, 1.86e-05, 0.00e+00, 3.80e-06, 0.00e+00, 8.51e-06,
        0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.26e-04, 0.00e+00, 6.92e-06, 0.00e+00, 0.00e+00,
    ),
)


"""
    solar_abundances(table = :angr)

Return XSPEC reference abundances as a named tuple from `H` through `Zn`.
Accepts a built-in table name or the path to an XSPEC-style text file with
30 numbers in atomic-number order. The values are linear number ratios to H.
"""
function solar_abundances(table::Symbol = :angr)
    haskey(XSPEC_ABUNDANCES, table) ||
        throw(ArgumentError("Unknown XSPEC abundance table $table. Choose one of $(keys(XSPEC_ABUNDANCES))"))
    return NamedTuple{ELEMENT_NAMES}(XSPEC_ABUNDANCES[table])
end

function solar_abundances(table::AbstractString)
    name = Symbol(table)
    haskey(XSPEC_ABUNDANCES, name) && return solar_abundances(name)
    isfile(table) || throw(ArgumentError("No abundance table or file named $table"))
    entries = split(read(table, String))
    length(entries) == 30 || throw(ArgumentError("An abundance file must contain 30 values from H through Zn"))
    values = parse.(Float64, entries)
    all(value -> isfinite(value) && value >= 0, values) ||
        throw(ArgumentError("Reference abundances must be finite and nonnegative"))
    return NamedTuple{ELEMENT_NAMES}(Tuple(values))
end

const APEC_METALS = (6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28)
const VAPEC_ELEMENTS = (2, APEC_METALS...)
const VAPEC_NAMES = map(z -> ELEMENT_NAMES[z], VAPEC_ELEMENTS)
const VAPEC_INDEX = ntuple(z -> something(findfirst(==(z), VAPEC_ELEMENTS), 0), 30)

struct ApecOptions
    ratios::NTuple{30,Float64}
    trace::Union{Float64,Int}
    thermal::Bool
end

function apec_options(abundance_table, trace_abundance, thermal_broadening::Bool)
    ratios = map(/, Tuple(solar_abundances(abundance_table)), XSPEC_ABUNDANCES.angr)
    return ApecOptions(ratios, trace_setting(trace_abundance), thermal_broadening)
end

function trace_setting(value::Real)
    isfinite(value) && value >= 0 || throw(ArgumentError("Trace abundance must be finite and nonnegative"))
    return Float64(value)
end

function trace_setting(value::Symbol)
    z = findfirst(==(value), ELEMENT_NAMES)
    (z == 1 || z in VAPEC_ELEMENTS) ||
        throw(ArgumentError("Trace abundance must refer to H, He, or a fitted metal, such as :Fe"))
    return z
end

function element_parameters(parameters, names)
    for name in keys(parameters)
        name in names || throw(ArgumentError("This APEC variant has no abundance parameter named $name"))
    end
    return map(name -> get(parameters, name, 1.0), names)
end

function tied_abundances(abundance, trace::Union{Float64,Int})
    trace_value = trace isa Int ? (trace in APEC_METALS ? abundance : 1.0) : trace
    return ntuple(30) do z
        z <= 2 && return 1.0
        return z in APEC_METALS ? abundance : trace_value
    end
end

function partial_abundances(values, trace::Union{Float64,Int})
    trace_value = trace isa Int ? (trace == 1 ? 1.0 : values[VAPEC_INDEX[trace]]) : trace
    return ntuple(30) do z
        z == 1 && return 1.0
        index = VAPEC_INDEX[z]
        return index == 0 ? trace_value : values[index]
    end
end
