# Run separately with LibXSPEC_jll available and the exact AtomDB prefix as ARGS[1].
using Plasma, LibXSPEC_jll, SHA, Test

const XSFUNCTIONS = LibXSPEC_jll.libXSFunctions
const XSUTIL = LibXSPEC_jll.libXSUtil

function set_model_string(name, value)
    ccall((:FPMSTR, XSUTIL), Cvoid, (Cstring, Cstring), name, value)
    return nothing
end

function reference_apec(edges; kT, abundance, redshift, norm)
    parameters = Float64[kT, abundance, redshift]
    counts = length(edges) - 1
    flux, error = fill(NaN, counts), zeros(counts)
    ccall((:C_apec, XSFUNCTIONS), Cvoid,
        (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cstring),
        edges, counts, parameters, 1, flux, error, "")
    return flux .* norm
end

length(ARGS) == 1 || error("Usage: julia test/xspec_parity.jl /path/to/apec_v3.0.9")
prefix = abspath(only(ARGS))
db = AtomDB(prefix)
@info "XSPEC APEC reference" libxspec_version = Base.pkgversion(LibXSPEC_jll) temperature_nodes = length(db.temperatures)
for path in (db.linefile, db.cocofile)
    @info "AtomDB emissivities" path sha256 = bytes2hex(open(sha256, path))
end

ccall((:FNINIT, XSUTIL), Cvoid, ())
status = Ref{Cint}(-1)
ccall((:FPSOLR, XSUTIL), Cvoid, (Cstring, Ref{Cint}), "angr", status)
status[] == 0 || error("Could not select the angr abundance table in XSPEC")
set_model_string("APECROOT", prefix)
set_model_string("APEC_TRACE_ABUND", "1.0")
set_model_string("APECVELOCITY", "0")
for key in ("APECTHERMAL", "APECNOLINES", "APECLOGINTERP", "APECLOGTINTERP", "APECUSENEI",
        "APECEEBREMSS", "APECDOECOR", "APECUSESPEX")
    set_model_string(key, "no")
end

broad = collect(range(0.1, 10.0; length = 2001))
fine = reverse(12.39841974 ./ collect(range(1.7, 10.5; length = 4401)))
temperatures = db.temperatures
middle = searchsortedfirst(temperatures, 1.0)
node = temperatures[middle]
between = (temperatures[middle - 1] + node) / 2

@testset "Bin-by-bin agreement with XSPEC" begin
    for (gridname, edges) in (("broad", broad), ("grating", fine)),
            kT in (temperatures[1] / 2, temperatures[1], node, between, 5.0,
                   temperatures[end], 2 * temperatures[end]),
            (abundance, redshift, norm) in ((1.0, 0.0, 1.0), (0.0, 0.0, 0.03), (0.2, 0.04, 2.0))
        @testset "$gridname, kT=$kT, abundance=$abundance, redshift=$redshift, norm=$norm" begin
            reference = reference_apec(edges; kT, abundance, redshift, norm)
            native = apec(edges, db; kT, abundance, redshift, norm)
            @test all(x -> isfinite(x) && x >= 0, reference)
            @test sum(reference) > 0
            @test sum(abs, native .- reference) <= 1e-5 * sum(reference)
            @test all(isapprox.(native, reference; rtol = 1e-4, atol = 1e-9 * maximum(reference)))
        end
    end
end
