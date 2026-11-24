struct PlasmaConstants
    kev_angstrom::Float64
    kev_per_kelvin::Float64
    light_speed_km_s::Float64
    kev_per_atomic_mass::Float64
    flow_particles_per_second::Float64
end

# Solar mass, year, and mean molecular weight used by the cooling-flow models.
const SOLAR_MASS_G = 1.98847e33
const JULIAN_YEAR_S = 31_557_600
const COOLING_MEAN_MOLECULAR_WEIGHT = 0.6

flow_particles_per_second(atomic_mass_g) =
    (5 / 2) * (SOLAR_MASS_G / JULIAN_YEAR_S) / (COOLING_MEAN_MOLECULAR_WEIGHT * atomic_mass_g)

const LEGACY_CONSTANTS = PlasmaConstants(
    12.39841974,
    1.380649e-23 / 1.602176634e-16,
    299792.458,
    931494.10242,
    flow_particles_per_second(1.66053906660e-24),
)

const CODATA2022_CONSTANTS = let
    c = CODATA2022.SpeedOfLightInVacuum
    h = CODATA2022.PlanckConstant
    k = CODATA2022.BoltzmannConstant
    m_u = CODATA2022.AtomicMassConstant
    PlasmaConstants(
        ustrip(u"keV*angstrom", h * c),
        ustrip(u"keV/K", k),
        ustrip(u"km/s", c),
        ustrip(u"keV", m_u * c^2),
        flow_particles_per_second(Float64(ustrip(u"g", m_u))),
    )
end

const DEFAULT_CONSTANTS_PRESET = :legacy

function plasma_constants(preset::Symbol)
    preset === :legacy && return LEGACY_CONSTANTS
    preset === :codata2022 && return CODATA2022_CONSTANTS
    throw(ArgumentError("Unknown physical constants preset $preset. Choose :legacy or :codata2022."))
end
