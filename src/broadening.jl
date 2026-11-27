# XSPEC's element masses in atomic mass units, ordered by atomic number.
const ATOMIC_MASSES = (
    1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994,
    18.9984032, 20.1797, 22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
    32.065, 35.4527, 39.948, 39.0983, 40.078, 44.955910, 47.867, 50.9415,
    51.9961, 54.938049, 55.845, 58.933200, 58.6934, 63.456, 65.39,
)
const GAUSSIAN_SIGMAS = 8

function line_width_fractions(kT, thermal, velocity, constants::PlasmaConstants)
    turbulent = (velocity / constants.light_speed_km_s)^2
    # The thermal variance uses the temperature of each emissivity slice.
    thermal_scale = thermal ? kT / constants.kev_per_atomic_mass : 0.0
    return map(mass -> sqrt(thermal_scale / mass + turbulent), ATOMIC_MASSES)
end

function add_broadened_lines!(output, edges, lines, coefficients, fractions, columns)
    energies = lines.energies
    radius_fraction = GAUSSIAN_SIGMAS * maximum(fractions)
    first_line = searchsortedfirst(energies, first(edges) / (1 + radius_fraction))
    last_line = radius_fraction < 1 ? searchsortedlast(energies, last(edges) / (1 - radius_fraction)) :
        length(energies)
    @inbounds for j in first_line:last_line
        center = energies[j]
        z = lines.elements[j]
        sigma = center * fractions[z]
        column = spectrum_column_index(columns, z)
        strength = coefficients[z] * lines.emissivities[j]
        if iszero(sigma)
            if first(edges) <= center <= last(edges)
                i = max(1, searchsortedfirst(edges, center) - 1)
                output[i, column] += strength
            end
            continue
        end
        radius = GAUSSIAN_SIGMAS * sigma
        first_bin = max(1, searchsortedlast(edges, center - radius))
        last_bin = min(size(output, 1), searchsortedfirst(edges, center + radius) - 1)
        first_bin > last_bin && continue
        scale = inv(sqrt(2) * sigma)
        previous = erf((edges[first_bin] - center) * scale)
        for i in first_bin:last_bin
            current = erf((edges[i + 1] - center) * scale)
            output[i, column] += (current - previous) / 2 * strength
            previous = current
        end
    end
    return nothing
end
