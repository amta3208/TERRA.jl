# -----------------------------------------------------------------------------
# Reactor CSV persistence
# -----------------------------------------------------------------------------

"""
$(SIGNATURES)

Save TERRA results to file.

# Arguments
- `results::ReactorResult`: Results to save
- `filename::String`: Output filename (CSV format)

# Returns
- `true` if save successful
"""
function save_results(results::ReactorResult, filename::String)
    try
        species_densities = species_density_matrix(results)
        temperatures = temperature_history(results)
        total_energy = total_energy_history(results)

        n_times = length(results.t)
        n_species = size(species_densities, 1)

        header = ["time", "total_energy", "T_trans", "T_electron", "T_vib"]
        for i in 1:n_species
            push!(header, "species_$(i)_density")
        end

        data = zeros(n_times, length(header))
        data[:, 1] = results.t
        data[:, 2] = total_energy
        data[:, 3] = temperatures.tt
        data[:, 4] = temperatures.te
        data[:, 5] = temperatures.tv

        for i in 1:n_species
            data[:, 5 + i] = species_densities[i, :]
        end

        open(filename, "w") do io
            println(io, join(header, ","))
            for row in eachrow(data)
                println(io, join(row, ","))
            end
        end

        @info "Results saved successfully" filename = filename
        return true
    catch e
        @error "Failed to save results" filename = filename exception = e
        return false
    end
end
