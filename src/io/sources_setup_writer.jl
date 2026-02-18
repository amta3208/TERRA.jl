"""
$(SIGNATURES)

Generate sources_setup.inp file from configuration.
"""
function generate_sources_setup_file(config::TERRAConfig, filepath::String)
    return generate_sources_setup_file(to_config(config), filepath)
end

function generate_sources_setup_file(config::Config, filepath::String)
    species_names = config.reactor.composition.species

    open(filepath, "w") do io
        println(io, "BEGIN SPECIES SOURCES")
        for species in species_names
            println(io, species)
        end
        println(io, "END SPECIES SOURCES")
        println(io)
        println(io, "BEGIN EXCITED STATE SOURCES")
        # For now, include common excited states for nitrogen
        # This should be made more general based on the species
        if "N" in species_names
            for i in 1:5
                println(io, "N($(i))")
            end
        end
        if "N2" in species_names
            for i in 1:9
                println(io, "N2($(i))")
            end
        end
        if "N2+" in species_names
            for i in 1:4
                println(io, "N2+($(i))")
            end
        end
        println(io, "END EXCITED STATE SOURCES")
    end
end
