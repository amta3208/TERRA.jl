"""
$(SIGNATURES)

Top-level configuration for refactored TERRA workflows.
"""
struct Config
    reactor::ReactorConfig
    models::ModelConfig
    sources::SourceTermsConfig
    numerics::NumericsConfig
    runtime::RuntimeConfig

    function Config(; reactor::ReactorConfig,
                    numerics::NumericsConfig,
                    sources::SourceTermsConfig = SourceTermsConfig(),
                    models::ModelConfig = ModelConfig(),
                    runtime::RuntimeConfig = RuntimeConfig())
        return new(reactor, models, sources, numerics, runtime)
    end
end

_coerce_residence_time_inlet_reactor(inlet::Config) = inlet.reactor
