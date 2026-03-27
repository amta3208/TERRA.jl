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
end

function Config(; reactor::ReactorConfig,
                numerics::NumericsConfig,
                sources::SourceTermsConfig = SourceTermsConfig(),
                models::ModelConfig = ModelConfig(),
                runtime::RuntimeConfig = RuntimeConfig())
    return Config(reactor, models, sources, numerics, runtime)
end
