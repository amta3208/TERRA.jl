"""
# TERRA Configuration Module

This module handles configuration management for TERRA simulations,
including parameter validation, default values, and input file generation.
"""

include("config/types_and_adapters.jl")
include("config/validation_core.jl")
include("io/input_generation.jl")
include("config/examples_and_species.jl")
include("config/validation_terra.jl")
include("config/unit_conversions.jl")
