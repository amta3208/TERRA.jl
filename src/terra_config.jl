"""
# TERRA Configuration Module

This module handles configuration management for TERRA simulations,
including parameter validation, default values, and input file generation.
"""

include("config/types.jl")
include("config/validation.jl")
include("io/prob_setup_writer.jl")
include("io/sources_setup_writer.jl")
include("io/tau_scaling_writer.jl")
include("io/input_generation.jl")
include("config/conversions.jl")
