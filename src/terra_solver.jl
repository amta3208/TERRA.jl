"""
# TERRA Solver Module

This module provides the high-level interface for running TERRA simulations
from Julia, hiding the complexity of the Fortran interface and providing
a clean, Julia-native API.
"""

include("solver/api_layout.jl")
include("solver/initial_state.jl")
include("solver/state_vector.jl")
include("solver/residence_time.jl")
include("solver/rhs.jl")
include("solver/integrate_0d.jl")
include("solver/driver.jl")
