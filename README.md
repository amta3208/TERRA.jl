# TERRA.jl – Julia Wrapper for the Multi-Temperature Collisional-Radiative Solver

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://amta3208.github.io/TERRA.jl/dev/)
[![Build Status](https://github.com/amta3208/TERRA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/amta3208/TERRA.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/amta3208/TERRA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/amta3208/TERRA.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

**TERRA (Toolkit for Excitation, Reactions, and Radiation with Applications)** is a general-purpose master-equation solver that can compute the evolution of a reacting ionized gas in thermodynamic and chemical nonequilibrium. It was developed in the Nonequilibrium Gas and Plasma Dynamics Laboratory at the University of Colorado Boulder and is currently maintained by Tim Aiken (timothy.aiken@colorado.edu). The native Fortran implementation can simulate either 0-D isochoric reactors or 1-D post-normal-shock flows. Users can choose between mode-level (multi-temperature) modeling and detailed state-to-state treatments for electronic and vibrational populations in arbitrary mixtures of atoms, molecules, ions, and electrons.

**TERRA.jl** is the Julia wrapper that interfaces with the TERRA Fortran API and is developed and maintained by Amin Taziny (amin.taziny@colorado.edu). It automates input generation, manages the Fortran runtime, exposes high-level Julia functions for integration and analysis, and optionally mirrors the native TERRA output layout to keep the established MATLAB/Tecplot tooling in play.
Within its current 0-D scope, it supports both closed-reactor integrations and open-flow 0-D reactors through a continuously stirred tank reactor (CSTR) model via an integrated convective term.

## Requirements

Before using the wrapper you must either build the native TERRA API library if you have access to the Fortran source code or reach out to Amin or Tim for a copy of the shared library file appropriate for your machine. At run time the Julia package expects:

- A compiled shared library (`libterra.so`, `libterra.dylib`, or `terra.dll`).
- An MPI runtime compatible with the library build.

Set the following environment variables prior to launching Julia:

| Variable | Purpose |
| --- | --- |
| `TERRA_LIB_PATH` | Absolute path to the TERRA shared library. |
| `JULIA_MPI_BINARY` *(optional)* | Override if you need a specific MPI binary; otherwise `MPI.jl` defaults are used. |

Install the wrapper in your Julia environment:

```julia
import Pkg
Pkg.add(url="https://github.com/amta3208/TERRA.jl")
```

## Quick Start – Molecular Nitrogen Example

The wrapper ships with a helper that reproduces a canonical nitrogen test case. This case runs an electronically-state-resolved simulation of $N_2$ plasma in a 0-D isochoric and adiabatic reactor with electron temperature $T_e=10$ eV, electron mole fraction $[e^-]=0.01$, and total number density $n_{tot}=10^{13}\;\text{cm}^{-3}$.

The example below uses the refactored nested config architecture with only required constructor inputs, relying on defaults for all optional settings:

```julia
using TERRA

composition = ReactorComposition(;
    species = ["N", "N2", "N+", "N2+", "E-"],
    mole_fractions = [1.0e-20, 0.9998, 1.0e-20, 0.0001, 0.0001],
    total_number_density = 1.0e13,  # 1/cm^3
)

thermal = ReactorThermalState(; Tt = 750.0, Tv = 750.0, Tee = 750.0, Te = 115000.0)
reactor = ReactorConfig(; composition = composition, thermal = thermal)

# Time inputs are in seconds in the Julia wrapper
time = TimeConfig(; dt = 5e-12, dt_output = 5e-6, duration = 1e-3)
numerics = NumericsConfig(; time = time)

config = Config(; reactor = reactor, numerics = numerics)
config = with_case_path(config, mktempdir())

initialize_terra(config)
try
    # Use closed-reactor behavior for this run
    results = solve_terra_0d(config; use_residence_time = false)
    @info "Final translational temperature (K)" results.temperatures.tt[end]
finally
    finalize_terra()
end
```

For convenience, this same case is also available out of the box:

```julia
using TERRA

results = nitrogen_10ev_example()
@info "Final translational temperature (K)" results.temperatures.tt[end]
```

The returned `SimulationResult` object contains the full time history of species densities, temperatures, and energy modes. Refer to the [package documentation](https://amta3208.github.io/TERRA.jl/stable/) for field descriptions and analysis utilities.

## Tools & MATLAB Post-Processing

The original TERRA toolchain is bundled under `tools/` and remains fully compatible when native-style outputs are enabled:

- `tools/postprocess.sh` – Shell wrapper that launches MATLAB and runs `postProcessTERRA.m`.
- `tools/matlab/` – MATLAB readers, plotting scripts, and the `DataViewerTERRA` app.
- `tools/run_terra.sh` – Convenience script for running the native executable (helpful for debugging outside of Julia).


Running the MATLAB post-processor requires adding your MATLAB root directory to your system `$PATH` variable. `MATLAB_ROOT_DIR` can be found by running `matlabroot` in the MATLAB Command Window. Once located, append it in your shell configuration:

```bash
export PATH="<MATLAB_ROOT_DIR>/bin:$PATH"
```

The entire `tools/matlab` directory, with subfolders, must be added to the MATLAB path before it can be used as intended. Once the command line executable for `matlab` is available following the above path addition, the required path additions to use everything in `tools/matlab` can be accomplished by navigating to the `tools/matlab` directory and running the path configuration shell file.

```bash
bash configure_matlab_path.sh
```

Once completed, MATLAB-ready results can be generated from a Julia run by:

1. Set native mirroring in runtime settings, for example with `config = with_runtime(config; write_native_outputs = true)`.
2. Run your simulation via `solve_terra_0d` (or `nitrogen_10ev_example` for the packaged reference case).
3. Execute the post-processing script, pointing it to the case directory:

   ```bash
   tools/postprocess.sh /path/to/case
   ```

   The script calls MATLAB in no-display mode and runs `postProcessTERRA.m`, which assembles a `.mat` summary of the case output for post-processing. This `.mat` file can then be opened within `tools/matlab/DataViewerTERRA.mlapp`, a provided MATLAB GUI application for quick plotting and post-processing of results.

The MATLAB tooling expects the native TERRA directory layout (`input/`, `output/`, `output/states/`, `output/sources/`). The wrapper guarantees this structure when native mirroring is enabled.
