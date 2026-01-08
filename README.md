# TERRA.jl – Julia Wrapper for the Multi-Temperature Collisional-Radiative Solver

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://amta3208.github.io/TERRA.jl/dev/)
[![Build Status](https://github.com/amta3208/TERRA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/amta3208/TERRA.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/amta3208/TERRA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/amta3208/TERRA.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

**TERRA (Toolkit for Excitation, Reactions, and Radiation with Applications)** is a general-purpose master-equation solver that can compute the evolution of a reacting ionized gas in thermodynamic and chemical nonequilibrium. It was developed in the Nonequilibrium Gas and Plasma Dynamics Laboratory at the University of Colorado Boulder and is currently maintained by Tim Aiken (timothy.aiken@colorado.edu). The native Fortran implementation can simulate either 0-D isochoric reactors or 1-D post-normal-shock flows. Users can choose between mode-level (multi-temperature) modeling and detailed state-to-state treatments for electronic and vibrational populations in arbitrary mixtures of atoms, molecules, ions, and electrons.

**TERRA.jl** is the Julia wrapper that interfaces with the TERRA Fortran API and is developed and maintained by Amin Taziny (amin.taziny@colorado.edu). It automates input generation, manages the Fortran runtime, exposes high-level Julia functions for integration and analysis, and optionally mirrors the native TERRA output layout to keep the established MATLAB/Tecplot tooling in play. Currently, only simulations of 0-D isochoric reactors are supported.

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

## Quick Start – Molcular Nitrogen Example

The wrapper ships with a helper that reproduces a canonical nitrogen test case. This example case runs an electronically-state resolved simulation of $N_2$ plasma in a 0-D isochoric and adiabatic reactor.
The initial conditions include an electron temperature of $T_e=10$ eV, electron mole fraction of $[e^-]=0.01$, and a total number density of $n_{tot}=10^{13}\;\text{cm}^{-3}$. This case can be ran simply via:

```julia
using TERRA

results = nitrogen_10ev_example()
@info "Final translational temperature" results.temperatures.tt[end]
```

To run the example in a specific case directory and mirror the native TERRA Tecplot outputs:

```julia
using TERRA

case_path = mktempdir()
config = nitrogen_10ev_config()

config = TERRAConfig(
    species = config.species,
    mole_fractions = config.mole_fractions,
    total_number_density = config.total_number_density,
    temperatures = config.temperatures,
    time_params = config.time_params,
    physics = config.physics,
    processes = config.processes,
    database_path = config.database_path,
    case_path = case_path,
    write_native_outputs = true,  # enable Tecplot-style outputs
)

initialize_terra(config)
results = solve_terra_0d(config)
finalize_terra()
```

The returned `TERRAResults` object contains the full time history of species densities, temperatures, and energy modes. Refer to the [package documentation](https://amta3208.github.io/TERRA.jl/stable/) for field descriptions and analysis utilities.

## Tools & MATLAB Post-Processing

The original TERRA toolchain is bundled under `tools/` and remains fully compatible when native-style outputs are enabled:

- `tools/postprocess.sh` – Shell wrapper that launches MATLAB and runs `postProcessTERRA.m`.
- `tools/matlab/` – MATLAB readers, plotting scripts, and the `DataViewerTERRA` app.
- `tools/run_terra.sh` – Convenience script for running the native executable (helpful for debugging outside of Julia).


Running the MATLAB post-processor requires the location of your MATLAB root directiry to be added to your system `$PATH` variable, where  `MATLAB_ROOT_DIR` can be found by executing `matlabroot` into the Command Window in a session of the MATLAB application. Once located, the path can be appended in the bash configuration file by executing the following command in your shell:

```bash
export PATH="<MATLAB_ROOT_DIR>/bin:$PATH"
```

The entire `tools/matlab` directory, with subfolders, must be added to the MATLAB path before it can be used as intended. Once the command line executable for `matlab` is available following the above path addition, the required path additions to use everything in `tools/matlab` can be accomplished by navigating to the `tools/matlab` directory and running the path configuration shell file.

```bash
bash configure_matlab_path.sh
```

Once completed, MATLAB-ready results can be generated from a Julia run by:

1. Set `write_native_outputs = true` in your `TERRAConfig` so the wrapper mirrors the native Tecplot files under `output/`.
2. Run your simulation via `integrate_0d_system`, `solve_terra_0d`, or a custom driver.
3. Execute the post-processing script, pointing it to the case directory:

   ```bash
   tools/postprocess.sh /path/to/case
   ```

   The script calls MATLAB in no-display mode and runs `postProcessTERRA.m`, which assembles a `.mat` summary of the case output for post-processing. This `.mat` file can then be opened within `tools/matlab/DataViewerTERRA.mlapp`, a provided MATLAB GUI application for quick plotting and post-processing of results.

The MATLAB tooling expects the native TERRA directory layout (`input/`, `output/`, `output/states/`, `output/sources/`). The wrapper guarantees this structure when native mirroring is enabled.
