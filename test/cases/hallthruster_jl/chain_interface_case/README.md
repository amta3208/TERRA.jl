# HallThruster Chain Interface Case

This fixture is a HallThruster-side source case for the TERRA chain-of-CSTR interface work.

Files:
1. `chain_interface_case.jl`
   Runs the HallThruster simulation and writes a HallThruster-native JSON output file.
2. `chain_interface_diagnostics.jl`
   Reuses the main case script to run the same simulation, print diagnostic summaries, and
   optionally display or save plots.
3. `output/`
   Default output directory for the HallThruster JSON artifact written by the main case.

Default generated source artifact:
1. `output/hallthruster_solution_avg.json`

This HallThruster JSON file is the source input for the later adapter script in
`tools/hallthruster_jl/export_chain_profile.jl`, which will reduce it to the TERRA-side
`chain_profile_v4.json` contract.

Example conversion command:
```bash
julia --project=tools/hallthruster_jl tools/hallthruster_jl/export_chain_profile.jl \
  output/hallthruster_solution_avg.json \
  /path/to/terra_case \
  --average_start_time=0.0005
```
