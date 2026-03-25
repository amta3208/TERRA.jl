# HallThruster.jl Tools

Standalone HallThruster.jl to TERRA.jl coupling scripts live here.

Current tool:
1. `export_chain_profile.jl`
   Converts a HallThruster JSON output file with an existing `output.average` block into the
   reduced TERRA chain-profile artifact at `<terra_case_path>/input/chain_profile_v4.json`.

Recommended invocation:
```bash
julia --project=tools/hallthruster_jl tools/hallthruster_jl/export_chain_profile.jl \
  <hallthruster_input.json> \
  <terra_case_path> \
  --average_start_time=<seconds> \
  [--ion_velocity_policy=<trim_to_first_positive|strict_positive>] \
  [--u_ion_floor=<float>] \
  [--min_consecutive_positive=<int>]
```

Notes:
1. Run `julia --project=tools/hallthruster_jl -e 'using Pkg; Pkg.instantiate()'` once before first use.
2. The HallThruster source file must already contain `output.average`.
3. For JSON inputs, `--average_start_time` is recorded as metadata only; the script does not recompute the average.
4. Default ion handling is `trim_to_first_positive` with `u_ion_floor=0.0` and `min_consecutive_positive=3`.
5. The exporter trims leading axial points until every exported ion species has reached a sustained positive-velocity region, then enforces the positivity contract on the exported profile.
6. All non-electron neutral and ion species present in HallThruster are exported automatically.
