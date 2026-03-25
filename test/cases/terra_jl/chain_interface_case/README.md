# TERRA Interface Fixture

Placeholder wrapper-side handoff target reserved for HallThruster coupling-interface tests.

For the current interface design, this fixture is not a full static TERRA case. Its purpose is to
hold the reduced JSON chain-profile artifact written by the HallThruster exporter.

For MVP, this fixture only needs:
1. `input/`
2. a generated `chain_profile_v4.json` artifact under `input/` when the exporter is refreshed

The Julia-side test code should create any runtime case/output directories separately when
exercising the actual TERRA chain solver.
