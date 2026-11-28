---
tags: [nd-production, sandreco, edep-sim, genie, debug, truth]
date: 2026-11-27
---

# Missing edepsim trajectories for primaries

## Error

run-sandreco stops partway through some files and produces no output:

```
[error] [g4_id_map.cpp:83] Interaction 74818: secondary track ids [754, 23957] are not contiguous (384 ids)
[warning] [runner.cpp:244] Caught exception 'ufw::exception' in context 1339: '...'
[error] [runner.cpp:249] Exception limit exceeded (1 out of 0 allowed)
terminate called after throwing an instance of 'ufw::exception'
```

Seen on 2026-11-27 in 6 files:

| production | run | context (progress) | interaction | message |
|---|---|---|---|---|
| nominal | 5 | 1339 (89%) | 74818 | secondary ids [754, 23957], 384 ids |
| horn_current_neg | 0 | 969 (64%) | 54854 | secondary ids [154, 2103], 79 ids |
| horn_current_neg | 4 | 164 (11%) | 9320 | secondary ids [221, 9472], 2625 ids |
| horn_current_neg | 41 | 844 (56%) | 47993 | secondary ids [136, 3605], 1121 ids |
| horn_current_neg | 49 | 120 (8%) | 6991 | secondary ids [684, 19821], 462 ids |
| horn_current_neg | 58 | 785 (52%) | 45226 | secondary ids [17, 1066], 104 ids |

In horn_current_neg these are the 5 missing files of `SANDRECO/0000000/` (95 of 100 present).

## Why sandreco stops

`make_g4_id_map` (`src/data/common/spill_maps/g4_id_map.cpp:64-87`) requires the track ids of each interaction to form two contiguous ranges: one for the primaries and one for the secondaries.

`convert2edepsim_spill_format` renumbers the ids spill-wide. The first N ids of each event go to the primary block, where N is the **number of vertex particles**; the rest go to the secondary block.

If a vertex particle has no trajectory, Geant4 never assigned it a track id, so all ids shift down. The first secondaries then land in the primary block, and the secondary ids of the interaction are no longer contiguous. Example, nominal run 5: the secondaries are {754} ∪ [23575, 23957]. 754 is a γ from a π⁰ decay that took the place of the missing muon.

sandreco reads the files correctly; the error message describes their content exactly. Its only fault is that it stops the whole file (exception limit 0) instead of skipping the interaction or the spill.

## Cause: traced back to edep-sim and GENIE

For each case the event was followed through every stage: OVERLAY → `EDEPSIM_SPILLS` (entry = interaction number) → hadd file (same entry) → original edep-sim file (hadd N = edep-sim files 10N…10N+9) → `DetSimPassThru/gRooTracker` in the same file. The vertex is identical at every stage, so spill build, hadd and conversion drop nothing.

| production / run | edep-sim file, entry (EventId) | reaction | vertex volume | roots / vertex particles | missing |
|---|---|---|---|---|---|
| nominal / 5 | `0000058`, 7771 (9093) | νμ CC DIS, Fe | `KLOEYokeEndcapDR` (iron) | 6 / 7 | **μ⁻** 2372 MeV (index 0) |
| hcn / 4 | `0000041`, 717 (832) | νμ CC DIS, Fe | `KLOEYokeBarrel` (iron) | 13 / 14 | **μ⁻** 8357 MeV (index 0) |
| hcn / 41 | `0000415`, 4889 (5603) | ν̄μ CC DIS charm, Pb | ECAL passive slab (lead) | 31 / 32 | **μ⁺** 7273 MeV (index 0) |
| hcn / 49 | `0000490`, 6991 (8027) | νμ CC DIS, Fe | `KLOEYokeBarrel` (iron) | 6 / 13 | **μ⁻** 3269 MeV + 6 slow protons (indices 0–6) |
| hcn / 58 | `0000585`, 1713 (1955) | νμ NC DIS, Fe | `KLOEYokeBarrel` (iron) | 2 / 6 | **ν, π⁻, π⁺, π⁰** (indices 0–3) |
| hcn / 0 | `0000006`, 3607 (4045) | νμ CC DIS charm, Fe | `KLOEYokeBarrel` (iron) | 2 / 3 | **diquark 2103** (index 2) |

(hcn = horn_current_neg)

In every case, in the original edep-sim file:
- trajectory ids are **contiguous** from 0 (e.g. 0–389 for nominal run 5): the missing particles never got a Geant4 track id;
- no hit segment has a contributor without a trajectory: they never deposited energy;
- in GENIE (`gRooTracker`) the missing particles are present as final-state particles (status 1).

### Two different causes

**1. Geant4/edep-sim drops the first vertex particles** (nominal run 5; hcn runs 4, 41, 49, 58).
- In GENIE these are normal final-state particles: muons, pions, protons, a neutrino. Example, nominal run 5: μ⁻ with E 2371.95 MeV, p (−486, −404, 2284) MeV, mass 105.66 MeV.
- The missing ones are always **the first N particles of the vertex list** (N = 1, 1, 1, 7, 4); the following ones are tracked normally.
- The vertex is always in **dense passive material**: the iron of the KLOE yoke, or the lead of the ECAL.
- Why Geant4 skips them is **not understood yet**.

**2. GENIE outputs a diquark as a final-state particle** (hcn run 0).
- `2103` (ud₁ diquark) with status 1 in a DIS charm event. A diquark should never be a final-state particle: this is a GENIE output problem.
- Geant4 cannot track a diquark, so it correctly has no trajectory.
- Here the missing particle is the last of the list, not the first.
- The energy differs: 999.12 MeV in GENIE, 635.09 MeV in the edep-sim vertex. edep-sim probably recomputes E from the momentum with a different mass for 2103 (not checked).

## Frequency

Scan of the **whole nominal edep-sim production** (999 files, 8 411 605 events), counting vertex particles vs root trajectories per event: **only one event** (nominal run 5). So about 1 in 8.4M, and only one sandreco file affected.

horn_current_neg: at least 5 of 100 files are affected. Each job stops at the first bad interaction, so a file may contain more. Not scanned yet.

Scan script: PyROOT on `EDepSimEvents`, only `Primaries*` and `Trajectories.ParentId` enabled; ~12 s per file, 8 in parallel.

## Tools and gotchas

- edep-sim v3.2: use the getters (`GetTrackId()`, `GetParentId()`, `GetPosition()`, …). The members are private. Load with `ROOT.gSystem.Load("libedepsim_io")`; the headers are not in the `sand-dev-cpu` container, so compiled ROOT macros fail.
- PyROOT crashes reading `gRooTracker` (2D arrays like `StdHepP4[StdHepN][4]`). Use uproot for it.
- The vertex volume: `EDepSimGeometry` in the same file, `FindNode(x, y, z)` + `GetPath()`.

## Next steps

- [ ] Decide how sandreco handles these interactions: skip the spill (`ufw-ignore-errors` in the ufw config?) or make `g4_id_map` tolerant (warn and skip the interaction). Then rerun nominal run 5 and hcn runs 0, 4, 41, 49, 58.
- [ ] Scan the whole horn_current_neg edep-sim production (and the other varied productions).
- [ ] Find out why Geant4 drops the first vertex particles: rerun edep-sim on one event (e.g. nominal EventId 9093, file 58) with verbose Geant4 output.
- [ ] Report the diquark in the final state to the GENIE side (hcn run 0, EventId 4045).
- [ ] `convert2edepsim_spill_format`: size the primary block on the number of root trajectories, not of vertex particles?

## Related
- [[ND Production]]
- [[SAND mc mismatch in run-cafmaker]]
- [[sandreco]]
