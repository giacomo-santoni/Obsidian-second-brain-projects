---
tags: [nd-production, cafmaker, sandreco, debug, truth, edep-sim, spill-build]
date: 2026-11-26
---

# SAND mc mismatch in run-cafmaker

## Error

run-cafmaker (ND_CAFMaker v5.3.1) aborts on some files:

```
terminate called after throwing an instance of 'std::runtime_error'
  what():  SAND: the mc of SAND entry 543 differs from sr.mc at nu[54]
```

Seen in job `23290759.1.1` (file 0000001, entry 543) and `23298861.62.62` (entry 1102, `nu[47]`).

It is not a branch/main problem. v5.3.1 added a check in `SANDRecoBranchFiller.cxx` (lines 159-170): the `mc` stored in the SANDRECO file must match the `sr.mc` that CAFMaker builds itself from GHEP + `EDEPSIM_SPILLS`. This covers ids, counts, pdg, G4ID, parent, momenta, and `start_pos`/`end_pos` (relative tolerance 1e-6). The SAND truth links are indices into the SAND file's own `mc`, so they are only valid if the two match.

## Cause

**A bug in sandreco v1.2.1.** In rare cases it gets the order of the trajectory points wrong.

1. `EDEPTree::fill_points` copies each point into a `std::map<subdetector, vector<point>>`. TRACKER points are added twice, under DRIFT and under STT.
2. `EDEPTrajectory::GetTrajectoryPointsVect()` joins the lists in subdetector order, then tries to restore the step order with `std::stable_sort` on the **absolute time** `T`.
3. The spill build (`overlaySinglesIntoSpillsSortedWithNuIntTime.cpp`, lines 450-487) writes `T = time in spill + spillN × 1.2e9 ns`. At ~1e12 ns a `double` resolves only ~1e-4 ns, so two points a few µm apart have **exactly the same T**. The sort can't separate them and the subdetector order wins.
4. `truth_filler` takes `front()`/`back()` as `start_pos`/`end_pos`. These are swapped, or wrong, for a few particles. CAFMaker uses `Points.front()`/`back()` in the Geant4 order, so the check fails.

The two mismatches in file 0000001 (in all other particles, ids, parents, momenta and counts are identical):

| spill | ixn | particle | SAND | CAF |
|---|---|---|---|---|
| 543 | 13006162, `nu[54]` | e⁻ G4ID 579, ~4 keV, 2 points | start = pt1, end = pt0 (**swapped**) | start = pt0, end = pt1 |
| 1194 | 17009499, `nu[35]` | e⁺ G4ID 1194, 7 points | start = pt1 | start = pt0 |

In both cases `t(pt0) == t(pt1)` exactly (`651600007866.98547`, `1432800006465.4956`). The trajectory points are identical in `EDEPSIM_SPILLS` and `OVERLAY`, so `convert2edepsim_spill_format` is **not** the cause.

The same bug was in other places too. The first test of the patch crashed on a `fast_reco_test` assertion at spill 540 (`ECAL cluster E 11.411 MeV isn't the deposit 11.410 MeV`). The ECAL reco decided "born inside the ECAL" by comparing times, and the test used the (now correct) truth `start_pos`. In production the two agreed only because both were wrong in the same way.

## Fix (sandreco)

Point order must come from the **index**, never from the time.

**Data structure (`edep_reader`)**
- `EDEPTrajectory::points_`: every point **once**, in Geant4 step order. `GetTrajectoryPointsVect()` returns it, with no sort and no copy.
- `trajectory_points_` (per subdetector) and `first_points_` / `last_points_` now hold **indices** into `points_` (`EDEPPointIndices = std::vector<std::size_t>`) instead of copies. Indices rather than pointers because the `EDEPTrajectory` copy constructor is `= default`: a pointer would point into the original object.
- New helpers: `GetPoint(i)`, `GetFirstPointInDetector(c)`, `IsBornInDetector(c)` (true if point 0 is in `c`).
- `IsTrajectorySaturated()` now uses `points_.size()`. Before, it also counted the STT/DRIFT duplicates.

**Users updated**

| where | before | after |
|---|---|---|
| `common/truth_filler` (start/end) | sort on `T` | `points_.front()` / `back()` |
| `fast/ecal` `born_inside_ecal` | compare first-point `T` | `IsBornInDetector(ECAL)` |
| `fast/tracker` `entry_point` (STT vs DRIFT) | smaller `T` | smaller index |
| `common/gluckstern_smearing` `tracker_entry_point` | smaller `T` | smaller index |
| `fast/grain`, `fast/ecal` entry point | copy of the point | `GetFirstPointInDetector` |
| `test/fast_reco_test` `find_entry_point`, born-in-ECAL | smaller `T` | smaller index |

**Left as is:** the hit sort in `EDEPTree.cpp:200` (`stable_sort` on `Start.T`). Hits are not duplicated, and within one detector they already come in edep-sim order. Adding the same offset can make two different times equal but can't swap them, and on ties `stable_sort` keeps the original order.

**Why a relative time doesn't help:** subtracting an offset in sandreco can't recover precision that the spill build already threw away. It would only work if the spill build saved times relative to the spill start, which changes the file convention that CAFMaker/TMS also use. The step index is exact anyway.

## Test

Everything is in `ND_Production_main/tests/debug_sand_mc/`. Production was not touched.

- `ND_CAFMaker/`: CAFMaker v5.3.1 with a field-by-field diff dump instead of the `throw`. This is how the 2 mismatches were found.
- `sandreco/`: sandreco v1.2.1 (commit `85dff1d`) + the patch, built in `sand-prod-cpu_v1.2.1.sif`.
- `run_full_test.sh`: run detached (`setsid nohup`). It runs patched sandreco on production input 0000001 (same `config_json.py` and parameters as production), then the **original** production CAFMaker v5.3.1 with the check enabled. The result goes to `status.txt`.
- `dump_parts.C`, `dump_traj.C`: macros that print the mismatched particles and their edep-sim points.

Result (2026-11-26):

```
END sandreco patched, exit=0          # all 1498 spills, no fast_reco_test assertion
END CAFMaker, exit=0
RESULT: CAFMaker check PASSED on all spills
```

- The `geoinfo: No object found` warnings: 36005, the same count as the production run. Unrelated to the patch.
- Patched SANDRECO vs production SANDRECO: 1001167223 vs 1001167835 bytes. The CAF output has the same `cafTree` size and branch sizes as the debug run.

## Gotchas met along the way

- The job uses the CAFMaker in `ND_Production/run-cafmaker/install/ND_CAFMaker` (tag v5.3.1), **not** `ND_Production_main/ND_CAFMaker`.
- In the debug CAF on disk, `nu.sec` is pruned after the SAND filler, so look particles up by interaction `id` and use the in-memory values.
- edep-sim v3.2: use the getters (`GetTrackId()`, `GetPosition()`, …). The members are private.
- To run patched sandreco, `ufw-ldpath` must list **only** the patched `lib64`. With `/usr/local/lib64` as well, ufw registers `caf_streamer` twice and aborts.
- Rebuilding/installing while a run is using the same libraries makes that run unreliable. Stop it first.

## Next steps

- [ ] Propose the patch to sandreco (baltig MR), together with the note on the time precision after the spill build.
- [ ] Rebuild `sand-prod-cpu` with the fix, then rerun run-sandreco and run-cafmaker on the files that failed.
- [ ] Retest file 62 (job `23298861.62.62`, entry 1102). Probably the same cause, not checked yet.
- [ ] Discuss with ND_Production whether the spill build should save times relative to the spill start.

## Related
- [[ND Production]]
- [[sandreco]]
