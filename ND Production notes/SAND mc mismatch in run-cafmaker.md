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

## Fix (sandreco, commit `5c4fb9a`)

Branch `hotfix/remove-sort-by-time-due-to-precision-error` on baltig, one commit on top of `85dff1d` (main): `5c4fb9a` "removed sort by time due to double precision error", 7 files, +89/−65.

Point order must come from the **Geant4 step order**, never from the time.

**1. Point order (`edep_reader`: `EDEPTree.cpp`, `EDEPTrajectory.h/.cpp`)**

| | `85dff1d` | `5c4fb9a` |
|---|---|---|
| `GetTrajectoryPointsVect()` | joins the per-subdetector lists, `stable_sort` on `T` | returns `ordered_points_` (new private member, filled in `fill_points`), Geant4 step order |
| STT+DRIFT points | twice | once |
| equal `T` | wrong/ambiguous order | correct |
| public interface | — | unchanged (same signature, still returns a copy) |

`trajectory_points_` and `first_points_` (per subdetector) are unchanged and still hold copies. `truth_filler` (`front()`/`back()` of `GetTrajectoryPointsVect()`) is correct now without changes.

**2. Tracker entry point** (`fast/tracker/reco.cpp`, `common/gluckstern_smearing/gluckstern_smearing.cpp`, test `find_entry_point`)

| `85dff1d` | `5c4fb9a` |
|---|---|
| first points of STT/DRIFT: smaller `T` wins | the fronts are candidates; the first one met along `GetTrajectoryPointsVect()` wins |
| test: smaller `T` over all points of the region | test: same logic as the reco |

**3. Born inside the ECAL** (`fast/ecal/reco.cpp` `born_inside_ecal`, `test_fast_reco.cpp`)

| | `85dff1d` | `5c4fb9a` |
|---|---|---|
| reco | ECAL first point not later (in `T`) than any other subdetector's | ECAL first point `==` first trajectory point |
| test | **any** ECAL point with `T == true.time` and float-cm position `== start_pos` | same criterion as the reco |

The old test criterion crashed run 92 (`fast_reco_test`, context 440): `ixn[25].part[0]: ECAL cluster var_position (14.6917, 22.4237, 18.4384) cm isn't the spread of its 1050 hits (14.6923, 22.4237, 18.4401) cm`.
- One hit differed between test and reco: 0.0021 MeV from an e⁻ (trj 9865, 5 generations below the π⁻ root) born just **outside** the ECAL.
- The reco correctly took it as not born in the ECAL. The test matched its **second** step, already in the ECAL, which had the same `T` (~5.28e11 ns) and the same float-cm position, so it added the hit to the π⁻ cluster.
- With 1050 hits, E and barycenter stayed within tolerance; the z spread did not (1.7e-3 cm vs ~1.2e-3 cm tolerance).

**4. ECAL cluster spread `var_position`** (reco and test, same formula)

| | `85dff1d` | `5c4fb9a` |
|---|---|---|
| sums | `e`, `e2 = Σde²` | `e`, `cross = Σ_{i<j} de_i de_j` (`cross += de * e` before `e += de`) |
| formula | `sqrt(neff/(neff−1) · Σw dx²/e)`, `neff = e²/e2` | `sqrt(e · Σw dx² / (2 cross))`, 0 if one hit or `cross <= 0` |
| result | — | identical (algebraic identity) |
| one dominant hit | `neff − 1` cancels, depends on the hit order | no subtraction, stable |

Not the cause of any crash: only a more stable rewrite of the same quantity.

**5. `var_position` assertion message (test):** now prints the reco values, the expected values and the number of hits.

**Left as is**
- The hit sort in `EDEPTree.cpp:205` (`stable_sort` on `Start.T`). Hits are not duplicated, and within one detector they already come in edep-sim order. Adding the same offset can make two different times equal but can't swap them, and on ties `stable_sort` keeps the original order.
- `IsTrajectorySaturated()` still counts the STT/DRIFT duplicates.

**Why a relative time doesn't help:** subtracting an offset in sandreco can't recover precision that the spill build already threw away. It would only work if the spill build saved times relative to the spill start, which changes the file convention that CAFMaker/TMS also use. The step order is exact anyway.

## Test

Everything is in `ND_Production_main/tests/debug_sand_mc/`.

- `ND_CAFMaker/`: CAFMaker v5.3.1 with a field-by-field diff dump instead of the `throw`. This is how the 2 mismatches were found.
- `run_full_test.sh`: runs patched sandreco on production input 0000001 (same `config_json.py` and parameters as production), then the **original** production CAFMaker v5.3.1 with the check enabled. The result goes to `status.txt`.
- `dump_parts.C`, `dump_traj.C`: macros that print the mismatched particles and their edep-sim points.

The 2026-11-26 run of `run_full_test.sh` (all 1498 spills, `CAFMaker check PASSED on all spills`) used an **earlier version** of the patch, not `5c4fb9a`.

Final version (`5c4fb9a`), 2026-11-27: production rerun of run-sandreco (cluster `23342671`) and run 92 (job `23340857`) complete with no `fast_reco_test` assertion. CAFMaker with the check on these outputs: not run yet.

## Other problems found on 2026-11-27 (not caused by the fix)

**Runs 0–3 (cluster `23342671`): output lost.** All ran to 100%, then `mv: cannot stat .../tmp/run-sandreco/...000000{0-3}.SANDRECO.root`.
- A second cluster, `23342805`, was submitted by mistake 2 min later and removed with `condor_rm` after ~7 s. In those seconds its `run_sandreco.sh` ran `rm -f "$SANDRECOfile"` on the **same temporary path** that the first jobs had just opened. ufwrun kept writing to the unlinked file, and at the end the file was gone.
- Run 4 survived: its duplicate started 1 s before the removal and never reached the `rm`.
- The files `0000000–0003` in `SANDRECO/0000000/` are from the previous night's production (old code).

**Run 5: stopped at context 1339 (`g4_id_map.cpp:83`, secondary track ids not contiguous).** A primary muon has no trajectory already in the edep-sim output. Moved to [[Missing edepsim trajectories for primaries]], together with the same problem in 5 horn_current_neg files.

## Gotchas met along the way

- The job uses the CAFMaker in `ND_Production/run-cafmaker/install/ND_CAFMaker` (tag v5.3.1), **not** `ND_Production_main/ND_CAFMaker`.
- In the debug CAF on disk, `nu.sec` is pruned after the SAND filler, so look particles up by interaction `id` and use the in-memory values.
- edep-sim v3.2: use the getters (`GetTrackId()`, `GetPosition()`, …). The members are private.
- To run patched sandreco, `ufw-ldpath` must list **only** the patched `lib64`. With `/usr/local/lib64` as well, ufw registers `caf_streamer` twice and aborts.
- Rebuilding/installing while a run is using the same libraries makes that run unreliable. Stop it first. Never use `install_sandreco.sh` with jobs running: it does `rm -rf build install_sandreco`; use `cmake --build build --target install`.
- The jobs load only `install_sandreco/lib64` (untracked), so `git checkout`/`stash`/`reset` in the source tree do not affect running jobs.

## Next steps

- [ ] Open the MR for `hotfix/remove-sort-by-time-due-to-precision-error` (`5c4fb9a`) into main, together with the note on the time precision after the spill build.
- [ ] Run CAFMaker v5.3.1 with the check on the `5c4fb9a` outputs.
- [ ] Rebuild `sand-prod-cpu` with the fix, then rerun run-sandreco and run-cafmaker on the files that failed.
- [ ] Rerun sandreco runs 0–3 (output lost to the duplicate submission).
- [ ] Retest file 62 (job `23298861.62.62`, entry 1102). Probably the same cause, not checked yet.
- [ ] Run 5 (missing edep-sim muon): see [[Missing edepsim trajectories for primaries]].
- [ ] Make the sandreco temporary output name unique per job.
- [ ] Discuss with ND_Production whether the spill build should save times relative to the spill start.

## Related
- [[ND Production]]
- [[Missing edepsim trajectories for primaries]]
- [[sandreco]]
