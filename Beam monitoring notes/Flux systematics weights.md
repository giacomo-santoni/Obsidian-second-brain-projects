---
tags: [beam-monitoring, sand, flux, systematics, duneanafluxtools, cafana]
date: 2026-11-23
---

# Flux systematics weights (duneanafluxtools)

How the flux systematic weights in the CAFs are computed, the bug that gave the SAND files the wrong weights, the fix, and how the reweighted HornCurrent −1σ compares with the varied-flux sample.

## TL;DR

- The weights are a **table lookup** (no interpolation): `w = 1 + σ · Δ(param, ν config, off-axis bin, E_true bin)`, where Δ is the fractional flux change at +1σ.
- In the SAND CAFs `mc.nu.baseline = NaN` and `beam.hornI = 0`. The old `caffluxweighter` therefore put every event in the **FD + RHC (antineutrino-mode)** tables, so the old weights were wrong. For numu, HornCurrent was ≈ 0 effect, which is why the old ratio plot was flat.
- Fix: new options `--nd --fhc` in `caffluxweighter` (branch `3-add-support-for-flatcafs-as-input-to-caffluxweighter`). Verified: **0 mismatches** over 1.1M weights.
- The HornCurrent −1σ effect is small (≤ ~3% at 3.5–5 GeV, −0.9% integrated). With 1 file it sits below the statistical error, and the varied sample shows a +1.6 ± 0.6% normalisation excess that reweighting can't produce. See [[#Reweighted −1σ vs varied sample]].

## How the weights are computed

Code: `duneanafluxtools/src/duneanafluxtools/FluxWeighter.cxx` and `CAFInterface.h`.

### Input tables

`duneanafluxtools/inputs/flux_variations_FD_and_PRISM_2023.root`. At runtime it is loaded from `$duneanafluxtools_ROOT/share/duneanafluxtools/inputs/`.

| Block | Path | Content |
|---|---|---|
| Focussing (45 params: HornCurrent, TargetDensity, …) | `FluxParameters/Focussing/<param>/` | `FD_<mode>_<sp>`: TH1 in E_ν. `ND_<mode>_<sp>/ND_<mode>_<sp>_<i>`: one TH1 per off-axis slice i (40 slices). `OffAxisTAxis`: edges −4, −3, …, 36 m |
| Hadron production (21 PCA comps) | `FluxParameters/HadronProduction/pca_N/` | Same structure, but each `ND_*` directory has its own `OffAxisTAxis` (−3 … 39 m, 42 slices); ND E range only 0–20 GeV |

- `<mode>` = `nu` (FHC) or `nubar` (RHC). `<sp>` = numu, nue, numubar, nuebar. There are also `_specrun_` / `_280kA_` variants, which the CAF tools never use.
- The energy binning **differs between off-axis slices**. For example `ND_nu_numu_4` has a 4.50–5.75 GeV bin, `_6` has 4.50–5.25.

### Inputs taken from each CAF interaction

| Input | Source |
|---|---|
| flavour | `nu.pdgorig` |
| energy | `nu.E` (true, GeV) |
| ND or FD | `nu.baseline < 1000` m ⚠️ NaN gives FD |
| FHC or RHC | `SR->beam.hornI > 0` ⚠️ 0 gives RHC |
| off-axis position | `nu.vtx.x * 0.01` (true vertex x, cm → m), ND only |

### Off-axis treatment (ND)

- The code uses `|vtx.x|` (`fabs`): only x, the sign is dropped, y is ignored. It assumes `vtx.x = 0` lies on the beam axis.
- Slice used: focussing `ND_..._<4 + floor(|x|)>`, hadron production `ND_..._<3 + floor(|x|)>`. The negative-x slices are never reached.
- `|x| ≥ 36 m` (focussing) or `≥ 39 m` (hadron production), or E outside the histogram range, gives an invalid bin and **w = 1**.
- No interpolation: the weight is a step function in 1 m × E-bin cells.
- In the SAND files the true vertices span |x| ≤ 3 m, roughly uniformly (43% of truth interactions have 2 ≤ |x| < 3 m), so slices `_4`, `_5`, `_6` are used. These do differ: for HornCurrent numu, the [0,1) m and [2,3) m slices differ by up to 0.027 on a shift of ≤ 0.045.

### Formula

- ND: `w = 1 + σ · content` (not clamped). FD: `w = max(0, 1 + σ · content)`.
- `caffluxweighter` stores `nu.flux_systs[idx].weights = {1, w(+1σ)}`, with header `vals = {0, 1}`.
- Header ids: focussing `3000 + i` at index `i` (0–44); hadron production `3500 + i` at index `45 + i` (45–65).
- In CAFAna (`FluxDialSyst::Shift`): `weight *= w[0] + σ·(w[1] − w[0])`, so −1σ gives `2 − w[1]`.

## The bug: SAND events got FD + RHC weights

File checked: `/pnfs/dune/persistent/users/gsantoni/beam_monitoring/data/fast_reco/nominal/sand_events_nominal_FHC.caffluxsyst.sand.full.sandstt1.0000001.CAF.flat.root` (the old one).

- `rec.mc.nu.baseline` = **NaN** for all 84,130 interactions. `NaN < 1000` is false, so `isND = false`.
- `rec.beam.hornI` = **0** for all 1,498 spills, so `is_neutrino_mode = false`.
- Every stored weight (8,335 of 8,335 with w ≠ 1, first 150 spills) matched `FluxParameters/Focussing/HornCurrent/FD_nubar_<sp>` bin by bin. None matched any ND table.
- Effect: the FD RHC table for numu HornCurrent is ≈ 0, so the ±1σ ratio was flat. The correct ND FHC table gives about +3 to +4.5% at 3–4.5 GeV.

Long-term fix: fill `baseline` and `beam.hornI` in the SAND CAFMaker.

## The fix: `--nd/--fd`, `--fhc/--rhc` in caffluxweighter

Repo `duneanafluxtools`, branch `3-add-support-for-flatcafs-as-input-to-caffluxweighter` (based on v0.9.1; writes both the structured and the `.flat.root` CAF). Not committed yet. The first attempt on `main` is in `git stash@{0}` ("add options for ND weights").

- `CAFInterface.h`: new overload `ana::GetFluxWeights(nu, is_neutrino_mode, isND)`. The old 2-argument version is unchanged.
- `app/caffluxweighter.cxx`: options `--nd|--fd`, `--fhc|--rhc` (default: old inference), unknown options rejected, an end-of-job summary (how many interactions got ND or FD tables), and warnings when `hornI <= 0` or `baseline = NaN` decided the tables.

### Build (needs the SL7 container)

```bash
apptainer exec -B /cvmfs,/exp,/nashome,/tmp,/pnfs /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest bash
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup cmake v3_22_2; setup gcc v12_1_0
setup root v6_28_12 -q e26:p3915:prof; setup duneanaobj v04_02_01 -q e26:prof
export CPATH=$SRPROXY_INC:$CPATH     # duneanaobj's CMake target doesn't export the SRProxy include path
cd duneanafluxtools/build && make caffluxweighter
```

- The cvmfs cmake needs SL7 libraries (`libssl.so.10`), so it doesn't work directly on EL9.
- On `main` the `--rntuple` code needs ROOT ≥ 6.34. It was guarded by a version check there; not relevant on this branch.

### Run

```bash
export duneanafluxtools_ROOT=<dir containing share/duneanafluxtools/inputs/flux_variations_FD_and_PRISM_2023.root>
caffluxweighter --nd --fhc input.CAF.root output.CAF.root   # also writes output.CAF.flat.root
```

### Verification

- New nominal file: `/pnfs/dune/persistent/users/gsantoni/beam_monitoring/sand_events_nominal_FHC.caffluxsyst.sand.full.sandstt1.0000001.CAF.flat.root`.
- All weights recomputed with the FluxWeighter library (ND + FHC): 16,793 ν × 66 params = 1,108,338 weights, **0 mismatches**, `w[0] = 1` everywhere.
- HornCurrent weights match `ND_nu_<sp>_<4+floor(|x|)>` at the expected slice.
- `globalTree` header: 66 params, all indices consistent with ids.

Sample of the per-neutrino comparison (HornCurrent, idx 13):

| pdg | E [GeV] | x [m] | table histogram | E bin | content | w[1] CAF | −1σ |
|---|---|---|---|---|---|---|---|
| 14 | 3.412 | 0.87 | `ND_nu_numu_4` | 3.25–3.50 | +0.02164 | 1.02164 | 0.97836 |
| 14 | 3.296 | 2.64 | `ND_nu_numu_6` | 3.25–3.50 | +0.01541 | 1.01541 | 0.98459 |
| 14 | 4.574 | 0.04 | `ND_nu_numu_4` | 4.50–5.75 | +0.02929 | 1.02929 | 0.97071 |
| 14 | 18.582 | −1.19 | `ND_nu_numu_5` | 17.25–119.5 | 0 | 1.00000 | 1.00000 |

## Analysis macro (`lblpwgtools/bm_analysis/my_analysis_syst.C`)

- HornCurrent index = 13 (id 3013 in `flux_systs.json`). Correct.
- `Multiverse::Hypercross(systs, 1)` universes: 0 = nominal, 1 = **−1σ**, 2 = +1σ. The macro plots universe 1, so it is the right one to compare with `horn_current_neg`.
- Lower pad changed to points only: universe/nominal as markers with no error bars (same events, so fully correlated), varied/nominal as markers with vertical stat errors (`E1 X0`).
- Open points:
  - `kNuMuCC_SAND_FV` is defined but **not applied** to the loaders.
  - `FluxParamIndex` (`id − 3000`) is only right for focussing params. For PCA it should be `id − 3500 + 45`, or better, look the index up by name in `globalTree`.
  - The legend `Form("Universe %u", i, universe_sigma[i])` ignores the sigma.
  - `ensemble` and the X/Z ensembles are filled but never used.

## Reweighted −1σ vs varied sample

numu CC, true E, no cuts, file `0000001` of each. The varied sample is scaled by POT (7.49e16 vs 7.495e16).

| E_true [GeV] | N nom | −1σ reweight / nom | varied / nom ± stat |
|---|---|---|---|
| 0–1 | 1122 | 1.002 | 1.061 ± 0.044 |
| 1–1.5 | 3165 | 1.000 | 1.051 ± 0.026 |
| 1.5–2 | 5641 | 1.002 | 0.985 ± 0.019 |
| 2–2.5 | 7450 | 0.999 | 1.008 ± 0.017 |
| 2.5–3 | 7755 | 0.991 | 1.017 ± 0.016 |
| 3–3.5 | 6283 | 0.983 | 1.019 ± 0.018 |
| 3.5–4 | 4371 | 0.976 | 0.984 ± 0.021 |
| 4–4.5 | 2718 | 0.972 | 1.055 ± 0.028 |
| 4.5–5 | 1825 | 0.973 | 1.001 ± 0.033 |
| 5–6 | 2298 | 0.979 | 0.969 ± 0.029 |
| 6–8 | 2394 | 0.991 | 1.090 ± 0.031 |
| 8–12 | 3913 | 0.996 | 1.012 ± 0.023 |
| 12–20 | 4528 | 0.997 | 1.032 ± 0.022 |
| **total** | 53463 | **0.991** | **1.016 ± 0.006** |

Interpretation:
- The expected −1σ effect (≤ 3%) is comparable to or below the stat error of one file (1.6–4.4% per bin). It is diluted further at reco level and without cuts.
- The varied sample is +1.6 ± 0.6% in total, the opposite direction to the expected −0.9%, with +9% at 6–8 GeV. χ² vs the reweighted prediction is ≈ 35 / 13 bins. It does not show the expected dip at 3.5–5 GeV, which a lower horn current should give.
- Possible causes to check:
  - the varied flux comes from an independent beam simulation with its own stat fluctuations;
  - a different definition of "1σ" (how many kA?) than the tables;
  - a different beam configuration from the "PRISM 2023" tables.

### Full statistics, reco level (2026-11-24)

Plot: `workspace/energy_reco_systs_weights_varied.png`, made by `my_analysis_syst.C`. Reco `Enu.calo`, no cuts applied. Inputs: nominal `CAF_fluxsyst.flat/0000000/*.root` (24 files) and varied `horn_current_neg/CAF.flat/*.root` (25 files), with the varied sample scaled to the nominal POT.

| Reco E [GeV] | varied / nom (red) | universe 1 (−1σ) / nom (blue) |
|---|---|---|
| 0.5–3 | ≈ 1.03–1.05 | ≈ 1.00 |
| 3.5–5 | ≈ 0.99–1.00 | ≈ 0.97–0.98 (dip) |
| > 6 | ≈ 1, scattered | ≈ 1 |

Conclusion: the two **are different**, and more statistics won't fix it.
- The stat error on red is ~0.5% near the peak, so a 3–4% gap over 5–6 bins below 3 GeV is many σ.
- The red error bar is also the right error for (red − blue): the universe is built from the same events as the nominal.
- The shapes differ too: red is an **excess at low E**, blue is a **deficit just above the peak**.
- Red/blue is roughly a flat +3–4% from 0 to 5 GeV, which looks more like a normalisation offset than a shape effect.

Likely causes, most suspicious first:
1. **The −1σ weight is a mirrored straight line.** With `vals = {0, 1}`, universe 1 gets `2·w[0] − w[1]`, i.e. the +1σ table extended out to −1σ. Nothing in the file gives the real −1σ response, and horn-current effects are often asymmetric. Checks: also plot universe 2 (+1σ), since the macro only draws universe 1, and fix the legend so it shows the sigma.
2. **Normalisation.** Print `sMuonRecoNuEnergyVaried.POT()` and `pot_n` next to the selected event counts. Check that the varied flux files count POT the same way as the nominal ones, and that "neg" really means −1σ in the units of the dial.
3. **A different physical change.** The varied flux might use a different horn current, or have other settings changed as well.

Quantitative check to add to the macro, after `r_us` is filled:
```cpp
double chi2 = 0; int ndf = 0;
for (int b = 1; b <= r_v->GetNbinsX(); ++b) {
  const double e = r_v->GetBinError(b);
  if (e <= 0 || h_n->GetBinContent(b) <= 0) continue;
  const double d = r_v->GetBinContent(b) - r_us[0]->GetBinContent(b);
  chi2 += d*d/(e*e); ++ndf;
}
std::cout << "chi2/ndf (varied vs universe 1) = " << chi2 << "/" << ndf << std::endl;
```
- Expect χ²/ndf ≫ 1, driven by the bins below 5 GeV.
- Then repeat with `h_v` scaled to the same integral as `h_n`. If χ² drops a lot, the problem is mostly normalisation (cause 2). If not, it is shape (cause 1 or 3).

## Next steps

- [x] Rerun `caffluxweighter --nd --fhc` on all nominal CAFs (24 files) and the varied ones.
- [x] Redo the reweight vs varied comparison on all files (done at reco level, 2026-11-24: still incompatible, see [[#Full statistics, reco level (2026-11-24)]]).
- [ ] Compute χ² varied vs −1σ, with and without normalising to the same integral.
- [ ] Plot universe 2 (+1σ) too; print the POT of both samples.
- [ ] Compare nominal and varied **flux histograms** directly at the SAND position, before GENIE and reco.
- [ ] Check the horn-current shift used for the varied flux vs the table's 1σ.
- [ ] Apply `kNuMuCC_SAND_FV` in the macro.
- [ ] Ask for `baseline` and `beam.hornI` to be filled in the SAND CAFMaker.
- [ ] Commit the caffluxweighter changes on the branch.

## Related
- [[Beam monitoring]]
- [[scp to pnfs Permission denied]]
