---
tags: [beam-monitoring, sand, flux, systematics, cafana, statistics]
date: 2026-11-25
---

# HornCurrent reweight vs varied sample

Does reweighting the nominal sample with the HornCurrent −1σ flux weights reproduce the `horn_current_neg` sample, which was generated with a different flux? If it does, the weights are a reliable way to model a horn-current shift in SAND, and can be trusted for beam monitoring.

How the weights are computed, and the earlier 1-file comparison: [[Flux systematics weights]].

## Setup

- Macro: `lblpwgtools/bm_analysis/my_analysis_syst.C`. Output: `workspace/energy_reco_systs_weights_varied.png`.
- Nominal: `fast_reco/nominal/CAF_fluxsyst.flat/0000000/*.root` (24 files).
- Varied: `fast_reco/varied/horn_current_neg/CAF.flat/*.root` (25 files). This is the only varied sample available (no +1σ).
- Variable: reco `Enu.calo`, 50 bins over 0–20 GeV, no cuts applied.
- Universes from `Multiverse::Hypercross(systs, 1)`: 0 = nominal, 1 = −1σ (blue), 2 = +1σ (green).
- Varied sample scaled to the nominal POT with `ToTH1(pot_n)`.

| | POT |
|---|---|
| Nominal | 1.73215e18 |
| Varied | 1.80045e18 |
| Ratio varied/nominal | 1.0394 |

## What the plot shows (all files)

| Reco E [GeV] | varied / nom (red) | −1σ / nom (blue) | +1σ / nom (green) |
|---|---|---|---|
| 0.5–3 | ≈ 1.03–1.05 | ≈ 1.00 | ≈ 1.00 |
| 3.5–5 | ≈ 0.99–1.00 | ≈ 0.97–0.98 (dip) | ≈ 1.02–1.03 |
| > 6 | ≈ 1, scattered | ≈ 1 | ≈ 1 |

- The stat error on red is ~0.5% near the peak, so the 3–5% gap below 3 GeV is many σ **on the event error bars**.
- Blue and green are exact mirror images around 1. That is by construction: weight = `1 ± Δ`. So this plot **cannot test linearity**.

## Possible explanations and how to test each

| # | Hypothesis | Test |
|---|---|---|
| 1 | **Beam-simulation fluctuations.** The varied flux comes from its own beam simulation, whose random fluctuations are shared by every event and are **not** in the event error bars. | Compare the nominal and varied **flux histograms** with their errors. Best: compare two nominal fluxes with different seeds (nominal A / nominal B). |
| 2 | **Normalisation (POT).** The two samples count POT differently. | Raw event ratio vs POT ratio (1.0394). χ² after area normalisation. |
| 3 | **A different "−1σ".** The varied sample isn't the same kA shift as the table's 1σ, or other beam settings changed too. | Ask the producer which horn current was used; compare with the 1σ that defines the flux tables. |
| 4 | **Non-linearity.** The −1σ weight is just the +1σ response mirrored. | Needs a +1σ (or ±2σ) varied sample, or a table at another shift (e.g. the `_280kA_` variants). Not possible with current samples. |

### Linearity test (when a +1σ varied sample exists)

Per bin, with `R±` = varied/nominal at ±1σ:

`A = (R₊ − 1) + (R₋ − 1)`

- Linear response: `A` consistent with 0.
- A coherent shape in `A`: there is a quadratic (non-linear) term.

### What "statistical fluctuation" means here

Nominal and varied are independent MC productions. Even if the horn current did nothing, varied/nominal would scatter around 1 by the statistical errors. The question is whether the scatter is larger than that.

There are two kinds of statistical error:
- **Event statistics** (number of GENIE events): in the error bars; shrinks with more files.
- **Beam-simulation statistics** (finite number of protons in the flux simulation): fully correlated across all events from that flux; **not** in the error bars, and doesn't shrink with more GENIE files. With ~1.8M events this may now dominate.

### Flux-level comparison: how many files

- One nominal flux file vs one varied flux file is enough, **if** each is the full flux its sample was generated from (use the merged one if the production merged several beam-sim runs), the histograms carry errors (`Sumw2`), and both are taken at the SAND position, numu FHC, same binning and POT normalisation.
- Reading the varied/nominal flux ratio:
  - error bars ≳ 1% per bin below 5 GeV: beam-sim fluctuations alone can explain the +3–5% excess;
  - tiny error bars and the excess already present: it is a real feature of the varied flux;
  - no excess at flux level: look downstream (GENIE, reco, POT accounting).
- Two CAF files (one nominal, one varied) are **not** enough: ~2–4% per bin, larger than the effect.

## χ² tests (added to the macro, 2026-11-25)

Printed after the POT lines:

```
POT ratio   varied/nominal = ...
Event ratio varied/nominal = ...   (raw counts, each sample at its own POT)
Events after POT scaling: varied/nominal = ...

Chi2 tests on varied/nominal (HornCurrent):
H0: no effect (ratio = 1)            chi2/ndf = ... / ...   p = ...
H1: -1 sigma reweight (universe 1)   chi2/ndf = ... / ...   p = ...
H0, shape only (area normalised)     chi2/ndf = ... / ...   p = ...
H1, shape only (area normalised)     chi2/ndf = ... / ...   p = ...
```

Details:
- Errors = stat errors of `r_v` (event statistics of both samples). Beam-sim fluctuations **not** included, so the p-values are a worst case for the weights.
- Bins with empty nominal are skipped.
- The shape-only H0 test scales varied to the **nominal** total. The shape-only H1 test scales it to the **−1σ universe** total, since the reweight itself lowers the total by ~1%. Both lose one degree of freedom (`ndf − 1`).
- `p = TMath::Prob(chi2, ndf)` (needs `#include "TMath.h"`).

Code:

```cpp
auto chi2_test = [&](const char* label, const TH1* ratio, const TH1* model, double scale, bool normalised) {
  double c = 0; int n = 0;
  for (int b = 1; b <= ratio->GetNbinsX(); ++b) {
    const double e = scale * ratio->GetBinError(b);
    if (e <= 0 || h_n->GetBinContent(b) <= 0) continue;
    const double m = model ? model->GetBinContent(b) : 1.;
    const double d = scale * ratio->GetBinContent(b) - m;
    c += d * d / (e * e); ++n;
  }
  if (normalised) --n;   // area normalisation uses up one degree of freedom
  std::cout << label << " chi2/ndf = " << c << " / " << n << "   p = " << TMath::Prob(c, n) << "\n";
};
const double area   = h_n->Integral()   / h_v->Integral();   // varied -> nominal area
const double area_u = hs[0]->Integral() / h_v->Integral();   // varied -> -1 sigma universe area
chi2_test("H0: no effect (ratio = 1)",          r_v, nullptr, 1.,     false);
chi2_test("H1: -1 sigma reweight (universe 1)", r_v, r_us[0], 1.,     false);
chi2_test("H0, shape only (area normalised)",   r_v, nullptr, area,   true);
chi2_test("H1, shape only (area normalised)",   r_v, r_us[0], area_u, true);
```

### How to read the results

- **Event ratio ≈ POT ratio (1.0394):** the POT in the files is consistent.
- **Full tests large, shape-only ≈ 1 per degree of freedom:** the disagreement is only normalisation.
- **Shape-only still large:** the shapes really differ. Beam-sim fluctuations (1) or a different horn-current shift (3) become the main suspects.
- **χ²(H1) < χ²(H0):** the −1σ reweight goes in the right direction, even if it doesn't match fully.

Possible outcomes:
- Normalisation or beam-sim fluctuations explain it: the weights aren't contradicted, but a higher-statistics beam simulation is needed to test them properly.
- The excess is real in the flux itself and the horn-current shift matches 1σ: the weights don't describe the variation well, which matters for beam monitoring.

## Macro fixes (2026-11-25)

- Legend: `SetHeader` keeps only one header, so the second call replaced the first; `AddEntry(Form(...))` with only a string looks up an object by name, giving empty entries. Now the event counts go in the labels: `Nominal (N)`, `Varied, horn current neg (N)`, `Universe i (±1σ, N)`.
- Universe colours chosen by hand: `kBlue + 1` (−1σ), `kGreen + 2` (+1σ), avoiding red (varied).

## Results (2026-11-25)

All files (24 nominal, 25 varied), reco `Enu.calo`, 50 bins 0–20 GeV, no cuts.

```
POT ratio   varied/nominal = 1.03943
Event ratio varied/nominal = 1.05608  (raw: 1.89231e+06 / 1.79182e+06)
Events after POT scaling: varied/nominal = 1.01602

H0: no effect (ratio = 1)                  chi2/ndf = 446.219  / 50   p = 5.31695e-65
H1: -1 sigma reweight (universe 1)         chi2/ndf = 665.478  / 50   p = 1.83422e-108
H0, shape only (area normalised)           chi2/ndf = 228.3    / 49   p = 5.94997e-25
H1, shape only (area normalised)           chi2/ndf = 131.915  / 49   p = 1.56102e-09
```

Universe ordering checked in `CAFAna/Core/Multiverse.cxx` (`Hypercross`): 0 = nominal, 1 = −1σ, 2 = +1σ. So `hs[0]` / `r_us[0]` = −1σ.

### Normalisation

| | varied / nominal |
|---|---|
| raw events | 1.056 |
| POT | 1.039 |
| events at equal POT | **1.016** |

- The varied histogram is scaled to the nominal POT, which is correct: the observable is compared at equal exposure. The 1.016 is **after** that scaling.
- At equal POT the varied sample has **+1.6%** more events. Event-stat error ≈ 0.1% (~1.8M events each), so ≈ 15σ. Same as the 1-file result (1.016 ± 0.006), so it is stable.
- The −1σ reweight predicts the opposite, about −1%. A weaker horn current should focus fewer pions, so fewer neutrinos, not more.
- So the total gap varied vs −1σ prediction is ~2.6%, as large as the horn-current effect itself.
- That's why the full H1 test (χ² 665) is **worse** than H0 (446): the reweight lowers the total while the varied sample raises it.

### Shape

- Shape only, the −1σ reweight improves χ² from 228 to 132 (Δχ² ≈ 96), so the reweight shape (dip at 3.5–5 GeV) goes in the right direction.
- Residual χ²/ndf ≈ 2.7 is still too much for event statistics alone.
- The tiny p-values only mean "rejected **if event statistics were the only error**". Beam-sim fluctuations are not included.

### Interpretation

- Event statistics are no longer the limit.
- The **normalisation excess is the main issue**. Most likely cause: the varied flux or its POT is normalised differently (per proton of its own beam sim, or a different POT per flux file), so "1 POT" in the varied files doesn't mean the same as in the nominal. POT scaling is only right if the stored POT matches the flux actually used.
- The leftover shape difference may be beam-sim fluctuations.

### Checks to do

- **Normalisation:** compare the **integrated flux per POT** at SAND in the nominal and varied flux files.
  - Varied already ~1.6% higher per POT: the excess comes from the flux normalisation, not GENIE or reco.
  - Equal per POT: look at how POT was accumulated in the varied CAFs (number of spills, POT per file).
- **Shape:** add an extra uncorrelated per-bin error `s` in quadrature and find the value giving χ²/ndf ≈ 1 for shape-only H1.
  - `s` ≲ 1%: plausible beam-sim fluctuation, so the shape is consistent with the reweight.
  - `s` of several %: the shapes really differ; check the horn current used for the varied sample.

## Next steps

- [x] Run the macro with the χ² tests and record the results above (2026-11-25).
- [ ] Compare the integrated flux per POT at SAND, nominal vs varied (explains the +1.6%?).
- [ ] Add the extra-error `s` scan to the macro for the shape-only H1 test.
- [ ] Compare nominal and varied flux histograms with errors; ideally nominal A vs nominal B.
- [ ] Find out which horn current (kA) was used for `horn_current_neg`, vs the table's 1σ.
- [ ] If possible, get a +1σ varied sample for the linearity test.

## Related
- [[Flux systematics weights]]
- [[Beam monitoring]]
