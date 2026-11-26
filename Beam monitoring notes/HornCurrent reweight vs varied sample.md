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

## Results without cuts (2026-11-25)

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

## Results with νμ CC + SAND FV selection (2026-11-25)

Same files and binning, with `kNuMuCC_SAND_FV` applied to both loaders: SAND particles, exactly 1 primary muon, and `inSANDFV` (= `inSANDInnerVol || inFrontECAL`). About 4% of the events survive (71k out of 1.8M).

```
POT ratio   varied/nominal = 1.03943
Event ratio varied/nominal = 1.02683  (raw: 72800 / 70898)
Events after POT scaling: varied/nominal = 0.987875

H0: no effect (ratio = 1)                  chi2/ndf = 63.8665  / 50   p = 0.0899272
H1: -1 sigma reweight (universe 1)         chi2/ndf = 52.4018  / 50   p = 0.380979
H0, shape only (area normalised)           chi2/ndf = 56.3967  / 49   p = 0.217937
H1, shape only (area normalised)           chi2/ndf = 51.8612  / 49   p = 0.362919
```

### Normalisation: matches the −1σ prediction

| | events at equal POT | ratio to nominal |
|---|---|---|
| Nominal | 70,898 | 1 |
| −1σ reweight (universe 1) | 70,199 | **0.990** |
| Varied (`horn_current_neg`) | 70,038 | **0.988 ± 0.005** |
| +1σ reweight (universe 2) | 71,597 | 1.010 |

- Error: `√(1/72800 + 1/70898)` ≈ 0.5%.
- Varied vs −1σ: 0.2% apart, well within errors. About 2.3σ from "no effect" and about 4σ from +1σ.
- The +1.6% excess seen without cuts is **gone**.

### Shape: consistent too

- H1 fits well: χ²/ndf ≈ 1 (52.4/50 full, 51.9/49 shape only, p ≈ 0.37).
- H1 preferred over H0 by Δχ² ≈ 11.5 (roughly 3σ). Most of it comes from the total; shape alone gives Δχ² ≈ 4.5.
- In the plot the red points follow the blue dip at 3.5–5 GeV (≈ 0.95–0.97); elsewhere flat within errors.

### Caveats

- **Less statistical power.** Errors ≈ 1.5% per bin at the peak, 5–10% above 6 GeV, for an effect ≤ 3%. "Consistent" means no disagreement at this precision; differences of ~1% would not be seen. Beam-sim fluctuations are probably small compared with these errors.
- **The no-cut excess lived in the rejected events**: interactions outside the SAND FV or not νμ CC with a reconstructed muon (mostly magnet, yoke, rock). A horn-current change can't produce an excess only outside SAND, so the varied production probably differs in something affecting only those events:
  - generation volume or geometry;
  - size or position of the flux window;
  - spatial (off-axis) extent of the flux.
- Check: run with the FV cut inverted (`!inSANDFV`) and plot vertex X and Z, nominal vs varied (the X/Z spectra are already filled in the macro but never plotted). An excess at large |x| or in a specific Z region points to geometry or flux-window differences.

### Conclusion

For the observable actually used (selected νμ CC in SAND FV), the HornCurrent −1σ weights **reproduce the varied sample within statistics**, in total and in shape. The no-cut excess outside the FV is a separate production issue, worth understanding but not affecting this analysis.

## Why χ² and not Kolmogorov–Smirnov

For this question χ² fits better. KS is still useful as a cross-check (see the end).

1. **KS compares shape only.** It compares the **normalised** cumulative distributions, so it ignores the total. Here half the information is in the normalisation (+1.6% without cuts; 0.988 vs 0.990 with selection). χ² sees both, and the shape-only version separates them (without cuts, H1 goes 665 → 132 once the normalisation is removed). ROOT's `TH1::KolmogorovTest` option `"N"` adds a separate normalisation comparison, but it is less transparent.
2. **KS is meant for unbinned, unweighted data.** Here we have:
   - binned histograms: KS p-values on binned data are too optimistic (too close to 1), as the ROOT docs say;
   - weighted events (POT scaling, −1σ weights): the KS p-value is not well defined without corrections.
   χ² on histograms uses the per-bin errors directly, which already account for the weights.
3. **KS assumes two independent samples.** The −1σ universe is made of the **same events as the nominal**, reweighted. KS of varied vs −1σ would treat the nominal statistics as independent, which they aren't. The χ² on ratios handles this explicitly (see below).
4. **χ² extends naturally:**
   - extra errors in quadrature (e.g. the beam-sim error `s`), or a covariance matrix for systematics;
   - Δχ² between H0 and H1 (e.g. 63.9 − 52.4 = 11.5) says which hypothesis is preferred and by how much; KS only gives "compatible or not" for each one separately;
   - the CAFAna fit for beam monitoring uses χ² / Poisson likelihood, so the numbers are directly comparable.

**Limits of χ²:**
- Needs enough events per bin (≳ 10–20) for Gaussian errors. With the selection the high-E bins have a few hundred, so fine. With fewer, use a Poisson likelihood (Baker–Cousins χ²).
- Depends somewhat on the binning; unbinned KS doesn't.

**When KS is useful:** as a binning-independent cross-check of the shape, on unbinned, unweighted events of varied vs nominal. KS is especially sensitive to an overall shift of the distribution, which is the kind of effect a horn-current change produces. Without weights it can only test H0 (varied vs nominal), not H1 (vs −1σ).

## Is it OK that the model comes from weighted MC?

Yes, as long as the model's own statistical error is in the χ². In this test it is, though not explicitly.

**The model is a MC prediction, which is standard.** H1 doesn't come from the sample being tested. It comes from the **nominal** sample, a separate MC production reweighted with the −1σ weights. Every analysis works this way. Two conditions:
1. The model must be **independent** of the tested sample. Nominal and varied are separate productions, so this holds.
2. The model's **statistical error** must be included, because the nominal has a finite number of events.

**Where the model's error sits.** Per bin the test computes

`d = r_v − r_u = h_v/h_n − hs₀/h_n = (h_v − hs₀) / h_n`

i.e. varied − (−1σ prediction), up to the common 1/h_n. Its true variance is `σ²(h_v) + σ²(hs₀)`:
- `σ²(h_v)`: statistics of the varied sample.
- `σ²(hs₀) = Σ wᵢ²`: statistics of the model, from the nominal events. Since the weights are close to 1 (within ±3%), this ≈ `σ²(h_n)`.

The error used is that of `r_v = h_v/h_n`. `Divide` propagates **both** histograms' errors: `σ(r_v) ≈ √(σ²(h_v) + σ²(h_n)) / h_n`. That is the right variance: the `σ²(h_n)` term **is** the model's statistical error.

**Why `r_u` has zero error** (`SetBinError(b, 0)`): the −1σ universe and the nominal are made of the same events, so their fluctuations almost cancel in the ratio. Giving `r_u` an error too would count the nominal statistics **twice** and make χ² too small.

The same holds for H0 (model = 1): `d = (h_v − h_n)/h_n`, variance `σ²(h_v) + σ²(h_n)`, again the error of `r_v`.

**Approximations:**
- Weights ≠ 1: `Σwᵢ²` differs from `Σ1` by a few % at most, negligible on χ². Exact version: bin error² = `(σ²(h_v) + Σwᵢ²) / h_n²`.
- The weights themselves have an uncertainty: the flux tables come from a finite-statistics beam simulation. Not included, same as the beam-sim fluctuations; makes the p-values a worst case.
- Shape-only tests: here one thing **does** come from the tested sample, the scale factor (from the varied total). That's why those tests use ndf − 1.

## Is the varied sample compatible with both H0 and H1? (2026-11-26)

With the selection over 0–20 GeV, both χ²/ndf are close to 1 (H0 63.9/50, p = 0.09; H1 52.4/50, p = 0.38), so neither is rejected by its own goodness of fit. That does **not** mean the data support both equally.

- **Dilution:** only ~10 of the 50 bins (≈ 3–6 GeV) are sensitive to the horn current. The other ~40 contribute ≈ 1 each to both χ², so both totals stay close to ndf. The goodness of fit over 50 bins is a weak test for this question.
- **The right quantity is Δχ² = χ²(H0) − χ²(H1)**, the likelihood ratio between the two hypotheses: `L(H1)/L(H0) = e^(Δχ²/2)`. Over 0–20 GeV: Δχ² = 11.5, so ≈ 300× more likely under H1, ≈ 3.4σ.

## χ² in a restricted energy range (2026-11-26)

Added to the macro: `chi2_range(label, model, scale, normalised)`, the same four tests but only on the bins whose centre is in `[emin, emax]`.
- Bins: `bmin = FindBin(emin + 1e-6)`, `bmax = FindBin(emax − 1e-6)`. With 0.4 GeV bins, 2–8 GeV = bins 6–20 (15 bins).
- Shape only: the scale factor uses the totals **inside** the range (`Integral(bmin, bmax)`), otherwise a normalisation difference outside the range would remain as an offset.
- Prints Δχ² (H0 − H1) and `√Δχ²` as an approximate number of σ, full and shape only.
- Compile bug met: without `return c;` at the end of the lambda its return type is `void`, so `const double c0 = chi2_range(...)` fails with "void value not ignored as it ought to be".

### Results (selection νμ CC + SAND FV)

0–10 GeV (bins 1–25):
```
H0: no effect (ratio = 1)                  chi2/ndf = 27.8399  / 25   p = 0.315278
H1: -1 sigma reweight (universe 1)         chi2/ndf = 17.5952  / 25   p = 0.859084
H0, shape only (area normalised)           chi2/ndf = 24.2921  / 24   p = 0.445001
H1, shape only (area normalised)           chi2/ndf = 17.6109  / 24   p = 0.821465
Delta chi2 (H0 - H1): full = 10.2447 (~3.20073 sigma),  shape only = 6.68118 (~2.5848 sigma)
```

2–8 GeV (bins 6–20):
```
H0: no effect (ratio = 1)                  chi2/ndf = 25.2816  / 15   p = 0.0462865
H1: -1 sigma reweight (universe 1)         chi2/ndf = 14.9028  / 15   p = 0.458444
H0, shape only (area normalised)           chi2/ndf = 22.4049  / 14   p = 0.0706696
H1, shape only (area normalised)           chi2/ndf = 14.5694  / 14   p = 0.408207
Delta chi2 (H0 - H1): full = 10.3789 (~3.22163 sigma),  shape only = 7.83544 (~2.79919 sigma)
```

| Range | H0 full | H1 full | Δχ² full | Δχ² shape only |
|---|---|---|---|---|
| 0–20 GeV (50 bins) | 63.9 / 50 (p 0.09) | 52.4 / 50 (p 0.38) | 11.5 (~3.4σ) | 4.5 (~2.1σ) |
| 0–10 GeV (25 bins) | 27.8 / 25 (p 0.32) | 17.6 / 25 (p 0.86) | 10.2 (~3.2σ) | 6.7 (~2.6σ) |
| 2–8 GeV (15 bins) | 25.3 / 15 (**p 0.046**) | 14.9 / 15 (p 0.46) | 10.4 (~3.2σ) | 7.8 (~2.8σ) |

### Interpretation

- **The full Δχ² is stable (≈ 10–11) in every range.** Cutting bins removes the same amount of χ² from H0 and H1: the removed bins are neutral and only add noise. All the discriminating information is inside 2–8 GeV.
- **The shape-only Δχ² grows as the range narrows** (4.5 → 6.7 → 7.8). The area normalisation is computed only where the effect is, so it stops absorbing the dip at 3.5–6 GeV.
- Of the ~10 units of full Δχ², about 8 come from the **shape** (the dip) and about 2 from the overall **normalisation** (−1%).
- **In 2–8 GeV, H0 is rejected at 5%** (p = 0.046, χ²/ndf 1.7), while **H1 fits perfectly** (χ²/ndf 0.99).
- H1 at 0–10 GeV has χ²/ndf ≈ 0.7: a bit low but normal for 25 bins (typical range ≈ 0.6–1.4).

### Conclusion

With the νμ CC + SAND FV selection:
1. The varied sample is **consistent with the −1σ reweight**, in total and in shape (χ²/ndf ≈ 1 in every range).
2. It is **not consistent with "no effect"** in the sensitive region (p ≈ 0.05), and the data prefer −1σ over H0 by **~3σ** in every range.
3. Most of the discriminating power is in the **shape** (dip at 3.5–6 GeV).

Reporting caveat: choosing the range after looking at the data is mild tuning. Here it is justified, because 2–8 GeV comes from where the flux tables put the effect, not from the varied data, and Δχ² is the same in all ranges. Report all three ranges anyway.

## Next: fit of the effect size μ

Instead of testing only μ = 0 (H0) and μ = 1 (H1), leave the size of the shift free:

`model(μ) = 1 + μ · (r_u − 1)`

The model is linear in μ, so the minimum has an exact solution (`x = r_u − 1`, `y = scale·r_v − 1`, `e` = error):

```
μ̂   = Σ(x·y/e²) / Σ(x²/e²)
σ_μ = 1 / √Σ(x²/e²)
```

Outputs:
- μ̂ ± σ_μ and χ²_min / (ndf − 1);
- significance of the effect: (μ̂ − 0)/σ_μ; compatibility with −1σ: (μ̂ − 1)/σ_μ;
- expected sensitivity 1/σ_μ (depends only on the errors, not on the varied data: how powerful the test is).
- Variants: full and shape only (scale recomputed at each μ, so a small scan), 0–20 and 2–8 GeV.
- Scan of χ²(μ) − χ²_min for μ ∈ [−1, 3], saved as `workspace/chi2_mu_scan.png`, with lines at Δχ² = 1, 4 and at μ = 0, 1.
- Rough expectation from the current numbers: μ̂ ≈ 0.8–1 ± ~0.3.
- Assumes a response linear in μ, the same assumption as the weights (`1 + σ·Δ`). Fine for μ ∈ [0, ~1.5].

## What the p-value means

**Definition:** `p = P(χ² ≥ χ²_observed | hypothesis true)`. It is what `TMath::Prob(χ², ndf)` returns.

- If the hypothesis is true and the errors are Gaussian and correct, χ² follows the χ² distribution with ndf degrees of freedom (mean = ndf; for 15 ndf, typical values ≈ 7–25). p is the area of its right tail beyond the observed value.
- **H0 in 2–8 GeV, χ² = 25.3/15, p = 0.046:** if the horn current did nothing, repeating the experiment with independent statistics, only 4.6% of the time would the disagreement be this large or larger. Just below the conventional 5% threshold.
- **H1 in 2–8 GeV, χ² = 14.9/15, p = 0.46:** if the −1σ weights are right, 46% of the time the disagreement would be at least this large. Entirely ordinary.

**What it is NOT:**
- **Not the probability that the hypothesis is true.** p = 0.46 doesn't mean "H1 is true at 46%". p assumes the hypothesis is true and measures how unusual the data are under it.
- **A large p doesn't prove the hypothesis,** it only says the data don't contradict it. That's why over 0–20 GeV both H0 and H1 had acceptable p.
- **It doesn't compare two hypotheses.** Each p tests one hypothesis alone; to choose between H0 and H1 use Δχ².

**Δχ² and "σ":**
- `L(H1)/L(H0) = e^(Δχ²/2)`. With Δχ² = 10.4: e^5.2 ≈ 180, i.e. the data are ~180× more likely under H1.
- "√Δχ² σ" (≈ 3.2σ) is an **approximation**, exact when H0 and H1 are two values of one parameter (μ = 0 and μ = 1) and χ² is a parabola in it. That's why the μ fit is the cleanest way: (μ̂ − 0)/σ_μ is a well-defined number of σ.

| σ (Gaussian, two-sided) | p |
|---|---|
| 1 | 0.32 |
| 2 | 0.046 |
| 3 | 0.0027 |
| 5 | 5.7 × 10⁻⁷ |

So H0's p = 0.046 in 2–8 GeV ≈ 2σ.

**Caveats for this case:**
1. **The errors must be right.** Here only event statistics; beam-sim fluctuations not included. If they matter, the true χ² are smaller and the p-values larger: H0 would be less excluded.
2. **Several tests on the same data.** 3 ranges × 4 tests. Picking the smallest p a posteriori underestimates the chance of finding "something" by accident. Small effect here (range motivated a priori, Δχ² stable), but report all of them.
3. **Enough events per bin** for Gaussian errors: thousands per bin in 2–8 GeV, so fine.

## Next steps

- [x] Run the macro with the χ² tests and record the results above (2026-11-25).
- [x] Rerun with the νμ CC + SAND FV selection (2026-11-25): consistent with −1σ.
- [x] χ² in restricted ranges 0–10 and 2–8 GeV (2026-11-26): H1 preferred by ~3σ in every range; H0 rejected at 5% in 2–8 GeV.
- [ ] Implement the μ fit (μ̂ ± σ_μ, χ²(μ) scan plot).
- [ ] Inverted FV (`!inSANDFV`): plot vertex X and Z, nominal vs varied, to locate the no-cut excess.
- [ ] Compare the integrated flux per POT at SAND, nominal vs varied (explains the +1.6%?).
- [ ] Add the extra-error `s` scan to the macro for the shape-only H1 test.
- [ ] Compare nominal and varied flux histograms with errors; ideally nominal A vs nominal B.
- [ ] Find out which horn current (kA) was used for `horn_current_neg`, vs the table's 1σ.
- [ ] If possible, get a +1σ varied sample for the linearity test.

## Related
- [[Flux systematics weights]]
- [[Beam monitoring]]
