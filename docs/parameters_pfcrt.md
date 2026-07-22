# Benchmark #1 — chloroquine / PfCRT: traceable parameters & method

Falsifiable test: **predict 7G8 chloroquine resistance from PfCRT transport,
calibrated only on Dd2, with 7G8's phenotype held out.** Every parameter below is
tagged by epistemic kind — `measured`, `derived`, `fitted`, `assumed` — and cited.

## Parameter table

| Parameter | Value | Uncertainty | Unit | Kind | Source |
|---|---|---|---|---|---|
| Digestive-vacuole pH | 5.2 | 5.0–5.6 | pH | measured | Kuhn 2007, *Cell Microbiol*, [10.1111/j.1462-5822.2006.00847.x](https://doi.org/10.1111/j.1462-5822.2006.00847.x) (pHluorin 5.18±0.05) |
| Parasite cytosol pH | 7.15 | ±0.07 | pH | measured | Kuhn 2007 |
| Chloroquine pKa₁ | 10.1 | 10.0–10.3 | pKa | measured | Warhurst (diprotic base; side-chain amine) |
| Chloroquine pKa₂ | 8.1 | 8.0–8.4 | pKa | measured | Warhurst (quinoline N) |
| Km (CQ, PfCRT‑Dd2) | 250 | ±40 | µM | measured | Summers 2014 *PNAS* [10.1073/pnas.1322965111](https://doi.org/10.1073/pnas.1322965111); confd [PMC6709616](https://pmc.ncbi.nlm.nih.gov/articles/PMC6709616/), [PMC4276893](https://pmc.ncbi.nlm.nih.gov/articles/PMC4276893/) |
| Vmax (CQ, PfCRT‑Dd2) | 33 | ±3 | pmol/oocyte/h | measured | Summers 2014; PMC6709616 |
| Km (CQ, PfCRT‑7G8) | 183 | 117–250 | µM | **assumed** | Summers 2014 *group* (exact Table‑1 value paywalled); high‑affinity group |
| Vmax (CQ, PfCRT‑7G8) | 20 | 10–33 | pmol/oocyte/h | **assumed** | Summers 2014 *group*; low‑capacity group; 7G8 is a transport‑vs‑resistance outlier |
| IC₅₀ CQ, isogenic WT (GC03) | 25 | ±1.1 | nM | measured | [PMC9067703](https://pmc.ncbi.nlm.nih.gov/articles/PMC9067703/) Table 1 (C2GC03) |
| IC₅₀ CQ, isogenic Dd2 | 132 | ±4.9 | nM | measured | PMC9067703 (C4Dd2) — **calibration** |
| IC₅₀ CQ, isogenic 7G8 | 84 | ±3.5 | nM | measured | PMC9067703 (C67G8) — **held‑out** |
| κ (lumped transport constant) | 32.4 | 22.7–44.2 | — | **fitted** | calibrated on Dd2 alone (this model) |
| fold(Dd2) = IC₅₀_Dd2/IC₅₀_WT | 5.28 | — | — | derived | from measured IC₅₀ |
| Ideal accumulation ratio (DV/cytosol) | 7.1×10³ | — | — | derived | Z(pH_DV)/Z(pH_cyt); *equilibrium limit, not experimental* |

Assay note: transport kinetics measured in *Xenopus* oocytes at pH 6.0; IC₅₀ from
isogenic transfectants differing only in the *pfcrt* allele (same assay), so
IC₅₀ differences are attributable to PfCRT.

## Model

First-order (specificity-constant) regime of the ion-trapping ODE (`exposure.py`),
valid when DV free concentration ≪ Km:

```
fold(strain) = IC50(strain)/IC50(WT) = 1 + κ · (Vmax/Km)_strain
```

`κ = α·cf_dv/(k_perm·f_dv)` lumps the oocyte→cell scaling α, passive permeability
k_perm, and pH/pKa speciation — all of which **cancel from the fold ratio** except
this one constant. `Vmax/Km` is the transporter's specificity constant.

## Calibration and held-out (the contract)

- **Fitted:** exactly one scalar, `κ`, from one datum (Dd2's fold). One parameter,
  one constraint → **identifiable**; we do not fit several free parameters to one IC₅₀.
- **Frozen:** `κ` uses only WT + Dd2. `predict_ic50()` takes no 7G8 IC₅₀ argument.
- **7G8 kinetics enter as input** → the held-out target is 7G8's **phenotype (IC₅₀)**.
- **Reproducing Dd2 is calibration, not validation.** Only 7G8 is the test.

## Result (Monte Carlo, n=50 000, seed=12345)

| Quantity | Value |
|---|---|
| Predicted fold(7G8) | 4.75× (90% CI 2.71–8.35) |
| Predicted IC₅₀(7G8) | **119 nM (90% CI 68–208)** |
| Observed IC₅₀(7G8) [held-out] | 84 ± 3.5 nM |
| Fold error (pred median / obs) | 1.41× (model over-predicts) |
| Observed inside predicted 90% CI? | **yes** |
| Dominant uncertainty (|Spearman ρ|) | Vmax_7G8 (0.75), Km_7G8 (0.47) — both **assumed** |

**Honest reading.** The frozen prediction is *consistent* with the held-out 7G8
IC₅₀ (84 nM sits inside the 68–208 nM interval), but the interval is wide and the
median over-predicts by 1.41×. The sensitivity analysis shows the result is
dominated by 7G8's **assumed** transport kinetics — so this is a *weak* pass, not a
confirmation. It also matches a known caveat: 7G8 is a documented outlier in the
transport↔resistance correlation (excluding it raises the literature R² from ~0.86
to ~0.98), so a transport-only model is expected to be imperfect for 7G8.

**Highest-value next data point:** the exact 7G8 CQ Km/Vmax (Summers 2014 Table 1).
Swapping the two `assumed` rows for `measured` values collapses the dominant
uncertainty and turns this into a sharp, strongly falsifiable test.
