# Benchmark #1 — chloroquine / PfCRT: traceable parameters & method

Falsifiable test: **predict 7G8 chloroquine resistance from PfCRT transport,
calibrated only on Dd2, with 7G8's phenotype held out.** Every parameter is tagged
by epistemic kind — `measured`, `derived`, `fitted`, `assumed` — and cited.

## Parameter table

| Parameter | Value | Uncertainty | Unit | Kind | Source |
|---|---|---|---|---|---|
| Digestive-vacuole pH | 5.2 | 5.0–5.6 | pH | measured | Kuhn 2007 *Cell Microbiol* [10.1111/j.1462-5822.2006.00847.x](https://doi.org/10.1111/j.1462-5822.2006.00847.x) (pHluorin 5.18±0.05) |
| Parasite cytosol pH | 7.15 | ±0.07 | pH | measured | Kuhn 2007 |
| Chloroquine pKa₁ | 10.1 | 10.0–10.3 | pKa | measured | Warhurst (diprotic base; side-chain amine) |
| Chloroquine pKa₂ | 8.1 | 8.0–8.4 | pKa | measured | Warhurst (quinoline N) |
| Km (CQ, PfCRT‑Dd2) | 232 | ±11 (n=15) | µM | measured | Summers 2014 *PNAS* [10.1073/pnas.1322965111](https://doi.org/10.1073/pnas.1322965111), **Table 1** (ANU open PDF, p.E1762) |
| Vmax (CQ, PfCRT‑Dd2) | 61 | ±6 (n=15) | pmol/oocyte/h | measured | Summers 2014 Table 1 |
| Km (CQ, PfCRT‑7G8) | 117 | ±6 (n=6) | µM | measured | Summers 2014 Table 1 — **same experiment as Dd2** |
| Vmax (CQ, PfCRT‑7G8) | 9 | ±1 (n=6) | pmol/oocyte/h | measured | Summers 2014 Table 1 — same experiment as Dd2 |
| IC₅₀ CQ, isogenic WT (GC03) | 25 | ±1.1 | nM | measured | [PMC9067703](https://pmc.ncbi.nlm.nih.gov/articles/PMC9067703/) Table 1 (C2GC03) |
| IC₅₀ CQ, isogenic Dd2 | 132 | ±4.9 | nM | measured | PMC9067703 (C4Dd2) — **calibration** |
| IC₅₀ CQ, isogenic 7G8 | 84 | ±3.5 | nM | measured | PMC9067703 (C67G8) — **held‑out** |
| κ (lumped transport constant) | 16.3 | 13.2–20.3 | — | **fitted** | calibrated on Dd2 alone (this model) |
| fold(Dd2) = IC₅₀_Dd2/IC₅₀_WT | 5.28 | — | — | derived | from measured IC₅₀ |
| Ideal accumulation ratio (DV/cytosol) | 7.1×10³ | — | — | derived | Z(pH_DV)/Z(pH_cyt); *equilibrium limit, not experimental* |

**Same-experiment rule (critical).** Vmax is an absolute rate whose scale varies
between assays (Dd2 Vmax is 61 pmol/oocyte/h in Summers 2014 Table 1 but 33 in the
2019 phosphomimetic study). Because the benchmark compares Vmax/Km *across* strains,
Dd2 and 7G8 are taken from the **same** Table 1. Assay: *Xenopus* oocytes, ~pH 6.0;
IC₅₀ from isogenic transfectants (pfcrt‑only difference), so IC₅₀ differences are
attributable to PfCRT.

## Model

First-order (specificity-constant) regime of the ion-trapping ODE (`exposure.py`),
valid when DV free concentration ≪ Km:

```
fold(strain) = IC50(strain)/IC50(WT) = 1 + κ · (Vmax/Km)_strain
```

`κ = α·cf_dv/(k_perm·f_dv)` lumps the oocyte→cell scaling α, passive permeability
k_perm and pH/pKa speciation — all of which **cancel from the fold ratio** except
this one constant. `Vmax/Km` is the transporter's specificity constant.

## Calibration and held-out (the contract)

- **Fitted:** exactly one scalar, `κ`, from one datum (Dd2's fold). **Structurally
  identifiable under this assumed linear model** — not a general practical-
  identifiability claim. We do not fit several free parameters to one IC₅₀.
- **Frozen:** `κ` uses only WT + Dd2. `predict_ic50()` takes no 7G8 IC₅₀ argument
  (enforced by a signature test).
- **7G8 kinetics enter as measured input** → the held-out target is 7G8's **phenotype (IC₅₀)**.
- **Reproducing Dd2 is calibration, not validation.** Only 7G8 is the test.

## Result (Monte Carlo, n=50 000, seed=12345)

| Quantity | Value |
|---|---|
| Predicted fold(7G8) | 2.25× (90% CI 1.93–2.68) |
| Predicted IC₅₀(7G8) | **56 nM (90% CI 48–66)** |
| Observed IC₅₀(7G8) [held-out] | 84 ± 3.5 nM |
| Fold error (pred median / obs) | 0.67× (model **under**-predicts ~1.5×) |
| Observed inside predicted 90% CI? | **NO** |
| Dominant uncertainty (\|Spearman ρ\|) | Vmax_7G8 (0.62), Vmax_Dd2 (0.55) |

**Honest reading — first falsifiable result.** With Dd2 and 7G8 kinetics measured
in the *same* experiment, the transport-only model calibrated on Dd2 predicts 56 nM
(90% CI 48–66), but the held-out 7G8 IC₅₀ is 84 nM — **outside** the interval. So a
CQ‑transport‑only account, tuned on Dd2, **under-predicts 7G8 resistance by ~1.5×**:
the transport-only hypothesis is **rejected for 7G8**. This is the benchmark working
as intended (it is falsifiable and it discriminates), and it reproduces a known
result — 7G8 is a documented transport‑vs‑resistance **outlier** (excluding it raises
the literature R² from ~0.86 to ~0.98 for Vmax). Interpretation: 7G8 carries CQ
resistance beyond what its (low) transport efficiency explains, pointing to
additional biology (fitness/host-factor contributions) to model in later slices.
The model does get the **rank** right (Dd2 > 7G8 > WT); it misses 7G8's magnitude.

**Not this benchmark:** the ~10⁴× ion-trapping ratio is the model's ideal
equilibrium limit from pH/pKa, not experimental cellular accumulation.
