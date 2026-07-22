# Benchmark #1b — PRE-REGISTRATION (nonlinear PfCRT transport → resistance index)

**Status: PRE-REGISTERED.** This document is committed and **tagged (`prereg-1b`)
BEFORE any Table S2 phenotype value is extracted.** The model form, the C\* handling,
the strain list and the pass/fail criteria below are **frozen**: they will not be
changed after seeing any result. The run happens **once**.

## Hypothesis under test
The saturable (Michaelis–Menten) PfCRT chloroquine-transport efficiency, calibrated
on Dd2 and frozen, predicts the chloroquine **resistance index (RI)** of held-out
PfCRT variants.

## Quantity — resistance index, not absolute IC50
`RI_s = IC50(CQR strain s) / IC50(CQS reference)`. Calibration **and** held-out use
the **same** phenotype set (Summers 2014 Table S2). We do **not** mix the isogenic-
line IC50 (132/25 nM, PMC9067703) with Table S2 indices. `RI_Dd2` from Table S2 is
the single calibration datum.

## Model (frozen form)
```
RI_s(C*) = 1 + (RI_Dd2 − 1) · [ Vmax_s /(Km_s + C*) ] / [ Vmax_Dd2 /(Km_Dd2 + C*) ]
```
- Eliminates the lumped constant `B` explicitly; one calibration datum (`RI_Dd2`).
- Limit `C*/Km → 0` recovers the linear #1a model
  `RI_s = 1 + (RI_Dd2 − 1)·(Vmax_s/Km_s)/(Vmax_Dd2/Km_Dd2)`.
- `Km_s`, `Vmax_s` from Summers 2014 **Table 1** (measured, same experiment):
  Dd2 232±11 / 61±6 · 7G8 117±6 / 9±1 · K1 293±10 / 79±17 · GB4 275±16 / 57±5 ·
  Ecu1110 191±17 / 36±4  (µM / pmol·oocyte⁻¹·h⁻¹).

## C\* is ASSUMED — scenario analysis, not a probability distribution
Summers estimates ~200–600 µM CQ in the digestive vacuole of resistant parasites
during treatment, near PfCRT saturation, but does **not** establish that this is the
*free* concentration producing 50% of the effect (total drug, heme binding and other
differences remain open in the paper). Therefore:
- `C*` is a **scenario** variable, `kind = assumed`.
- Sweep `C* ∈ {0 (linear limit), 200, 400, 600} µM`; report across the **whole** range.
- We do **not** sample `C*` from a distribution and do not pick a "most likely" value.
- Monte-Carlo uncertainty is over the **measured** parameters only (Km, Vmax, RI
  SEMs), run **separately at each C\* scenario**.

## Strains
- **Calibration:** Dd2 (one datum, `RI_Dd2` from Table S2).
- **Held-out (frozen, predicted simultaneously):** K1, GB4, Ecu1110, **and re-test 7G8**.

## Table S2 audit (to complete AT extraction, before running)
For each strain, document in the results file: resistant IC50; the **CQS reference**
used in the denominator; the assay and original source; and whether the RIs are
mutually comparable. **"Same table" ≠ "same experiment".** These are whole strains,
not isogenic lines — so a success supports **predictive adequacy on this benchmark**,
NOT causal proof that PfCRT alone explains resistance.

## Pass/fail criteria (frozen)
Per strain, **PASS iff all three**:
1. observed RI lies within the predicted **90% CI**;
2. median fold-error `max(pred/obs, obs/pred) ≤ 1.5×`;
3. predicted interval width `upper/lower ≤ 2×` (guards against a "giant interval" pass).

Benchmark **ADEQUATE iff**:
- **≥ 3 of 4** held-out strains PASS; **and**
- **no** individual median fold-error exceeds **2×**.

Report additionally (diagnostics, not pass/fail): `RMSE(log2 RI)` and geometric-mean
fold-error. All of the above are reported **per C\* scenario**.

## Saturation conclusion (separate, frozen criterion)
> The linearization was the problem **iff** #1b, relative to #1a: (a) **reduces**
> held-out `RMSE(log2 RI)`; **and** (b) converts **7G8 from outside to inside** its
> predicted 90% CI at some `C* ∈ [200,600]`.
> Otherwise: **saturation does not explain the 7G8 failure**, and the candidate shifts
> to κ non-transferability across systems, or the transport→phenotype hypothesis itself.

## Execution
Extract Table S2 phenotypes (with the audit above), then run **once**. No parameter,
model choice, or criterion is altered after seeing results. Slice 2 (enzymatic) is
held until this run reports.
