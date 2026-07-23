# Benchmark #1b — RESULT: NOT EVALUABLE (phenotype heterogeneity)

**Pre-registered outcome (frozen).** The confirmatory Benchmark #1b (`prereg-1b`) is
declared **NOT EVALUABLE**. The held-out phenotypes (Table S2 CQ resistance indices)
fail the pre-registered comparability audit, which is an explicit **stop rule** in the
protocol ("Se forem incompatíveis, declarar o benchmark 'não avaliável', sem adaptar
os dados"). No PASS/FAIL is assigned to the model.

This is not a failure of the method — it is the pre-registration working: it prevented
a spurious PASS/FAIL from being read off non-comparable data.

## Data extracted (frozen) — Table S2

Source: **Summers PhD thesis (ANU), Table 6.3**, p.162, which reproduces PNAS Table S2
(the main paper compiled its Fig. 6 correlation from these values). CQ resistance index
`RI = IC50(CQR strain) / IC50(CQS strain)`.

| Strain | RI (mean ± SEM) | n studies |
|---|---:|---:|
| Dd2 | 11.13 ± 1.40 | 8 |
| K1 | 11.79 ± 2.01 | 4 |
| GB4 | 8.56 ± 1.70 | 3 |
| Ecu1110 | 6.50 ± 0.89 | 2 |
| 7G8 | 6.63 ± 0.60 | 9 |
| Ph1 | 5.04 | 1 |

Per-study values (sources): Dd2 — Fidock 2000 10.69, Mu 2003 14.50, Johnson 2004 11.47,
Sidhu 2007 5.58, Sá 2009 13.59, Yuan 2009 17.42, Sanchez 2011 7.51, Ch'ng 2013 8.28.
7G8 — Mu 2003 7.92, Chen 2003 5.68, Johnson 2004 8.24, Sidhu 2007 4.71, Lehane 2008 7.43,
Sá 2009 8.11, Yuan 2009 8.71, Sanchez 2011 4.17, Ch'ng 2013 4.73. (K1/GB4/Ecu1110/Ph1
likewise pooled across the same study set.)

## Audit → why NOT EVALUABLE

The table's own note states each RI uses the CQS strain measured **in that same study**;
when a study has several CQS lines, their mean is used; then results from **different
studies** are aggregated. Therefore:
- there is **no single, common CQS denominator**;
- there is **no common assay or laboratory** (per-study IC50 protocols differ);
- the per-study spread is large (e.g. Dd2 5.58–17.42; 7G8 4.17–8.71).

This is exactly the "heterogeneous compilation" the pre-registered stop rule guards
against. A confirmatory PASS/FAIL on these numbers would be an artifact of pooling
non-comparable assays, so none is assigned.

Note also: these are **whole strains**, not isogenic lines — so even if they were
comparable, a success would support *predictive adequacy*, never causal proof that
PfCRT alone explains resistance.

## C\* audit — `derived`/`assumed`, not `measured`

The 200–600 µM digestive-vacuole [CQ] used for the saturation scenario is **not a direct
measurement**. The thesis derives it (external pH 7.3, DV pH 5.1, pKa 8.1/10.2, ~60%
plasma binding, ~10× lower accumulation in CQR):

| Derived quantity | Day 3 | Day 7 |
|---|---:|---:|
| Mean plasma CQ | 564 ± 54 nM | 246 ± 28 nM |
| Predicted DV, CQS | 4 890 ± 469 µM | 2 132 ± 244 µM |
| Predicted DV, CQR | 489 ± 47 µM | 213 ± 24 µM |

So the CQR DV range (~213–489 µM) motivates the `C* ∈ {200,400,600}` scenario sweep, and
`C*` is correctly classified **derived/assumed**, treated as scenarios (not sampled).

## Exploratory analysis (clearly separated — NO pass/fail)

Run `examples/benchmark_1b_exploratory.py`. Model calibrated on the Table-S2 `RI_Dd2`
(11.13), predicting the held-outs. **Hypothesis-generating only**; the confirmatory
verdict remains NOT EVALUABLE.

| Strain | pred RI (C*=0) | pred RI (C*=400) | obs RI | fold err (C*=0 → 400) |
|---|---:|---:|---:|---|
| K1 | 11.22 | 12.79 | 11.79 | 1.05× → 1.08× |
| GB4 | 8.94 | 9.82 | 8.56 | 1.04× → 1.15× |
| Ecu1110 | 8.21 | 7.35 | 6.50 | 1.26× → 1.13× |
| 7G8 | 3.94 | 2.81 | 6.63 | 1.68× → 2.36× |

Exploratory observations (not verdicts): the transport model tracks K1 and GB4 closely,
mildly over-predicts Ecu1110, and **under-predicts 7G8** — and adding saturation (larger
`C*`) makes 7G8 **worse**, not better. So, tentatively, the earlier "maybe it was just
the linearization" idea is **not** supported here: 7G8 stays anomalous in both regimes.
But because the RIs are non-comparable, this cannot be promoted to a finding.

## Conclusion

- **Confirmatory #1b:** NOT EVALUABLE (phenotype heterogeneity). Frozen.
- **Next honest step:** a real #1b needs a phenotype set measured with a **single common
  CQS reference and assay** across Dd2, K1, GB4, Ecu1110, 7G8 (ideally isogenic lines).
  Until such data exist, the multi-strain transport→RI hypothesis stays untested.
- The `RMSE`/saturation criteria from the pre-registration are **not** applied, by rule.

Credit: Table S2 / Table 6.3 and the C\* derivation located by the maintainer in the
Summers PhD thesis (ANU open repository).
