# Benchmark #2b — RESULT: static docking does NOT detect variant resistance (FAIL)

**Pre-registered outcome (frozen).** Executed once under `exec-manifest-2b` (which
implements `prereg-2b`). Structural gate **passed** for both tools, so the rankings are
**interpretable**; the cross-variant ranking **FAILS** for both GNINA and AutoDock-Vina-
GPU-2.1, in both template builds. Conclusion, phrased as pre-registered:

> The docking **score** does **not** order the seven experimentally-derived pyrimethamine
> `Ki` across PfDHFR variants — i.e. **static docking does not detect the resistance
> caused by these mutations**.

This is a clean, publishable **negative** result. Both tools and both builds are reported
(no cherry-picking).

## Structural gate (cognate re-dock, symmetry RMSD via obrms) — PASS
| Control | GNINA RMSD (Å) | Vina RMSD (Å) |
|---|---:|---:|
| 3QGT (WT) | 0.56 | 0.38 |
| 1J3J (C59R+S108N) | 0.31 | 0.30 |
| 3QG2 (quadruple) | 0.33 | 0.26 |

All ≤ 2.0 Å → both tools reproduce the cognate pyrimethamine pose; the ranking is not a
preparation artifact.

## Rankings — Spearman ρ(predicted strength, pKi), n = 7
PASS iff ρ ≥ 0.786 (two-sided exact p<0.05 at n=7); clean pass = **both** builds.

| Scorer (strength field) | wt-build ρ | quad-build ρ | Verdict |
|---|---:|---:|---|
| **GNINA** (CNNaffinity) | +0.393 (p=0.38) | −0.071 (p=0.88) | **FAIL** |
| **Vina-GPU-2.1** (−kcal/mol) | +0.579 (p=0.17) | −0.116 (p=0.81) | **FAIL** |

pKi target (from `mdae.pfdhfr_data`, Ki 1.5→859 nM): WT 8.82, S108N 7.89,
N51I+S108N 7.43, C59R+S108N 7.14, N51I+C59R+S108N 6.92, C59R+S108N+I164L 6.42,
quadruple 6.07.

## Interpretation (careful, not over-claimed)
- Both tools correctly rank **WT as the strongest binder** in the wt-build (right
  direction), but the **six mutants are near-flat** (Vina 8.6–8.7; GNINA 6.2–6.5): the
  ~570× resistance gradient in `Ki` is **not resolved** by the score.
- In the quad-build (WT reached by four reversions), the correlation **vanishes/inverts**.
- So a **static** docking score — empirical (Vina) or CNN (GNINA) — does not capture the
  affinity shifts these resistance mutations produce for pyrimethamine at PfDHFR. This is
  a known limitation of scoring functions for close point mutants; it motivates
  **free-energy / MD** methods (e.g. alchemical ΔΔG) to capture mutation effects — the
  static score is insufficient.
- The gate rules out "the pose is wrong": the failure is in the **scoring/ranking**, not
  the docking geometry.

## Consequence for #2c (per pre-registration)
`#2c` (freeze a score→pKi calibration and predict an excluded variant) is **not attempted
on this basis**: there is no justification that score→pKi transfers when the score does
not even order the variants. A meaningful #2c would need a scoring method that first
passes #2b (e.g. free-energy).

## Reproducibility
- Manifest: `exec-manifest-2b` (binaries+hashes, chain A, NADPH kept, CP6 N1 monocation,
  PyMOL mutagenesis+sculpt, box, seeds 42, CNNaffinity/−kcal, obrms gate).
- Pipeline scripts + score CSVs are committed alongside this doc; heavy structures/poses
  live under `results/benchmark_2b/` (gitignored).
- Note captured during execution: `obabel -p 7` returned the **neutral** pyrimethamine;
  the frozen manifest specifies the **N1 monocation**, which was therefore built
  explicitly (RDKit) rather than accepting the tool default.
