# Benchmark #2b — PRE-REGISTRATION: does docking detect variant-caused resistance?

**Status: PRE-REGISTERED.** Committed and tagged (`prereg-2b`) **before the first
docking score**. Structures, model-building protocol, metric, and pass/fail criteria
below are **frozen**; they will not change after seeing any score. Runs happen once.

## Question (correctly scoped)
Not "docking predicts affinity of any molecule". We have **one ligand (pyrimethamine)
against seven PfDHFR variants** → this tests **sensitivity to resistance mutations**:

> Does the docking score **order** the seven experimentally-derived pyrimethamine `Ki`
> across PfDHFR variants?

The conclusion will be phrased exactly that way — **not** "docking computed a binding
energy". The `Ki` are derived from a functional competitive assay (Sirawaraporn 1997);
the target is their `pKi`, not a thermodynamic ΔG.

## Target (frozen)
`pKi = −log10(Ki [M])` for the 7 variants in `mdae.pfdhfr_data.PANEL` (Ki 1.5→859 nM).

## Structures (frozen)
Pyrimethamine + NADPH co-crystals, 2.30 Å (RCSB):
- **3QGT** — WT.
- **1J3J** — C59R+S108N (double).
- **3QG2** — N51I+C59R+S108N+I164L (quadruple).
- **1J3K excluded** from the primary benchmark (it holds **WR99210**, not pyrimethamine).

## Protocol (frozen)

1. **Structural gate (cognate re-dock).** Both GNINA and Vina must reproduce the
   crystallographic **pyrimethamine** pose in WT (3QGT), double (1J3J) and quadruple
   (3QG2) with **heavy-atom RMSD ≤ 2.0 Å**. If a tool fails the gate, its cross-variant
   ranking is declared **"non-interpretable"** (no PASS/FAIL for that tool).

2. **Comparable models (no crystal-to-crystal confound).** The 7 variants are built from
   **one common template** with **identical** preparation: protonation at pH 7.0, NADPH
   present, identical box/search settings, identical relaxation. We do **not** use a
   different crystal per variant in the primary ranking (that would mix the mutation
   with crystal-to-crystal differences).

3. **Template sensitivity (two builds).** Repeat the whole pipeline built from:
   - the **WT** template (3QGT); and
   - the **quadruple** template (3QG2), introducing the reversions.
   The experimental WT / double / quadruple structures serve as **external controls** of
   the modeled variants.

4. **Primary metric.** Spearman `ρ(predicted binding strength, pKi)`, n = 7. Sign
   convention fixed now: "binding strength" is defined so a stronger predicted binder is
   **larger** (GNINA: CNNaffinity; Vina: −affinity[kcal/mol]), so a working model gives
   **positive** ρ. **PASS iff ρ ≥ 0.786** — the smallest discrete Spearman value reaching
   two-sided exact **p < 0.05** at n = 7.

5. **No cherry-picking.** **GNINA is primary, Vina is secondary**; each receives its own
   PASS / FAIL / non-interpretable, and **both are reported** regardless of outcome. No
   "report only what worked".

6. **Both builds must agree for a clean PASS.** A tool PASSES only if ρ ≥ 0.786 in **both**
   template builds (WT-based and quadruple-based); passing one build only is reported as
   **template-dependent** (weaker).

## Interpretation (frozen)
- **PASS:** the static docking score orders the variant `Ki` → concrete justification to
  proceed to a frozen score→`pKi` calibration (#2c).
- **FAIL:** static docking does **not** capture variant resistance for this target → an
  honest, publishable negative result; #2c is not attempted on this basis.
- **Non-interpretable (gate fail):** the tool cannot even reproduce the cognate pose, so
  its ranking says nothing.

## Scope note for #2c (recorded now, not changed later)
With this panel (one ligand × 7 variants), a future #2c can predict an **excluded
variant**, but **not an excluded ligand** — we only have pyrimethamine. Demonstrating
**cross-molecule transfer** will require a matrix of several inhibitors × variants
(future data, to be sourced with the same comparability audit as #1b/#2).

## Execution
Docking tools are present (`/home/girias/Documents/Pfpk2/tools/{gnina,Vina-GPU-2.1}`).
Run once after this file is tagged. Report GNINA and Vina, both template builds, the
structural-gate RMSDs, and the two Spearman ρ per tool — pass or fail.
