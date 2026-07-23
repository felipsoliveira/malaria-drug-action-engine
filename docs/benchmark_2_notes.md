# Slice 2 — PfDHFR + pyrimethamine: status & plan (Layer C)

Opening Layer C (molecular function) on a target with strong ground truth: crystal
structures (WT and mutants), a single-system kinetic panel, and a clear structural
resistance mechanism. PfPK2 comes later, once the architecture is validated here.

## Data incorporated (measured, one experimental system)

`mdae.pfdhfr_data.PANEL` — Sirawaraporn et al. 1997 PNAS, Table 3 (PMC19755): `kcat`,
`Km(DHF)`, `Ki(pyrimethamine)` for the resistance series (two evolutionary routes
converge on the quadruple mutant). Assay: [DHF]=100 µM, [NADPH]=100 µM, pH 7.0, 25 °C.

**Provenance caveats (frozen in the data module):**
- The Table-3 `Ki` were **derived from IC50** via Cheng-Prusoff `IC50 = Ki(1+[S]/Km)`.
  Reconstructing IC50 from these `Ki` is therefore an **algebraic round-trip** — a
  units/implementation/provenance check, **not** a biological validation.
- **Do not merge** with Sandefur et al. 2007 (PMC2020854): different construct; WT `Km`
  13 → 44.4 µM and `kcat` 88 → 1.24 /s. The two are not interchangeable as one panel.

## Done in this commit

- `enzyme.py` — Layer-C competitive-inhibition model (`v(I)`, residual activity,
  Cheng-Prusoff IC50, inverse). Synthetic-validated (5/5).
- Panel incorporated as `measured`/`derived` with provenance; round-trip **self-
  consistency** test (4/4). Explicitly **not** a confirmatory benchmark.
- Illustration only: the same `Ki` maps to different functional IC50 as `Km` varies
  (IC50/Ki spans 5–18×) — the Km-dependence Layer C encodes.

## What a FALSIFIABLE #2 still needs

- **#2a (binding→function, falsifiable):** freeze `Ki`, `Km`; predict **independent**
  data — an activity/IC50 curve at a **different** [DHF], or a **directly measured**
  Kd/Ki plus a **functional** IC50. The Cheng-Prusoff round-trip does not count.
- **#2b (docking ranking):** gnina/Vina rank these `Ki` as the experimental target,
  evaluated **separately** (rank correlation), **without** turning score → ΔG/Kd.
- **#2c (end-to-end held-out, pre-registered):** calibrate score→pKi on several
  compounds/variants, **freeze**, then predict an excluded ligand or mutant.

## Status

Layer-C machinery + a comparable measured panel are in place. No confirmatory claim is
made yet. #2b can start now (docking → Ki ranking); #2a/#2c need the independent
readout above before any PASS/FAIL.
