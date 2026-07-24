# Benchmark #2b — EXECUTION MANIFEST (frozen before the first PfDHFR score)

Feasibility is GREEN (2026-07-24, after reboot): RTX 5090, driver 595.84, kernel
595.84 (matched); `nvidia-smi` OK. Smoke tests on **bundled/synthetic** inputs (not
PfDHFR): GNINA CNNaffinity ran (synthetic pair) and AutoDock-Vina-GPU-2.1 ran the
bundled 2bm2 example on the GPU. This manifest freezes every operational choice so
there is no freedom after seeing scores. It does not alter the scientific
pre-registration (`prereg-2b`). Tag: `exec-manifest-2b`.

## Binaries (versions + SHA256)
- **GNINA** v1.3.2 (master f23dd2b, built 2025-07-29); binary
  `714ef2928a22c7b20680ccfeade3c2a652a4e452a1f64343da541434528ed9cd`; run via `gnina_wrap.sh`.
- **Vina = AutoDock-Vina-GPU-2.1** (confirmed exactly this, not CPU Vina); binary
  `9d27a4e88f3dd597297168d74943e7f8a151d7986ffb118f6b30b9bd349826b1`; run via `vina_gpu_wrap.sh`.

## Structures (RCSB, SHA256) — chain A used
- **3QGT** (WT) `df2701ca377a39eeff82e96ce50a352599f97d14972ea07bb3ea131f910646a1`
- **1J3J** (C59R+S108N, double) `2faf10870d0130979625b141e90f6916627f16410691ace1fc7e53ce89b0de23`
- **3QG2** (quadruple) `d90a5bcd4ce8fbb96f2a8d8ff469df0dae23bb1072f2c9545db13b6d998b9597`
- Each has CP6 (pyrimethamine) + NDP (NADPH) in chains A and B. **Use chain A**; drop
  chain B, waters, UMP, and the TS domain. Cognate ligand = **CP6**.

## Receptor preparation (identical for all)
- Keep chain A protein + **NDP (NADPH) retained in the rigid receptor**; remove waters/UMP/chain B.
- Protonate at **pH 7.0** with `reduce`; GNINA takes the PDB receptor directly; for Vina,
  convert to a **rigid** PDBQT with `obabel -xr` (Gasteiger). NADPH kept in both.

## Ligand — pyrimethamine (CP6) at pH 7
- Reference pose = crystal **CP6, chain A**, per template.
- Docking input = CP6 re-protonated to the **monocation** (2,4-diaminopyrimidine **N1
  protonated** — the canonical DHFR-binding cationic form that salt-bridges Asp54),
  Gasteiger charges; identical ligand file for every variant.

## Mutant models — two builds (PyMOL)
- Tool: **PyMOL mutagenesis wizard**, rotamer = **lowest strain**; then **local
  relaxation of mutated residues + neighbours within 5 Å via PyMOL `sculpt` (100
  cycles), backbone restrained**. Same protocol for all 14 models.
- **WT build:** template **3QGT** chain A; WT = as-is; introduce mutations for
  S108N / N51I+S108N / C59R+S108N / N51I+C59R+S108N / C59R+S108N+I164L / quadruple.
- **Quadruple build:** template **3QG2** chain A; quadruple = as-is; introduce
  reversions to reach each variant.
- Experimental **WT (3QGT), double (1J3J), quadruple (3QG2)** = external controls of the models.

## Docking box (same within a build)
- Centred on the crystallographic **CP6 centroid** of the template (chain A); cubic
  **22.5 Å** per side. GNINA `--autobox_ligand <CP6> --autobox_add 4`; Vina center =
  CP6 centroid, size 22.5³.

## Search settings (fixed seeds)
- GNINA: `--cnn_scoring rescore --exhaustiveness 16 --num_modes 9 --seed 42`.
- Vina-GPU-2.1: `thread 8000 --seed 42`, default search depth, `num_modes 9`.

## Score fields (frozen)
- **GNINA primary = CNNaffinity** (predicted pKd-like; higher = stronger).
- **Vina = −affinity [kcal/mol]** (higher = stronger).
- Per variant: the **top-ranked pose's** value.

## Structural gate (cognate re-dock)
- Dock CP6 into its own crystal (3QGT, 1J3J, 3QG2); top pose; **symmetry/automorphism-aware
  heavy-atom RMSD to crystal CP6 via `obrms`**. **PASS iff RMSD ≤ 2.0 Å** for that structure.
- A tool that clears all three controls → its cross-variant ranking is interpretable;
  any control failing → that tool's ranking is **non-interpretable** (no PASS/FAIL).

## Primary metric & pass rule (frozen)
- Spearman **ρ(strength, pKi)**, n = 7; `pKi = −log10(Ki[M])` from `mdae.pfdhfr_data`:
  WT 8.82, S108N 7.89, N51I+S108N 7.43, C59R+S108N 7.14, N51I+C59R+S108N 6.92,
  C59R+S108N+I164L 6.42, quadruple 6.07.
- **PASS iff ρ ≥ 0.786** (two-sided exact p<0.05 at n=7); **clean pass = both builds**.
- **GNINA primary, Vina secondary**; each reported with its own PASS / FAIL /
  non-interpretable. **No cherry-picking; both tools + both builds always reported.**

## Determinism & outputs
- Fixed seeds; single scripted execution. **All poses/logs preserved even if the gate
  fails.** Structures/intermediates live under `results/benchmark_2b/` (gitignored);
  only the summary/result doc is committed.

## Scope note for #2c (unchanged)
One ligand × 7 variants → a future #2c can predict an excluded **variant**, not an
excluded **ligand**; cross-molecule transfer needs a multi-inhibitor × variant matrix.
