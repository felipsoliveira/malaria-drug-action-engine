# Benchmark #2b — EXECUTION BLOCKED (GPU driver/library mismatch)

**Feasibility check 2026-07-24. Verdict: BLOCKED before the first score.** No PfDHFR
docking was performed and no score was generated. The pre-registration `prereg-2b` is
untouched and remains valid; the one-shot resumes once the GPU is fixed.

## The blocker (hard)
NVIDIA userspace/kernel versions mismatch:
- kernel module (NVRM): **595.71.05**
- NVML library: **595.84**

`nvidia-smi` → *"Failed to initialize NVML: Driver/library version mismatch"*. The
CUDA/GPU runtime is unavailable, which blocks **both** pre-registered scorers:
- **GNINA CNNaffinity** (the frozen *primary* score field) — needs the GPU CNN;
- **AutoDock-Vina-GPU-2.1** (the frozen *secondary*) — needs the GPU (OpenCL).

Per the pre-registered rule ("faltar … CUDA → para antes do primeiro score e registra o
bloqueio"), execution stops here. We do **not** substitute CPU-only empirical scoring:
that would change the frozen primary field (operational freedom the protocol forbids).

**To unblock:** repair the driver/library mismatch (reboot into, or reinstall, the
NVIDIA userspace libraries matching kernel **595.71.05**), then re-verify `nvidia-smi`
and a GNINA CNN smoke test on a *bundled* example before resuming.

## Everything else is ready (inventory)
| Component | Status |
|---|---|
| GNINA binary | **v1.3.2** (master f23dd2b, built 2025-07-29); loads via `gnina_wrap.sh` (libcudnn.so.9 in `miniforge3/lib`) |
| Vina | `AutoDock-Vina-GPU-2-1` binary + `libOpenCL.so` present (env `ai_md_cuda`) — confirmed = Vina-GPU-2.1 |
| Prep: obabel, reduce, pymol, rdkit 2026.03.3 | present |
| Prep: meeko / MGLTools (`prepare_receptor4.py`) | **absent** — substitutable (obabel → PDBQT; reduce/obabel → protonation; PyMOL → mutations) |
| GPU runtime | **BLOCKED** (driver mismatch) |

## Facts captured now (head-start for the future execution manifest)
Structures (RCSB, downloaded, SHA256):
- **3QGT** (WT): `df2701ca377a39eeff82e96ce50a352599f97d14972ea07bb3ea131f910646a1`
- **1J3J** (C59R+S108N): `2faf10870d0130979625b141e90f6916627f16410691ace1fc7e53ce89b0de23`
- **3QG2** (quadruple): `d90a5bcd4ce8fbb96f2a8d8ff469df0dae23bb1072f2c9545db13b6d998b9597`

Cognate ligand in all three = **CP6** (pyrimethamine); cofactor **NDP** (NADPH — keep);
**UMP** present (thymidylate-synthase domain of the bifunctional PfDHFR-TS). The RMSD
gate will target **CP6**.

## Not done (by rule)
- No execution manifest frozen (that happens only when feasible).
- No score, no ranking, no RMSD gate applied. `prereg-2b` stays the sole frozen artifact.
- When the GPU is repaired: re-run this feasibility check, then freeze the full execution
  manifest (binary hashes, chain, NADPH handling, CP6 protonation at pH 7, mutation tool/
  rotamer/relax, box/seed/exhaustiveness/poses, primary GNINA field, RMSD symmetry mapping)
  **before** the first score, and run the scripted one-shot.
