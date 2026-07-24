"""One-shot #2b — prep the 3 crystal controls per the frozen manifest (exec-manifest-2b).
Chain A only; keep NADPH (NDP); drop waters, UMP (TS domain), chain B. Save receptor
(protein+NDP, no CP6) and the crystal CP6 reference pose."""
from pymol import cmd
import os

base = "/home/girias/projects/malaria_body/results/benchmark_2b"
struct, outdir = f"{base}/structures", f"{base}/gate"
os.makedirs(outdir, exist_ok=True)

for pid in ["3QGT", "1J3J", "3QG2"]:
    cmd.reinitialize()
    cmd.load(f"{struct}/{pid}.pdb", "st")
    cmd.remove("solvent")
    cmd.remove("not chain A")     # keep chain A only
    cmd.remove("resn UMP")        # drop TS-domain ligand
    cmd.create("rec", "st and (polymer or resn NDP)")   # receptor = protein + NADPH
    cmd.save(f"{outdir}/{pid}_rec.pdb", "rec")
    cmd.create("ref", "st and resn CP6")                # crystal CP6 reference pose
    cmd.save(f"{outdir}/{pid}_cp6_ref.sdf", "ref")
    cmd.save(f"{outdir}/{pid}_cp6_ref.pdb", "ref")
    print(f"{pid}: rec_atoms={cmd.count_atoms('rec')} cp6_atoms={cmd.count_atoms('ref')} "
          f"ndp={cmd.count_atoms('rec and resn NDP')}")
