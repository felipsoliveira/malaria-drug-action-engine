"""One-shot #2b — build the 14 variant receptors per the frozen manifest.
Two builds (from WT 3QGT and from quadruple 3QG2). PyMOL mutagenesis (lowest-strain
rotamer) + sculpt relaxation of mutated sidechains + 5A neighbours, backbone & NDP fixed."""
from pymol import cmd
import os

base = "/home/girias/projects/malaria_body/results/benchmark_2b"
gate, outdir = f"{base}/gate", f"{base}/models"
os.makedirs(outdir, exist_ok=True)

VARIANTS = {
 "WT":                    {51:"ASN",59:"CYS",108:"SER",164:"ILE"},
 "S108N":                 {51:"ASN",59:"CYS",108:"ASN",164:"ILE"},
 "N51I_S108N":            {51:"ILE",59:"CYS",108:"ASN",164:"ILE"},
 "C59R_S108N":            {51:"ASN",59:"ARG",108:"ASN",164:"ILE"},
 "N51I_C59R_S108N":       {51:"ILE",59:"ARG",108:"ASN",164:"ILE"},
 "C59R_S108N_I164L":      {51:"ASN",59:"ARG",108:"ASN",164:"LEU"},
 "N51I_C59R_S108N_I164L": {51:"ILE",59:"ARG",108:"ASN",164:"LEU"},
}
TEMPLATES = {"wt": f"{gate}/3QGT_rec.pdb", "quad": f"{gate}/3QG2_rec.pdb"}


def current_resn(resi):
    r = []
    cmd.iterate(f"chain A and resi {resi} and name CA", "r.append(resn)", space={"r": r})
    return r[0] if r else None


def mutate(resi, target):
    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()
    w = cmd.get_wizard()
    w.set_mode(target)
    w.do_select(f"chain A and resi {resi}")
    try:
        bumps = list(w.bump_scores)
        best = bumps.index(min(bumps)) + 1 if bumps else 1
    except Exception:
        best = 1
    cmd.frame(best)
    w.apply()
    cmd.set_wizard()


for build, tmpl in TEMPLATES.items():
    for vname, targets in VARIANTS.items():
        cmd.reinitialize()
        cmd.load(tmpl, "mol")
        muts = []
        for resi, tgt in targets.items():
            cur = current_resn(resi)
            if cur != tgt:
                mutate(resi, tgt)
                muts.append(f"{cur}{resi}{tgt}")
        if muts:
            resis = "+".join(str(r) for r in targets)
            cmd.flag("fix", "all", "set")
            cmd.flag("fix", f"(byres (chain A and resi {resis} around 5) or (chain A and resi {resis})) and sidechain", "clear")
            cmd.sculpt_activate("mol")
            cmd.sculpt_iterate("mol", cycles=100)
            cmd.sculpt_deactivate("mol")
            cmd.flag("fix", "all", "clear")
        out = f"{outdir}/{build}_{vname}.pdb"
        cmd.save(out, "mol")
        ok = all(current_resn(r) == t for r, t in targets.items())
        ndp = cmd.count_atoms("resn NDP")
        print(f"{build}/{vname}: muts={muts or 'none'} verify={'OK' if ok else 'FAIL'} ndp={ndp}")
