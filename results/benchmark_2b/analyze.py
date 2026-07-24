"""One-shot #2b — final analysis. Spearman rho(strength, pKi), n=7, per tool per build.
PASS iff rho>=0.786 (two-sided exact p<0.05 at n=7); clean pass = both builds."""
import gzip, glob, os
from rdkit import Chem
from scipy.stats import spearmanr

pKi = {"WT":8.824,"S108N":7.886,"N51I_S108N":7.432,"C59R_S108N":7.143,
       "N51I_C59R_S108N":6.921,"C59R_S108N_I164L":6.417,"N51I_C59R_S108N_I164L":6.066}
order = list(pKi.keys())

def gnina_cnnaff(path):
    for m in Chem.ForwardSDMolSupplier(gzip.open(path)):
        if m is not None:
            return float(m.GetProp("CNNaffinity"))   # pose 1 = top-ranked
    return None

gnina = {"wt":{}, "quad":{}}
for path in glob.glob("models/*_gnina_out.sdf.gz"):
    name = os.path.basename(path).replace("_gnina_out.sdf.gz","")
    build = name.split("_")[0]; variant = name[len(build)+1:]
    gnina[build][variant] = gnina_cnnaff(path)

vina = {"wt":{}, "quad":{}}
for line in open("scores_vina.csv").read().splitlines()[1:]:
    build,variant,aff,strength = line.split(",")
    vina[build][variant] = float(strength)

def report(tool, build, scores):
    s=[scores[v] for v in order]; p=[pKi[v] for v in order]
    rho,pval = spearmanr(s,p)
    print(f"\n{tool} / {build}-build:  Spearman rho = {rho:+.3f} (p={pval:.3f}) -> "
          f"{'PASS' if rho>=0.786 else 'FAIL'}")
    print(f"  {'variant':24s} {'strength':>9s} {'pKi':>6s}")
    for v in order:
        print(f"  {v:24s} {scores[v]:+9.3f} {pKi[v]:6.3f}")
    return rho

print("="*64); print("Benchmark #2b RESULT — does docking order the variant Ki?"); print("="*64)
gw=report("GNINA(CNNaffinity)","wt",gnina["wt"]); gq=report("GNINA(CNNaffinity)","quad",gnina["quad"])
vw=report("Vina(-kcal)","wt",vina["wt"]); vq=report("Vina(-kcal)","quad",vina["quad"])
print("\n"+"-"*64)
print(f"GNINA clean pass (both builds rho>=0.786)? {gw>=0.786 and gq>=0.786}")
print(f"Vina  clean pass (both builds rho>=0.786)? {vw>=0.786 and vq>=0.786}")
