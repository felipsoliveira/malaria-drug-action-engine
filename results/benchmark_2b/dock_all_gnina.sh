#!/bin/bash
# One-shot #2b — dock pyrimethamine (N1 monocation) into all 14 variant models with GNINA.
cd /home/girias/projects/malaria_body/results/benchmark_2b
G=/home/girias/Documents/Pfpk2/tools/gnina/gnina_wrap.sh
for f in models/*.pdb; do
  name=$(basename "$f" .pdb); build=${name%%_*}
  if [ "$build" = wt ]; then ref=gate/3QGT_cp6_ref.pdb; else ref=gate/3QG2_cp6_ref.pdb; fi
  "$G" -r "$f" -l gate/cp6_lig_pH7.sdf --autobox_ligand "$ref" --autobox_add 4 \
    --exhaustiveness 16 --num_modes 9 --seed 42 --cnn_scoring rescore \
    -o "models/${name}_gnina_out.sdf.gz" > "models/${name}_gnina.log" 2>&1
  echo "done $name"
done
echo "ALL_GNINA_DONE"
