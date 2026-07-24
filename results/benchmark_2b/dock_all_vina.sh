#!/bin/bash
# One-shot #2b — dock pyrimethamine into all 14 variant models with AutoDock-Vina-GPU-2.1.
cd /home/girias/projects/malaria_body/results/benchmark_2b
V=/home/girias/Documents/Pfpk2/tools/Vina-GPU-2.1/vina_gpu_wrap.sh
echo "build,variant,vina_affinity,strength" > scores_vina.csv
for f in models/*.pdb; do
  name=$(basename "$f" .pdb); build=${name%%_*}; variant=${name#*_}
  if [ "$build" = wt ]; then cx=28.035; cy=5.436; cz=58.756; else cx=28.290; cy=5.201; cz=58.803; fi
  obabel "$f" -O "vina/${name}_rec.pdbqt" -xr 2>/dev/null
  mkdir -p "vina/${name}_lig" "vina/${name}_out"; cp vina/cp6_lig.pdbqt "vina/${name}_lig/"
  cat > "vina/${name}_cfg.txt" <<EOF
receptor = $(pwd)/vina/${name}_rec.pdbqt
ligand_directory = $(pwd)/vina/${name}_lig
output_directory = $(pwd)/vina/${name}_out
center_x = $cx
center_y = $cy
center_z = $cz
size_x = 22.5
size_y = 22.5
size_z = 22.5
thread = 8000
seed = 42
EOF
  "$V" --config "vina/${name}_cfg.txt" > "vina/${name}_vina.log" 2>&1
  aff=$(grep -m1 "VINA RESULT" "vina/${name}_out/cp6_lig_out.pdbqt" 2>/dev/null | awk '{print $4}')
  strength=$(python3 -c "print(round(-1*float('$aff'),3))" 2>/dev/null)
  echo "$build,$variant,$aff,$strength" >> scores_vina.csv
  echo "$name: affinity=$aff strength=$strength"
done
echo "ALL_VINA_DONE"
