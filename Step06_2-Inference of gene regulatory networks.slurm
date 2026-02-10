#!/bin/bash
#SBATCH --job-name=pyscenic 
#SBATCH --output=pyscenic_%j.out  
#SBATCH --error=pyscenic_%j.err   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24                
#SBATCH --mem=400G                       
#SBATCH --time=196:00:00                 
WORK_DIR=~/space/pyscenic/
cd $WORK_DIR

# Files below can be downloaded from https://resources.aertslab.org/cistarget/
TFS=~/genomeprogram/scenic/allTFs_hg38.txt
RANKINGS=~/genomeprogram/scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
MOTIFS=~/genomeprogram/scenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

N_WORKERS=24

LOOM_FILES=($(ls *.loom 2>/dev/null))

if [ ${#LOOM_FILES[@]} -eq 0 ]; then
  echo "Error: No *.loom files found in $WORK_DIR"
  exit 1
fi

echo "Found ${#LOOM_FILES[@]} loom files. Starting pySCENIC pipeline..."
echo "Files: ${LOOM_FILES[@]}"

for input_loom in "${LOOM_FILES[@]}"; do
  base_name=$(basename "$input_loom" .loom)
  cancer_dir="$WORK_DIR/$base_name"

  mkdir -p "$cancer_dir"
  cd "$cancer_dir"

  cp "$WORK_DIR/$input_loom" ./input.loom

  # Step 1: GRN (GRNboost2)
  pyscenic grn \
    --num_workers $N_WORKERS \
    --output adj.tsv \
    --method grnboost2 \
    input.loom $TFS

  # Step 2: cisTarget
  pyscenic ctx \
    adj.tsv \
    $RANKINGS \
    --annotations_fname $MOTIFS \
    --expression_mtx_fname input.loom \
    --mode "dask_multiprocessing" \
    --output regulons.csv \
    --num_workers $N_WORKERS \
    --mask_dropouts

  # Step 3: AUCell
  pyscenic aucell \
    input.loom \
    regulons.csv \
    --output "${base_name}_aucell.loom" \
    --num_workers $N_WORKERS

  echo "$input_loom completed! Results in $cancer_dir"
  cd "$WORK_DIR"
done

