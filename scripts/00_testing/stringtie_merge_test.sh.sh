#!/bin/bash

#SBATCH --job-name=merge_bams
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4

set -uo pipefail
wd="/hpc/pmc_vanheesch/projects/Damon/Neuroblastoma_neoantigens/RNA_nbl_transcriptome_assembly"

source ${wd}/scripts/NBL_full.config

gtfmergefile="${outdir}/stringtie_test/gtflist.txt"

echo "`date` running Stringtie --merge"

# Run stringtie merge
apptainer exec -B /hpc:/hpc ${container_dir}/stringtie-2.1.5.sif stringtie --version
apptainer exec -B /hpc:/hpc ${container_dir}/stringtie-2.1.5.sif stringtie ${gtfmergefile} \
  --merge \
  -G ${reference_gtf} \
  -f 0.05 \
  -m 200 \
  -T 1 \
  -o "${outdir}/stringtie_test/${merged_gtf_basename}/${merged_gtf_basename}_m200_T1_f005_2.gtf" \
  -p 4

echo "`date` finished Stringtie Merge"
