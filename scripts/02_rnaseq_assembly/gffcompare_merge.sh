#!/bin/bash

set -uo pipefail

# Get sample IDs
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Get list of GTF files
> "${outdir}/stringtie/gtflist.txt"  # Clear or create the file
for id in "${sample_ids[@]}"; do
  find "${outdir}/stringtie/${id}" -type f -name "${id}.gtf" >> "${outdir}/stringtie/gtflist.txt"
done

echo "`date` running gffcompare"
mkdir -p "${outdir}/gffcompare/${merged_gtf_basename}"

# Run stringtie merge
apptainer exec -B /hpc:/hpc ${container_dir}/gffcompare-0.12.6.sif gffcompare \
  -V \
  -r "${reference_gtf}" \
  -s ${masked_fasta} \
  -o "${outdir}/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}" \
  -i "${outdir}/stringtie/gtflist.txt"    