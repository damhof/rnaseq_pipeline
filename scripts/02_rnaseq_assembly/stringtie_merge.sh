#!/bin/bash

set -uo pipefail

# Get sample IDs
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Get gtfmergefile
gtfmergefile="${outdir}/stringtie/gtfmergefile.txt"

# # Check whether script needs to run
# if [[ -f "${outdir}/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf" ]]; then
#   echo "Merged GTF already present"
#   exit 0
# fi

# # Create output dirs
# mkdir -p "${outdir}/gffcompare/${merged_gtf_basename}/"
# mkdir -p "${outdir}/stringtie/${merged_gtf_basename}/"
# cd "${outdir}/stringtie/${merged_gtf_basename}/"

# echo "`date` running Stringtie --merge"

# # Run stringtie merge
# apptainer exec -B /hpc:/hpc ${container_dir}/stringtie-2.1.5.sif stringtie --version
# apptainer exec -B /hpc:/hpc ${container_dir}/stringtie-2.1.5.sif stringtie ${gtfmergefile} \
#   --merge \
#   -G ${reference_gtf} \
#   -f 0.05 \
#   -m 50 \
#   -T 0 \
#   -o "${outdir}/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf" \
#   -p 4

# echo "`date` finished Stringtie Merge"

# echo "`date` running GFFcompare for annotation"

# # Create output dirs
# cd "${outdir}/gffcompare/${merged_gtf_basename}/"

# # Run GFFcompare to annotate novel assembly .GTF
# apptainer exec -B /hpc:/hpc ${container_dir}/gffcompare-0.12.6.sif gffcompare --version
# apptainer exec -B /hpc:/hpc ${container_dir}/gffcompare-0.12.6.sif gffcompare \
#   -V \
#   -M \
#   -r ${reference_gtf} \
#   -s ${masked_fasta} \
#   -o "${merged_gtf_basename}_gffcompare_new" \
#   "${outdir}/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf"

echo "`date` running GFFcompare for transcript occurence"

# Run GFFcompare to merge stringtie outputs and generate .tracking file
# I'm testing this step to skip stringtie --merge to see if gffcompare direct merging works better
# Also test later with -D mode, as I'm still not sure what it does. Also try -C, -A, and -X.

> "${outdir}/gffcompare/gtflist.txt"  # Clear or create the file

for id in "${sample_ids[@]}"; do
  find "${outdir}/stringtie/${id}" -type f -name "${id}.gtf" >> "${outdir}/gffcompare/gtflist.txt"
done

apptainer exec -B /hpc:/hpc ${container_dir}/gffcompare-0.12.6.sif gffcompare \
  -V \
  -M \
  -r ${reference_gtf} \
  -o "${merged_gtf_basename}_gffcompare_merged" \
  -i "${outdir}/gffcompare/gtflist.txt"

echo "`date` finished running GFFcompare"