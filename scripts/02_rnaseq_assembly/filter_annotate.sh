#!/bin/bash

# Check whether script needs to run
if [[ -f "${outdir}/customannotation/${merged_gtf_basename}_novel_filtered.gtf" ]]; then
  echo "Merged filtered annotated GTF already present"
  exit 0
fi

# Create output dirs
mkdir -p "${outdir}/customannotation"

echo "`date` Filter and annotate novel GTF"

# Process and filter novel GTF
apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/r_rna_filter-4.1.2.sif Rscript "${scriptdir}/02_rnaseq_assembly/filter_annotate_v2.R" \
  "${scriptdir}/02_rnaseq_assembly/filter_annotate_functions.R" \
  "${merged_basename}" \
  "${reference_gtf}" \
  "${refseq_gtf}" \
  "${outdir}/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}.annotated.gtf" \
  "${outdir}/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}_combined.combined.gtf" \
  "${outdir}/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}.tracking" \
  "${outdir}/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}_combined.tracking" \
  "${min_occurrence}" \
  "${outdir}/customannotation/${merged_gtf_basename}_novel_filtered_4_stringtie.gtf"

echo "`date` Finished GTF filtering"