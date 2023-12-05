#!/bin/bash

# run_rnaseq_transcriptome_assembly.sh
#
# Short description of pipeline script
#
#

set -uo pipefail

function usage() {
    cat <<EOF
SYNOPSIS
  run_rnaseq_transcriptome_assembly.sh [-c <config file>] [-h]
DESCRIPTION
  1. Assemble transcriptomes for each sample with STRINGTIE, using BAM files as input
  2. Merge GTF files per sample group with gffcompare
  3. Annotate, filter, and fix merged GTF with custom R scripts
  4. Create custom annotation with custom merged GTF for ribo-seq analysis
  5. Generate Salmon index
  6. Quantify reads with salmon quant
  7. Run MultiQC on new data
OPTIONS
  -c, --config <file>    Configuration file to use
  -h, --help             Display this help message
AUTHOR
  Jip van Dinter, MSc
  Damon Hofman, MSc
EOF
}

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -c|--config)
    CONFIG="$2"
    shift
    shift
    ;;
    -h|--help)
    usage
    exit
    ;;
    "")
    echo "Error: no option provided"
    usage
    exit 1
    ;;
    *)
    echo "Unknown option: $key"
    usage
    exit 1
    ;;
esac
done

# Check that configuration file is provided
if [[ -z ${CONFIG+x} ]]; then 
    echo "Error: no configuration file provided"
    usage
    exit 1
fi

# Load configuration variables
source $CONFIG

# Load general functions
source ${scriptdir}/general_functions.sh

# Find samples
echo "$(date '+%Y-%m-%d %H:%M:%S') Finding samples..."
get_samples $project_data_folder $data_folder $paired_end

# Create a unique prefix for the names for this run of the pipeline. 
# This makes sure that runs can be identified
run_id=$(uuidgen | tr '-' ' ' | awk '{print $1}')

# Create output directories
mkdir -p ${project_folder}/log/${run_id}/{gffcompare,stringtie,salmon_quant}

################################################################################

# Step 1: Transcript assembly with stringtie
stringtie_jobid=()
stringtie_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.stringtie \
  --output=${project_folder}/log/${run_id}/stringtie/%A_%a.out \
  --export=ALL\
  ${scriptdir}/02_rnaseq_assembly/stringtie.sh
))
info "stringtie jobid: ${stringtie_jobid[@]}"

# Step 2: Merge all sample GTF files into single GTF file using gffcompare with reference GTF file as guide
gffcompare_merge_jobid=()
gffcompare_merge_jobid+=($(sbatch --parsable \
  --mem=48G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --job-name=${run_id}.gffcompare_merge \
  --output=${project_folder}/log/${run_id}/%A_gffcompare_merge.out \
  --dependency=afterany:${stringtie_jobid} \
  --export=ALL \
  ${scriptdir}/02_rnaseq_assembly/gffcompare_merge.sh
))
info "Gffcompare merge jobid: ${gffcompare_merge_jobid[@]}"

# Step 3: Fix and filter the merged gffcompare output file
filter_annotate_jobid=()
filter_annotate_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=2 \
  --time=24:00:00 \
  --job-name=${run_id}.filter_annotate \
  --output=${project_folder}/log/${run_id}/%A_filter_annotate.out \
  --export=ALL \
  ${scriptdir}/02_rnaseq_assembly/filter_annotate.sh
))
info "Filter and annotation jobid: ${filter_annotate_jobid[@]}"

# Step 4: Create custom annotation file for RiboseQC and ORFquant based on filtered custom GTF
custom_annotation_jobid=()
if [[ ${create_annotation} =~ "TRUE" ]]; then
  if [[ $(find ${wd}/ -name '*.Rannot' | wc -l) -eq 0 ]]; then
    custom_annotation_jobid+=($(sbatch --parsable \
      --mem=10G \
      --cpus-per-task=2 \
      --gres=tmpspace:50G \
      --time=24:00:00 \
      --job-name=${run_id}.custom_annotation \
      --output=${project_folder}/log/${run_id}/%A_custom_annotation.out \
      --export=ALL \
      ${scriptdir}/02_rnaseq_assembly/custom_annotation.sh
    ))
    info "Custom annotation jobid: ${custom_annotation_jobid[@]}"
  else
  echo "Annotation file already present"
  fi
else
  echo "Creation of annotation not specified"
fi