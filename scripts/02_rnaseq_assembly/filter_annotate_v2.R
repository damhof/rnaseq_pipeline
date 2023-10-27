# Input requirements: =========================================================
# 1. gtf_annot_file: generated during stringtie_merge.sh step, using gffcompare with stringtie --merge output GTF as input and ensembl ref GTF as reference
# 2. gtf_novel_file: generated during stringtie_merge.sh step, using gffcompare with individual sample StringTie GTFs as input and novel annotated GTF (input 1.) as reference
# 3. orig_tracking_file: generated in same way as gtf_annot_file
# 4. combined_tracking_file: generated in same way as gtf_novel_file
# 5. min_occurrence: minimum number of samples in which a transcript must be found to be included in the final GTF
# 6. outfile: output file name
#
# Output: =====================================================================
# 1. outfile: custom annotated GTF

# Load required packages
message(paste(Sys.time(), "Loading required libraries ..."), sep = "\t")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
})

# Global arguments
args <- commandArgs(trailingOnly = TRUE)

functions_file = args[1]
merged_basename = args[2]
gtf_ref_loc = args[3]
gtf_refseq_basename = args[4]
gtf_annot_file = args[5]
gtf_novel_file = args[6]
orig_tracking_file = args[7]
combined_tracking_file = args[8]
min_occurrence = args[9]
outfile = args[10]

min_occurrence <- as.numeric(min_occurrence)

# Load filtering functions
source(functions_file)


# Load files
message(paste(Sys.time(), "Loading files ..."), sep = "\t")

## Load annotated novel GTF file
gtf_annot <- rtracklayer::import(gtf_annot_file)
gtf_annot <- data.table::as.data.table(gtf_annot)

## Load combined GTF file (generated with gffcompare, using individual sample StringTie GTFs as input and novel annotated GTF as reference)
gtf_novel_GR <- rtracklayer::import(gtf_novel_file)

## Load original tracking file
orig_tracking <- data.table::fread(orig_tracking_file, header = F, select = c(1:5)) %>%
  tidyr::separate(V3, into = c("ref_gene_id_tracking", "ref_transcript_id"), sep = "\\|") %>%
  tidyr::separate(V5, into = c("query_gene_id", "query_transcript_id", "query_num_exons", "query_FPKM", "query_TPM", "query_cov", "query_len"), sep = "\\|")
colnames(orig_tracking) <- c("TCONS", "xloc", "ref_gene_id_orig", "ref_transcript_id", "class_code_orig", "query_gene_id", "query_transcript_id", "query_num_exons", "query_FPKM", "query_TPM", "query_cov", "query_len")
orig_tracking$query_gene_id <- sub("^q1:", "", orig_tracking$query_gene_id)

### Merge orig_tracking_file with annotated GTF (generated with gffcompare using stringtie --merge GTF as input and ref GTF as reference)
orig_tracking_annot <- orig_tracking %>%
  left_join(gtf_annot[, c("gene_name", "ref_gene_id", "cmp_ref", "cmp_ref_gene", "transcript_id")], by = c("query_transcript_id" = "transcript_id"))

### Remove redundant lines from orig_tracking_annot
orig_tracking_annot <- orig_tracking_annot %>%
  group_by(TCONS) %>%
  mutate(
    gene_name = first(gene_name[gene_name != "NA"]),
    ref_gene_id = first(ref_gene_id[ref_gene_id != "NA"]),
    cmp_ref = first(cmp_ref[cmp_ref != "NA"]),
    cmp_ref_gene = first(cmp_ref_gene[cmp_ref_gene != "NA"]),
  ) %>%
  ungroup() %>%
  distinct()

## Load combined gffcompare tracking file
combined_tracking <- data.table::fread(combined_tracking_file, header = F, select = c(1:5)) %>%
  tidyr::separate(V3, into = c("combined_gene_id", "ref_transcript_id"), sep = "\\|") %>%
  tidyr::separate(V5, into = c("query_gene_id", "query_transcript_id", "query_num_exons", "query_FPKM", "query_TPM", "query_cov", "query_len"), sep = "\\|")
colnames(combined_tracking) <- c("TCONS", "xloc", "combined_gene_id", "combined_transcript_id", "class_code_combined", "query_gene_id", "query_transcript_id", "query_num_exons", "query_FPKM", "query_TPM", "query_cov", "query_len")

## Merge combined_tracking file with orig_tracking file to annotate combined GTF
combined_tracking_annotated <- combined_tracking[, c("TCONS", "combined_gene_id", "combined_transcript_id", "class_code_combined")] %>%
  left_join(orig_tracking_annot[, c("ref_gene_id", "ref_transcript_id", "query_gene_id", "query_transcript_id", "class_code_orig", "gene_name", "cmp_ref", "cmp_ref_gene")], na_matches = "never", by = c("combined_transcript_id" = "query_transcript_id"))

## Load reference GTF
gtf_reference <- rtracklayer::import.gff("/hpc/pmc_vanheesch/shared_resources/GENOMES/Homo_sapiens.GRCh38/102/annotation/Homo_sapiens.GRCh38.102.gtf", colnames = c(
  "type",
  "source",
  "gene_id",
  "gene_name",
  "gene_biotype",
  "transcript_id",
  "transcript_name"
))
gtf_ref_df <- as.data.frame(gtf_reference)


# Filtering and annotation #####################################################
message(paste(Sys.time(), "Removing transcripts on scaffolds ..."), sep = "\t")
gtf_novel_GR <- gtf_novel_GR[seqnames(gtf_novel_GR) %in% c(1:22, "X", "Y"), ]

message(paste(Sys.time(), "Removing unstranded transcripts ..."), sep = "\t")
gtf_novel_GR <- gtf_novel_GR[strand(gtf_novel_GR) != "*", ]

gtf_novel <- data.table::as.data.table(gtf_novel_GR)
gtf_novel$num_samples <- as.numeric(gtf_novel$num_samples)

message(paste0(Sys.time(), "\t", "Filtering transcripts found in fewer than ", min_occurrence, " samples ..."))
transcripts_keep <- subset(gtf_novel, type == "transcript" & num_samples >= min_occurrence)$transcript_id
gtf_novel <- gtf_novel[which(gtf_novel$transcript_id %in% transcripts_keep),]

## Set exon cmp_ref IDs and transcript_ids same as parent transcript
gtf_novel <- gtf_novel %>%
  group_by(transcript_id) %>%
  mutate(
    gene_name = first(gene_name[gene_name != "NA"]),
    cmp_ref = first(cmp_ref[cmp_ref != "NA"]),
    cmp_ref_gene = first(cmp_ref_gene[cmp_ref_gene != "NA"]),
    oId = first(oId[oId != "NA"])
  ) %>%
  ungroup()

## Merge gtf_novel with combined tracking files
gtf_novel <- gtf_novel %>%
  dplyr::left_join(combined_tracking_annotated, suffix = c(".gtf", ".annot"), by = c("transcript_id" = "TCONS"))

## Remove unwanted transcript classes
transcripts_discard <- gtf_novel[which(gtf_novel$type == "transcript" &
                                                 gtf_novel$class_code_orig %in% c("=", "c", "j", "m", "n", "e", "r", "s")
  ), ]
gtf_novel <-
  gtf_novel[which(!(
    gtf_novel$transcript_id %in% transcripts_discard$transcript_id
  )), ]

# Set unique transcript IDs for each TCONS ID (gffcompare finds different TCONSs within same annotated parent transcript isoforms, resulting in multiple different isoforms with same MSTRG.X.X ID if not fixed)
setDT(gtf_novel)
setorder(gtf_novel, combined_gene_id, transcript_id)
gtf_novel[, counter := cumsum(!duplicated(transcript_id)), by = combined_gene_id]
gtf_novel[, transcript_id_new := paste0(combined_gene_id, ".", counter)]
gtf_novel[, counter := NULL]

## Annotate gene IDs and transcript IDs
gtf_novel <- gtf_novel %>%
  group_by(transcript_id) %>%
  mutate(
    ref_gene_id = ref_gene_id,
    gene_id = ifelse(!is.na(ref_gene_id), ref_gene_id, combined_gene_id),
    ref_ranscript_id = ref_transcript_id,
    transcript_id = combined_transcript_id,
    gene_name = ifelse(!is.na(gene_name.annot), gene_name.annot, gene_name.gtf),
    class_code = first(class_code[class_code != "NA"]),
    class_code_orig = first(class_code_orig[class_code_orig != "NA"]),
    # gene_name = ifelse(is.na(gene_name), gene_id, gene_name)
    ) %>%
  ungroup()
gtf_novel <- subset(gtf_novel, class_code != "u")  # Need to find out what to do with these guys, as some of them occur in quite a lot of samples

## Select relevant columns
gtf_novel <- gtf_novel %>%
  mutate(class_code = class_code_orig) %>%
  dplyr::select(c(seqnames, start, end, width, strand, source, type, score, phase, gene_id, gene_name, transcript_id, ref_gene_id, class_code, tss_id, num_samples, exon_number, cmp_ref.annot, cmp_ref_gene.annot))
colnames(gtf_novel) <- c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "phase", "gene_id", "gene_name", "transcript_id", "ref_gene_id", "class_code", "tss_id", "num_samples", "exon_number", "cmp_ref", "cmp_ref_gene")

## Flagging RefSeq transcripts
message(paste(Sys.time(), "Flagging XR transcript overlap ..."), sep = "\t")
gtf_novel <- annotate_overlap(gtf = gtf_novel,
                              gtf_refseq_basename = gtf_refseq_basename,
                              x_name = "xr")

message(paste(Sys.time(), "Flagging NR transcript overlap ..."), sep = "\t")
gtf_novel <- annotate_overlap(gtf = gtf_novel,
                              gtf_refseq_basename = gtf_refseq_basename,
                              x_name = "nr")

gtf_novel <- data.table::as.data.table(gtf_novel)

## Add biotype to custom annotation based on reference ID 
message(paste(Sys.time(), "Adding biotype to StringTie transcripts ..."),
        sep = "\t")
gtf_novel$gene_biotype <-
  gtf_reference$gene_biotype[match(gtf_novel$ref_gene_id,
                                   gtf_reference$gene_id)]
gtf_novel[which(is.na(gtf_novel$gene_biotype)), "gene_biotype"] <- "stringtie"

## Check ref tx overlap 
message(paste(Sys.time(), "Checking reference transcript overlap ..."),
        sep = "\t")
ref_overlap_txs <- suppressWarnings(check_tx_overlap(gtf = gtf_novel,
                                    gtf_reference = gtf_reference))

## Annotate same-stranded i class
message(paste(Sys.time(), "Remove same sense i class transcripts ..."),
        sep = "\t")
same_strand_i_txs <- suppressWarnings(filter_i_class(gtf_df = gtf_novel,
                                    reference_granges = gtf_reference))

## Filter for final transcript subset
gtf_novel <- subset(
  gtf_novel,!(gtf_novel$transcript_id %in% same_strand_i_txs) &
    !(gtf_novel$transcript_id %in% ref_overlap_txs)
)

## Rename new stringtie genes
message(paste(Sys.time(), "Rename novel genes ..."), sep = "\t")
gtf_novel_stringtie <-
  rename_stringtie_transcripts(gtf_novel_df = gtf_novel)
gtf_novel_GR <-
  GenomicRanges::makeGRangesFromDataFrame(gtf_novel_stringtie, keep.extra.columns = T)

# Merge novel GTF with hg38 GTF
gtf_novel_merged <- c(gtf_reference, gtf_novel_GR)

# Export GTF
message(paste(Sys.time(), "Exporting custom gtf ... "), sep = "\t")
rtracklayer::export.gff(object = gtf_novel_merged, con = outfile)
message(paste(Sys.time(), "Exporting custom gtf ... Done!"), sep = "\t")
