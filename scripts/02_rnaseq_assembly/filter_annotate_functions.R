# GTF filtering functions 

count_mono_exonics <- function(gtf) {
  # Input:
  # gtf_novel = Granges object containing novel transcripts
  #
  # Finds transcripts with a single exon and removes those from the object
  
  # Remove transcripts that are found by StringTie with only single exon
  t <- as.data.frame(table(gtf$transcript_id))
  colnames(t) <- c("transcript_id",
                   "exon_count")
  
  # Each transcript has a transcript row and exon row
  t$exon_count <- t$exon_count - 1
  gtf_count <- left_join(as.data.frame(gtf),
                         t,
                         by = "transcript_id")
  
  mono_transcripts <- subset(gtf_count,
                             exon_count < 2 &
                               is.na(ref_gene_id) &
                               source == "StringTie")
  
  return(mono_transcripts)
  
}

###########################################################################

calculate_tx_occurence <-
  function(transcript_tracking_loc,
           tracking_file,
           min_occurence = 1,
           min_tpm = 1,
           gtf) {
    # Input:
    # tracking_file = table output from GFFcompare that list occurences of each
    #                 transfrag in each sample
    # gtf_novel = Granges output from GFFcompare GTF
    # min_occurence = integer (default = 1) how many samples should share a transfrag
    # min_tpm = integer (default = 1) how much coverage each sample should have to count as covered
    #
    # Calculate the occurence in the pool of samples, subsets transcripts on
    # minimal occurence parameter. Output occurence file in directory of tracking file.
    # for further filtering.
    
    colnames(tracking_file) <-
      c("transfrag_id", "locus_id", "ref_gene_id", "class_code")
    tracking_file <- subset(tracking_file,
                            class_code == "=")
    
    tracking_genes <- data.frame(samples = tracking_file[, 1:4])
    tracking_genes$rowid <- 1:nrow(tracking_genes)
    
    tracking_samples <-
      data.frame(samples = tracking_file[, 5:length(tracking_file)])
    # tracking_samples$rowid <- 1:nrow(tracking_samples)
    
    library(stringr)

    tracking_samples_tpm <- str_extract(tracking_samples, "(?<=\\|)[^|]+$")
    
    tracking_samples_tpm <- data.frame(tracking_samples_tpm)
    tracking_samples_tpm <-
      dplyr::mutate_all(tracking_samples_tpm, function(x)
        as.numeric(as.character(x)))
    
    # Count transfrag occurrence
    occurrence_df <-
      data.frame(ifelse(tracking_samples_tpm > min_tpm, 1, 0))
    occurrence_df[is.na(occurrence_df)] <- 0
    occurrence_df <- data.frame(rowSums(occurrence_df))
    colnames(occurrence_df) <- "occurrence"
    sum(occurrence_df > 2)
    
    # Subset transcripts on their occurence
    tracking_genes <- cbind(tracking_genes, occurrence_df)
    tracking_genes <-
      tracking_genes[tracking_genes$occurrence >= min_occurence, ]
    tracking_genes$transcript_id <-
      gsub(".*\\|", "", tracking_genes$samples.ref_gene_id)
    
    # Grab colnames for sample matrix
    sample_ids <- readLines(paste(wd, "documentation/sample_ids.txt", sep = "/"))
    
    # Subset sample matrix and write to tracking file directory
    tracking_file <-
      tracking_file[tracking_file$ref_gene_id %in% tracking_genes$samples.ref_gene_id, ]
    
    colnames(tracking_file) <-
      c("transfrag_id",
        "locus_id",
        "ref_gene_id",
        "class_code",
        sample_ids)
    write.table(
      tracking_file,
      file = gsub(
        ".tracking",
        "sample_matrix.txt",
        transcript_tracking_loc
      ),
      quote = F,
      row.names = F,
      sep = "\t"
    )
    
    gtf_keep <- subset(gtf, gtf$transcript_id %in%
                         tracking_genes$transcript_id)
    
    
    return(gtf_keep)
  }

filter_tx_occurrence <- function(combined_gtf_loc, min_occurence, gtf){
  
  combined_gtf <- rtracklayer::import(combined_gtf_loc)
  combined_gtf_filtered <- subset(combined_gtf, num_samples > 2)  

  # for each gtf$transcript_id transcript structure that has the same id as a combined_gtf_filtered$cmp_ref transcript structure, set the gtf structure to be the same as the combined_gtf_filtered$cmp_ref structure
  for (i in 1:nrow(gtf)){
    if (gtf$transcript_id[i] %in% combined_gtf_filtered$cmp_ref){
      gtf[i, "transcript_id"] <- combined_gtf_filtered$cmp_ref[which(combined_gtf_filtered$cmp_ref == gtf$transcript_id[i])]
    }
  }


  gtf_keep <- subset(gtf, gtf$transcript_id %in% combined_gtf_filtered$cmp_ref)


  return(gtf_keep)
}

###########################################################################

check_tx_overlap <- function(gtf, gtf_reference) {
  # gtf_novel_processed = input from previous function
  # gtf_reference = ensembl reference GTF Grange
  #
  # This function looks for exonic overlap between query (novel tx exons)
  # and subject (ref tx exons). Novel transcripts that overlap with multiple
  # reference genes are flagged for removal.
  
  # These genes should have no reference overlap
  gtf_no_overlap <- gtf[gtf$class_code %in% c("u", "x", "i", "y"), ]
  no_overlap_tx <- gtf_no_overlap$transcript_id
  
  # Create Granges with reference transcript exons using Ensembl canonical transcript
  gene_tx_refs <-
    gtf_reference[which(gsub(".*-", "", gtf_reference$transcript_name) == "201"), ]$transcript_id
  gene_tx_refs <- unique(gene_tx_refs)
  
  gene_ref_gtf <-
    gtf_reference[which(gtf_reference$type == "exon" &
                          gtf_reference$transcript_id %in% gene_tx_refs), ]
  
  # Grab novel exons
  gtf_novel_exons <-
    GenomicRanges::makeGRangesFromDataFrame(gtf[which(gtf$type == "exon"), ],
                                            keep.extra.columns = T)
  
  gtf_no_overlap_exons <-
    GenomicRanges::makeGRangesFromDataFrame(gtf[which(gtf$type == "exon" &
                                                        gtf$transcript_id %in% no_overlap_tx), ],
                                            keep.extra.columns = T)
  
  # Calculate overlapping exons
  no_overlap <-
    GenomicRanges::findOverlaps(query = gtf_no_overlap_exons,
                                subject = gene_ref_gtf,
                                type = "any")
  if (!(length(no_overlap) == 0)) {
    print("WARNING: Overlap in unannotated genes detected")
    print(length(no_overlap))
  }
  print(length(no_overlap))
  
  overlap <- GenomicRanges::findOverlaps(query = gtf_novel_exons,
                                         subject = gene_ref_gtf,
                                         type = "any")
  overlap_df <- data.frame(gtf_novel_exons[from(overlap),
                                           c("transcript_id",
                                             "type",
                                             "class_code",
                                             "cmp_ref",
                                             "ref_gene_id")],
                           gene_ref_gtf[to(overlap), c("gene_id", "gene_name")])
  
  # Calculate overlapping ref genes
  split_by_tx <- split(overlap_df, overlap_df$transcript_id)
  
  check <- bind_rows(lapply(split_by_tx, function(x) {
    # Count per transcript the number of unique gene IDs with overlapping
    # exons. Label transcripts based on this.
    
    single_check = ifelse(length(unique(x$gene_id)) > 1, "multi_gene", "single_gene")
    df = data.frame(transcript_id = unique(x$transcript_id), single_check)
    
    return(df)
    
  }))
  
  check_no_overlap <- unique(subset(check,
                                    check$single_check == "multi_gene")$transcript_id)
  
  return(check_no_overlap)
}

###########################################################################

filter_i_class <- function(gtf_df, reference_granges) {
  # Input:
  # reference_granges = GRanges object of the reference GTF
  # gtf_novel_df = data frame or data table containing new transcripts
  #
  # Searches for same-strand overlap with known genes for i class
  # transcripts.
  tx_i <- subset(gtf_df, type == "transcript" & class_code == "i")$transcript_id
    # gtf_df[which(gtf_df$type == "transcript" &
    #                gtf_df$class_code == "i")]$transcript_id
  gtf_i <- subset(gtf_df, transcript_id %in% tx_i)
  # Convert the entire gtf_i data frame to a GRanges object once
  gtf_i_granges <- makeGRangesFromDataFrame(gtf_i, keep.extra.columns = TRUE)

  # Function to check overlaps
  check_overlap <- function(tx_id) {
    # Subset GRanges object for the current transcript_id
    x <- gtf_i_granges[gtf_i_granges$transcript_id == tx_id]
  
    # Find overlaps
    overlap <- findOverlaps(query = x, subject = reference_granges, type = "any")
  
    # Check if there is any overlap
    if (length(overlap) > 0) {
      return(tx_id)
    }
  }

# Get unique transcript IDs
unique_tx_ids <- unique(gtf_i$transcript_id)

# Apply the function to each transcript_id and unlist results
i_no_pass <- unlist(lapply(unique_tx_ids, check_overlap))

# Return the results
return(i_no_pass)
  
}

###########################################################################

annotate_overlap <- function(gtf, gtf_refseq_basename, x_name) {
  # Input:
  # gtf_novel = Previous Granges output
  # gtf_refseq_basename = name of custom refseq gtf that holds XR and NR transcripts respectively
  # x_name = character vector of the name of the new column and the name of
  #          the RefSeq GFF that contains the transcripts
  #
  # Annotates transcripts found in gtf_novel with any overlap in the X
  # granges object. Currently used to annotate XR_### and NR_### transcripts
  # as these might not be applicable for neo-antigen detection.
  
  x <-
    rtracklayer::import(paste(gtf_refseq_basename, x_name, "gff", sep = "."))
  
  # Convert refseq seqnames to ensembl seqnames
  x <-
    x[seqnames(x) %in% c(
      "NC_000001.11",
      "NC_000002.12",
      "NC_000003.12",
      "NC_000004.12",
      "NC_000005.10",
      "NC_000006.12",
      "NC_000007.14",
      "NC_000008.11",
      "NC_000009.12",
      "NC_000010.11",
      "NC_000011.10",
      "NC_000012.12",
      "NC_000013.11",
      "NC_000014.9",
      "NC_000015.10",
      "NC_000016.10",
      "NC_000017.11",
      "NC_000018.10",
      "NC_000019.10",
      "NC_000020.11",
      "NC_000021.9",
      "NC_000022.11",
      "NC_000023.11",
      "NC_000024.10"
    )]
  seqlevels(x) <- as.character(unique(seqnames(x)))
  x <- GenomeInfoDb::renameSeqlevels(x, c(1:22, "X", "Y"))
  x <- subset(x, x$type == "exon")
  
  gtf_novel_gr <-
    GenomicRanges::makeGRangesFromDataFrame(gtf, keep.extra.columns = T)
  
  gtf_novel_exons <-
    subset(gtf_novel_gr, gtf_novel_gr$type == "exon")
  
  x_overlap <-
    GenomicRanges::findOverlaps(query = x,
                                subject = gtf_novel_exons,
                                type = "any")
  
  novel_x_hits <-
    as.data.frame(unique(gtf_novel_exons[subjectHits(x_overlap)]))
  novel_x_hits <- unique(novel_x_hits[, "transcript_id"])
  elementMetadata(gtf_novel_gr)[[paste0(x_name, "_overlap")]] <-
    ifelse(gtf_novel_gr$transcript_id %in% novel_x_hits,
           paste0(x_name, "_hit"),
           "none")
  return(gtf_novel_gr)
}

###########################################################################

rename_stringtie_transcripts <- function(gtf_novel_df) {
  # Input:
  # gtf_novel_df = data frame or data table containing new transcripts
  #
  # Renames genes without a reference gene in a similar fashion as
  # we name our ORFs using chromosome, start and stop position. Linked
  # transcripts will have the same gene name
  
  gtf_by_gene <- split(gtf_novel_df, gtf_novel_df$gene_id)
  
  gtf_renamed <-
    suppressWarnings(bind_rows(lapply(gtf_by_gene, function(x) {
      if (is.na(x$ref_gene_id)[1]) {
        chr <- unique(x$seqnames)
        start <- min(x$start)
        end <- max(x$end)
        x$gene_name <- paste0(chr, ":", start, "-", end)
      }
      
      return(x)
      
    })))
}