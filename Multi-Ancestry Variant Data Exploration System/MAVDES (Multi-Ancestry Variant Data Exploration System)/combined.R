### 0. Set up


suppressPackageStartupMessages(library(optparse))
#Argument parser with help
option_list <- list(
  make_option(c("-g","--gene"), type="character", default="GENE",
              help = "Gene basename (default: %default)"),
  
  make_option(c("-i","--input_dir"), type="character", default=".",
              help = paste(
                "Input directory containing basic required files.",
                "Expected files:",
                "\n  ├── <GENE>.xlsx                  # gene locus coordinates",
                "\n  SampleArea.xlsx             -- mapping of samples to regions",
                "├── SampleArea.xlsx              # sample-to-region mapping",
                "├── snv/                          # SNV subdirectory",
                "│     ├── <GENE>_snv_vep_output.vcf",
                "├── indel/                          # INDEL subdirectory",
                "│     ├── <GENE>_indel_vep_output.vcf",
                "├── sv/                          # SV subdirectory",
                "│     ├── <GENE>_sv_vep_output.vcf",
                "",
                "If the structure differs, adjust input paths accordingly.",
                sep = "\n"
              )),
  
  make_option(c("-o","--output_dir"), type="character", default="combined_results",
              help="Output folder path"),
  
  make_option(c("-s","--step"), type="character", default="1,2,3",
              help = paste(
                "Pipeline steps to run (comma-separated).",
                "Available steps:",
                "1 = Variability score", 
                "2 = Comparison with evolution metrics", 
                "3 = Haplotypes", 
                "*Step 2 must be run together step 1",
                "Example: --step 1,2,3",
                "(default: %default)",
                sep = "\n"
              )),
  
  make_option(c("--msa_dir"), type="character", default="msa",
              help=paste(
                "Directory containing MSA and evolutionary matrices for Step 2.",
                "Expected files:",
                "\n  ├── *.fas or *.fasta  -- one or multiple MSA files starting with <GENE> (e.g., GENE_A.fas, GENE_B.fas)",
                "\n  ├── distance/      -- subfolder with *.txt distance matrices (matching file names with MSA files)",
                "\n  ├── dn-ds/         -- subfolder with *.txt dn-ds matrices (matching file names with MSA files)",
                "Files will be detected dynamically based on the gene basename provided.",
                "",
                "FASTA naming rules:",
                "  >GENE_Hsa        -- human sequence (must end with '_Hsa')",
                "  >GENE_species    -- non-human sequences (any species tag not ending with '_Hsa')",
                "  The human gene prefix (e.g. “GENE1” from “GENE1-a_Hsa”) is used to detect related sequences.",
                "  Any sequence whose name begins with this prefix followed by “_”, “-”, “−” or nothing (e.g. GENE1_Mmu, GENE1_2) is accepted.",
                "  Names starting with a different prefix (e.g. GENE10_Mmu) are excluded.",
                "",
                "Matrices must contain matching sequence names and be square NxN tables.",
                "Pairwise distance matrix can be fully filled or triangular (upper or lower)",
                "In dN-dS matrix upper triangular part contains dN–dS estimates and lower triangular part contains the corresponding p-values",
                sep="\n")),
  
  #OPTIONAL ARGUMENTS
  make_option("--variant_missing_threshold", type = "double", default = 0.5,
              help = "Threshold for removing variants with missing data proportion above this value (default: %default)"),
  
  make_option("--sample_missing_threshold", type = "double", default = 0.7,
              help = "Threshold for removing samples with missing data proportion above this value (default: %default)"),
  
  make_option("--use_length_factor", type = "logical", default = TRUE,
              help = "Whether variability scores should be boosted by variant length (default: %default)"),
  
  make_option("--normalize_by_samples", type = "logical", default = TRUE,
              help = "Whether contributions should be normalized by available sample count (default: %default)"),
  
  make_option("--detailed", type = "logical", default = FALSE,
              help = "If TRUE, produce detailed regression summaries for gene conservation analysis (default: %default)")
  
)

opt_parser <- OptionParser(option_list = option_list)

#Parse arguments
args <- if (interactive()) character(0) else commandArgs(trailingOnly = TRUE)
opt <- parse_args(opt_parser, args = args)

#Extract variables
gene_name <- opt$gene
input_dir <- opt$input_dir
output_dir <- opt$output_dir
steps <- as.integer(unlist(strsplit(opt$step, ",")))
msa_dir <- opt$msa_dir
variant_missing_threshold <- opt$variant_missing_threshold
sample_missing_threshold  <- opt$sample_missing_threshold
use_length_factor      <- opt$use_length_factor
normalize_by_samples   <- opt$normalize_by_samples
detailed <- opt$detailed


#Print summary
cat("Gene:", gene_name, "\n")
cat("Input dir:", input_dir, "\n")
cat("Output dir:", output_dir, "\n")
cat("Steps to run:", paste(steps, collapse=", "), "\n")
cat("MSA dir:", msa_dir, "\n")
cat("Variant missing threshold:", variant_missing_threshold, "\n")
cat("Sample missing threshold:", sample_missing_threshold, "\n")
cat("Use length factor:", use_length_factor, "\n")
cat("Normalize by samples:", normalize_by_samples, "\n")
cat("Detailed conservation analysis:", detailed, "\n")


### Libraries


suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(tibble)
  library(GenomicRanges)
  library(stringr)
  library(tidyr)
  library(forcats)
  library(ggrepel)
  library(reshape2)
  library(pheatmap)
  library(readxl)
  library(tidyverse)
  library(ggplot2)
  library(Biostrings)
  library(reshape2)
  library(ggrepel)
  library(MASS)
  library(purrr)
  library(haplo.stats)
  library(readr)
  library(data.table)
  library(RColorBrewer)
})

#dplyr preferences
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("desc", "dplyr")

#base preferences
conflicted::conflict_prefer("intersect", "base")
conflicted::conflict_prefer("setdiff", "base")

#data.table preferences
conflicted::conflict_prefer("melt", "data.table")
conflicted::conflict_prefer("dcast", "data.table")

#stringr preference
conflicted::conflict_prefer("fixed", "stringr")


### Helping functions


#Safely extract ALT alleles
#Returns comma-separated string, or "." if missing
extract_alt_safely <- function(alt_sequences) {
  sapply(alt_sequences, function(x) {
    if (length(x) > 0) paste(as.character(x), collapse = ",") else "."
  })
}

#Safely extract REF alleles as character
extract_ref_safely <- function(ref_sequences) {
  as.character(ref_sequences)
}

#Safely extract CSQ (VEP consequence) field
#Handles missing or empty entries
csq_safe <- function(csq_entry) {
  if (is.null(csq_entry) || length(csq_entry) == 0) return(NA)
  sapply(csq_entry, function(x) {
    if (length(x) == 0) return(NA)
    paste(as.character(x), collapse = ",")
  }, USE.NAMES = FALSE)
}

#Categorize VEP consequences into human-readable categories
#Uses severity hierarchy from HIGH -> MODERATE -> LOW -> MODIFIER
categorize_vep_consequences <- function(consequence_vector) {
  if (is.null(consequence_vector) || length(consequence_vector) == 0) {
    return("Unknown")
  }
  
  #Split multi-consequence strings like "a&b&c" into individual terms
  consequences <- unique(unlist(strsplit(consequence_vector, "&")))
  
  #HIGH impact
  if ("transcript_ablation" %in% consequences) return("Transcript ablation")
  if ("splice_acceptor_variant" %in% consequences) return("Splice acceptor")
  if ("splice_donor_variant" %in% consequences) return("Splice donor")
  if ("stop_gained" %in% consequences) return("Stop gained")
  if ("frameshift_variant" %in% consequences) return("Frameshift")
  if ("stop_lost" %in% consequences) return("Stop lost")
  if ("start_lost" %in% consequences) return("Start lost")
  if ("transcript_amplification" %in% consequences) return("Transcript amplification")
  if ("feature_elongation" %in% consequences) return("Feature elongation")
  if ("feature_truncation" %in% consequences) return("Feature truncation")
  
  #MODERATE impact
  if ("inframe_insertion" %in% consequences) return("Inframe insertion")
  if ("inframe_deletion" %in% consequences) return("Inframe deletion")
  if ("missense_variant" %in% consequences) return("Missense")
  if ("protein_altering_variant" %in% consequences) return("Protein altering")
  
  #LOW impact
  if ("splice_polypyrimidine_tract_variant" %in% consequences ||
      "splice_region_variant" %in% consequences ||
      "splice_donor_5th_base_variant" %in% consequences ||
      "splice_donor_region_variant" %in% consequences) return("Splice region")
  if ("incomplete_terminal_codon_variant" %in% consequences) return("Incomplete terminal codon")
  if ("start_retained_variant" %in% consequences) return("Start retained")
  if ("stop_retained_variant" %in% consequences) return("Stop retained")
  if ("synonymous_variant" %in% consequences) return("Synonymous")
  
  #MODIFIER impact
  if ("coding_sequence_variant" %in% consequences) return("Coding sequence")
  if ("mature_miRNA_variant" %in% consequences) return("Mature miRNA variant")
  if ("3_prime_UTR_variant" %in% consequences ||
      "5_prime_UTR_variant" %in% consequences) return("UTR")
  if ("non_coding_transcript_exon_variant" %in% consequences) return("Non coding transcript exon")
  if ("intron_variant" %in% consequences) return("Intron")
  if ("NMD_transcript_variant" %in% consequences) return("NMD transcript")
  if ("non_coding_transcript_variant" %in% consequences) return("Non coding transcript")
  if ("coding_transcript_variant" %in% consequences) return("Coding transcript")
  if ("upstream_gene_variant" %in% consequences) return("Upstream gene")
  if ("downstream_gene_variant" %in% consequences) return("Downstream gene")
  if ("TFBS_ablation" %in% consequences ||
      "TFBS_amplification" %in% consequences ||
      "TF_binding_site_variant" %in% consequences) return("TF binding site")
  if ("regulatory_region_ablation" %in% consequences || 
      "regulatory_region_variant" %in% consequences ||
      "regulatory_region_amplification" %in% consequences) return("Regulatory region")
  if ("intergenic_variant" %in% consequences) return("Intergenic")
  if ("sequence_variant" %in% consequences) return("Sequence")
  
  #Fallback category if none match
  return("Other")
}

## Additional information

#Load sample-region information
load_region_information <- function(sample_area_file) {
  #Read sample metadata
  sampleInfo <- read_xlsx(sample_area_file)
  
  #Get unique population/region names
  regions <- unique(sampleInfo$population)
  
  #Assign colors (Set1 palette, repeated if needed)
  n <- length(regions)
  colors <- brewer.pal(min(max(n,3), 9), "Set1")
  if(n > 9){
    colors <- rep(colors, length.out = n)
  }
  palette <- setNames(colors, regions)
  
  #Create a list of samples per region
  sampleArea <- lapply(regions, function(r) {
    sampleInfo %>% filter(population == r) %>% pull(sample)
  })
  names(sampleArea) <- regions
  
  return(list(
    sampleArea = sampleArea,
    region_colors = palette
  ))
}

#Load gene locus information
load_gene_information <- function(gene_list_file) {
  #Read gene table from xlsx
  gene_row_table <- read_xlsx(gene_list_file)
  
  #Select relevant columns: gene name, chromosome, start, end
  gene_locus <- gene_row_table %>% dplyr::select(5,2,3,4)
  
  return(gene_locus)
}

#Variant type hierarchy
variant_hierarchy <- c(
  #HIGH impact
  "Transcript ablation", "Frameshift", "Stop gained", "Stop lost",
  "Splice acceptor", "Splice donor", "Start lost",
  "Transcript amplification", "Feature truncation", "Feature elongation",
  
  #MODERATE impact
  "Missense", "Inframe insertion", "Inframe deletion", "Protein altering",
  
  #LOW impact
  "Splice region", "Incomplete terminal codon", "Start retained",
  "Stop retained", "Synonymous",
  
  #MODIFIER impact
  "UTR", "Non coding transcript exon", "Non coding transcript", 
  "Coding transcript", "Coding sequence", "Intron", "NMD transcript",
  "Upstream gene", "Downstream gene", "TF binding site",
  "Regulatory region", "Intergenic", "Sequence"
)

##snv

#Main function to analyze VEP-annotated SNV variants from a VCF
#Returns a list with:
#- results: data.frame with samples as columns, variants as rows
#- categories: vector of VEP categories for each variant
analyze_vep_snv_variants <- function(vcf_file) {
  
  cat("Reading VCF file...\n")
  vcf <- readVcf(vcf_file)
  
  sample_names <- samples(header(vcf))
  cat("Found", length(sample_names), "samples:", paste(sample_names, collapse = ", "), "\n")
  
  #Filter only SNVs
  snv_mask <- isSNV(vcf)
  vcf_snv <- vcf[snv_mask]
  cat("Found", length(vcf_snv), "SNV variants\n")
  
  if (length(vcf_snv) == 0) stop("No SNV variants found in VCF")
  
  #Create variant identifiers
  variant_positions <- paste0(
    as.character(seqnames(vcf_snv)), ":",
    start(vcf_snv), "_",
    extract_ref_safely(ref(vcf_snv)), ">",
    extract_alt_safely(alt(vcf_snv))
  )
  
  #Read genotype and CSQ fields
  geno_data <- geno(vcf_snv)
  csq_data <- csq_safe(info(vcf_snv)$CSQ)
  
  #Convert genotype to data.table in long format
  gt_df <- as.data.frame(geno_data$GT, stringsAsFactors = FALSE)
  colnames(gt_df) <- sample_names
  gt_df$Variant <- variant_positions
  geno_dt <- as.data.table(gt_df)
  
  geno_long <- data.table::melt(
    geno_dt,
    id.vars = "Variant",
    variable.name = "Sample",
    value.name = "GT"
  )
  
  #Keep only mutated genotypes
  geno_long[, GT := as.character(GT)]
  geno_long <- geno_long[GT %in% c("0|1","1|0","1|1","0/1","1/0","./1","1/.",".|1","1|.")]
  
  #Merge with CSQ annotations
  csq_df <- data.table(Variant = variant_positions, CSQ = csq_data)
  geno_long <- merge(geno_long, csq_df, by = "Variant", all.x = TRUE)
  
  #Map CSQ field to human-readable categories
  extract_category <- function(csq_entry) {
    if (is.na(csq_entry)) return("Variant (no CSQ)")
    cons_split <- unlist(strsplit(csq_entry, ","))
    cons_categories <- sapply(cons_split, function(x) {
      parts <- unlist(strsplit(x, "\\|"))
      if (length(parts) >= 2) parts[2] else NA
    })
    categorize_vep_consequences(cons_categories)
  }
  
  geno_long[, Category := vapply(CSQ, extract_category, character(1))]
  
  #Convert back to wide format: Variants x Samples
  wide_dt <- dcast(geno_long, Variant ~ Sample, value.var = "Category", fill = NA)
  wide_dt[, Variant := as.character(Variant)]
  
  full_dt <- data.table(Variant = as.character(variant_positions))
  result_dt_full <- merge(full_dt, wide_dt, by = "Variant", all.x = TRUE, sort = FALSE)
  
  #Ensure all samples are columns
  missing_samples <- setdiff(sample_names, colnames(result_dt_full))
  if (length(missing_samples) > 0) {
    result_dt_full[, (missing_samples) := NA_character_]
  }
  setcolorder(result_dt_full, c("Variant", sample_names))
  
  result_df <- as.data.frame(result_dt_full, stringsAsFactors = FALSE)
  colnames(result_df)[1] <- "Position"
  
  #Vector of overall categories per variant
  all_categories <- sapply(seq_along(csq_data), function(i) {
    if (!is.na(csq_data[i])) {
      cons_split <- unlist(strsplit(csq_data[i], ","))
      cons_categories <- sapply(cons_split, function(x) {
        parts <- unlist(strsplit(x, "\\|"))
        if (length(parts) >= 2) parts[2] else NA
      })
      categorize_vep_consequences(cons_categories)
    } else {
      "No CSQ"
    }
  })
  
  return(list(
    results = result_df,
    categories = all_categories
  ))
}

process_snv_vcf <- function(
    vcf_file, 
    locus,
    missing_rate_threshold = 0.2,
    sample_missing_threshold = 0.2
) {
  # Read VCF
  vcf <- readVcf(vcf_file)
  sample_names <- samples(header(vcf))
  n_samples_total <- length(sample_names)
  n_variants_total <- length(SummarizedExperiment::rowRanges(vcf))
  message("Initial: ", n_samples_total, " samples | ", n_variants_total, " variants")
  
  # Extract genotypes
  geno_df <- as.data.frame(geno(vcf)$GT, stringsAsFactors = FALSE)
  geno_df$Variant <- paste0(
    seqnames(vcf), ":", start(vcf), "_",
    as.character(ref(vcf)), ">",
    sapply(alt(vcf), function(x) as.character(x[1]))
  )
  
  # Long-format genotype table
  geno_long <- data.table::melt(
    as.data.table(geno_df),
    id.vars = "Variant",
    variable.name = "Sample",
    value.name = "GT"
  )
  data.table::set(geno_long, j = "GT", value = as.character(geno_long$GT))
  setDT(geno_long)
  
  # Missing GT patterns
  missing_gt <- c(".", "./.", ".|.", ".|0", "0|.", "0/.", "./0")
  geno_long[, missing := GT %in% missing_gt]
  
  # Remove samples with high missing rate
  sample_missing <- geno_long[, .(missing_rate = mean(missing)), by = Sample]
  bad_samples <- sample_missing[missing_rate > sample_missing_threshold, Sample]
  
  if (length(bad_samples) > 0) {
    message("Removing ", length(bad_samples),
            " samples with missing_rate > ", sample_missing_threshold)
    geno_long <- geno_long[!Sample %in% bad_samples]
  }
  
  n_samples_after <- length(unique(geno_long$Sample))
  
  # Allele decomposition
  geno_long[, c("allele1", "allele2") := tstrsplit(GT, "[/|]")]
  geno_long[, allele1 := as.numeric(allele1)]
  geno_long[, allele2 := as.numeric(allele2)]
  geno_long[is.na(allele1), allele1 := 0]
  geno_long[is.na(allele2), allele2 := 0]
  geno_long[, allele_count := allele1 + allele2]
  
  # Load VEP results
  result <- analyze_vep_snv_variants(vcf_file)$result
  
  result_long <- result %>%
    pivot_longer(
      cols = -Position,
      names_to = "Sample",
      values_to = "Variant_Type"
    ) %>%
    filter(!is.na(Variant_Type))
  
  result_dt <- as.data.table(result_long)
  setnames(result_dt, "Position", "Variant")
  
  # Merge VCF + VEP
  snv_dt <- merge(geno_long, result_dt, 
                  by = c("Variant", "Sample"), all.x = TRUE)
  
  # Filter variants by missing rate
  variant_missing_rate <- snv_dt[, .(missing_rate = mean(missing)), by = Variant]
  
  variants_keep <- variant_missing_rate[
    missing_rate <= missing_rate_threshold, 
    Variant
  ]
  
  snv_dt <- snv_dt[Variant %in% variants_keep]
  
  # Remove REF homozygotes and partially missing calls
  remove_gt <- c("0|0","0/0",".","./.","0|.","0|0",".|0","./0",".|.")
  snv_dt <- snv_dt[!GT %in% remove_gt]
  
  message("Filtering summary")
  message("Samples: ", n_samples_total, " → ", n_samples_after)
  message("Variants: ", n_variants_total,
          " → ", length(variants_keep),
          " (", nrow(variant_missing_rate) - length(variants_keep),
          " removed after missing filter)")
  
  # Gene overlap using foverlaps
  locus_dt <- as.data.table(locus)
  locus_dt[, chr := paste0("chr", Chr)]
  locus_dt[, start := Start]
  locus_dt[, end := End]
  locus_dt[, gene := gene_name]
  locus_dt <- locus_dt[, .(chr, start, end, gene)]
  
  snv_dt[, chr := paste0("chr", sub(":.*", "", Variant))]
  snv_dt[, pos := as.numeric(str_extract(Variant, "(?<=:)[0-9]+"))]
  snv_dt[, start := pos]
  snv_dt[, end := pos]
  
  setkey(snv_dt, chr, start, end)
  setkey(locus_dt, chr, start, end)
  
  snv_with_gene <- foverlaps(snv_dt, locus_dt, 
                             type = "within", nomatch = 0L)
  
  # Final table
  snv_final <- snv_with_gene[, .(
    gene,
    Sample,
    Variant,
    Variant_Type,
    GT,
    allele1,
    allele2,
    allele_count,
    missing
  )]
  
  return(list(
    snv_final = snv_final,
    n_samples_after = n_samples_after
  ))
}

## indel

#Main function to analyze VEP-annotated INDEL variants from a VCF
#Returns a list with:
#- results: data.frame with samples as columns, variants as rows
#- categories: vector of VEP categories for each variant
analyze_vep_indel_variants <- function(vcf_file) {
  
  cat("Reading VCF file...\n")
  vcf <- readVcf(vcf_file)
  
  sample_names <- samples(header(vcf))
  cat("Found", length(sample_names), "samples:", paste(sample_names, collapse = ", "), "\n")
  
  #Filter to only INDELs (insertions or deletions)
  indel_mask <- isIndel(vcf)
  vcf_indel <- vcf[indel_mask]
  cat("Found", length(vcf_indel), "INDEL variants\n")
  
  if (length(vcf_indel) == 0) stop("No INDEL variants found in VCF")
  
  info_df <- info(vcf_indel)
  
  ##-Extract SVLEN-
  if ("SVLEN" %in% colnames(info_df)) {
    svlen_info <- as.numeric(info_df$SVLEN)
  } else {
    info_strings <- tryCatch({
      as.character(info(vcf_indel)$INFO)
    }, error = function(e) {
      rep(NA_character_, length(vcf_indel))
    })
    if (!is.null(info_strings) && length(info_strings) == length(vcf_indel)) {
      svlen_info <- as.numeric(sub(".*SVLEN=([-0-9]+).*", "\\1", info_strings))
    } else {
      svlen_info <- rep(NA_real_, length(vcf_indel))
    }
  }
  svlen_info[is.na(svlen_info)] <- 0
  
  ##-Extract SVTYPE-
  if ("SVTYPE" %in% colnames(info_df)) {
    svtype_info <- as.character(info_df$SVTYPE)
  } else {
    #try to parse from raw INFO field
    info_strings <- tryCatch({
      as.character(info(vcf_indel)$INFO)
    }, error = function(e) rep(NA_character_, length(vcf_indel)))
    
    if (!is.null(info_strings) && length(info_strings) == length(vcf_indel)) {
      svtype_info <- sub(".*SVTYPE=([^;]+).*", "\\1", info_strings)
      svtype_info[grepl("^\\*|^\\.", svtype_info) | svtype_info == info_strings] <- NA
    } else {
      svtype_info <- rep(NA_character_, length(vcf_indel))
    }
  }
  svtype_info[is.na(svtype_info)] <- "UNK"
  
  #Create variant identifiers
  variant_positions <- paste0(
    as.character(seqnames(vcf_indel)), ":",
    start(vcf_indel), "_",
    extract_ref_safely(ref(vcf_indel)), ">",
    extract_alt_safely(alt(vcf_indel))
  )
  
  geno_data <- geno(vcf_indel)
  csq_data <- csq_safe(info(vcf_indel)$CSQ)
  
  gt_df <- as.data.frame(geno_data$GT, stringsAsFactors = FALSE)
  colnames(gt_df) <- sample_names
  gt_df$Variant <- variant_positions
  
  geno_dt <- as.data.table(gt_df)
  
  geno_long <- data.table::melt(
    geno_dt,
    id.vars = "Variant",
    variable.name = "Sample",
    value.name = "GT"
  )
  
  geno_long[, GT := as.character(GT)]
  geno_long <- geno_long[!GT %in% c("0|0", ".|.", ".|0", "0|.")]
  geno_long <- geno_long[GT %in% c("0|1", "1|0", "1|1", "0/1", "1/0", "1/1", ".|1", "1|.", "./1", "1/.")]
  
  csq_df <- data.table(Variant = variant_positions, CSQ = csq_data)
  geno_long <- merge(geno_long, csq_df, by = "Variant", all.x = TRUE)
  
  extract_category <- function(csq_entry) {
    if (is.na(csq_entry)) return("Variant (no CSQ)")
    cons_split <- unlist(strsplit(csq_entry, ","))
    cons_categories <- sapply(cons_split, function(x) {
      parts <- unlist(strsplit(x, "\\|"))
      if (length(parts) >= 2) parts[2] else NA
    })
    categorize_vep_consequences(cons_categories)
  }
  
  geno_long[, Category := vapply(CSQ, extract_category, character(1))]
  
  wide_dt <- dcast(geno_long, Variant ~ Sample, value.var = "Category", fill = NA)
  wide_dt[, Variant := as.character(Variant)]
  
  full_dt <- data.table(
    Variant = as.character(variant_positions),
    SVLEN = svlen_info,
    SVTYPE = svtype_info
  )
  result_dt_full <- merge(full_dt, wide_dt, by = "Variant", all.x = TRUE, sort = FALSE)
  
  missing_samples <- setdiff(sample_names, colnames(result_dt_full))
  if (length(missing_samples) > 0) {
    result_dt_full[, (missing_samples) := NA_character_]
  }
  
  setcolorder(result_dt_full, c("Variant", "SVLEN", "SVTYPE", sample_names))
  
  result_df <- as.data.frame(result_dt_full, stringsAsFactors = FALSE)
  colnames(result_df)[1] <- "Position"
  
  all_categories <- sapply(seq_along(csq_data), function(i) {
    if (!is.na(csq_data[i])) {
      cons_split <- unlist(strsplit(csq_data[i], ","))
      cons_categories <- sapply(cons_split, function(x) {
        parts <- unlist(strsplit(x, "\\|"))
        if (length(parts) >= 2) parts[2] else NA
      })
      categorize_vep_consequences(cons_categories)
    } else {
      "No CSQ"
    }
  })
  
  cat("Analysis completed.\n")
  
  return(list(
    results = result_df,
    categories = all_categories
  ))
}

process_indel_vcf <- function(
    vcf_file,
    locus,
    missing_rate_threshold = 0.2,
    sample_missing_threshold = 0.2
) {
  # Read VCF
  vcf <- readVcf(vcf_file)
  sample_names <- samples(header(vcf))
  n_samples_total <- length(sample_names)
  n_variants_total <- length(SummarizedExperiment::rowRanges(vcf))
  message("Initial: ", n_samples_total, " samples | ", n_variants_total, " variants")
  
  # Filter to INDEL variants only
  indel_mask <- isIndel(vcf)
  vcf_indel <- vcf[indel_mask]
  if (length(vcf_indel) == 0) stop("No INDEL variants found")
  
  # Extract SVLEN and SVTYPE from INFO
  info_df <- info(vcf_indel)
  svlen <- if ("SVLEN" %in% colnames(info_df)) as.numeric(info_df$SVLEN) else rep(0, length(vcf_indel))
  svtype <- if ("SVTYPE" %in% colnames(info_df)) as.character(info_df$SVTYPE) else rep("UNK", length(vcf_indel))
  svlen[is.na(svlen)] <- 0
  svtype[is.na(svtype)] <- "UNK"
  
  # Variant IDs
  variant_ids <- paste0(
    as.character(seqnames(vcf_indel)), ":",
    start(vcf_indel), "_",
    extract_ref_safely(ref(vcf_indel)), ">",
    extract_alt_safely(alt(vcf_indel))
  )
  
  # Genotypes → long format
  gt_df <- as.data.frame(geno(vcf_indel)$GT, stringsAsFactors = FALSE)
  colnames(gt_df) <- sample_names
  gt_df$Variant <- variant_ids
  
  geno_long <- melt(as.data.table(gt_df),
                    id.vars = "Variant",
                    variable.name = "Sample",
                    value.name = "GT")
  setDT(geno_long)
  geno_long[, GT := as.character(GT)]
  
  # Missing genotypes
  missing_gt <- c(".", "./.", ".|.", ".|0", "0|.", "0/.", "./0")
  geno_long[, missing := GT %in% missing_gt]
  
  # Remove samples with high missing rate
  sample_missing <- geno_long[, .(missing_rate = mean(missing)), by = Sample]
  bad_samples <- sample_missing[missing_rate > sample_missing_threshold, Sample]
  
  if (length(bad_samples) > 0) {
    message("Removing ", length(bad_samples),
            " samples with missing_rate > ", sample_missing_threshold)
    geno_long <- geno_long[!Sample %in% bad_samples]
  }
  
  n_samples_after <- length(unique(geno_long$Sample))
  
  # Allele decomposition
  geno_long[, c("allele1", "allele2") := tstrsplit(GT, "[/|]")]
  geno_long[, allele1 := as.numeric(allele1)]
  geno_long[, allele2 := as.numeric(allele2)]
  geno_long[is.na(allele1), allele1 := 0]
  geno_long[is.na(allele2), allele2 := 0]
  geno_long[, allele_count := allele1 + allele2]
  
  # Recompute missing flag cleanly
  geno_long[, missing := GT %in% missing_gt]
  
  # Variant missing filter
  variant_missing_rate <- geno_long[, .(missing_rate = mean(missing)), by = Variant]
  variants_keep <- variant_missing_rate[missing_rate <= missing_rate_threshold, Variant]
  
  n_variants_before_filter <- nrow(variant_missing_rate)
  n_variants_after_filter <- length(variants_keep)
  
  geno_long <- geno_long[Variant %in% variants_keep]
  
  # Remove REF-homozygous / partially missing calls
  no_variant_gt <- c("0|0", "0/0", ".", "./.", "0|.", ".|0", "0/.", "./0", ".|.")
  geno_long <- geno_long[!GT %in% no_variant_gt]
  
  message("Filtering summary")
  message("Samples: ", n_samples_total, " → ", n_samples_after)
  message("Variants: ", n_variants_total, " → ", n_variants_after_filter,
          " (", n_variants_before_filter - n_variants_after_filter, " removed after missing filter)")
  
  # Load VEP annotation
  result <- analyze_vep_indel_variants(vcf_file)$result
  
  result_long <- result %>%
    pivot_longer(
      cols = -c(Position, SVLEN, SVTYPE),
      names_to = "Sample",
      values_to = "Variant_Type"
    ) %>%
    filter(!is.na(Variant_Type))
  
  result_dt <- as.data.table(result_long)
  setnames(result_dt, "Position", "Variant")
  
  # Merge VCF genotypes + VEP annotation
  geno_long <- merge(geno_long, result_dt,
                     by = c("Variant", "Sample"), all.x = TRUE)
  
  # Add SVLEN / SVTYPE
  geno_long[, SVLEN := svlen[match(Variant, variant_ids)]]
  geno_long[, SVTYPE := svtype[match(Variant, variant_ids)]]
  
  geno_long[is.na(SVLEN), SVLEN := 0]
  
  # Gene locus table
  locus_dt <- as.data.table(locus)
  locus_dt[, chr := paste0("chr", Chr)]
  locus_dt[, start := Start]
  locus_dt[, end := End]
  locus_dt[, gene := gene_name]
  locus_dt <- locus_dt[, .(chr, start, end, gene)]
  
  # Variant positions
  geno_long[, chr := paste0("chr", sub(":.*", "", Variant))]
  geno_long[, pos := as.numeric(str_extract(Variant, "(?<=:)[0-9]+"))]
  
  # Define event ranges
  geno_del <- geno_long[SVTYPE == "DEL"]
  geno_del[, start := pos]
  geno_del[, end := pos + abs(SVLEN) - 1]
  
  geno_ins <- geno_long[SVTYPE == "INS"]
  geno_ins[, start := pos]
  geno_ins[, end := pos]
  
  # Overlap with locus
  setkey(locus_dt, chr, start, end)
  setkey(geno_del, chr, start, end)
  setkey(geno_ins, chr, start, end)
  
  del_with_gene <- foverlaps(geno_del, locus_dt, type = "any", nomatch = 0L)
  ins_with_gene <- foverlaps(geno_ins, locus_dt, type = "within", nomatch = 0L)
  
  indel_with_gene <- rbindlist(list(del_with_gene, ins_with_gene))
  
  # Final output table
  indel_final <- indel_with_gene[, .(
    gene,
    Sample,
    Variant,
    Variant_Type,
    SVTYPE,
    SVLEN,
    GT,
    allele1,
    allele2,
    allele_count,
    missing
  )]
  
  return(list(
    indel_final = indel_final,
    n_samples_after = n_samples_after
  ))
}

## sv

#Main function to analyze VEP-annotated SV variants from a VCF
#Returns a list with:
#- results: data.frame with samples as columns, variants as rows
#- categories: vector of VEP categories for each variant
analyze_vep_sv_variants <- function(vcf_file) {
  
  cat("Reading VCF file...\n")
  vcf_sv <- readVcf(vcf_file)
  
  sample_names <- samples(header(vcf_sv))
  cat("Found", length(sample_names), "samples:", paste(sample_names, collapse = ", "), "\n")
  
  cat("Found", length(vcf_sv), "SV variants\n")
  if (length(vcf_sv) == 0) stop("No SV variants found in VCF")
  
  info_df <- info(vcf_sv)
  
  ##-Extract SVLEN-
  if ("SVLEN" %in% colnames(info_df)) {
    svlen_info <- as.numeric(info_df$SVLEN)
  } else {
    info_strings <- tryCatch({
      as.character(info(vcf_sv)$INFO)
    }, error = function(e) {
      rep(NA_character_, length(vcf_sv))
    })
    if (!is.null(info_strings) && length(info_strings) == length(vcf_sv)) {
      svlen_info <- as.numeric(sub(".*SVLEN=([-0-9]+).*", "\\1", info_strings))
    } else {
      svlen_info <- rep(NA_real_, length(vcf_sv))
    }
  }
  svlen_info[is.na(svlen_info)] <- 0
  
  ##-Extract SVTYPE-
  if ("SVTYPE" %in% colnames(info_df)) {
    svtype_info <- as.character(info_df$SVTYPE)
  } else {
    #try to parse from raw INFO field
    info_strings <- tryCatch({
      as.character(info(vcf_sv)$INFO)
    }, error = function(e) rep(NA_character_, length(vcf_sv)))
    
    if (!is.null(info_strings) && length(info_strings) == length(vcf_sv)) {
      svtype_info <- sub(".*SVTYPE=([^;]+).*", "\\1", info_strings)
      svtype_info[grepl("^\\*|^\\.", svtype_info) | svtype_info == info_strings] <- NA
    } else {
      svtype_info <- rep(NA_character_, length(vcf_sv))
    }
  }
  svtype_info[is.na(svtype_info)] <- "UNK"
  
  #Create variant identifiers
  variant_positions <- paste0(
    as.character(seqnames(vcf_sv)), ":",
    start(vcf_sv), "_",
    extract_ref_safely(ref(vcf_sv)), ">",
    extract_alt_safely(alt(vcf_sv))
  )
  
  geno_data <- geno(vcf_sv)
  csq_data <- csq_safe(info(vcf_sv)$CSQ)
  
  gt_df <- as.data.frame(geno_data$GT, stringsAsFactors = FALSE)
  colnames(gt_df) <- sample_names
  gt_df$Variant <- variant_positions
  
  geno_dt <- as.data.table(gt_df)
  
  geno_long <- data.table::melt(
    geno_dt,
    id.vars = "Variant",
    variable.name = "Sample",
    value.name = "GT"
  )
  
  geno_long[, GT := as.character(GT)]
  geno_long <- geno_long[!GT %in% c("0|0", ".|.", ".|0", "0|.")]
  geno_long <- geno_long[GT %in% c("0|1", "1|0", "1|1", "0/1", "1/0", "1/1", ".|1", "1|.", "./1", "1/.")]
  
  csq_df <- data.table(Variant = variant_positions, CSQ = csq_data)
  geno_long <- merge(geno_long, csq_df, by = "Variant", all.x = TRUE)
  
  extract_category <- function(csq_entry) {
    if (is.na(csq_entry)) return("Variant (no CSQ)")
    cons_split <- unlist(strsplit(csq_entry, ","))
    cons_categories <- sapply(cons_split, function(x) {
      parts <- unlist(strsplit(x, "\\|"))
      if (length(parts) >= 2) parts[2] else NA
    })
    categorize_vep_consequences(cons_categories)
  }
  
  geno_long[, Category := vapply(CSQ, extract_category, character(1))]
  
  wide_dt <- dcast(geno_long, Variant ~ Sample, value.var = "Category", fill = NA)
  wide_dt[, Variant := as.character(Variant)]
  
  full_dt <- data.table(
    Variant = as.character(variant_positions),
    SVLEN = svlen_info,
    SVTYPE = svtype_info
  )
  result_dt_full <- merge(full_dt, wide_dt, by = "Variant", all.x = TRUE, sort = FALSE)
  
  missing_samples <- setdiff(sample_names, colnames(result_dt_full))
  if (length(missing_samples) > 0) {
    result_dt_full[, (missing_samples) := NA_character_]
  }
  
  setcolorder(result_dt_full, c("Variant", "SVLEN", "SVTYPE", sample_names))
  
  result_df <- as.data.frame(result_dt_full, stringsAsFactors = FALSE)
  colnames(result_df)[1] <- "Position"
  
  all_categories <- sapply(seq_along(csq_data), function(i) {
    if (!is.na(csq_data[i])) {
      cons_split <- unlist(strsplit(csq_data[i], ","))
      cons_categories <- sapply(cons_split, function(x) {
        parts <- unlist(strsplit(x, "\\|"))
        if (length(parts) >= 2) parts[2] else NA
      })
      categorize_vep_consequences(cons_categories)
    } else {
      "No CSQ"
    }
  })
  
  cat("Analysis completed.\n")
  
  return(list(
    results = result_df,
    categories = all_categories
  ))
}

process_sv_vcf <- function(
    vcf_file,
    locus,
    missing_rate_threshold = 0.2,
    sample_missing_threshold = 0.2
) {
  # Read VCF
  vcf <- readVcf(vcf_file)
  sample_names <- samples(header(vcf))
  n_samples_total <- length(sample_names)
  n_variants_total <- length(SummarizedExperiment::rowRanges(vcf))
  message("Initial: ", n_samples_total, " samples | ", n_variants_total, " variants")
  
  # Select structural variants (Indel-based filter — consistent with your code)
  sv_mask <- isIndel(vcf)
  vcf_sv <- vcf[sv_mask]
  if (length(vcf_sv) == 0) stop("No SV variants found")
  
  # Extract SVLEN / SVTYPE
  info_df <- info(vcf_sv)
  svlen <- if ("SVLEN" %in% colnames(info_df)) as.numeric(info_df$SVLEN) else rep(0, length(vcf_sv))
  svtype <- if ("SVTYPE" %in% colnames(info_df)) as.character(info_df$SVTYPE) else rep("UNK", length(vcf_sv))
  svlen[is.na(svlen)] <- 0
  svtype[is.na(svtype)] <- "UNK"
  
  # Variant ID construction
  variant_ids <- paste0(
    as.character(seqnames(vcf_sv)), ":",
    start(vcf_sv), "_",
    extract_ref_safely(ref(vcf_sv)), ">",
    extract_alt_safely(alt(vcf_sv))
  )
  
  # Genotypes → long format
  gt_df <- as.data.frame(geno(vcf_sv)$GT, stringsAsFactors = FALSE)
  colnames(gt_df) <- sample_names
  gt_df$Variant <- variant_ids
  
  geno_long <- melt(as.data.table(gt_df),
                    id.vars = "Variant",
                    variable.name = "Sample",
                    value.name = "GT")
  setDT(geno_long)
  geno_long[, GT := as.character(GT)]
  
  # Missing genotype definitions
  missing_gt <- c(".", "./.", ".|.", ".|0", "0|.", "0/.", "./0")
  geno_long[, missing := GT %in% missing_gt]
  
  # Remove samples with high missing rate
  sample_missing <- geno_long[, .(missing_rate = mean(missing)), by = Sample]
  bad_samples <- sample_missing[missing_rate > sample_missing_threshold, Sample]
  
  if (length(bad_samples) > 0) {
    message("Removing ", length(bad_samples),
            " samples with missing_rate > ", sample_missing_threshold)
    geno_long <- geno_long[!Sample %in% bad_samples]
  }
  
  n_samples_after <- length(unique(geno_long$Sample))
  
  # Allele decomposition
  geno_long[, c("allele1", "allele2") := tstrsplit(GT, "[/|]")]
  geno_long[, allele1 := suppressWarnings(as.numeric(allele1))]
  geno_long[, allele2 := suppressWarnings(as.numeric(allele2))]
  
  geno_long[is.na(allele1), allele1 := 0]
  geno_long[is.na(allele2), allele2 := 0]
  
  geno_long[, allele_count := allele1 + allele2]
  
  # Recompute missing
  geno_long[, missing := GT %in% missing_gt]
  
  # Variant-level missing filter
  variant_missing_rate <- geno_long[, .(missing_rate = mean(missing)), by = Variant]
  variants_keep <- variant_missing_rate[missing_rate <= missing_rate_threshold, Variant]
  
  n_variants_before_filter <- nrow(variant_missing_rate)
  n_variants_after_filter <- length(variants_keep)
  
  geno_long <- geno_long[Variant %in% variants_keep]
  
  # Remove REF-homozygous + partially missing calls
  no_variant_gt <- c("0|0", "0/0", ".", "./.", "0|.", ".|0", "0/.", "./0", ".|.")
  geno_long <- geno_long[!GT %in% no_variant_gt]
  
  message("Filtering summary")
  message("Samples: ", n_samples_total, " → ", n_samples_after)
  message("Variants: ", n_variants_total, " → ", n_variants_after_filter,
          " (", n_variants_before_filter - n_variants_after_filter, " removed after missing filter)")
  
  # Add VEP annotation
  result <- analyze_vep_sv_variants(vcf_file)$result
  
  result_long <- result %>%
    pivot_longer(cols = -c(Position, SVLEN, SVTYPE),
                 names_to = "Sample",
                 values_to = "Variant_Type") %>%
    filter(!is.na(Variant_Type))
  
  result_dt <- as.data.table(result_long)
  setnames(result_dt, "Position", "Variant")
  
  # Merge genotype + VEP
  geno_long <- merge(geno_long, result_dt,
                     by = c("Variant", "Sample"), all.x = TRUE)
  
  # Add SVLEN/SVTYPE explicitly
  geno_long[, SVLEN := svlen[match(Variant, variant_ids)]]
  geno_long[, SVTYPE := svtype[match(Variant, variant_ids)]]
  
  geno_long[is.na(SVLEN), SVLEN := 0]
  
  # Locus table
  locus_dt <- as.data.table(locus)
  locus_dt[, chr := paste0("chr", Chr)]
  locus_dt[, start := Start]
  locus_dt[, end := End]
  locus_dt[, gene := gene_name]
  locus_dt <- locus_dt[, .(chr, start, end, gene)]
  
  # Variant positions
  geno_long[, chr := paste0("chr", sub(":.*", "", Variant))]
  geno_long[, pos := as.numeric(str_extract(Variant, "(?<=:)[0-9]+"))]
  
  # Event ranges
  geno_del_inv <- geno_long[SVTYPE %in% c("DEL", "INV")]
  geno_del_inv[, start := pos]
  geno_del_inv[, end := pos + abs(SVLEN) - 1]
  
  geno_ins <- geno_long[SVTYPE == "INS"]
  geno_ins[, start := pos]
  geno_ins[, end := pos]
  
  # Overlap with locus
  setkey(locus_dt, chr, start, end)
  setkey(geno_del_inv, chr, start, end)
  setkey(geno_ins, chr, start, end)
  
  del_inv_with_gene <- foverlaps(geno_del_inv, locus_dt,
                                 type = "any", nomatch = 0L)
  ins_with_gene <- foverlaps(geno_ins, locus_dt,
                             type = "within", nomatch = 0L)
  
  sv_with_gene <- rbindlist(list(del_inv_with_gene, ins_with_gene))
  
  # Final output
  sv_final <- sv_with_gene[, .(
    gene,
    Sample,
    Variant,
    Variant_Type,
    SVTYPE,
    SVLEN,
    GT,
    allele1,
    allele2,
    allele_count,
    missing
  )]
  
  return(list(
    sv_final = sv_final,
    n_samples_after = n_samples_after
  ))
}

## 1. variability score

# Calculate gene-level conservation / mutational burden scores
calculate_gene_conservation_score <- function(
    snv_final,          # SNV table: gene, Sample, Variant, Variant_Type, GT, allele_count, missing
    indel_final,        # INDEL table with same structure as SNV table
    sv_final,           # SV table with additional SVTYPE, SVLEN columns
    gene_locus,         # Gene coordinates: gene_name (or gene), Start, End
    total_samples_snv,  # Number of samples contributing SNV data
    total_samples_indel,# Number of samples contributing INDEL data
    total_samples_sv,   # Number of samples contributing SV data
    use_length_factor = TRUE,      # Whether SV scores should be boosted by SV length
    normalize_by_samples = TRUE   # Whether contributions should be normalized by available sample count
) {
  
  # Convert to data.table
  snv_final   <- as.data.table(snv_final)
  indel_final <- as.data.table(indel_final)
  sv_final    <- as.data.table(sv_final)
  gene_locus  <- as.data.table(gene_locus)
  
  # Variant type → weight mapping
  variant_weights <- data.table(
    Variant_Type = c(
      "Transcript ablation", "Splice acceptor", "Splice donor",
      "Stop gained", "Frameshift", "Stop lost", "Start lost",
      "Transcript amplification", "Feature elongation", "Feature truncation",
      "Missense", "Inframe insertion", "Inframe deletion", "Protein altering",
      "Splice region", "Incomplete terminal codon",
      "Start retained", "Stop retained", "Synonymous",
      "Coding sequence", "UTR", "Mature miRNA variant",
      "Non coding transcript exon", "Intron", "NMD transcript",
      "Non coding transcript", "Coding transcript",
      "Upstream gene", "Downstream gene",
      "TF binding site", "Regulatory region",
      "Intergenic", "Sequence", "Other"
    ),
    weight = c(
      3.0, 2.6, 2.6,
      2.5, 2.5, 2.2, 2.2,
      2.0, 2.0, 2.0,
      1.8, 1.6, 1.6, 1.6,
      1.5, 1.0,
      0.8, 0.8, 0.5,
      0.4, 0.3, 0.3,
      0.25, 0.2, 0.25,
      0.25, 0.25,
      0.15, 0.15,
      0.3, 0.25,
      0.1, 0.1, 0.1
    )
  )
  
  message("Detected variant types:")
  message("  SNV:   ", paste(unique(snv_final$Variant_Type), collapse = ", "))
  message("  INDEL: ", paste(unique(indel_final$Variant_Type), collapse = ", "))
  message("  SV:    ", paste(unique(sv_final$Variant_Type), collapse = ", "))
  
  # SNVs do not have SVTYPE / SVLEN → harmonize columns
  snv_final[, SVTYPE := NA_character_]
  snv_final[, SVLEN  := NA_character_]
  
  common_cols <- c(
    "gene", "Sample", "Variant", "Variant_Type",
    "SVTYPE", "SVLEN", "GT", "allele_count", "missing"
  )
  
  # Combine all variant types into a single table
  combined <- rbindlist(list(
    snv_final[, ..common_cols][, variant_group := "SNV"],
    indel_final[, ..common_cols][, variant_group := "INDEL"],
    sv_final[, ..common_cols][, variant_group := "SV"]
  ), fill = TRUE)
  
  # Add variant weights
  combined <- merge(combined, variant_weights, by = "Variant_Type", all.x = TRUE)
  
  # SV length scoring factor
  if (use_length_factor) {
    combined[, SVLEN_num := as.numeric(SVLEN)]
    combined[, length_factor :=
               ifelse(!is.na(SVLEN_num) & SVLEN_num != 0,
                      1 + pmin(log10(abs(SVLEN_num) + 1), 2),
                      1)]
  } else {
    combined[, length_factor := 1]
  }
  
  # Compute allele_count contribution
  combined[, allele_count_num := as.numeric(allele_count)]
  combined[is.na(allele_count_num), allele_count_num := 0]
  
  combined[, contribution := weight * allele_count_num * length_factor]
  
  # Normalize by sample count
  if (normalize_by_samples) {
    combined[variant_group == "SNV",
             contribution := contribution / total_samples_snv * 100]
    combined[variant_group == "INDEL",
             contribution := contribution / total_samples_indel * 100]
    combined[variant_group == "SV",
             contribution := contribution / total_samples_sv * 100]
  }
  
  # Aggregate at gene level
  gene_scores <- combined[, .(
    total_contribution = sum(contribution, na.rm = TRUE),
    n_variants = uniqueN(Variant),
    n_samples  = uniqueN(Sample)
  ), by = gene]
  
  # Gene locus normalisation
  setnames(gene_locus, "gene_name", "gene", skip_absent = TRUE)
  
  gene_scores_full <- merge(
    gene_locus[, .(gene, Start, End)],
    gene_scores,
    by = "gene",
    all.x = TRUE
  )
  
  gene_scores_full[is.na(total_contribution), total_contribution := 0]
  gene_scores_full[is.na(n_variants),         n_variants := 0]
  gene_scores_full[is.na(n_samples),          n_samples := 0]
  
  # Normalize by gene length
  gene_scores_full[, total_contribution_len :=
                     total_contribution / (End - Start)]
  
  gene_scores_full <- gene_scores_full[order(total_contribution_len)]
  
  message("Processed genes: ", nrow(gene_scores_full))
  rng <- range(gene_scores_full$total_contribution_len, na.rm = TRUE)
  message(sprintf(
    "Normalized contribution range: [%.4f, %.4f]",
    rng[1], rng[2]
  ))
  
  return(gene_scores_full)
}

## visualization

plot_gene_variation_scores <- function(gene_scores_full, gene_name, output_dir = "part2/graphs/score/") {
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 1. Histogram — distribution of variation scores per base (log1p scale)
  p1 <- ggplot(gene_scores_full, aes(x = log1p(total_contribution_len))) +
    geom_histogram(bins = 30, fill = "#5a78c2", color = "white", alpha = 0.9) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(color = "gray20")
    ) +
    labs(
      title = paste("Distribution of Variation Score per bp in", gene_name),
      x = "log(1 + Variation score per bp)",
      y = "Number of genes"
    )
  
  ggsave(file.path(output_dir, paste0(gene_name, "_distribution_log.png")),
         p1, width = 8, height = 5, dpi = 300)
  message("→ Histogram of variation scores saved and plotted.")
  
  # 2. Barplot — top 20 most variable genes (log1p scale)
  top_genes <- gene_scores_full[order(-total_contribution_len)]
  top_genes <- head(top_genes, 20)
  
  p2 <- ggplot(top_genes, aes(x = reorder(gene, total_contribution_len), y = log1p(total_contribution_len))) +
    geom_col(fill = "#e34a33", alpha = 0.9, width = 0.7) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(size = 9, color = "gray20"),
      axis.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("Top 20 Most Variable Genes in", gene_name),
      x = "Gene",
      y = "log(1 + Variation score per bp)"
    )
  
  ggsave(file.path(output_dir, paste0(gene_name, "_top20_log.png")),
         p2, width = 8, height = 6, dpi = 300)
  message("→ Barplot of top 20 genes saved and plotted.")
  
  # 3. Scatter — total variation vs. gene length (log1p scale)
  gene_scores_full[, gene_length := End - Start]
  
  p3 <- ggplot(gene_scores_full, aes(x = log1p(gene_length), y = log1p(total_contribution))) +
    geom_point(alpha = 0.7, size = 2.5, color = "#41ab5d") +
    geom_smooth(method = "lm", se = FALSE, color = "gray40", linetype = "dashed") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray95"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(color = "gray20")
    ) +
    labs(
      title = paste("Total Variation vs Gene Length in", gene_name),
      x = "log(1 + Gene length [bp])",
      y = "log(1 + Total variation score)"
    )
  
  ggsave(file.path(output_dir, paste0(gene_name, "_length_correlation_log.png")),
         p3, width = 8, height = 5.5, dpi = 300)
  message("→ Scatter plot of variation vs gene length saved and plotted.")
  
  # 4. Line plot — ranked gene variation score (log1p scale)
  gene_scores_full <- gene_scores_full[order(total_contribution_len)]
  gene_scores_full[, gene := factor(gene, levels = gene)]
  
  n_genes <- nrow(gene_scores_full)
  base_width <- 8
  width_per_gene <- 0.15
  plot_width <- max(base_width, n_genes * width_per_gene)
  
  p4 <- ggplot(gene_scores_full, aes(
    x = reorder(gene, total_contribution_len),
    y = log1p(total_contribution_len)
  )) +
    geom_segment(
      aes(xend = gene, y = 0, yend = log1p(total_contribution_len)),
      color = "#5a78c2", linewidth = 0.8, alpha = 0.7
    ) +
    geom_point(color = "#e34a33", size = 2.5, alpha = 0.9) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 8, color = "gray20"),
      axis.text.y = element_text(size = 9, color = "gray20"),
      axis.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
      plot.margin = margin(20, 20, 40, 20)
    ) +
    labs(
      title = paste("Gene Variation Score per Base Pair (", gene_name, " locus)", sep = ""),
      x = "Gene (ranked by normalized score)",
      y = "log(1 + Variation score per bp)"
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(gene_name, "_log.png")),
    plot = p4,
    width = plot_width,
    height = 6,
    dpi = 300,
    limitsize = FALSE
  )
  message("→ Ranked gene variation score plot saved and plotted.")
}


### 2. add evolutionary results


## map (distance & dn-ds)

# Load distance and dN-dS matrices for all segments
read_matrices <- function(path_dist_dir, path_dnds_dir, segments) {
  dist_list <- list()
  dnds_list <- list()
  
  for (seg in segments) {
    # Distance matrix
    if (has_distance) {
      dist_mat <- as.matrix(read.csv(
        file.path(path_dist_dir, paste0(seg, ".txt")),
        row.names = 1, na.strings = "NA", check.names = FALSE
      ))
      dist_long <- reshape2::melt(dist_mat, na.rm = TRUE)
      colnames(dist_long) <- c("Pair1", "Pair2", "Distance")
      dist_long$Distance <- suppressWarnings(as.numeric(dist_long$Distance))
      dist_long <- dist_long[!is.na(dist_long$Distance), , drop = FALSE]
      dist_list[[seg]] <- dist_long
    }
    
    # dN-dS and p-value matrix
    if (has_dnds) {
      full_matrix <- as.matrix(read.csv(
        file.path(path_dnds_dir, paste0(seg, ".txt")),
        row.names = 1, na.strings = "NA", check.names = FALSE
      ))
      
      dnds_matrix <- full_matrix
      dnds_matrix[lower.tri(dnds_matrix, diag = TRUE)] <- NA
      dnds_matrix[lower.tri(dnds_matrix)] <- t(dnds_matrix)[lower.tri(dnds_matrix)]
      dnds_matrix[is.na(dnds_matrix)] <- 0
      
      pvalue_matrix <- full_matrix
      pvalue_matrix[upper.tri(pvalue_matrix, diag = TRUE)] <- NA
      pvalue_matrix[upper.tri(pvalue_matrix)] <- t(pvalue_matrix)[upper.tri(pvalue_matrix)]
      pvalue_matrix[is.na(pvalue_matrix)] <- 1
      
      dnds_long <- reshape2::melt(dnds_matrix, na.rm = TRUE)
      colnames(dnds_long) <- c("Pair1", "Pair2", "dN_dS")
      dnds_long$dN_dS <- suppressWarnings(as.numeric(dnds_long$dN_dS))
      
      pvalue_long <- reshape2::melt(pvalue_matrix, na.rm = TRUE)
      colnames(pvalue_long) <- c("Pair1", "Pair2", "p")
      pvalue_long$p <- suppressWarnings(as.numeric(pvalue_long$p))
      
      dnds_data <- merge(dnds_long, pvalue_long, by = c("Pair1", "Pair2"))
      dnds_data$p_adj <- p.adjust(dnds_data$p, method = "BH")
      
      dnds_list[[seg]] <- dnds_data
    }
  }
  
  dist_all <- if (has_distance) do.call(rbind, dist_list) else NULL
  dnds_all <- if (has_dnds) do.call(rbind, dnds_list) else NULL
  
  return(list(dist = dist_all, dnds = dnds_all))
}

# Extract all human gene names and their core identifiers
get_human_genes_and_cores <- function(msa_files) {
  human_full <- character(0)
  
  for (msa in msa_files) {
    aln <- readAAStringSet(msa)
    idx <- grep("_Hsa$", names(aln))
    if (length(idx) > 0) human_full <- c(human_full, names(aln)[idx])
  }
  
  human_full <- unique(human_full)
  cores <- str_extract(sub("_Hsa$", "", human_full), "^[A-Za-z]+[0-9]+")
  data.frame(human_full = human_full, core = cores, stringsAsFactors = FALSE)
}

# Find non-human genes matching a given core
match_nonhuman_by_core <- function(core, all_genes) {
  pattern <- paste0("\\b", core, "\\b")
  matches <- all_genes[str_detect(all_genes, regex(pattern))]
  matches[!str_detect(matches, "_Hsa")]
}

# Compare each human gene to related genes based on core naming
compare_genes_core <- function(msa_files, dnds_data = NULL, dist_data = NULL, locus = NULL) {
  
  has_dnds <- !is.null(dnds_data)
  has_dist <- !is.null(dist_data)
  
  # Extract all sequence names from MSA files
  all_genes <- unlist(lapply(msa_files, function(f) {
    lines <- readLines(f)
    gsub("^>", "", lines[grepl("^>", lines)])
  }))
  
  human_genes <- grep("_Hsa$", all_genes, value = TRUE)
  
  #restrict to genes present in locus
  if (!is.null(locus)) {
    base_names <- gsub("_Hsa$", "", human_genes)
    invalid <- base_names[!base_names %in% locus$gene]
    if (length(invalid) > 0) {
      message("Removing genes outside locus: ", paste(invalid, collapse = ", "))
      human_genes <- human_genes[!gsub("_Hsa$", "", human_genes) %in% invalid]
    }
  }
  
  results <- list()
  
  for (gene_name in human_genes) {
    core_name <- sub("([-−_].*)?$", "", gene_name)
    pattern <- paste0("^", core_name, "([_−-]|$)")
    
    related <- all_genes[
      grepl(pattern, all_genes) &
        !grepl("_Hsa$", all_genes)
    ]
    
    # --- safe initialization ---
    mean_dnds <- NA; median_dnds <- NA
    mean_distance <- NA; median_distance <- NA
    
    if (length(related) > 0) {
      
      if (has_dnds) {
        dnds_pairs <- dnds_data %>%
          filter(
            p_adj < 0.05 &
              ((Pair1 == gene_name & Pair2 %in% related) |
                 (Pair2 == gene_name & Pair1 %in% related))
          )
        if (nrow(dnds_pairs) > 0) {
          mean_dnds <- mean(dnds_pairs$dN_dS, na.rm = TRUE)
          median_dnds <- median(dnds_pairs$dN_dS, na.rm = TRUE)
        }
      }
      
      if (has_dist) {
        dist_pairs <- dist_data %>%
          filter(
            (Pair1 == gene_name & Pair2 %in% related) |
              (Pair2 == gene_name & Pair1 %in% related)
          )
        if (nrow(dist_pairs) > 0) {
          mean_distance <- mean(dist_pairs$Distance, na.rm = TRUE)
          median_distance <- median(dist_pairs$Distance, na.rm = TRUE)
        }
      }
    }
    
    results[[gene_name]] <- data.frame(
      human_gene = gsub("_Hsa$", "", gene_name),
      mean_distance = mean_distance,
      median_distance = median_distance,
      mean_dnds = mean_dnds,
      median_dnds = median_dnds,
      n_compared = length(related),
      matched_genes = if(length(related)>0) paste(related, collapse = ";") else NA
    )
  }
  
  do.call(rbind, results)
}

## distribution and correlation

analyze_gene_conservation <- function(
    complete_data, 
    gene_name = "GENE",
    output_dir = "part2/graphs/compare",
    detailed = FALSE
) {
  # Create output directory if needed
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  plots <- list()
  labels <- list()
  
  # Histograms
  plots[[length(plots)+1]] <- ggplot(complete_data, aes(x = total_contribution_len)) + 
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(title = paste(gene_name, "gene conservation score distribution"),
         x = "Score", y = "Gene count")
  labels[[length(labels)+1]] <- "score"
  
  if (any(!is.na(complete_data$mean_dnds))) {
    plots[[length(plots)+1]] <- ggplot(complete_data, aes(x = mean_dnds)) +
      geom_histogram(bins = 30, fill = "tomato", color = "black") +
      theme_minimal() +
      labs(title = paste(gene_name, "mean dN-dS distribution"),
           x = "Mean dN-dS", y = "Gene count")
    labels[[length(labels)+1]] <- "dnds"
  }
  
  if (any(!is.na(complete_data$mean_distance))) {
    plots[[length(plots)+1]] <- ggplot(complete_data, aes(x = mean_distance)) +
      geom_histogram(bins = 30, fill = "khaki", color = "black") +
      theme_minimal() +
      labs(title = paste(gene_name, "mean distance distribution"),
           x = "Mean distance", y = "Gene count")
    labels[[length(labels)+1]] <- "dist"
  }
  
  # Save histograms
  for (i in seq_along(plots)) {
    ggsave(file.path(output_dir, paste0(gene_name, "_hist_", labels[[i]], ".png")),
           plots[[i]], width=6, height=4, dpi=300)
  }
  
  # Spearman correlations
  cor_score_dnds <- if (all(is.na(complete_data$mean_dnds))) NA else
    cor(complete_data$total_contribution_len, complete_data$mean_dnds, use = "complete.obs", method = "spearman")
  
  cor_score_dist <- if (all(is.na(complete_data$mean_distance))) NA else
    cor(complete_data$total_contribution_len, complete_data$mean_distance, use = "complete.obs", method = "spearman")
  
  corr_summary <- data.frame(cor_score_dnds = cor_score_dnds,
                             cor_score_dist = cor_score_dist)
  
  
  message("Correlation summary:")
  print(corr_summary)
  
  # Optional regression models
  if (detailed) {
    message("Detailed regression model summaries:")
    
    model_interaction <- lm(total_contribution_len ~ mean_dnds * mean_distance, data = complete_data)
    model_dnds <- lm(mean_dnds ~ total_contribution_len, data = complete_data)
    model_dist <- lm(mean_distance ~ total_contribution_len, data = complete_data)
    
    print(summary(model_interaction))
    print(summary(model_dnds))
    print(summary(model_dist))
  }
  
  return(list(
    data = complete_data,
    correlations = corr_summary
  ))
}

##gene comparison

# Prepare data for conservation analysis
prepare_conservation_data <- function(complete_data) {
  complete_data <- complete_data %>%
    mutate(
      z_score = scale(total_contribution_len),
      z_dnds = scale(mean_dnds),
      z_dist = scale(mean_distance)
    )
  
  complete_data <- complete_data %>%
    mutate(
      # Categorize genes based on quantiles
      conservation_level = case_when(
        total_contribution_len >= quantile(total_contribution_len, 0.66, na.rm = TRUE) ~ "var",
        total_contribution_len <= quantile(total_contribution_len, 0.33, na.rm = TRUE) ~ "cons",
        TRUE ~ "mid"
      ),
      dnds_level = case_when(
        mean_dnds >= quantile(mean_dnds, 0.66, na.rm = TRUE) ~ "var",
        mean_dnds <= quantile(mean_dnds, 0.33, na.rm = TRUE) ~ "cons",
        TRUE ~ "mid"
      ),
      dist_level = case_when(
        mean_distance >= quantile(mean_distance, 0.66, na.rm = TRUE) ~ "var",
        mean_distance <= quantile(mean_distance, 0.33, na.rm = TRUE) ~ "cons",
        TRUE ~ "mid"
      ),
      # Squeezed z-scores for plotting
      z_score_squeezed = atan(z_score) / (pi / 2),
      z_dnds_squeezed = atan(z_dnds) / (pi / 2),
      z_dist_squeezed = atan(z_dist) / (pi / 2)
    )
  
  return(complete_data)
}

# Plot human gene variability vs evolutionary dN-dS
plot_dnds_comparison <- function(complete_data, gene_name, output_dir = "part2/graphs/compare/") {
  if (all(is.na(complete_data$mean_dnds))) {
    message("No dN-dS data available — skipping dN-dS comparison plot.")
    return(invisible(NULL))
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  complete_data <- prepare_conservation_data(complete_data)
  complete_data <- complete_data %>%
    mutate(pattern = paste0("human_", conservation_level, "_evo_", dnds_level))
  
  p <- ggplot(complete_data, aes(x = z_dnds_squeezed, y = z_score_squeezed, color = pattern)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_text_repel(aes(label = gene), size = 2.5, max.overlaps = 50) +
    theme_minimal() +
    labs(
      x = "Evolutionary variability (squeezed z_dN-dS)",
      y = "Human population variability (squeezed z-score)",
      title = paste("Comparison of gene conservation patterns in", gene_name, "(squeezed scale)")
    )
  
  ggsave(file.path(output_dir, paste0(gene_name, "_dnds.png")),
         p, width = 10, height = 8, dpi = 300, limitsize = FALSE)
}

# Plot human gene variability vs evolutionary distance
plot_distance_comparison <- function(complete_data, gene_name, output_dir = "part2/graphs/compare/") {
  if (all(is.na(complete_data$mean_distance))) {
    message("No distance data available — skipping distance comparison plot.")
    return(invisible(NULL))
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  complete_data <- prepare_conservation_data(complete_data)
  complete_data <- complete_data %>%
    mutate(pattern = paste0("human_", conservation_level, "_evo_", dist_level))
  
  p <- ggplot(complete_data, aes(x = z_dist_squeezed, y = z_score_squeezed, color = pattern)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_text_repel(aes(label = gene), size = 2.5, max.overlaps = 50) +
    theme_minimal() +
    labs(
      x = "Evolutionary variability (squeezed z-distance)",
      y = "Human population variability (squeezed z-score)",
      title = paste("Comparison of gene conservation patterns in", gene_name, "(distance-based, squeezed scale)")
    )
  
  ggsave(file.path(output_dir, paste0(gene_name, "_dist.png")),
         p, width = 10, height = 8, dpi = 300, limitsize = FALSE)
}

# Combined comparison: human variability vs both dN-dS and distance
plot_combined_conservation_comparison <- function(complete_data, gene_name, output_dir = "part2/graphs/compare/") {
  if (all(is.na(complete_data$mean_dnds)) || 
      all(is.na(complete_data$mean_distance))) {
    message("Both metrics required for combined comparison — skipping combined plot.")
    return(invisible(NULL))
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  complete_data <- prepare_conservation_data(complete_data)
  
  plot_data <- complete_data %>%
    mutate(
      metric1_x = z_dnds_squeezed,
      metric2_x = z_dist_squeezed
    ) %>%
    select(gene, z_score_squeezed, metric1_x, metric2_x, conservation_level) %>%
    tidyr::pivot_longer(
      cols = c(metric1_x, metric2_x),
      names_to = "metric",
      values_to = "x_value"
    ) %>%
    mutate(
      metric = recode(metric,
                      "metric1_x" = "Evolutionary variability (dN-dS)",
                      "metric2_x" = "Evolutionary variability (mean distance)")
    )
  
  p <- ggplot(plot_data, aes(x = x_value, y = z_score_squeezed, color = conservation_level)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_text_repel(aes(label = gene), size = 2.2, max.overlaps = 50) +
    facet_wrap(~metric, scales = "fixed") +
    scale_color_manual(values = c("cons" = "steelblue", "mid" = "grey60", "var" = "tomato")) +
    theme_minimal(base_size = 11) +
    labs(
      x = "Evolutionary variability (squeezed z-values)",
      y = "Human population variability (squeezed z-score)",
      color = "Human-level\nvariability",
      title = paste("Comparison of gene conservation patterns across metrics in", gene_name, "(squeezed scale)")
    ) +
    theme(strip.text = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, paste0(gene_name, "_both.png")),
         p, width = 10, height = 8, dpi = 300, limitsize = FALSE)
}


###3. haplotype


run_haplotype_analysis <- function(snv_final, indel_final, sv_final, sample_region,
                                   output_file = "haplo_results_partial.csv") {
  # Combine all variant types into a single table
  all_variants <- bind_rows(
    snv_final %>% mutate(Source = "SNV"),
    indel_final %>% mutate(Source = "INDEL"),
    sv_final %>% mutate(Source = "SV")
  ) %>%
    group_by(gene, Source) %>%
    mutate(Variant_idx = dense_rank(Variant),
           Variant = paste0(gene, "_", Source, "_V", Variant_idx)) %>%
    ungroup() %>%
    select(-any_of(c("Variant_Type", "SVTYPE", "SVLEN", "Variant_idx")))
  
  # Map each sample to its region
  sample_region_map <- data.frame(
    Sample = unlist(sample_region),
    Region = rep(names(sample_region), times = lengths(sample_region)),
    stringsAsFactors = FALSE
  )
  
  all_samples <- unique(all_variants$Sample)
  
  # Create complete Sample × Variant grid
  complete_variants <- all_variants %>%
    distinct(gene, Variant) %>%
    tidyr::expand(gene, Variant, Sample = all_samples)
  
  haplo_input <- complete_variants %>%
    left_join(
      all_variants %>% select(gene, Variant, Sample, allele1, allele2, allele_count),
      by = c("gene", "Variant", "Sample")
    ) %>%
    mutate(
      allele1 = ifelse(is.na(allele1), 0, allele1),
      allele2 = ifelse(is.na(allele2), 0, allele2),
      allele_count = ifelse(is.na(allele_count), 0, allele_count)
    ) %>%
    left_join(sample_region_map, by = "Sample")
  
  message("Data preparation complete. Starting per-gene haplotype analysis...")
  gc()
  
  # Analyze haplotypes for a single gene
  analyze_gene_haplotypes <- function(gene_name, data) {
    message("\nAnalyzing gene: ", gene_name)
    
    geno_gene <- data %>%
      filter(gene == gene_name) %>%
      select(Sample, Variant, allele1, allele2)
    
    # Complete grid for alleles
    geno_gene_complete <- tidyr::expand_grid(
      Sample = unique(geno_gene$Sample),
      Variant = unique(geno_gene$Variant)
    ) %>%
      left_join(geno_gene, by = c("Sample", "Variant")) %>%
      mutate(
        allele1 = ifelse(is.na(allele1), 0, allele1) + 1,
        allele2 = ifelse(is.na(allele2), 0, allele2) + 1
      )
    
    # Convert to matrix for haplo.em
    geno_mat <- geno_gene_complete %>%
      pivot_wider(
        names_from = Variant,
        values_from = c(allele1, allele2),
        names_glue = "{Variant}_{.value}",
        values_fill = 0
      ) %>%
      select(-Sample) %>%
      as.matrix()
    
    locus_names <- sub("_(allele1|allele2)$", "", colnames(geno_mat))
    locus_table <- table(locus_names)
    complete_loci <- names(locus_table[locus_table == 2])
    
    if (length(complete_loci) < 1) {
      message("→ No complete loci, skipping gene")
      return(NULL)
    }
    
    geno_mat <- geno_mat[, locus_names %in% complete_loci, drop = FALSE]
    
    # Keep only loci with variability
    keep_loci <- sapply(complete_loci, function(loc) {
      cols <- grep(paste0("^", loc, "_allele[12]$"), colnames(geno_mat), value = TRUE)
      any(apply(geno_mat[, cols, drop = FALSE], 2, \(x) length(unique(x)) > 1))
    })
    final_loci <- complete_loci[keep_loci]
    
    message("Variable loci kept: ", paste(final_loci, collapse = ", "))
    
    # Handle case with only one variable locus
    if (length(final_loci) < 2) {
      message("→ Single variable locus – computing haplotypes directly")
      loc <- final_loci[1]
      cols <- grep(paste0("^", loc, "_allele[12]$"), colnames(geno_mat), value = TRUE)
      alleles <- as.vector(geno_mat[, cols, drop = FALSE])
      alleles <- alleles[alleles > 0]  # remove missing
      
      if (length(alleles) == 0) {
        message("→ No alleles left, locus is monomorphic")
        return(data.frame(
          gene = gene_name,
          haplotype_id = "SINGLE_1",
          variant = loc,
          allele = 1,
          freq = 1.0
        ))
      }
      
      unique_alleles <- sort(unique(alleles))
      haplo_freq <- sapply(unique_alleles, function(a) mean(alleles == a))
      
      return(data.frame(
        gene = rep(gene_name, length(unique_alleles)),
        haplotype_id = paste0("SINGLE_", seq_along(unique_alleles)),
        variant = loc,
        allele = unique_alleles,
        freq = haplo_freq
      ))
    }
    
    # Run EM haplotype estimation
    keep_cols <- unlist(lapply(final_loci, function(loc) {
      grep(paste0("^", loc, "_allele[12]$"), colnames(geno_mat), value = TRUE)
    }))
    geno_mat <- geno_mat[, keep_cols, drop = FALSE]
    
    haplo_res <- tryCatch({
      haplo.em(geno = geno_mat, locus.label = final_loci)
    }, error = function(e) {
      message("haplo.em failed for ", gene_name, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(haplo_res)) return(NULL)
    
    hap_df <- as.data.frame(haplo_res$haplotype)
    hap_df$gene <- gene_name
    hap_df$freq <- haplo_res$hap.prob
    hap_df <- hap_df %>%
      mutate(haplotype_id = paste0("HAP_", seq_len(n()))) %>%
      pivot_longer(cols = all_of(final_loci), names_to = "variant", values_to = "allele") %>%
      select(gene, haplotype_id, variant, allele, freq)
    
    return(hap_df)
  }
  
  # Save results incrementally
  if (file.exists(output_file)) file.remove(output_file)
  genes <- unique(haplo_input$gene)
  first_write <- TRUE
  
  for (g in genes) {
    res <- analyze_gene_haplotypes(g, haplo_input)
    if (!is.null(res)) {
      if (first_write) {
        readr::write_csv(res, output_file, append = FALSE)
        first_write <- FALSE
      } else {
        readr::write_csv(res, output_file, append = TRUE, col_names = FALSE)
      }
    }
    rm(res); gc()
  }
  
  haplo_results <- readr::read_csv(output_file, show_col_types = FALSE)
  
  return(list(
    haplo_input = haplo_input,
    haplo_results = haplo_results
  ))
}

##graphs

analyze_haplotype_variation <- function(haplo_input, haplo_results, locus,
                                        gene_name, output_dir = "part2/graphs/haplotypes") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 1. Compute average allele count per gene and region
  hap_freq_region <- haplo_input %>%
    group_by(gene, Region) %>%
    summarise(mean_allele_count = mean(allele_count, na.rm = TRUE), .groups = "drop")
  
  # Identify top 16 genes with highest variance across regions
  gene_var <- hap_freq_region %>%
    group_by(gene) %>%
    summarise(var_across_regions = var(mean_allele_count, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(var_across_regions)) %>%
    slice_head(n = 16)
  
  # Plot: allele count for top variable genes
  p1 <- hap_freq_region %>%
    left_join(gene_var, by = "gene") %>%
    filter(gene %in% gene_var$gene) %>%
    ggplot(aes(x = Region, y = mean_allele_count, color = Region)) +
    geom_point(size = 2) +
    facet_wrap(~fct_reorder(gene, -var_across_regions), scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(title = paste("Average allele count for top variable", gene_name, "genes"),
         y = "Average allele count", x = NULL) +
    theme(legend.position = "top",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(file.path(output_dir, paste0(gene_name, "_top_var.png")), p1, width = 10, height = 6, dpi = 300)
  
  # 2. Plot allele variant levels per gene across regions
  p2 <- hap_freq_region %>%
    ggplot(aes(x = gene, y = mean_allele_count, color = Region)) +
    geom_point(position = position_dodge(width = 0.5), size = 3, alpha = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Allele variant levels between regions in", gene_name),
         y = "Average allele count", x = "Gene") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9),
          legend.position = "top",
          panel.grid.minor = element_blank())
  
  ggsave(file.path(output_dir, paste0(gene_name, "_allele_variants.png")), p2, width = 10, height = 5, dpi = 300)
  
  # 3. Global comparison of allele counts across regions using boxplots
  top_genes_by_region <- hap_freq_region %>%
    group_by(Region) %>%
    slice_max(order_by = mean_allele_count, n = 5, with_ties = FALSE)
  
  p3 <- hap_freq_region %>%
    ggplot(aes(x = Region, y = mean_allele_count, fill = Region)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
    geom_text_repel(data = top_genes_by_region,
                    aes(label = gene),
                    size = 2, color = "black",
                    box.padding = 0.15, max.overlaps = 30,
                    segment.color = "grey50") +
    theme_minimal(base_size = 12) +
    labs(title = paste("Allele count between regions in", gene_name),
         y = "Average allele count", x = NULL) +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, paste0(gene_name, "_global_region_boxplot.png")), p3, width = 7, height = 5, dpi = 300)
  
  # 4. Summarize haplotypes per gene and normalize by gene length
  locus <- locus %>%
    mutate(length_bp = End - Start + 1)
  
  haplo_summary <- haplo_results %>%
    group_by(gene, haplotype_id) %>%
    summarise(freq = mean(freq), .groups = "drop") %>%
    group_by(gene) %>%
    summarise(
      n_haplotypes = n(),
      total_freq = sum(freq),
      top_freq = max(freq),
      diversity = -sum(freq * log(freq)),
      .groups = "drop"
    ) %>%
    left_join(locus %>% select(gene_name, length_bp), by = c("gene" = "gene_name")) %>%
    mutate(n_haplotypes_per_kb = n_haplotypes / (length_bp / 1000)) %>%
    arrange(desc(n_haplotypes_per_kb))
  
  plot_width <- max(6, 3 + nrow(haplo_summary) * 0.3)
  
  p4 <- ggplot(haplo_summary,
               aes(x = reorder(gene, -n_haplotypes_per_kb),
                   y = n_haplotypes_per_kb,
                   fill = diversity)) +
    geom_col() +
    scale_fill_viridis_c(option = "plasma") +
    theme_minimal(base_size = 13) +
    labs(title = paste("Normalized number of haplotypes per gene in", gene_name),
         x = "Gene",
         y = "Number of haplotypes per kb",
         fill = "Diversity\n(Shannon)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, paste0(gene_name, "_haplotype_summary.png")),
         p4, width = plot_width, height = 6, dpi = 300)
  
  # 5. Correlation heatmap of allele counts between regions
  region_matrix <- hap_freq_region %>%
    pivot_wider(names_from = Region, values_from = mean_allele_count) %>%
    filter(if_any(-gene, ~ !is.na(.)))
  
  region_corr <- region_matrix %>%
    select(-gene) %>%
    cor(use = "pairwise.complete.obs", method = "spearman")
  
  heatmap_path <- file.path(output_dir, paste0(gene_name, "_correlation.png"))
  
  # Generate correlation heatmap
  breaks <- unique(seq(min(region_corr, na.rm=TRUE),
                       max(region_corr, na.rm=TRUE), length.out = 100))
  if(length(breaks) < 2) breaks <- c(0,1)
  
  png(heatmap_path, width = 1500, height = 1200, res = 200)
  
  pheat <- pheatmap(
    region_corr,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("white", "orange", "red"))(length(breaks)-1),
    breaks = breaks,
    main = paste("Correlation of allele count between regions in", gene_name)
  )
  
  dev.off()
  
  return(list(
    hap_freq_region = hap_freq_region,
    haplo_summary = haplo_summary,
    region_corr = region_corr
  ))
}

## Per region

# Function to compute haplotypes for a single gene within a dataset
analyze_gene_haplotypes <- function(gene_name, region_data) {
  message("\n------ Analyzing gene: ", gene_name, " ------")
  
  # Extract alleles for the given gene
  geno_gene <- region_data %>%
    filter(gene == gene_name) %>%
    select(Sample, Variant, allele1, allele2)
  
  # Complete grid of Sample × Variant, fill missing alleles with 0
  geno_gene_complete <- tidyr::expand_grid(
    Sample = unique(geno_gene$Sample),
    Variant = unique(geno_gene$Variant)
  ) %>%
    left_join(geno_gene, by = c("Sample", "Variant")) %>%
    mutate(
      allele1 = ifelse(is.na(allele1), 0, allele1) + 1,
      allele2 = ifelse(is.na(allele2), 0, allele2) + 1
    )
  
  # Pivot to wide matrix format for haplotype EM
  geno_mat <- geno_gene_complete %>%
    pivot_wider(
      names_from = Variant,
      values_from = c(allele1, allele2),
      names_glue = "{Variant}_{.value}",
      values_fill = 0
    ) %>%
    select(-Sample) %>%
    as.matrix()
  
  # Identify loci with both alleles present
  locus_names <- sub("_(allele1|allele2)$", "", colnames(geno_mat))
  locus_table <- table(locus_names)
  complete_loci <- names(locus_table[locus_table == 2])
  
  # No complete loci → monomorphic
  if (length(complete_loci) < 1) {
    message("→ No complete loci → monomorphic SINGLE_1 haplotype")
    return(data.frame(
      gene = gene_name,
      haplotype_id = "SINGLE_1",
      variant = NA,
      allele = 1,
      freq = 1
    ))
  }
  
  # Keep only complete loci
  geno_mat <- geno_mat[, locus_names %in% complete_loci, drop = FALSE]
  
  # Filter loci that show variability
  keep_loci <- sapply(complete_loci, function(loc) {
    cols <- grep(paste0("^", loc, "_allele[12]$"), colnames(geno_mat), value = TRUE)
    any(apply(geno_mat[, cols, drop = FALSE], 2, function(x) length(unique(x)) > 1))
  })
  final_loci <- complete_loci[keep_loci]
  
  message("Loci kept after variability filter: ", 
          ifelse(length(final_loci) > 0, paste(final_loci, collapse = ", "), "None"))
  
  # No variable loci → monomorphic
  if (length(final_loci) < 1) {
    message("→ No variable loci → monomorphic SINGLE_1 haplotype")
    return(data.frame(
      gene = gene_name,
      haplotype_id = "SINGLE_1",
      variant = NA,
      allele = 1,
      freq = 1
    ))
  }
  
  # Single variable locus → compute haplotype frequencies directly
  if (length(final_loci) == 1) {
    message("→ Single variable locus – computing haplotypes")
    loc <- final_loci[1]
    cols <- grep(paste0("^", loc, "_allele[12]$"), colnames(geno_mat), value = TRUE)
    alleles <- as.vector(geno_mat[, cols, drop = FALSE])
    alleles <- alleles[alleles > 0]  # remove missing
    
    if (length(alleles) == 0) {
      message("→ No non-zero alleles – monomorphic")
      return(data.frame(
        gene = gene_name,
        haplotype_id = "SINGLE_1",
        variant = loc,
        allele = 1,
        freq = 1
      ))
    }
    
    # Compute allele frequencies
    unique_alleles <- sort(unique(alleles))
    haplo_freq <- sapply(unique_alleles, function(a) mean(alleles == a))
    
    return(data.frame(
      gene = rep(gene_name, length(unique_alleles)),
      haplotype_id = paste0("SINGLE_", seq_along(unique_alleles)),
      variant = loc,
      allele = unique_alleles,
      freq = haplo_freq
    ))
  }
  
  # Multiple variable loci → run EM haplotype estimation
  keep_cols <- unlist(lapply(final_loci, function(loc) {
    grep(paste0("^", loc, "_allele[12]$"), colnames(geno_mat), value = TRUE)
  }))
  geno_mat <- geno_mat[, keep_cols, drop = FALSE]
  
  haplo_res <- tryCatch({
    haplo.em(geno = geno_mat, locus.label = final_loci)
  }, error = function(e) {
    message("haplo.em FAILED for ", gene_name, ": ", e$message)
    return(NULL)
  })
  
  # EM failure → fallback monomorphic
  if (is.null(haplo_res)) {
    message("→ haplo.em failed → fallback SINGLE_1 haplotype")
    return(data.frame(
      gene = gene_name,
      haplotype_id = "SINGLE_1",
      variant = NA,
      allele = 1,
      freq = 1
    ))
  }
  
  # Format EM results
  hap_df <- as.data.frame(haplo_res$haplotype)
  hap_df$gene <- gene_name
  hap_df$freq <- haplo_res$hap.prob
  
  hap_df <- hap_df %>%
    mutate(haplotype_id = paste0("HAP_", seq_len(n()))) %>%
    pivot_longer(cols = all_of(final_loci), names_to = "variant", values_to = "allele") %>%
    select(gene, haplotype_id, variant, allele, freq)
  
  message("→ haplotypes computed: ", nrow(hap_df), " rows")
  return(hap_df)
}


# Function to run haplotype analysis for all genes, separated by population regions
run_haplotype_analysis_by_region <- function(snv_final, indel_final, sv_final, sample_region,
                                             output_file = "haplo_results_by_region.csv") {
  # Combine all variant types into a single table and assign unique Variant IDs
  all_variants <- bind_rows(
    snv_final %>% mutate(Source = "SNV"),
    indel_final %>% mutate(Source = "INDEL"),
    sv_final %>% mutate(Source = "SV")
  ) %>%
    group_by(gene, Source) %>%
    mutate(Variant_idx = dense_rank(Variant),
           Variant = paste0(gene, "_", Source, "_V", Variant_idx)) %>%
    ungroup() %>%
    select(-any_of(c("Variant_Type", "SVTYPE", "SVLEN", "Variant_idx")))
  
  # Map each sample to its population region
  sample_region_map <- data.frame(
    Sample = unlist(sample_region),
    Region = rep(names(sample_region), times = lengths(sample_region)),
    stringsAsFactors = FALSE
  )
  
  # Prepare haplo_input per region (Sample × Variant × Gene)
  haplo_input <- lapply(unique(sample_region_map$Region), function(region_name) {
    region_samples <- sample_region_map$Sample[sample_region_map$Region == region_name]
    region_data <- all_variants %>% filter(Sample %in% region_samples)
    
    complete_variants <- region_data %>%
      distinct(gene, Variant) %>%
      tidyr::expand(gene, Variant, Sample = region_samples)
    
    haplo_input_region <- complete_variants %>%
      left_join(region_data %>% select(gene, Variant, Sample, allele1, allele2, allele_count),
                by = c("gene", "Variant", "Sample")) %>%
      mutate(
        allele1 = ifelse(is.na(allele1), 0, allele1),
        allele2 = ifelse(is.na(allele2), 0, allele2),
        allele_count = ifelse(is.na(allele_count), 0, allele_count),
        Region = region_name
      )
    
    return(haplo_input_region)
  }) %>% bind_rows()
  
  message("Input prepared. Running haplotype EM per region...")
  gc()
  
  # Helper: analyze all genes in a single region
  process_region <- function(region_name) {
    message("\n=== Processing region: ", region_name, " ===")
    region_data <- haplo_input %>% filter(Region == region_name)
    
    region_genes <- unique(region_data$gene)
    region_results <- list()
    
    for (g in region_genes) {
      # Skip genes with excessive memory requirements
      res <- tryCatch({
        analyze_gene_haplotypes(g, region_data)
      }, error = function(e) {
        message("Error for gene ", g, ": ", e$message)
        return(NULL)
      })
      
      if (!is.null(res)) {
        res$Region <- region_name
        region_results[[g]] <- res
      }
      
      rm(res); gc()
    }
    
    if (length(region_results) > 0) {
      return(do.call(rbind, region_results))
    } else {
      return(NULL)
    }
  }
  
  # Run haplotype analysis for each region
  final_results <- lapply(unique(haplo_input$Region), process_region)
  final_results <- final_results[!sapply(final_results, is.null)]  # remove NULLs
  haplo_results_region <- do.call(rbind, final_results)
  
  # Save results to CSV
  readr::write_csv(haplo_results_region, output_file)
  
  return(list(haplo_input = haplo_input, haplo_results_region = haplo_results_region))
}

## graph per region

visualize_haplo_results_by_region <- function(haplo_results_region, haplo_input, locus, gene_name, output_dir) {
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  message(">>> Generating haplotype summary plot per region...")
  
  # Ensure key columns are character for joins and plotting
  haplo_results_region <- as_tibble(haplo_results_region) %>%
    mutate(Region = as.character(Region), gene = as.character(gene))
  haplo_input <- as_tibble(haplo_input) %>%
    mutate(Region = as.character(Region), gene = as.character(gene))
  locus <- as_tibble(locus)
  
  # Count number of samples per gene per region
  samples_per_gene_region <- haplo_input %>%
    group_by(Region, gene) %>%
    summarise(n_samples_gene = n_distinct(Sample), .groups = "drop")
  
  # Compute haplotype summary per gene and region
  # - sum frequencies per haplotype
  # - normalize frequencies
  # - compute Shannon diversity and number of haplotypes
  haplo_region_summary <- haplo_results_region %>%
    group_by(gene, Region, haplotype_id) %>%
    summarise(freq = sum(freq, na.rm = TRUE), .groups = "drop_last") %>%
    mutate(freq_norm = freq / sum(freq, na.rm = TRUE)) %>%
    summarise(
      n_haplotypes = n(),
      diversity = -sum(freq_norm * log(freq_norm + 1e-12)),  # avoid log(0)
      .groups = "drop"
    ) %>%
    left_join(samples_per_gene_region, by = c("Region", "gene")) %>%
    mutate(
      n_samples_gene = replace_na(n_samples_gene, 0),
      n_haplotypes_per_sample = ifelse(n_samples_gene > 0,
                                       n_haplotypes / n_samples_gene,
                                       NA_real_)
    )
  
  # Determine plotting width based on number of genes
  n_genes <- length(unique(haplo_region_summary$gene))
  n_regions <- length(unique(haplo_region_summary$Region))
  plot_width_reg <- max(8, 3 + n_genes * 0.25)
  
  # Generate bar plot of haplotypes per sample per region
  p_region <- ggplot(haplo_region_summary,
                     aes(
                       x = reorder(gene, -ifelse(is.na(n_haplotypes_per_sample), 0, n_haplotypes_per_sample)),
                       y = n_haplotypes_per_sample,
                       fill = diversity
                     )) +
    geom_col() +
    facet_wrap(~Region, scales = "fixed", nrow = n_regions) +  #y-scale across regions
    scale_fill_viridis_c(option = "plasma") +
    coord_cartesian(ylim = c(0, 1)) +  # y-axis limited 0–1
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Haplotype count per sample per region in", gene_name),
      x = "Gene",
      y = "Number of haplotypes per sample",
      fill = "Diversity\n(Shannon)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot to file
  out_file <- file.path(output_dir, paste0(gene_name, "_haplotype_region_summary_v2.png"))
  ggsave(out_file, p_region, width = plot_width_reg, height = 6, dpi = 300)
  
  return(
    summary_region = haplo_region_summary  # per-region haplotype summary
    )
}


### run


## Data processing

#1. Load data
locus <- load_gene_information(file.path(input_dir, paste0(gene_name,".xlsx")))
res <- load_region_information(file.path(input_dir, "SampleArea.xlsx"))

sample_region <- res$sampleArea
region_colors <- res$region_colors

#2. Process variants
snv_final <- indel_final <- sv_final <- NULL

snv_res <- process_snv_vcf(file.path(input_dir,"snv", paste0(gene_name,"_snv_vep_output.vcf")), locus, variant_missing_threshold, sample_missing_threshold)
snv_final <- snv_res$snv_final
total_samples_snv <- snv_res$n_samples_after

indel_res <- process_indel_vcf(file.path(input_dir,"indel",paste0(gene_name,"_indel_vep_output.vcf")), locus, variant_missing_threshold, sample_missing_threshold)
indel_final <- indel_res$indel_final
total_samples_indel <- indel_res$n_samples_after

sv_res <- process_sv_vcf(file.path(input_dir, "sv", paste0(gene_name, "_sv_vep_output.vcf")), locus, variant_missing_threshold, sample_missing_threshold)
sv_final <- sv_res$sv_final
total_samples_sv <- sv_res$n_samples_after

if (1 %in% steps) {
  message("Running Step 1: Variability score")
  #Gene conservation scores
  cat("SNV samples:", total_samples_snv, "\n")
  cat("INDEL samples:", total_samples_indel, "\n")
  cat("SV samples:", total_samples_sv, "\n")
  
  gene_scores <- calculate_gene_conservation_score(
    snv_final = snv_final,
    indel_final = indel_final,
    sv_final = sv_final,
    gene_locus = locus,
    total_samples_snv = total_samples_snv,
    total_samples_indel = total_samples_indel,
    total_samples_sv = total_samples_sv,
    use_length_factor     = use_length_factor,
    normalize_by_samples  = normalize_by_samples
  )
  
  score_dir <- file.path(output_dir,"score")
  if (!dir.exists(score_dir)) dir.create(score_dir, recursive=TRUE)
  
  plot_gene_variation_scores(gene_scores, gene_name, score_dir)
  message("Step 1 finished.")
}

if (2 %in% steps) {
  if (!exists("gene_scores")) stop("Step 1 must be run before evolutionary comparison")
  message("Running Step 2: Comparison with evolutionary metrics")
  
  #1. Define base directory for MSA files
  msa_dir <- opt$msa_dir
  
  #2. Find all MSA files for the gene (case-insensitive)
  all_files <- list.files(
    msa_dir,
    pattern = "\\.(fas|fasta)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  msa_files <- all_files[grepl(paste0("^", gene_name), basename(all_files), ignore.case = TRUE)]
  
  if (length(msa_files) == 0) stop("No MSA files found for gene: ", gene_name, " in ", msa_dir)
  
  #3. Extract segments (basename without extension)
  segments <- tools::file_path_sans_ext(basename(msa_files))
  
  cat("Found MSA files:\n"); print(msa_files)
  cat("Segments:\n"); print(segments)
  
  #4. Directories with distance and dN/dS matrices
  dist_dir <- file.path(msa_dir, "distance")
  dnds_dir <- file.path(msa_dir, "dn-ds")
  
  has_distance <- dir.exists(dist_dir)
  has_dnds     <- dir.exists(dnds_dir)
  
  if (!has_distance && !has_dnds)
    stop("Neither distance nor dN/dS directories are present.")
  
  #5. Read matrices
  matrices <- read_matrices(dist_dir, dnds_dir, segments)
  
  #6. Compare genes
  results <- compare_genes_core(msa_files, dnds_data = matrices$dnds, dist_data = matrices$dist)
  
  #7. Merge with gene scores
  merged_data <- gene_scores %>% left_join(results, by = c("gene" = "human_gene"))
  
  #8. Create output directory
  compare_dir <- file.path(output_dir, "compare")
  if (!dir.exists(compare_dir)) dir.create(compare_dir, recursive = TRUE)
  
  #9. Plot comparisons
  analyze_gene_conservation(
    complete_data = merged_data,
    gene_name = gene_name,
    output_dir = compare_dir,
    detailed = detailed
  )
  
  plot_dnds_comparison(merged_data, gene_name, compare_dir)
  plot_distance_comparison(merged_data, gene_name, compare_dir)
  plot_combined_conservation_comparison(merged_data, gene_name, compare_dir)
  message("Step 2 finished.")
}

if (3 %in% steps) {
  message("Running Step 3: Haplotypes")
  haplo_dir <- file.path(output_dir, "haplotypes")
  if (!dir.exists(haplo_dir)) dir.create(haplo_dir, recursive = TRUE)
  temp <- file.path(output_dir, "temp")
  if (!dir.exists(temp)) dir.create(temp, recursive = TRUE)
  
  # Full analysis (all regions merged)
  haplo_out <- run_haplotype_analysis(
    snv_final, indel_final, sv_final, sample_region,
    output_file = file.path(temp,
                            paste0("haplo_results_partial_", gene_name, ".csv"))
  )
  
  analyze_haplotype_variation(
    haplo_out$haplo_input,
    haplo_out$haplo_results,
    locus,
    gene_name,
    haplo_dir
  )
}

if (3 %in% steps) {
  # By-region analysis
  haplo_out_region <- run_haplotype_analysis_by_region(
    snv_final,
    indel_final,
    sv_final,
    sample_region,
    output_file = file.path(output_dir, "temp",
                            paste0("haplo_by_region_", gene_name, ".csv"))
  )
  
  visualize_haplo_results_by_region(
    haplo_results_region = haplo_out_region$haplo_results_region,
    haplo_input = haplo_out_region$haplo_input,
    locus = locus,
    gene_name = gene_name,
    output_dir = file.path(output_dir, "haplotypes")
  )
  message("Step 3 finished.")
}

