### 0. Set up


suppressPackageStartupMessages(library(optparse))
#Argument parser with help
option_list <- list(
  make_option(c("-g","--gene"), type="character", default="GENE",
              help = "Gene basename (default: %default)"),
  
  make_option(c("-i","--input_dir"), type="character", default=".",
              help = paste(
                "Input directory containing all required files.",
                "Expected files:",
                "\n  ├── <GENE>.xlsx                  # gene locus coordinates",
                "├── SampleArea.xlsx              # sample-to-region mapping",
                "├── indel/                          # INDEL subdirectory",
                "│     ├── <GENE>_indel_vep_output.vcf",
                "│     └── <GENE>_indel_raw.vcf",
                "",
                "If the structure differs, adjust input paths accordingly.",
                sep = "\n"
              )),
  
  make_option(c("-o","--output_dir"), type="character", default="indel_results",
              help = "Output directory for generated plots (default: %default)"),
  
  make_option(c("-s","--step"), type="character", default="1,2,3,4,5,6",
              help = paste(
                "Pipeline steps to run (comma-separated).",
                "Available steps:",
                "  1 = Variant burden heatmap",
                "  2 = INDEL bubble plot",
                "  3 = PTV vs non-PTV bubble comparison",
                "  4 = Summary plots (all consequences)",
                "  5 = Summary plots (non-modifier only)",
                "  6 = Variant track map",
                "Example: --step 1,3,5",
                "(default: %default)",
                sep = "\n"
              )),
  
  #OPTIONAL ARGUMENTS
  make_option("--min_size", type="integer", default=1,
              help = "Minimum bubble size for bubble plots (default: %default)"),
  
  make_option("--max_size", type="integer", default=10,
              help = "Maximum bubble size for bubble plots (default: %default)"),
  
  make_option("--point_alpha", type="double", default=0.4,
              help="Alpha transparency for INDEL tracks (default: %default)"),
  
  make_option("--indel_colors", type="character", default=NULL,
              help="Comma-separated list of INDEL colors (default: %default)")
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
min_size <- opt$min_size
max_size <- opt$max_size
point_alpha <- opt$point_alpha

if (!is.null(opt$indel_colors)) {
  # Parse "DEL=#e41a1c,INS=#377eb8, INV"=#4daf4a
  indel_colors <- sapply(
    strsplit(opt$indel_colors, ",")[[1]],
    function(x) {
      kv <- strsplit(x, "=")[[1]]
      setNames(kv[2], kv[1])
    }
  )
} else {
  # default palette
  indel_colors <- c("DEL"="#e41a1c", "INS"="#377eb8", "INV"="#4daf4a")
}

#Print summary
cat("Gene:", gene_name, "\n")
cat("Input dir:", input_dir, "\n")
cat("Output dir:", output_dir, "\n")
cat("Steps to run:", paste(steps, collapse=", "), "\n")
cat("min_size:", min_size, "\n")
cat("max_size:", max_size, "\n")
cat("point_alpha:", point_alpha, "\n")
cat("indel_colors:", indel_colors, "\n")


### Libraries


suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(viridis)
  library(tidyverse)
  library(readxl)
  library(data.table)
  library(viridis)
  library(reshape2)
  library(GenomicRanges)
  library(VariantAnnotation)
  library(dplyr)
  library(tibble)
  library(RColorBrewer)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(purrr)
  library(cowplot)
  library(ggnewscale)
  library(stringr)
  library(patchwork)
  library(ggtext)
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


### Functions


##Annotate

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
  
  ##Extract SVLEN
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

#Same hierarchy but excluding MODIFIER variants
variant_hierarchy_no_modifier <- variant_hierarchy[1:19]

#Variant colors
variant_colors <- c(
  #HIGH impact
  "Transcript ablation"       = "#99000D",
  "Splice acceptor"           = "#D7301F",
  "Splice donor"              = "#EF6548",
  "Stop gained"               = "#FB6A4A",
  "Frameshift"                = "#FC9272",
  "Stop lost"                 = "#FCBBA1",
  "Start lost"                = "#FEE0D2",
  
  #MODERATE impact
  "Inframe insertion"         = "#DD8D1F",
  "Inframe deletion"          = "#E39F3A",
  "Missense"                  = "#6A51A3",
  "Protein altering"          = "#9E9AC8",
  
  #LOW impact
  "Splice region"             = "#FDBF6F",
  "Incomplete terminal codon" = "#FDD49E",
  "Start retained"            = "#CAB2D6",
  "Stop retained"             = "#BDBDBD",
  "Synonymous"                = "#BC80BD",
  
  #MODIFIER impact
  "Coding sequence"           = "#80CDC1",
  "Mature miRNA variant"      = "#35978F",
  "UTR"                       = "#8DD3C7",
  "Non coding transcript exon"= "#A6CEE3",
  "Intron"                    = "#B2DF8A",
  "NMD transcript"            = "#63A598",
  "Non coding transcript"     = "#89C5DF",
  "Coding transcript"         = "#A1D9E8",
  "Upstream gene"             = "#33A02C",
  "Downstream gene"           = "#66A61E",
  "TF binding site"           = "#41B6C4",
  "Regulatory region"         = "#74C476",
  "Intergenic"                = "#FFFF99",
  "Sequence"                  = "#F7F7F7",
  
  #fallback
  "Other"                     = "#CCCCCC"
)

#Generate color palettes for each variant (white to variant color gradient)
variant_palettes <- lapply(variant_colors, function(col) colorRampPalette(c("#ffffff", col))(100))

## 1. Heatmap (log1p counts)

run_indel_heatmap <- function(
    gene_loc,        # data.frame with gene locations: Chr, Start, End, gene_name, Type
    indel_data,      # data.frame with INDELs: CHROM, POS, REF, ALT, INFO, SVTYPE, SVLEN, sample columns...
    sample_region,   # list of samples per region, e.g., list(EUR = c(...), AMR = c(...))
    region_colors,   # vector of colors for regions, names = region
    output_dir,      # output directory
    gene_name = NULL # optional, used for plot title
) {
  if(is.null(gene_name)) gene_name <- unique(gene_loc$Type)
  message("Running INDEL heatmap for gene: ", gene_name)
  
  # 1. Identify sample columns
  fixed_cols <- c("CHROM","POS","REF","ALT","INFO","SVTYPE","SVLEN")
  sample_cols <- setdiff(colnames(indel_data), fixed_cols)
  indel_data <- indel_data[, c(fixed_cols, sample_cols)]
  
  # 2. Genotype converter
  convert_genotype <- function(gt) {
    if (gt %in% c("0|0", ".|.", "0/0", "./.")) return(0)
    if (gt %in% c("0|1","1|0","0/1","1/0",".|1","1|.")) return(1)
    if (gt %in% c("1|1","1/1")) return(2)
    return(0)
  }
  
  # 3. Assign INDELs to genes using GRanges
  gene_gr <- GRanges(
    seqnames = paste0("chr", gene_loc$Chr),
    ranges = IRanges(start = gene_loc$Start, end = gene_loc$End),
    gene_id = gene_loc$gene_name
  )
  
  # Compute start/end positions for INDEL
  indel_data <- indel_data %>%
    mutate(
      CHROM = ifelse(grepl("^chr", CHROM), CHROM, paste0("chr", CHROM)),
      indel_start = POS,
      indel_end = ifelse(SVTYPE == "INS", POS, POS + abs(SVLEN)),
      indel_start = pmin(indel_start, indel_end),
      indel_end = pmax(indel_start, indel_end),
      indel_id = seq_len(nrow(indel_data))
    )
  
  # GRanges for DEL/INV and INS
  gr_del_inv <- GRanges(
    seqnames = indel_data$CHROM[indel_data$SVTYPE %in% c("DEL","INV")],
    ranges = IRanges(
      start = indel_data$indel_start[indel_data$SVTYPE %in% c("DEL","INV")],
      end   = indel_data$indel_end[indel_data$SVTYPE %in% c("DEL","INV")]
    ),
    indel_id = indel_data$indel_id[indel_data$SVTYPE %in% c("DEL","INV")]
  )
  
  gr_ins <- GRanges(
    seqnames = indel_data$CHROM[indel_data$SVTYPE == "INS"],
    ranges = IRanges(
      start = indel_data$indel_start[indel_data$SVTYPE == "INS"],
      end   = indel_data$indel_end[indel_data$SVTYPE == "INS"]
    ),
    indel_id = indel_data$indel_id[indel_data$SVTYPE == "INS"]
  )
  
  # 4. Find overlaps
  hits_del_inv <- findOverlaps(gr_del_inv, gene_gr)
  hits_ins <- findOverlaps(gr_ins, gene_gr, type="within")
  
  matches <- bind_rows(
    data.frame(indel_id = mcols(gr_del_inv)$indel_id[queryHits(hits_del_inv)],
               gene_id = mcols(gene_gr)$gene_id[subjectHits(hits_del_inv)], stringsAsFactors=FALSE),
    data.frame(indel_id = mcols(gr_ins)$indel_id[queryHits(hits_ins)],
               gene_id = mcols(gene_gr)$gene_id[subjectHits(hits_ins)], stringsAsFactors=FALSE)
  ) %>% distinct()
  
  if(nrow(matches) == 0) {
    stop("No INDELs overlap genes with the chosen logic.")
  }
  
  # 5. Join back to indel_data
  indel_matched <- indel_data %>% inner_join(matches, by="indel_id")
  
  # Keep only relevant columns
  keep_cols <- c("indel_id","CHROM","POS","SVTYPE","SVLEN","gene_id", sample_cols)
  indel_matched <- indel_matched[, keep_cols]
  
  # 6. Melt to long format and convert genotypes
  long_dt <- data.table::melt(as.data.table(indel_matched),
                              id.vars = c("indel_id","CHROM","POS","SVTYPE","SVLEN","gene_id"),
                              measure.vars = sample_cols,
                              variable.name = "Sample",
                              value.name = "GT")
  long_dt[, mut_count := vapply(GT, convert_genotype, numeric(1))]
  
  # Aggregate mutations per gene and sample
  result_dt <- long_dt[, .(mutations = sum(mut_count)), by = .(gene_id, Sample)]
  
  # 7. Assign populations
  sample_population <- unlist(lapply(names(sample_region), function(r) {
    setNames(rep(r, length(sample_region[[r]])), sample_region[[r]])
  }))
  result_dt$Population <- sample_population[result_dt$Sample]
  
  # Order samples and genes
  result_dt$Sample <- factor(result_dt$Sample, levels = unlist(sample_region))
  gene_order <- gene_loc %>% arrange(Chr, Start) %>% pull(gene_name)
  result_dt$gene_id <- factor(result_dt$gene_id, levels = gene_order)
  
  sample_regions <- data.frame(
    Sample = unlist(sample_region),
    Population = rep(names(sample_region), lengths(sample_region))
  )
  
  # 8. Plot heatmap
  p <- ggplot(result_dt, aes(x = gene_id, y = Sample, fill = mutations)) +
    geom_tile(color = NA) +
    scale_fill_viridis_c(option = "plasma", name = "Mutations") +
    theme_minimal() +
    labs(title = paste("INDEL Heatmap for", gene_name),
         x = "Gene", y = "Sample") +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          axis.text.y = element_text(size=7),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(size=14, face="bold", hjust=0.5))
  
  # Color sample names according to population
  sample_order <- levels(result_dt$Sample)
  sample_colors <- region_colors[sample_regions$Population][match(sample_order, sample_regions$Sample)]
  p <- p + theme(axis.text.y = element_markdown(size=7, color = sample_colors))
  
  # Add hidden points for region legend
  p <- p + geom_point(data=sample_regions,
                      aes(x=1, y=Sample, color=Population),
                      inherit.aes = FALSE, size=0, alpha=0) +
    scale_color_manual(values = region_colors, name = "Population",
                       guide = guide_legend(override.aes=list(size=3, alpha=1)))
  
  # 9. Save plot
  n_genes <- length(unique(result_dt$gene_id))
  n_samples <- length(unique(result_dt$Sample))
  width_dynamic <- max(6, n_genes*0.20)
  height_dynamic <- max(4, n_samples*0.13)
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  outfile <- file.path(output_dir, paste0("indel_heatmap_", gene_name, ".png"))
  ggsave(filename = outfile, plot = p, width=width_dynamic, height=height_dynamic, dpi=300, 
         limitsize = FALSE)
  message("INDEL heatmap saved to: ", outfile)
  
  return(p)
}

## 2. Bubble map visualization

plot_indel_bubble <- function(
    indel_result,       #data.frame with INDEL results, must contain columns: Sample, Position, Variant_Type
    sample_region,    #list of samples per region, e.g., list(EUR = c(...), AMR = c(...))
    gene_locus,       #data.frame with gene coordinates, columns: gene, chromosome, start, end
    indel_type,         #vector of variant types to include, e.g., c("Missense", "Stop gained")
    color_palette,    #vector of colors for gradient fill
    min_size = 1,     #minimum bubble size
    max_size = 10,    #maximum bubble size
    title = "INDEL Distribution by Region and Gene",
    output_file = NULL
) {
  #1. Prepare gene annotations
  colnames(gene_locus) <- c("gene","chromosome","start","end")
  gene_locus <- gene_locus %>%
    mutate(
      chromosome = as.character(chromosome),
      start = as.numeric(start),
      end = as.numeric(end)
    )
  
  #2. Preprocess variant table (filter by type, extract coordinates)
  indel_long <- indel_result %>%
    filter(Variant_Type %in% indel_type) %>%
    #split Position into chromosome and pos_info (works for "chr:start" or "chr:start-end")
    separate(Position, into = c("chromosome","pos_info"), sep = ":", remove = FALSE, fill = "right") %>%
    mutate(
      #handle possible "start-end" or just "start"
      position = as.numeric(str_extract(pos_info, "^[0-9]+")),
      SVLEN_raw = as.numeric(SVLEN),
      SVLEN = ifelse(is.na(SVLEN_raw), 0, SVLEN_raw),
      SVLEN = abs(SVLEN),               #<- IMPORTANT: use abs for length
      chromosome = as.character(chromosome),
      SVTYPE = as.character(SVTYPE),
      Sample = as.character(Sample),
      #compute indel_start / indel_end: for INS it is a point; for DEL/INV use position -> position + abs(SVLEN)
      indel_start = case_when(
        SVTYPE %in% c("DEL","INV") ~ position,
        SVTYPE == "INS" ~ position,
        TRUE ~ position
      ),
      indel_end = case_when(
        SVTYPE %in% c("DEL","INV") ~ position + SVLEN,
        SVTYPE == "INS" ~ position,
        TRUE ~ position + SVLEN
      ),
      #ensure ordering
      indel_start = pmin(indel_start, indel_end),
      indel_end   = pmax(indel_start, indel_end)
    )
  
  if (nrow(indel_long) == 0) {
    message(
      "No variants found for INDEL type: ",
      paste(indel_type, collapse = ", ")
    )
    return(NULL)
  }
  
  indel_long$chromosome <- sub("^chr", "", indel_long$chromosome)
  gene_locus$chromosome <- sub("^chr", "", gene_locus$chromosome)
  
  #3. Match variants to genes
  gr_genes <- GRanges(
    seqnames = gene_locus$chromosome,
    ranges = IRanges(start = gene_locus$start, end = gene_locus$end),
    gene = gene_locus$gene
  )
  
  #DEL/INV ranges
  idx_del_inv <- which(indel_long$SVTYPE %in% c("DEL","INV"))
  gr_del_inv <- GRanges(
    seqnames = indel_long$chromosome[idx_del_inv],
    ranges = IRanges(
      start = indel_long$indel_start[idx_del_inv],
      end   = indel_long$indel_end[idx_del_inv]
    ),
    Sample = indel_long$Sample[idx_del_inv],
    SVTYPE = indel_long$SVTYPE[idx_del_inv],
    Variant_Type = indel_long$Variant_Type[idx_del_inv]
  )
  
  #INS ranges (point)
  idx_ins <- which(indel_long$SVTYPE == "INS")
  gr_ins <- GRanges(
    seqnames = indel_long$chromosome[idx_ins],
    ranges = IRanges(
      start = indel_long$indel_start[idx_ins],
      end   = indel_long$indel_start[idx_ins]
    ),
    Sample = indel_long$Sample[idx_ins],
    SVTYPE = indel_long$SVTYPE[idx_ins],
    Variant_Type = indel_long$Variant_Type[idx_ins]
  )
  
  
  #4. Find overlaps
  
  hits_del_inv <- findOverlaps(gr_del_inv, gr_genes)
  hits_ins <- findOverlaps(gr_ins, gr_genes, type="within")
  
  if(length(hits_del_inv) + length(hits_ins) == 0){
    message(
      "No variants match gene coordinates for INDEL type: ",
      paste(indel_type, collapse = ", ")
    )
    return(NULL)
  }
  
  df_del_inv <- if (length(hits_del_inv) > 0) {
    data.frame(
      as.data.frame(mcols(gr_del_inv))[queryHits(hits_del_inv), , drop = FALSE],
      gene = mcols(gr_genes)$gene[subjectHits(hits_del_inv)],
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        indel_start = start(gr_del_inv)[queryHits(hits_del_inv)],
        indel_end   = end(gr_del_inv)[queryHits(hits_del_inv)]
      )
  } else tibble::tibble()
  
  df_ins <- if (length(hits_ins) > 0) {
    data.frame(
      as.data.frame(mcols(gr_ins))[queryHits(hits_ins), , drop = FALSE],
      gene = mcols(gr_genes)$gene[subjectHits(hits_ins)],
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        indel_start = start(gr_ins)[queryHits(hits_ins)],
        indel_end   = end(gr_ins)[queryHits(hits_ins)]
      )
  } else tibble::tibble()
  
  matched <- bind_rows(df_del_inv, df_ins)
  
  if (nrow(matched) == 0) {
    message(
      "No variants match gene coordinates for INDEL type: ",
      paste(indel_type, collapse = ", ")
    )
    return(NULL)
  }
  
  #4. Map samples to regions
  sample_region_map <- data.frame(
    Sample = unlist(sample_region),
    Region = rep(names(sample_region), times = lengths(sample_region)),
    stringsAsFactors = FALSE
  )
  matched <- matched %>% left_join(sample_region_map, by = "Sample")
  
  #5. Summarize variants per gene and per region
  region_total_samples <- lengths(sample_region)
  
  summary_data <- matched %>%
    mutate(Region = as.character(Region)) %>%
    group_by(gene, Region, Variant_Type) %>%
    summarise(
      mutated_samples = n_distinct(Sample),
      total_mutations = n(),
      .groups = "drop"
    ) %>%
    mutate(
      region_total = as.numeric(region_total_samples[Region]),
      region_total = ifelse(is.na(region_total) | region_total == 0, NA, region_total),
      proportion = mutated_samples / region_total,
      mutations_per_sample = total_mutations / region_total
    ) %>%
    filter(!is.na(region_total))
  
  summary_data <- summary_data %>% filter(proportion > 0)
  
  if (nrow(summary_data) == 0) {
    cat(
      "No data available for plotting for variant type: ",
      paste(indel_type, collapse = ", "),
      "\n"
    )
    return(NULL)
  }
  
  #6. Prepare data for plotting
  n_genes <- length(unique(summary_data$gene))
  n_regions <- length(unique(summary_data$Region))
  
  plot_width <- max(5, 3 + n_genes * 0.3)
  plot_height <- max(4, 2.2 + n_regions * 0.4)
  
  #7. Generate bubble plot
  p <- ggplot(summary_data, aes(x = gene, y = Region)) +
    geom_point(
      aes(size = proportion, fill = mutations_per_sample),
      shape = 21, color = "black", alpha = 0.8
    ) +
    scale_size_continuous(
      name = "Proportion of samples",
      range = c(min_size, max_size),
      labels = scales::percent_format(accuracy = 1),
      guide = guide_legend(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        ncol = 2
      )
    ) +
    scale_fill_gradientn(
      name = "Mutations per sample",
      colours = color_palette,
      values = scales::rescale(seq(0, 1, length.out = length(color_palette)))
    ) +
    labs(title = title, x = "Gene", y = "Region") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
  
  #8. Save to file (optional)
  if (!is.null(output_file) && nzchar(output_file)) {
    ggsave(output_file, p, width = plot_width, height = plot_height,
           dpi = 300, limitsize = FALSE)
    cat("Plot saved to:", output_file, "\n")
  }
  
  return(p)
}

## 3. PTV - bubble track

plot_indel_bubble_dual <- function(
    indel_result, 
    sample_region, 
    gene_locus,
    min_size = 1,
    max_size = 10,
    title = "INDEL Mutation Distribution by Region and Gene",
    output_file = NULL
) {
  #1. Prepare gene info
  colnames(gene_locus) <- c("gene", "chromosome", "start", "end")
  gene_locus <- gene_locus %>%
    mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      gene = as.character(gene),
      chromosome = as.character(chromosome)
    )
  
  #2. Process INDELs
  indel_long <- indel_result %>%
    separate(Position, into = c("chromosome", "position_info"), sep = ":", remove = FALSE) %>%
    mutate(
      position = as.numeric(str_extract(position_info, "^[0-9]+")),
      chromosome = as.character(chromosome),
      SVLEN = as.numeric(SVLEN),
      SVLEN = ifelse(is.na(SVLEN), 0, SVLEN),
      indel_start = position,
      indel_end = case_when(
        SVLEN > 0 ~ position + SVLEN - 1,
        SVLEN < 0 ~ position + abs(SVLEN) - 1,
        TRUE ~ position
      ),
      PTV_status = ifelse(
        Variant_Type %in% c("Frameshift", "Transcript ablation", "Stop gained", 
                            "Stop lost", "Start lost"),
        "PTV",
        "non-PTV"
      ),
      # optionally remove haplotype suffix for matching sample_region
      Sample = sub("-h[12]$", "", Sample)
    )
  
  # Debug: liczba wierszy
  message("Total INDELs processed: ", nrow(indel_long))
  
  #3. Match INDELs to genes
  matched_variants <- map_dfr(1:nrow(gene_locus), function(i) {
    gene_row <- gene_locus[i, ]
    indel_long %>%
      filter(
        chromosome == gene_row$chromosome,
        indel_start <= gene_row$end & indel_end >= gene_row$start
      ) %>%
      mutate(gene = gene_row$gene)
  })
  
  message("Variants matched to genes: ", nrow(matched_variants))
  
  if(nrow(matched_variants) == 0) {
    message("No INDELs overlap genes. Check coordinates or SVLEN.")
    return(NULL)
  }
  
  #4. Map samples to regions
  sample_region_map <- data.frame(
    Sample = unlist(sample_region),
    Region = rep(names(sample_region), times = lengths(sample_region))
  )
  
  matched_variants <- matched_variants %>%
    left_join(sample_region_map, by = "Sample") %>%
    filter(!is.na(Region))
  
  message("Variants after sample_region mapping: ", nrow(matched_variants))
  if(nrow(matched_variants) == 0) {
    message("No variants match sample_region. Check Sample names or haplotype suffix.")
    return(NULL)
  }
  
  #5. Summarize per gene × region × PTV_status
  region_total_samples <- lengths(sample_region)
  
  summary_data <- matched_variants %>%
    group_by(gene, Region, PTV_status) %>%
    summarise(
      mutated_samples = n_distinct(Sample),
      total_mutations = n(),
      .groups = "drop"
    ) %>%
    mutate(
      region_total = region_total_samples[Region],
      proportion = mutated_samples / region_total,
      mutations_per_sample = total_mutations / region_total
    )
  
  #6. Fill missing combinations
  all_combinations <- expand.grid(
    gene = unique(summary_data$gene),
    Region = unique(summary_data$Region),
    PTV_status = c("PTV", "non-PTV")
  )
  
  summary_data <- all_combinations %>%
    left_join(summary_data, by = c("gene", "Region", "PTV_status")) %>%
    mutate(
      mutated_samples = ifelse(is.na(mutated_samples), 0, mutated_samples),
      total_mutations = ifelse(is.na(total_mutations), 0, total_mutations),
      proportion = ifelse(is.na(proportion), 0, proportion),
      mutations_per_sample = ifelse(is.na(mutations_per_sample), 0, mutations_per_sample),
      region_total = ifelse(is.na(region_total), 1, region_total)
    )
  
  #7. Prepare dual-track Y positions
  regions <- unique(summary_data$Region)
  y_positions <- data.frame(
    Region = rep(regions, each = 2),
    PTV_status = rep(c("PTV", "non-PTV"), times = length(regions)),
    y_pos = 1:(length(regions) * 2)
  )
  
  summary_data <- summary_data %>% left_join(y_positions, by = c("Region", "PTV_status"))
  
  # Filter zero rows
  plot_data <- summary_data %>% filter(proportion > 0)
  if(nrow(plot_data) == 0) {
    message("No INDEL variants to plot after filtering zero proportion.")
    return(NULL)
  }
  
  #8. Plot
  p <- ggplot(plot_data, aes(x = gene, y = y_pos)) +
    geom_point(aes(size = proportion, fill = mutations_per_sample),
               shape = 21, color = "black", alpha = 0.8) +
    scale_fill_viridis_c(name = "Mutations per sample", option = "magma") +
    scale_size_continuous(name = "Proportion of samples",
                          range = c(min_size, max_size),
                          labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(
      breaks = y_positions$y_pos,
      labels = paste0(y_positions$Region, " - ", y_positions$PTV_status),
      expand = expansion(mult = 0.1)
    ) +
    scale_x_discrete(expand = expansion(mult = 0.01)) +
    labs(y = NULL, x = "Gene", title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
  
  #9. Save plot if requested
  n_genes <- length(unique(plot_data$gene))
  n_regions <- length(unique(plot_data$Region))
  
  plot_width <- max(6, 3 + n_genes * 0.3)
  
  # Każdy region -> 2 wiersze (PTV i non-PTV)
  plot_height <- max(4, 2.2 + n_regions * 0.6)
  
  if (!is.null(output_file)) {
    ggsave(output_file, p,
           width = plot_width,
           height = plot_height,
           dpi = 300, limitsize = FALSE)
    message("Plot saved to: ", output_file)
  }
  
  
  return(p)
}

## 4. Summary plots

generate_indel_summary_plots <- function(
    indel_result,
    result_long,
    sample_region,
    region_colors,
    variant_colors,
    gene_name,
    output_prefix = NULL,
    variant_types = NULL
) {
  
  message("Preparing INDEL summary plots")
  
  # 1. Prepare input data
  plot_data <- result_long %>%
    filter(!is.na(Variant_Type)) %>%
    mutate(
      chromosome = sub(":.*", "", Position),
      position   = as.numeric(sub(".*:(\\d+)_.*", "\\1", Position))
    )
  
  if (!is.null(variant_types)) {
    plot_data <- plot_data %>% filter(Variant_Type %in% variant_types)
  }
  
  # sample_region = list(region → vector(sample))
  sample_to_region <- stack(sample_region)
  colnames(sample_to_region) <- c("Sample", "Region")
  
  plot_data <- plot_data %>% left_join(sample_to_region, by = "Sample")
  
  # order samples
  ordered_samples <- unlist(sample_region)
  plot_data$Sample <- factor(plot_data$Sample, levels = rev(ordered_samples))
  plot_data$Sample <- droplevels(plot_data$Sample)
  
  # matching colors
  type_colors <- variant_colors[names(variant_colors) %in% unique(plot_data$Variant_Type)]
  
  # pretty labels
  plot_data$position_label <- paste0(
    plot_data$chromosome, ":",
    format(plot_data$position, big.mark = ",")
  )
  
  plot_data <- plot_data %>% arrange(chromosome, position)
  plot_data$Position <- factor(plot_data$Position, levels = unique(plot_data$Position))
  
  # region colors for y-axis
  sample_regions <- plot_data %>% distinct(Sample, Region)
  y_label_colors <- region_colors[
    as.character(
      sample_regions$Region[
        match(levels(plot_data$Sample), sample_regions$Sample)
      ]
    )
  ]
  
  # 2. Main INDEL variant track plot
  p_main <- ggplot(plot_data, aes(x = Position, y = Sample, fill = Variant_Type)) +
    geom_tile(color = NA) +
    scale_fill_manual(values = type_colors, name = "Variant Type") +
    labs(
      title = paste(gene_name, "INDEL Variant Track"),
      x = "Genomic Position",
      y = "Sample",
      caption = paste("Total variants:", nrow(plot_data))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_markdown(size = 9, color = y_label_colors),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid = element_blank()
    )
  
  #add invisible points to create Region legend in the main plot
  p_main <- p_main +
    geom_point(
      data = sample_regions,
      aes(x = 1, y = Sample, color = Region),
      inherit.aes = FALSE,
      size = 0, alpha = 0
    ) +
    scale_color_manual(
      values = region_colors,
      name = "Region",
      guide = guide_legend(override.aes = list(size = 4, alpha = 1))
    )
  
  # region boundaries
  region_boundaries <- plot_data %>%
    group_by(Region) %>%
    summarise(
      min_idx = min(as.numeric(Sample)),
      max_idx = max(as.numeric(Sample)),
      .groups = "drop"
    ) %>%
    arrange(desc(max_idx))
  
  for (i in 1:(nrow(region_boundaries) - 1)) {
    p_main <- p_main +
      geom_hline(
        yintercept = region_boundaries$min_idx[i] - 0.5,
        color = "black", linetype = "dashed", linewidth = 0.5
      )
  }
  
  # position density adjustment
  n_pos <- length(unique(plot_data$Position))
  all_labels <- plot_data %>% distinct(Position, position_label)
  
  if (n_pos > 50) {
    idx <- seq(1, n_pos, length.out = 50)
    p_main <- p_main +
      scale_x_discrete(
        breaks = all_labels$Position[idx],
        labels = all_labels$position_label[idx]
      )
  } else {
    p_main <- p_main +
      scale_x_discrete(
        labels = setNames(all_labels$position_label, all_labels$Position)
      )
  }
  
  # 3. Sample summary plot
  sample_summary_data <- plot_data %>%
    count(Sample, Variant_Type, Region)
  
  x_levels <- unique(sample_summary_data$Sample)
  
  sample_regions_x <- plot_data %>%
    distinct(Sample, Region) %>%
    filter(Sample %in% x_levels) %>%
    mutate(Sample = factor(Sample, levels = x_levels)) %>%
    arrange(Sample)
  
  x_label_colors <- region_colors[sample_regions_x$Region]
  
  p_samples <- ggplot(sample_summary_data,
                      aes(x = Sample, y = n, fill = Variant_Type)) +
    geom_col() +
    scale_fill_manual(values = type_colors, name = "Variant Type") +
    labs(
      title = paste(gene_name, "— INDEL Count by Sample"),
      x = "Sample",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1, color = x_label_colors)
    ) +
    geom_point(
      data = sample_regions_x,
      aes(x = Sample, y = 0, color = Region),
      inherit.aes = FALSE, alpha = 0
    ) +
    scale_color_manual(
      values = region_colors,
      name = "Region",
      guide = guide_legend(override.aes = list(size = 4, alpha = 1))
    )
  
  # 4. Type distribution plot
  p_types <- plot_data %>%
    count(Variant_Type) %>%
    ggplot(aes(x = reorder(Variant_Type, n), y = n, fill = Variant_Type)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.4, size = 3) +
    scale_fill_manual(values = type_colors) +
    labs(
      title = paste(gene_name, "— INDEL Type Distribution"),
      x = "Variant Type",
      y = "Count"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # 5. region-level plot
  region_plot_data <- plot_data %>%
    group_by(Region, Sample, Variant_Type) %>%
    summarise(n_mut = n(), .groups = "drop") %>%
    group_by(Region, Variant_Type) %>%
    summarise(mean_mut = mean(n_mut), .groups = "drop")
  
  p_region <- ggplot(region_plot_data,
                     aes(x = Region, y = mean_mut, fill = Variant_Type)) +
    geom_col() +
    scale_fill_manual(values = type_colors) +
    labs(
      title = paste(gene_name, "— Average INDELs per Sample by Region"),
      x = "Region",
      y = "Mean INDELs per sample"
    ) +
    theme_minimal()
  
  # 6. Save outputs
  n_samples  <- length(unique(plot_data$Sample))
  n_positions <- length(unique(plot_data$Position))
  n_types <- length(unique(plot_data$Variant_Type))
  
  w_main <- max(10, min(30, n_positions * 0.05))
  h_main <- max(6,  min(25, n_samples * 0.18))
  
  h_samples <- max(5, min(20, n_samples * 0.14))
  w_samples <- max(10, min(20, n_samples * 0.2))
  
  w_types <- max(6, min(16, n_types * 0.5))
  h_types <- 6
  
  w_regions <- max(6, min(16, n_types * 0.5))
  h_regions <- 6
  
  if (!is.null(output_prefix)) {
    ggsave(paste0(output_prefix, "_main.png"),    p_main,    width = w_main,    height = h_main,  dpi = 300, limitsize = FALSE)
    ggsave(paste0(output_prefix, "_samples.png"), p_samples, width = w_samples, height = h_samples, dpi = 300, limitsize = FALSE)
    ggsave(paste0(output_prefix, "_types.png"),   p_types,   width = w_types,   height = h_types, dpi = 300, limitsize = FALSE)
    ggsave(paste0(output_prefix, "_regions.png"), p_region, width = w_regions, height = h_regions, dpi = 300, limitsize = FALSE)
  }
  
  return(list(
    main_plot = p_main,
    sample_plot = p_samples,
    type_plot = p_types,
    region_plot = p_region,
    plot_data = plot_data
  ))
}

## 5. Track map

create_indel_tracks <- function(
    indel_data,              # wide table (VCF-like): CHROM, POS, ..., INFO, sample columns...
    gene_df,                 # gene coordinates data.frame with columns Start, End, gene_name (or gene_name col)
    gene_name,               # gene name string for title/filename
    x_range,                 # numeric vector c(xmin, xmax)
    region_colors,           # named vector of colors, names = region names (order defines region order)
    output_dir,              # where to save SVG
    sample_region = NULL,    # optional: list(region -> c(samples)), if NULL try to infer Region column after pivot
    sample_region_map = NULL,# optional: named character vector sample -> region (alternative to sample_region)
    indel_colors = c("DEL"="#e41a1c", "INS"="#377eb8", "INV" = "#4daf4a"),
    point_alpha = 0.4
) {
  
  message("Preparing SV track map")
  
  # --- 0. sanity: make sure required cols exist ---
  if (!all(c("CHROM","POS") %in% colnames(indel_data))) {
    stop("indel_data must contain at least columns: CHROM and POS")
  }
  
  # --- 1. Ensure SVTYPE and SVLEN exist (caller may have computed them) ---
  if (!"SVTYPE" %in% colnames(indel_data)) {
    indel_data$SVTYPE <- sub(".*SVTYPE=([^;]+);.*", "\\1", indel_data$INFO)
  }
  if (!"SVLEN" %in% colnames(indel_data)) {
    indel_data$SVLEN <- as.numeric(sub(".*SVLEN=([-0-9]+).*", "\\1", indel_data$INFO))
    indel_data$SVLEN[is.na(indel_data$SVLEN)] <- 0
  }
  
  # --- 2. Identify sample columns (everything except fixed VCF fields) ---
  fixed_cols <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SVTYPE","SVLEN")
  sample_cols <- setdiff(colnames(indel_data), fixed_cols)
  if (length(sample_cols) == 0) stop("No sample columns found in indel_data")
  
  # --- 3. Pivot to long format (vcf_long) ---
  vcf_long <- indel_data %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "Sample",
      values_to = "Genotype"
    ) %>%
    mutate(
      Genotype = as.character(Genotype),
      Has_Variant = !Genotype %in% c("0|0","0/0",".|.", "./.", ".|0","0|.","./0","0/.",".", NA),
      Variant_Info = case_when(
        Genotype %in% c("0|0","0/0",".|.", "./.", ".|0","0|.","./0","0/.",".") ~ "Wild_Type",
        Genotype %in% c("0|1","1|0","0/1","1/0","1|.","./1",".|1","1/.") ~ "Heterozygous",
        Genotype %in% c("1|1","1/1") ~ "Homozygous",
        TRUE ~ "Unknown"
      )
    ) %>%
    filter(Has_Variant)
  
  # attach SVTYPE/SVLEN columns to vcf_long if missing (they come from indel_data)
  if (!"SVTYPE" %in% colnames(vcf_long)) vcf_long$SVTYPE <- indel_data$SVTYPE[match(vcf_long$POS, indel_data$POS)]
  if (!"SVLEN" %in% colnames(vcf_long)) vcf_long$SVLEN <- indel_data$SVLEN[match(vcf_long$POS, indel_data$POS)]
  
  # --- 4. Map samples -> region ---
  if (!is.null(sample_region_map)) {
    # sample_region_map: named vector sample->region
    vcf_long$Region <- sample_region_map[vcf_long$Sample]
  } else if (!is.null(sample_region)) {
    # sample_region: list(region -> c(samples))
    sample_region_map_df <- data.frame(
      Sample = unlist(sample_region),
      Region = rep(names(sample_region), times = lengths(sample_region)),
      stringsAsFactors = FALSE
    )
    vcf_long <- vcf_long %>% left_join(sample_region_map_df, by = "Sample")
  } else if ("Region" %in% colnames(vcf_long)) {
    # already present, do nothing
    vcf_long$Region <- as.character(vcf_long$Region)
  } else {
    # no region info -> set NA and warn
    vcf_long$Region <- NA_character_
    warning("No sample_region or sample_region_map provided and vcf_long has no Region column; y-axis coloring by region will be skipped.")
  }
  
  # Keep only records within x_range (optional but helpful)
  if (!is.null(x_range) && length(x_range) == 2) {
    vcf_long <- vcf_long %>% filter(POS >= x_range[1] & POS <= x_range[2])
  }
  
  # --- 5. Build sample_regions table and ordering by region_colors order ---
  region_order <- names(region_colors)
  sample_regions <- vcf_long %>%
    select(Sample, Region) %>%
    distinct() %>%
    arrange(factor(Region, levels = region_order))
  
  # make Sample factor in that order
  vcf_long$Sample <- factor(vcf_long$Sample, levels = sample_regions$Sample)
  # Region factor in region_order
  vcf_long$Region <- factor(vcf_long$Region, levels = region_order)
  
  total_samples <- length(levels(vcf_long$Sample))
  
  # --- 6. Gene track boundaries ---
  gene_track_ymin <- -0.5
  gene_track_ymax <- 0
  
  # --- 7. Plot ---
  p <- ggplot() +
    geom_rect(
      data = gene_df,
      aes(xmin = Start, xmax = End, ymin = gene_track_ymin, ymax = gene_track_ymax),
      fill = "#9467bd", alpha = 0.8
    ) +
    geom_text(
      data = gene_df,
      aes(x = (Start + End)/2, y = gene_track_ymin - 0.05, label = gene_name),
      angle = 45, hjust = 1, vjust = 1, size = 1.2
    ) +
    geom_segment(
      data = vcf_long %>%
        mutate(
          x_start = POS,
          x_end = ifelse(SVTYPE == "INS", POS, POS + abs(SVLEN)),
          xmin = pmin(x_start, x_end),
          xmax = pmax(x_start, x_end)
        ),
      aes(x = xmin, xend = xmax, y = Sample, yend = Sample, color = SVTYPE),
      alpha = point_alpha,
      linewidth = 1,
      lineend = "round"
    ) +
    geom_point(
      data = sample_regions,
      aes(x = 1, y = Sample, shape = Region, fill = Region),
      size = 3, alpha = 0
    ) +
    scale_color_manual(values = indel_colors, name = "SV Type") +
    scale_fill_manual(
      values = region_colors,
      name = "Region",
      guide = guide_legend(override.aes = list(alpha = 1, size = 3))
    ) +
    scale_shape_manual(values = 21:25, name = "Region") +
    scale_x_continuous(labels = scales::comma_format(), limits = x_range) +
    scale_y_discrete(
      limits = levels(vcf_long$Sample),
      expand = expansion(add = c(5,0.3))
    ) +
    labs(
      x = "Locus (bp)",
      y = "Samples",
      title = paste0(gene_name, " — INDEL track map")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 6, hjust = 0.5),
      axis.text.y = element_markdown(
        size = 4,
        color = {
          # if Region exists and colours available, map; else default black
          col_map <- rep("black", length(levels(vcf_long$Sample)))
          if (!all(is.na(sample_regions$Region))) {
            map_idx <- match(levels(vcf_long$Sample), sample_regions$Sample)
            col_map <- region_colors[ as.character(sample_regions$Region[map_idx]) ]
            col_map[is.na(col_map)] <- "black"
          }
          col_map
        }
      ),
      axis.text.x = element_text(size = 4),
      axis.title.x = element_text(size = 4, face = "bold"),
      axis.title.y = element_text(size = 4, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 4),
      legend.text = element_text(size = 4),
      legend.spacing.y = unit(0.1, "cm"),
      legend.box.spacing = unit(0.1, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.1),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.1),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )
  
  # --- 8. Size and save ---
  plot_height <- 1 + 0.05 * total_samples
  plot_width  <- 10
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  out_file <- file.path(output_dir, paste0(gene_name, "_INDEL_track.svg"))
  
  ggsave(
    filename = out_file,
    plot = p,
    width = plot_width,
    height = plot_height,
    device = "svg"
  )
  
  message("Saved INDEL track to: ", out_file)
  
  return(list(plot = p, width = plot_width, height = plot_height, sample_order = sample_regions))
}


### Run


## Data processing

# 1. Load gene and region metadata
locus <- load_gene_information(file.path(input_dir, paste0(gene_name,".xlsx")))
res <- load_region_information(file.path(input_dir, "SampleArea.xlsx"))

sample_region <- res$sampleArea
region_colors <- res$region_colors

# 2. Load and preprocess INDEL VEP results (only if required by 'steps')
if (2 %in% steps || 3 %in% steps || 4 %in% steps || 5 %in% steps) {
  message("Loading and processing INDEL VEP results for: ", gene_name)
  
  df <- analyze_vep_indel_variants(
    file.path(input_dir, "indel", paste0(gene_name, "_indel_vep_output.vcf"))
  )
  
  result <- df$result
  all_categories <- df$categories
  
  # Convert wide matrix (Position + samples) to long format
  result_long <- result %>%
    pivot_longer(
      cols      = -c(Position, SVLEN, SVTYPE),
      names_to  = "Sample",
      values_to = "Variant_Type"
    ) %>%
    filter(!is.na(Variant_Type))
  
  # Report missing positions if any
  missing_positions <- sum(is.na(result$Position))
  if (missing_positions > 0) {
    message("Warning: ", missing_positions, " INDELs have missing genomic positions.")
  }
  
  # Count variant categories across all samples
  category_counts <- as.data.frame(table(all_categories)) %>%
    arrange(desc(Freq))
  
  print(category_counts)
}

#Step 1: Variant burden heatmap
if (1 %in% steps) {
  message("Running Step 1: Variant burden heatmap")
  
  indel_data <- fread(file.path(input_dir, "indel", paste0(gene_name,"_indel_raw.vcf")), data.table = FALSE)
  
  expected_cols <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  colnames(indel_data)[1:length(expected_cols)] <- expected_cols
  indel_data <- indel_data %>% select(-ID, -QUAL, -FILTER, -FORMAT)
  
  #extract SVTYPE, SVLEN
  indel_data$SVTYPE <- sub(".*SVTYPE=([^;]+);.*", "\\1", indel_data$INFO)
  indel_data$SVLEN  <- as.numeric(sub(".*SVLEN=([-0-9]+).*", "\\1", indel_data$INFO))
  indel_data$SVLEN[is.na(indel_data$SVLEN)] <- 0
  
  p <- run_indel_heatmap(
    gene_loc      = locus,
    indel_data    = indel_data,
    sample_region = sample_region,
    region_colors = region_colors,
    output_dir    = output_dir,
    gene_name     = gene_name
  )
  message("Step 1 finished.")
}


#Step 2: Bubble plots
if (2 %in% steps) {
  
  message("Running Step 2: Bubble map plots")
  
  plots_list <- list()
  
  for (var in names(variant_palettes)) {
    p <- plot_indel_bubble(
      indel_result  = result_long,
      sample_region = sample_region,
      gene_locus    = locus,
      indel_type    = var,
      color_palette = variant_palettes[[var]],
      min_size      = min_size,
      max_size      = max_size,
      title         = paste(var, "INDEL in", gene_name),
      output_file   = file.path(
        output_dir,
        paste0(gsub(" ", "_", var), "_INDEL_in_", gene_name, ".png")
      )
    )
    
    plots_list[[var]] <- p
  }
  
  message("Step 2 finished.")
}

#Step 3: Dual-bubble plot
if (3 %in% steps) {
  
  message("Running Step 3: Dual-bubble INDEL plot")
  
  output_file_dual <- file.path(
    output_dir,
    paste0("INDEL_Dual_Bubble_",gene_name,".png")
  )
  
  p_indel <- plot_indel_bubble_dual(
    indel_result  = result_long,
    sample_region = sample_region,
    gene_locus    = locus,
    min_size      = min_size,
    max_size      = max_size,
    title         = paste(gene_name, "INDEL Dual-Bubble (PTV vs non-PTV)"),
    output_file   = output_file_dual
  )
  
  message("Step 3 finished.")
}


#Step 4&5: Summary plots
if (4 %in% steps) {
  message("Running Step 4: Summary plots (all consequences)")
  summary_plots <- generate_indel_summary_plots(
    indel_result   = result,
    result_long    = result_long,
    sample_region  = sample_region,
    region_colors  = region_colors,
    variant_colors = variant_colors,
    gene_name      = gene_name,
    output_prefix  = file.path(output_dir, paste0(gene_name, "_indel_summary"))
  )
  message("Step 4 finished.")
}

if (5 %in% steps) {
  message("Running Step 5: Summary plots without modifier consequences")
  summary_plots_focus <- generate_indel_summary_plots(
    indel_result   = result,
    result_long    = result_long,
    sample_region  = sample_region,
    region_colors  = region_colors,
    variant_colors = variant_colors,
    gene_name      = gene_name,
    variant_types  = variant_hierarchy_no_modifier,
    output_prefix  = file.path(output_dir, paste0(gene_name, "_indel_summary_focus"))
  )
  message("Step 5 finished.")
}

#Step 6: Track plot
if (6 %in% steps) {
  message("Running Step 6: Track map")
  indel_data <- fread(file.path(input_dir, "indel", paste0(gene_name,"_indel_raw.vcf")), data.table = FALSE)
  
  expected_cols <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  colnames(indel_data)[1:length(expected_cols)] <- expected_cols
  indel_data <- indel_data %>% select(-ID, -QUAL, -FILTER, -FORMAT)
  
  indel_data$SVTYPE <- sub(".*SVTYPE=([^;]+);.*", "\\1", indel_data$INFO)
  indel_data$SVLEN  <- as.numeric(sub(".*SVLEN=([-0-9]+).*", "\\1", indel_data$INFO))
  indel_data$SVLEN[is.na(indel_data$SVLEN)] <- 0
  
  indel_tracks <- create_indel_tracks(
    indel_data     = indel_data,
    gene_df        = locus,
    gene_name      = gene_name,
    x_range        = c(min(locus$Start), max(locus$End)),
    region_colors  = region_colors,
    sample_region  = sample_region,
    output_dir     = output_dir,
    indel_colors   = indel_colors,
    point_alpha    = point_alpha 
  )
  message("Step 6 finished.")
}