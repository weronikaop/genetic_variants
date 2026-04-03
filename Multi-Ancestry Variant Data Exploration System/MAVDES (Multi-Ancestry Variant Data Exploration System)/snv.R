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
                "├── snv/                          # SNV subdirectory",
                "│     ├── <GENE>_snv_vep_output.vcf",
                "│     └── <GENE>_snv_raw.vcf",
                "",
                "If the structure differs, adjust input paths accordingly.",
                sep = "\n"
              )),
  
  make_option(c("-o","--output_dir"), type="character", default="snv_results",
              help = "Output directory for generated plots (default: %default)"),
  
  make_option(c("-s","--step"), type="character", default="1,2,3,4,5",
              help = paste(
                "Pipeline steps to run (comma-separated).",
                "Available steps:",
                "  1 = Variant burden heatmap",
                "  2 = SNV bubble plot",
                "  3 = Synonymous vs nonsynonymous bubble comparison",
                "  4 = Summary plots (all consequences)",
                "  5 = Summary plots (non-modifier only)",
                "Example: --step 1,3,5",
                "(default: %default)",
                sep = "\n"
              )),
  
  # Optional parameters
  make_option("--min_size", type="integer", default=1,
              help = "Minimum bubble size for bubble plots (default: %default)"),
  
  make_option("--max_size", type="integer", default=10,
              help = "Maximum bubble size for bubble plots (default: %default)")
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

#Print summary
cat("Gene:", gene_name, "\n")
cat("Input dir:", input_dir, "\n")
cat("Output dir:", output_dir, "\n")
cat("Steps to run:", paste(steps, collapse=", "), "\n")
cat("min_size:", min_size, "\n")
cat("max_size:", max_size, "\n")


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


## Annotate

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

run_snv_heatmap <- function(
    gene_loc,           #data.frame with gene locations, columns: Chr, Start, End, gene_name, Type
    snv_data,           #data.frame with SNVs, columns: CHROM, POS, REF, ALT, INFO, sample columns ...
    sample_region,      #list of samples per region, e.g., list(EUR = c(...), AMR = c(...))
    region_colors,      #vector of colors for regions, names = region
    output_dir,         #output directory
    gene_name = NULL    #gene name (optional, used for plot title)
) {
  
  if(is.null(gene_name)) gene_name <- unique(gene_loc$Type)
  message("Running SNV heatmap for gene: ", gene_name)
  
  #1. Identify sample columns
  fixed_cols <- c("CHROM", "POS", "REF", "ALT", "INFO")
  sample_cols <- setdiff(colnames(snv_data), fixed_cols)
  snv_data <- snv_data[, c(fixed_cols, sample_cols)]
  
  #2. Convert genotype to variant count
  convert_genotype <- function(gt) {
    if (gt %in% c("0|0", ".|.")) return(0)
    else if (gt %in% c("0|1", "1|0", ".|1", "1|.")) return(1)
    else if (gt == "1|1") return(2)
    else return(0)
  }
  
  #3. Assign SNVs to genes using genomic ranges
  gene_gr <- GRanges(
    seqnames = paste0("chr", gene_loc$Chr),
    ranges = IRanges(start = gene_loc$Start, end = gene_loc$End),
    gene_id = gene_loc$gene_name
  )
  snv_gr <- GRanges(
    seqnames = snv_data$CHROM,
    ranges = IRanges(start = snv_data$POS, end = snv_data$POS)
  )
  overlaps <- findOverlaps(snv_gr, gene_gr)
  snv_data$gene_id <- NA
  snv_data$gene_id[queryHits(overlaps)] <- gene_gr$gene_id[subjectHits(overlaps)]
  
  valid_snvs <- snv_data[!is.na(snv_data$gene_id), c("gene_id", sample_cols)]
  
  #4. Count variants per sample and gene
  long_dt <- data.table::melt(as.data.table(valid_snvs),
                              id.vars = "gene_id",
                              variable.name = "Sample",
                              value.name = "GT")
  long_dt[, variant_count := sapply(GT, convert_genotype)]
  
  result_dt <- long_dt[, .(variants = sum(variant_count)), by = .(gene_id, Sample)]
  result_matrix_long <- result_dt
  result_matrix_long$variants_log <- log1p(result_matrix_long$variants)
  
  #5. Assign population/region to samples
  sample_population <- unlist(lapply(names(sample_region), function(r) {
    setNames(rep(r, length(sample_region[[r]])), sample_region[[r]])
  }))
  result_matrix_long$Population <- sample_population[result_matrix_long$Sample]
  
  #6. Set the order of samples and genes
  result_matrix_long$Sample <- factor(result_matrix_long$Sample,
                                      levels = unlist(sample_region))
  
  gene_order <- gene_loc %>% arrange(Chr, Start) %>% pull(gene_name)
  result_matrix_long$gene_id <- factor(result_matrix_long$gene_id, levels = gene_order)
  
  sample_regions <- data.frame(
    Sample = unlist(sample_region),
    Population = rep(names(sample_region), lengths(sample_region))
  )
  
  #7. Create heatmap
  p <- ggplot(result_matrix_long, aes(x = gene_id, y = Sample, fill = variants_log)) +
    geom_tile(color = NA) +
    scale_fill_viridis(option="plasma", name="log1p(Variants)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          axis.text.y = element_text(size=7),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(size=14, face="bold", hjust=0.5)) +
    labs(title = paste("SNV Heatmap for", gene_name, "(log1p variants)"),
         x = "Gene", y = "Sample")
  
  #Color sample names according to region
  sample_order <- levels(result_matrix_long$Sample)
  sample_colors <- region_colors[sample_regions$Population][match(sample_order, sample_regions$Sample)]
  
  p <- p + scale_y_discrete(labels = function(x) paste0(x)) +
    theme(axis.text.y = element_markdown(size = 7, colour = sample_colors))
  
  #Add hidden points for region legend
  p <- p + geom_point(data = sample_regions,
                      aes(x = 1, y = Sample, color = Population),
                      inherit.aes = FALSE, size = 0, alpha = 0) +
    scale_color_manual(values = region_colors, name = "Population",
                       guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  #Create y-position index for samples in the order used in heatmap
  sample_regions$idx <- match(sample_regions$Sample, sample_order)
  
  #Compute boundaries for dashed lines
  region_boundaries <- sample_regions %>%
    group_by(Population) %>%
    summarise(min_idx = min(idx),
              max_idx = max(idx),
              .groups = "drop") %>%
    arrange(desc(max_idx))
  
  #Add dashed lines
  for (i in 1:(nrow(region_boundaries)-1)) {
    boundary_y <- region_boundaries$min_idx[i] - 0.5
    p <- p + geom_hline(yintercept = boundary_y, 
                        color="black", linetype="dashed", linewidth=0.5)
  }
  
  #Add region labels on the right
  region_labels <- region_boundaries %>% mutate(y_pos = (min_idx + max_idx)/2)
  p <- p + annotate("text", x = Inf, y = region_labels$y_pos,
                    label = region_labels$Population,
                    hjust = -0.1, vjust=0.5, size=3,
                    color="gray50", fontface="bold")
  
  #8. Save plot to file
  n_genes <- length(unique(result_matrix_long$gene_id))
  n_samples <- length(unique(result_matrix_long$Sample))
  width_dynamic  <- max(6, n_genes * 0.20)
  height_dynamic <- max(4, n_samples * 0.13)
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  outfile <- file.path(output_dir, paste0("snv_heatmap_", gene_name, ".png"))
  ggsave(
    filename = outfile,
    plot = p,
    width = width_dynamic,
    height = height_dynamic,
    dpi = 300, 
    limitsize = FALSE
  )
  message("SNV heatmap saved to: ", outfile)
  
  return(p)
}

## 2. Bubble map visualization

plot_snv_bubble <- function(
    snv_result,       #data.frame with SNV results, must contain columns: Sample, Position, Variant_Type
    sample_region,    #list of samples per region, e.g., list(EUR = c(...), AMR = c(...))
    gene_locus,       #data.frame with gene coordinates, columns: gene, chromosome, start, end
    snv_type,         #vector of variant types to include, e.g., c("Missense", "Stop gained")
    color_palette,    #vector of colors for gradient fill
    min_size = 1,     #minimum bubble size
    max_size = 10,    #maximum bubble size
    title = "SNV Distribution by Region and Gene",
    output_file = NULL
) {
  #1. Prepare gene annotations
  colnames(gene_locus) <- c("gene", "chromosome", "start", "end")
  gene_locus <- gene_locus %>%
    mutate(across(c(start, end), as.numeric),
           gene = as.character(gene),
           chromosome = as.character(chromosome))
  
  #2. Preprocess variant table (filter by type, extract coordinates)
  snv_long <- snv_result %>%
    filter(Variant_Type %in% snv_type) %>%
    separate(Position, into = c("chromosome", "position_info"), sep = ":", remove = FALSE) %>%
    separate(position_info, into = c("position", "ref_alt"), sep = "_", extra = "merge") %>%
    mutate(position = as.numeric(position),
           chromosome = as.character(chromosome))
  
  #3. Match variants to genes
  matched_variants <- map_dfr(1:nrow(gene_locus), function(i) {
    gene_row <- gene_locus[i, ]
    snv_long %>%
      filter(chromosome == gene_row$chromosome,
             position >= gene_row$start,
             position <= gene_row$end) %>%
      mutate(gene = gene_row$gene)
  })
  
  #4. Map samples to regions
  sample_region_map <- data.frame(
    Sample = unlist(sample_region),
    Region = rep(names(sample_region), times = lengths(sample_region))
  )
  matched_variants <- matched_variants %>% left_join(sample_region_map, by = "Sample")
  
  #5. Summarize variants per gene and per region
  region_total_samples <- lengths(sample_region)
  
  summary_data <- matched_variants %>%
    group_by(gene, Region, Variant_Type) %>%
    summarise(
      mutated_samples = n_distinct(Sample),
      total_mutations = n(),
      .groups = "drop"
    ) %>%
    mutate(
      region_total = region_total_samples[Region],
      proportion = mutated_samples / region_total,
      mutations_per_sample = total_mutations / region_total
    ) %>%
    filter(!is.na(Region))
  
  #6. Prepare data for plotting
  plot_data <- summary_data %>% filter(proportion > 0)
  
  if (nrow(plot_data) == 0) {
    cat(
      "No data available for plotting for variant type: ",
      paste(snv_type, collapse = ", "),
      "\n"
    )
    return(NULL)
  }
  
  n_genes <- length(unique(plot_data$gene))
  plot_width <- max(5, 3 + n_genes * 0.3)
  
  #7. Generate bubble plot
  p <- ggplot(plot_data, aes(x = gene, y = Region)) +
    geom_point(
      aes(size = proportion, fill = mutations_per_sample),
      shape = 21,
      color = "black",
      alpha = 0.8
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
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = plot_width, height = 4, dpi = 300, limitsize = FALSE)
    cat("Plot saved to:", output_file, "\n")
  }
  
  return(p)
}

## 3. Bubble map visualization for groups

plot_snv_bubble_dual <- function(
    snv_result,       #data.frame with SNV results, must contain columns: Sample, Position, Variant_Type
    sample_region,    #list of samples per region, e.g., list(EUR = c(...), AMR = c(...))
    gene_locus,       #data.frame with gene coordinates, columns: gene, chromosome, start, end
    min_size = 1,     #minimum bubble size
    max_size = 10,    #maximum bubble size
    title = "SNV Mutation Distribution by Region and Gene",
    output_file = NULL
) {
  
  #1. Prepare gene information
  colnames(gene_locus) <- c("gene", "chromosome", "start", "end")
  gene_locus <- gene_locus %>%
    mutate(start = as.numeric(start),
           end = as.numeric(end), 
           gene = as.character(gene),
           chromosome = as.character(chromosome))
  
  #2. Process SNV data: split Position, assign variant category (synonymous vs nonsynonymous)
  snv_long <- snv_result %>%
    separate(Position, into = c("chromosome", "position_info"), sep=":", remove=FALSE) %>%
    separate(position_info, into = c("position", "ref_alt"), sep="_", extra="merge") %>%
    mutate(
      position = as.numeric(position),
      chromosome = as.character(chromosome)
    ) %>%
    mutate(
      Variant_Category = case_when(
        Variant_Type %in% c("Missense", "Stop gained", "Stop lost", "Start lost") ~ "nonsynonymous",
        Variant_Type %in% c("Synonymous", "Stop retained") ~ "synonymous",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Variant_Category))
  
  #3. Match SNVs to genes based on genomic coordinates
  matched_variants <- map_dfr(1:nrow(gene_locus), function(i) {
    gene_row <- gene_locus[i, ]
    snv_long %>%
      filter(chromosome==gene_row$chromosome,
             position>=gene_row$start,
             position<=gene_row$end) %>%
      mutate(gene=gene_row$gene)
  })
  
  #4. Map samples to regions
  sample_region_map <- data.frame(
    Sample = unlist(sample_region),
    Region = rep(names(sample_region), times = lengths(sample_region))
  )
  matched_variants <- matched_variants %>% 
    left_join(sample_region_map, by="Sample") %>%
    filter(!is.na(Region))
  
  #5. Summarize mutations per gene, region, and variant category
  region_total_samples <- lengths(sample_region)
  summary_data <- matched_variants %>%
    group_by(gene, Region, Variant_Category) %>%
    summarise(
      mutated_samples = n_distinct(Sample),
      total_mutations = n(),
      .groups="drop"
    ) %>%
    mutate(
      region_total = region_total_samples[Region],
      proportion = mutated_samples / region_total,        #proportion of samples with mutation
      mutations_per_sample = total_mutations / region_total #normalized mutations per sample
    )
  
  #6. Fill missing gene-region-category combinations with zeros
  all_combinations <- expand.grid(
    gene = unique(summary_data$gene),
    Region = unique(summary_data$Region),
    Variant_Category = c("nonsynonymous","synonymous")
  )
  summary_data <- all_combinations %>%
    left_join(summary_data, by=c("gene","Region","Variant_Category")) %>%
    mutate(
      mutated_samples = ifelse(is.na(mutated_samples),0,mutated_samples),
      total_mutations = ifelse(is.na(total_mutations),0,total_mutations),
      proportion = ifelse(is.na(proportion),0,proportion),
      mutations_per_sample = ifelse(is.na(mutations_per_sample),0,mutations_per_sample),
      region_total = ifelse(is.na(region_total),1,region_total)
    )
  
  #7. Assign y positions for dual track (nonsynonymous/synonymous)
  regions <- unique(summary_data$Region)
  
  if (length(regions) == 0) {
    message("No regions available for plotting — check sample_region mapping.")
    return(NULL)
  }
  
  y_positions <- data.frame(
    Region = rep(regions, each=2),
    Variant_Category = rep(c("nonsynonymous","synonymous"), times=length(regions)),
    y_pos = 1:(length(regions)*2)
  )
  summary_data <- summary_data %>% left_join(y_positions, by=c("Region","Variant_Category"))
  
  plot_data <- summary_data %>% filter(proportion > 0)
  
  if (nrow(plot_data) > 0) {
    
    #8. Create dual-track bubble plot
    p <- ggplot(plot_data, aes(x = gene, y = y_pos)) +
      geom_point(aes(size = proportion, fill = mutations_per_sample),
                 shape = 21, color = "black", alpha = 0.8) +
      scale_fill_viridis_c(
        name = "Mutations per sample",
        option = "plasma"
      ) +
      scale_size_continuous(
        name = "Proportion of samples",
        range = c(min_size, max_size),
        labels = scales::percent_format(accuracy = 1)
      ) +
      scale_y_continuous(
        breaks = y_positions$y_pos,
        labels = paste0(
          y_positions$Region, " - ",
          ifelse(y_positions$Variant_Category=="nonsynonymous","Non-syn","Syn")
        ),
        expand = expansion(mult = 0.1)
      ) +
      scale_x_discrete(expand = expansion(mult = 0.01)) +
      labs(y=NULL, x="Gene", title = title) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle=45, hjust=1, size=10),
        axis.text.y = element_text(size=10),
        axis.title = element_text(face="bold"),
        plot.title = element_text(face="bold", hjust=0.5),
        panel.grid.major = element_line(color="gray90", linewidth=0.2),
        panel.grid.minor = element_blank()
      )
    
    #9. Adjust plot width dynamically
    n_genes <- length(unique(plot_data$gene))
    plot_width <- max(6, 3 + n_genes*0.3)
    
    #10. Save plot to file if specified
    if(!is.null(output_file)) {
      ggsave(output_file, p, width=plot_width, height=6, dpi=300, limitsize = FALSE)
      cat("Plot saved to:", output_file, "\n")
    }
    
  } else {
    #No data to plot
    cat("Plot not created: no data to display\n")
    p <- NULL
  }
  
  return(p)
}

## 4. Summary plots

generate_snv_summary_plots <- function(
    snv_result,        #data.frame with SNV summary statistics (not directly used but kept for interface consistency)
    result_long,       #long-format SNV dataset; must contain Sample, Variant_Type, Position
    sample_region,     #named list of samples per region, e.g. list(EUR = c(...), AFR = c(...))
    region_colors,     #named vector of colors for regions
    variant_colors,    #named vector of colors for variant types
    gene_name,         #gene name for plot titles
    output_prefix = NULL,  #output prefix for saved files (optional)
    variant_types = NULL   #optional filter of variant types to include
) {
  
  message("Generating SNV summary plots for gene: ", gene_name)
  
  #1. Prepare input data
  #   - Extract chromosome/position from Position field
  #   - Filter by selected variant types (if provided)
  #   - Assign region labels to samples
  
  plot_data <- result_long %>%
    filter(!is.na(Variant_Type)) %>%
    mutate(
      chromosome = sub(":.*", "", Position),
      position   = as.numeric(sub(".*:(\\d+)_.*", "\\1", Position))
    )
  
  #Filter variant types if specified
  if (!is.null(variant_types)) {
    plot_data <- plot_data %>% filter(Variant_Type %in% variant_types)
  }
  
  #Map samples to regions
  sample_to_region <- stack(sample_region)
  colnames(sample_to_region) <- c("Sample", "Region")
  
  plot_data <- plot_data %>%
    left_join(sample_to_region, by = "Sample")
  
  #Order samples according to sample_region list (top-down)
  ordered_samples <- unlist(sample_region)
  plot_data$Sample <- factor(plot_data$Sample, levels = rev(ordered_samples))
  plot_data$Sample <- droplevels(plot_data$Sample)
  
  #Restrict color palette to used types
  type_colors <- variant_colors[names(variant_colors) %in% unique(plot_data$Variant_Type)]
  
  #Pretty labels for genomic positions
  plot_data$position_label <- paste0(
    plot_data$chromosome, ":",
    format(plot_data$position, big.mark = ",")
  )
  
  #Order genomic positions in natural order
  plot_data <- plot_data %>% arrange(chromosome, position)
  plot_data$Position <- factor(plot_data$Position, levels = unique(plot_data$Position))
  
  #2. SNV track plot
  #   Heatmap-like tile plot showing samples (y-axis) vs SNV positions (x-axis)
  
  sample_regions <- plot_data %>% distinct(Sample, Region)
  
  #Color sample labels by region
  y_label_colors <- region_colors[
    as.character(
      sample_regions$Region[
        match(levels(plot_data$Sample), sample_regions$Sample)
      ]
    )
  ]
  
  p_main <- ggplot(plot_data, aes(x = Position, y = Sample, fill = Variant_Type)) +
    geom_tile(color = NA) +
    scale_fill_manual(values = type_colors, name = "Variant Type") +
    labs(
      title   = paste(gene_name, "SNV Variant Track"),
      x       = "Genomic Position",
      y       = "Sample",
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
  
  #Add invisible points to create region legend
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
  
  #Add dashed separators between regions
  region_boundaries <- plot_data %>%
    group_by(Region) %>%
    summarise(
      min_idx = min(as.numeric(Sample)),
      max_idx = max(as.numeric(Sample)),
      .groups = "drop"
    ) %>%
    arrange(desc(max_idx))
  
  for (i in 1:(nrow(region_boundaries) - 1)) {
    boundary_y <- region_boundaries$min_idx[i] - 0.5
    p_main <- p_main +
      geom_hline(yintercept = boundary_y,
                 color = "black", linetype = "dashed", linewidth = 0.5)
  }
  
  #Control density of genomic position labels
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
  
  #3. Per-sample summary plot: counts of variant types per sample
  
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
      title = paste(gene_name, "— Variant Count by Sample"),
      x = "Sample",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1,
                                     color = x_label_colors)
    ) +
    #Invisible points for region legend
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
  
  #4. Variant-type distribution plot
  
  p_types <- plot_data %>%
    count(Variant_Type) %>%
    ggplot(aes(x = reorder(Variant_Type, n), y = n, fill = Variant_Type)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.4, size = 3) +
    scale_fill_manual(values = type_colors) +
    labs(
      title = paste(gene_name, "— Variant Type Distribution"),
      x = "Variant Type",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  #5. Region-level summary plot: mean variants per sample per region
  
  region_plot_data <- plot_data %>%
    group_by(Region, Sample, Variant_Type) %>%
    summarise(n_mut = n(), .groups = "drop") %>%
    group_by(Region, Variant_Type) %>%
    summarise(mean_mut = mean(n_mut), .groups = "drop")
  
  region_plot <- ggplot(region_plot_data,
                        aes(x = Region, y = mean_mut, fill = Variant_Type)) +
    geom_col() +
    scale_fill_manual(values = type_colors) +
    labs(
      title = paste(gene_name, "— Average Variants per Sample by Region"),
      x = "Region",
      y = "Mean mutations per sample"
    ) +
    theme_minimal()
  
  #6. Save plots (adaptive size)
  
  n_samples   <- length(unique(plot_data$Sample))
  n_positions <- length(unique(plot_data$Position))
  n_types     <- length(unique(plot_data$Variant_Type))
  
  w_main <- max(10, min(30, n_positions * 0.05))
  h_main <- max(6,  min(25, n_samples * 0.18))
  
  h_samples <- max(5, min(20, n_samples * 0.14))
  w_samples <- max(10, min(20, n_samples * 0.2))
  
  w_types <- max(6, min(16, n_types * 0.5))
  h_types <- 6
  
  w_regions <- w_types
  h_regions <- 6
  
  if (!is.null(output_prefix)) {
    ggsave(paste0(output_prefix, "_main.png"),    p_main,       width = w_main,    height = h_main,    dpi = 300, 
           limitsize = FALSE)
    ggsave(paste0(output_prefix, "_samples.png"), p_samples,    width = w_samples, height = h_samples, dpi = 300, 
           limitsize = FALSE)
    ggsave(paste0(output_prefix, "_types.png"),   p_types,      width = w_types,   height = h_types,   dpi = 300, 
           limitsize = FALSE)
    ggsave(paste0(output_prefix, "_regions.png"), region_plot,  width = w_regions, height = h_regions, dpi = 300, 
           limitsize = FALSE)
  }
  
  return(list(
    main_plot    = p_main,
    sample_plot  = p_samples,
    type_plot    = p_types,
    region_plot  = region_plot,
    plot_data    = plot_data
  ))
}


### Run


## Data processing

# 1. Load gene and region metadata
locus <- load_gene_information(
  file.path(input_dir, paste0(gene_name, ".xlsx"))
)

res <- load_region_information(
  file.path(input_dir, "SampleArea.xlsx")
)

sample_region <- res$sampleArea
region_colors <- res$region_colors

# 2. Load and preprocess SNV VEP results (only if required by 'steps')
if (any(c(2, 3, 4, 5) %in% steps)) {
  
  message("Loading and processing SNV VEP results for: ", gene_name)
  
  df <- analyze_vep_snv_variants(
    file.path(input_dir, "snv", paste0(gene_name, "_snv_vep_output.vcf"))
  )
  
  result <- df$result
  all_categories <- df$categories
  
  # Convert wide matrix (Position + samples) to long format
  result_long <- result %>%
    pivot_longer(
      cols      = -Position,
      names_to  = "Sample",
      values_to = "Variant_Type"
    ) %>%
    filter(!is.na(Variant_Type))
  
  # Report missing positions if any
  missing_positions <- sum(is.na(result$Position))
  if (missing_positions > 0) {
    message("Warning: ", missing_positions, " SNVs have missing genomic positions.")
  }
  
  # Count variant categories across all samples
  category_counts <- as.data.frame(table(all_categories)) %>%
    arrange(desc(Freq))
  
  print(category_counts)
}

#Step 1: Variant burden heatmap
if (1 %in% steps) {
  message("Running Step 1: Variant burden heatmap")
  snv_data <- fread(file.path(input_dir, "snv", paste0(gene_name,"_snv_raw.vcf")), data.table = FALSE)
  
  expected_cols <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  colnames(snv_data)[1:length(expected_cols)] <- expected_cols
  snv_data <- snv_data %>% select(-ID, -QUAL, -FILTER, -FORMAT)
  
  p <- run_snv_heatmap(
    gene_loc      = locus,
    snv_data      = snv_data,
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
    p <- plot_snv_bubble(
      snv_result    = result_long,
      sample_region = sample_region,
      gene_locus    = locus,
      snv_type      = var,
      color_palette = variant_palettes[[var]],
      min_size      = min_size,
      max_size      = max_size,
      title         = paste(var, "SNV in", gene_name),
      output_file   = file.path(
        output_dir,
        paste0(gsub(" ", "_", var), "_SNV_in_", gene_name, ".png")
      )
    )
    
    plots_list[[var]] <- p
  }
  
  message("Step 2 finished.")
}

#Step 3: Dual-bubble plot
if (3 %in% steps) {
  
  message("Running Step 3: Dual-bubble SNV plot")
  
  output_file_dual <- file.path(
    output_dir,
    paste0("SNV_Dual_Bubble_",gene_name,".png")
  )
  
  p_dual <- plot_snv_bubble_dual(
    snv_result    = result_long,
    sample_region = sample_region,
    gene_locus    = locus,
    min_size      = min_size,
    max_size      = max_size,
    title         = paste(gene_name, "SNV Dual-Bubble"),
    output_file   = output_file_dual
  )
  
  message("Step 3 finished.")
}


#Step 4&5: Summary plots
if (4 %in% steps) {
  message("Running Step 4: Summary plots (all consequences)")
  summary_plots <- generate_snv_summary_plots(
    snv_result     = result,
    result_long    = result_long,
    sample_region  = sample_region,
    region_colors  = region_colors,
    variant_colors = variant_colors,
    gene_name      = gene_name,
    output_prefix  = file.path(output_dir, paste0(gene_name, "_summary"))
  )
  message("Step 4 finished.")
}

if (5 %in% steps) {
  message("Running Step 5: Summary plots without modifier consequences")
  summary_plots_focus <- generate_snv_summary_plots(
    snv_result     = result,
    result_long    = result_long,
    sample_region  = sample_region,
    region_colors  = region_colors,
    variant_colors = variant_colors,
    gene_name      = gene_name,
    variant_types  = variant_hierarchy_no_modifier,
    output_prefix  = file.path(output_dir, paste0(gene_name, "_summary_focus"))
  )
  message("Step 5 finished.")
}