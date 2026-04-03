# step5-visualization.R
# ================================================
# Evolutionary Pattern Visualization Script
# ================================================

# ---- Load required libraries ----
library(ggplot2)
library(ggdendro)
library(igraph)
library(gridExtra)
library(stringr)
library(dplyr)
library(RColorBrewer)

# Avoid automatic Rplots.pdf creation
pdf(NULL)

# ---- Parse command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("请提供输出目录路径，例如: Rscript step5-visualization.R output/GENENAME")
}
output_dir <- args[1]
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Custom colors and parameters ----
gtb_palette <- c("#ff9797","#ffb887","#ffdd50","#c2e76e","#9be4ff",
                 "#34c6ff","#cbb2de","#ffc2dc","#a8f0de","#e8d9bc")
get_colors <- function(n) colorRampPalette(gtb_palette)(max(10,n))
base_color <- "#ff9797"

title_size <- 12
axis_title_size <- 8
axis_text_size <- 6
legend_title_size <- 8
legend_text_size <- 6
dendro_label_size <- 3.5

# ---- Read summary file ----
summary_files <- list.files(output_dir, pattern="_summary\\.txt$", full.names=TRUE)
if(length(summary_files) == 0){
  stop(paste0("Summary file not found: ", file.path(output_dir, "*_summary.txt")))
}
summary_file <- summary_files[1]
file_content <- readLines(summary_file, encoding="UTF-8")

# ---- Extract cluster data ----
cluster_start <- grep("=== 分析基因: .+ ===", file_content)
cluster_end <- grep("==================================================", file_content)
cluster_end <- cluster_end[cluster_end > cluster_start][1]
cluster_lines <- file_content[cluster_start:cluster_end]
cluster_pattern <- "簇 (\\d+):"
cluster_indices <- grep(cluster_pattern, cluster_lines)

cluster_data <- data.frame()
if(length(cluster_indices) > 0){
  for(i in cluster_indices){
    cluster_num <- as.numeric(str_extract(cluster_lines[i], "\\d+"))
    prefix_line <- cluster_lines[i+1]
    prefix <- str_extract(prefix_line, "(?<=公共前缀: ')[^']+")
    prefix_length <- as.numeric(str_extract(prefix_line, "(?<=长度: )\\d+"))
    species_line <- cluster_lines[i+2]
    species <- str_split(str_extract(species_line, "(?<=包含物种: ).+"), ", ")[[1]]
    size_line <- cluster_lines[i+3]
    size <- as.numeric(str_extract(size_line, "(?<=簇大小: )\\d+"))
    dist_line <- cluster_lines[i+4]
    avg_dist <- as.numeric(str_extract(dist_line, "(?<=平均加权汉明距离: )\\d+\\.\\d+"))
    std_dev <- as.numeric(str_extract(dist_line, "(?<=± )\\d+\\.\\d+"))
    
    cluster_data <- rbind(cluster_data, data.frame(
      Cluster=cluster_num, Prefix=prefix, PrefixLength=prefix_length, 
      Size=size, AvgDistance=avg_dist, StdDev=std_dev, 
      Species=paste(species, collapse=", "), stringsAsFactors=FALSE))
  }
}

# ---- Extract HGT data ----
hgt_start <- grep("水平基因转移\\(HGT\\)检测:", file_content)
hgt_end <- grep("加速进化检测", file_content)
hgt_data <- data.frame()
if(length(hgt_start) > 0 && length(hgt_end) > 0 && hgt_end > hgt_start){
  hgt_lines <- file_content[hgt_start:(hgt_end-1)]
  hgt_indices <- grep("HGT候选:", hgt_lines)
  if(length(hgt_indices) > 0){
    for(i in hgt_indices){
      candidate <- str_extract(hgt_lines[i], "(?<=HGT候选: )[^\\s]+")
      topo_line <- hgt_lines[i+1]
      topo_neighbor <- str_extract(topo_line, "(?<=拓扑近邻: )[^\\s]+")
      topo_distance <- as.numeric(str_extract(topo_line, "(?<=距离: )\\d+"))
      feature_line <- hgt_lines[i+2]
      feature_neighbor <- str_extract(feature_line, "(?<=特征近邻: )[^\\s]+")
      feature_distance <- as.numeric(str_extract(feature_line, "(?<=距离: )\\d+"))
      feature_diff <- as.numeric(str_extract(feature_line, "(?<=特征差异: )\\d+"))
      
      hgt_data <- rbind(hgt_data, data.frame(
        Candidate=candidate, TopoNeighbor=topo_neighbor, FeatureNeighbor=feature_neighbor,
        TopoDistance=topo_distance, FeatureDistance=feature_distance, FeatureDiff=feature_diff,
        stringsAsFactors=FALSE))
    }
  }
}

# ---- Extract convergent evolution data ----
conv_start <- grep("趋同进化检测:", file_content)
conv_end <- grep("分析完成!", file_content)
conv_lines <- character(0)
if(length(conv_start) > 0 && length(conv_end) > 0 && conv_end > conv_start){
  conv_lines <- file_content[conv_start:(conv_end-1)]
}
conv_pattern <- "[\\w-]+ 和 [\\w-]+: 拓扑距离=\\d+, 特征相似度=\\d+\\.\\d+"
conv_matches <- str_extract(conv_lines, conv_pattern)
conv_matches <- conv_matches[!is.na(conv_matches)]
conv_data <- data.frame()
if(length(conv_matches) > 0){
  for(match in conv_matches){
    pair <- str_extract(match, "^.+?(?=:)")
    topo_dist <- as.numeric(str_extract(match, "(?<=拓扑距离=)\\d+"))
    feature_sim <- as.numeric(str_extract(match, "(?<=特征相似度=)\\d+\\.\\d+"))
    conv_data <- rbind(conv_data, data.frame(Pair=pair, TopoDistance=topo_dist, FeatureSimilarity=feature_sim, stringsAsFactors=FALSE))
  }
}

# ---- Build species-prefix mapping ----
species_groups <- list()
if(nrow(cluster_data) > 0){
  for(i in 1:nrow(cluster_data)){
    prefix <- cluster_data$Prefix[i]
    species <- str_split(cluster_data$Species[i], ", ")[[1]]
    species_groups[[prefix]] <- species
  }
}
species_list <- unique(unlist(species_groups))
if(length(species_list)==0) species_list <- character(0)
prefix_mapping <- data.frame(Species=species_list, ShortPrefix=NA, stringsAsFactors=FALSE)
if(length(species_groups)>0){
  for(prefix in names(species_groups)){
    idx <- prefix_mapping$Species %in% species_groups[[prefix]]
    prefix_mapping$ShortPrefix[idx] <- ifelse(is.na(prefix_mapping$ShortPrefix[idx]), prefix,
                                              paste(prefix_mapping$ShortPrefix[idx], prefix, sep=","))
  }
}

# ---- Cluster scatter plot ----
if(nrow(cluster_data) > 0){
  prefix_levels <- sort(unique(cluster_data$PrefixLength))
  colors_for_prefixlen <- setNames(get_colors(length(prefix_levels))[seq_along(prefix_levels)], as.character(prefix_levels))
  p1 <- ggplot(cluster_data %>% filter(PrefixLength>1), 
               aes(x=Size, y=AvgDistance, color=as.factor(PrefixLength))) +
    geom_point(size=2, alpha=0.8) +
    geom_errorbar(aes(ymin=AvgDistance-StdDev, ymax=AvgDistance+StdDev), width=0.2, linewidth=0.3, alpha=0.5) +
    scale_color_manual(values=colors_for_prefixlen) +
    labs(title="Clustering Statistics (prefix length > 1)",
         x="Cluster Size", y="Average Weighted Hamming Distance",
         color="Shared prefix length\n(tree depth)") +
    theme_minimal() +
    theme(plot.title=element_text(size=title_size, face="bold"),
          axis.title=element_text(size=axis_title_size),
          axis.text=element_text(size=axis_text_size),
          legend.title=element_text(size=legend_title_size),
          legend.text=element_text(size=legend_text_size))
}else{
  p1 <- ggplot() + ggtitle("No cluster data") + theme_minimal()
}

# ---- Dendrogram ----
if(length(species_list) > 0){
  species_prefix_matrix <- matrix(0, nrow=length(species_list), ncol=length(names(species_groups)))
  rownames(species_prefix_matrix) <- species_list
  colnames(species_prefix_matrix) <- names(species_groups)
  for(prefix in names(species_groups)){
    species_prefix_matrix[species_groups[[prefix]], prefix] <- 1
  }
  dist_matrix <- dist(species_prefix_matrix, method="binary")
  hc <- hclust(dist_matrix)
  dendro <- as.dendrogram(hc)
  dd <- dendro_data(dendro)
  label_data <- dd$labels
  label_data$MainPrefix <- sapply(label_data$label, function(species){
    prefixes <- prefix_mapping$ShortPrefix[prefix_mapping$Species==species]
    if(length(prefixes)==0 || is.na(prefixes)) return(NA)
    prefix_list <- unlist(strsplit(prefixes, ",")) %>% str_trim()
    depth2 <- prefix_list[nchar(prefix_list)==3]
    if(length(depth2)>0) return(depth2[1])
    longer <- prefix_list[nchar(prefix_list)>1]
    if(length(longer)>0) return(longer[which.min(nchar(longer))])
    prefix_list[1]
  })
  unique_prefixes <- na.omit(unique(label_data$MainPrefix))
  prefix_colors <- setNames(get_colors(length(unique_prefixes))[seq_along(unique_prefixes)], unique_prefixes)
  
  segs <- dd$segments
  p2 <- ggplot(segs) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend), linewidth=0.5) +
    geom_text(data=label_data, aes(x=x, y=y, label=label, color=MainPrefix), hjust=0, size=dendro_label_size, nudge_y=-0.1) +
    scale_color_manual(values=prefix_colors, na.value="gray40") +
    coord_flip() + scale_y_reverse(expand=expansion(mult=c(0.2,0.1))) +
    labs(title="Species Clustering Dendrogram Based on Prefix Clusters", x="", y="Distance", color="Main Prefix") +
    theme_minimal() +
    theme(plot.title=element_text(size=title_size, face="bold"),
          axis.title=element_text(size=axis_title_size),
          axis.text=element_text(size=axis_text_size),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
}else{
  p2 <- ggplot() + ggtitle("No species/prefix data") + theme_minimal()
}

# ---- Convergent evolution plot ----
if(nrow(conv_data)>0){
  p3 <- ggplot(conv_data, aes(x=FeatureSimilarity, y=factor(TopoDistance))) +
    geom_point(alpha=0.7, color=base_color, position=position_jitter(width=0.01)) +
    labs(title="Convergent Evolution Detection", x="Feature Similarity", y="Topological Distance") +
    theme_minimal() +
    theme(plot.title=element_text(size=title_size, face="bold"),
          axis.title=element_text(size=axis_title_size),
          axis.text=element_text(size=axis_text_size))
}else{
  p3 <- ggplot() + annotate("text", x=0.5, y=0.5, label="No convergent evolution data available", size=5, hjust=0.5) + theme_void()
}

# ---- Save CSV files ----
if(nrow(cluster_data)>0) write.csv(cluster_data, file.path(output_dir,"cluster_data.csv"), row.names=FALSE)
if(nrow(hgt_data)>0) write.csv(hgt_data, file.path(output_dir,"hgt_data.csv"), row.names=FALSE)
if(nrow(conv_data)>0) write.csv(conv_data, file.path(output_dir,"convergence_data.csv"), row.names=FALSE)

# ---- Save plots ----
plot_height <- max(8, nrow(cluster_data)/2)  # dendrogram height scale
ggsave(file.path(output_dir,"cluster_analysis.png"), p1, width=10, height=8, dpi=300)
ggsave(file.path(output_dir,"dendrogram.png"), p2, width=10, height=plot_height, dpi=300, limitsize=FALSE)
ggsave(file.path(output_dir,"convergence_analysis.png"), p3, width=10, height=10, dpi=300)

# ---- HGT network ----
png(file.path(output_dir,"hgt_network.png"), width=2000, height=2000, res=300)
if(nrow(hgt_data)>0){
  all_nodes <- unique(c(hgt_data$Candidate,hgt_data$TopoNeighbor,hgt_data$FeatureNeighbor))
  node_types <- ifelse(all_nodes %in% hgt_data$Candidate,"Candidate","Neighbor")
  hgt_nodes <- data.frame(id=all_nodes, type=node_types, stringsAsFactors=FALSE)
  hgt_edges <- data.frame(
    from=c(hgt_data$Candidate,hgt_data$Candidate),
    to=c(hgt_data$TopoNeighbor,hgt_data$FeatureNeighbor),
    type=rep(c("Topological","Feature"), each=nrow(hgt_data)),
    weight=c(hgt_data$TopoDistance,hgt_data$FeatureDistance),
    stringsAsFactors=FALSE
  )
  hgt_graph <- graph_from_data_frame(hgt_edges, vertices=hgt_nodes)
  lay <- layout_nicely(hgt_graph)
  plot(hgt_graph, layout=lay, edge.arrow.size=0.5, vertex.color=base_color,
       vertex.size=10, vertex.frame.color="gray", vertex.label.color="black",
       vertex.label.cex=0.3, edge.label=E(hgt_graph)$weight,
       edge.label.cex=0.7, edge.label.color="darkred",
       main="Horizontal Gene Transfer (HGT) Detection", cex.main=1.2)
}else{
  plot.new()
  text(0.5,0.5,"No horizontal gene transfer detected", cex=1.2)
}
dev.off()

cat("All plots and CSV files saved to:", output_dir, "\n")