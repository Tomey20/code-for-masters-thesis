############ Run in terminal first
# sudo apt update
# sudo apt install libgsl27 libgsl-dev
# sudo apt-get install libhdf5-dev libsz2
# sudo apt-get install libcurl4-openssl-dev libxml2-dev libssl-dev libbz2-dev liblzma-dev
# sudo apt-get install libjpeg-dev

# load packages
.libPaths("/work/Home/Data_for_masters_project/R packages")
#myPaths <- .libPaths()

#install.packages("devtools")
# install.packages('foreach')
# install.packages('doParallel')
# install.packages('glmnet')
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# install.packages("Seurat")# install.packages("Seurat", version = "4")
# install.packages("Signac")
# install.packages("dplyr")
# install.packages("patchwork")
# install.packages("readr")
# BiocManager::install("glmGamPoi")
# a
# #install.packages("bigstatsr")
# BiocManager::install("DirichletMultinomial")
# BiocManager::install("monaLisa")
# BiocManager::install("motifmatchr")
# BiocManager::install("edgeR")
# #BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# install.packages("hdf5r")
# BiocManager::install("AnnotationHub")
# BiocManager::install("ensembldb")
# BiocManager::install("biovizBase")
# BiocManager::install("DESeq2")
# BiocManager::install("SummarizedExperiment", force = TRUE)
# install.packages("XML")
# install.packages("R.utils")

library(GenomicRanges)
library(Biostrings)
library(jpeg)
library(ShortRead)
library(EDASeq)
library(Biobase)
library(R.utils)
library(XML)
library(SparseArray)
library(biovizBase)
library(AnnotationHub)
library(ensembldb)
library(dplyr)
library(Seurat)
library(patchwork)
library(glmGamPoi)
library(Signac)
library(readr)
library(DirichletMultinomial)
#library(hdf5r)
# BiocManager::install("S4Vectors", force = TRUE)
# BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38", force = TRUE)
# BiocManager::install("BSgenome", force = TRUE)
# BiocManager::install("motifmatchr") # this one had a non zero exit status
# install.packages("motifmatchr")
# BiocManager::install("edgeR")
# install.packages("edgeR")
# install.packages("doParallel")
# BiocManager::install("glmnet", force = TRUE)
# BiocManager::install("Matrix")
# install.packages("Matrix")
# BiocManager::install("GenomicRanges", force = TRUE)
# BiocManager::install("BiocParallel", force = TRUE)
# install.packages("foreach")
# install.packages("parallel")
#load packages
#library(BiocManager)
library(S4Vectors)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)#NCBI.HG38)#GRCh38)
library(monaLisa)
library(edgeR)
library(BSgenome)
library(motifmatchr)
library(glmnet)
library(Matrix)
library(BiocParallel)
library(foreach)
library(parallel)
library(doParallel)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(DESeq2)
library(devtools)
###########
# Day 4 and 7 analysis.
###########
# load data
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
# subset data into day 4 and day 7
metadata4days <- SHARE_seq_subsample_data@meta.data[SHARE_seq_subsample_data@meta.data$timept == "day 4",]
metadata7days <- SHARE_seq_subsample_data@meta.data[SHARE_seq_subsample_data@meta.data$timept == "day 7",]
# UMAP of seurat clusters seperated by timepoint
UMAP_days <- DimPlot(SHARE_seq_subsample_data, reduction = "wnn.umap", split.by = "timept", group.by = "seurat_clusters") + ggtitle("WNN UMAP of Seurat Clusters Separated by Timepoint") +
  labs(color = "Clusters") +
  theme_bw(base_size = 16) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold", colour = "black"),   
    legend.position = "right",  legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"))

library(ggplot2)
library(dplyr)

# Add a timepoint column to each data frame
metadata4days$Timepoint <- "Day 4"
metadata7days$Timepoint <- "Day 7"

# Combine all datasets
combined_metadata <- bind_rows(metadata4days, metadata7days)#, metadata)
combined_metadata$Timepoint <- factor(combined_metadata$Timepoint, levels = c("Day 7", "Day 4"))# "Day 4 and 7"

# Create the combined histogram
cells_per_cluster <- ggplot() +
  # Add stacked bars for Day 4 and Day 7
  geom_bar(data = combined_metadata %>% filter(Timepoint %in% c("Day 7", "Day 4")),
           aes(x = seurat_clusters, fill = Timepoint),
           stat = "count", color = "black", position = "stack", width = 0.8) +
  labs(title = "Number of Cells per Seurat Cluster by Timepoint",
       x = "Cluster",
       y = "Number of Cells") +
  theme_minimal() +
  theme_bw(base_size = 16) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold title
    axis.title = element_text(face = "bold"),  # Bold axis titles (X and Y)
    axis.text = element_text(face = "bold", colour = "black"),    # Bold axis tick labels
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  ) +
  scale_fill_manual(values = c("#1F77B4", "#E15759")) # c("skyblue", "lightcoral")


# subset data into only cells that have been overexpressed (OE)
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")

# Define the values to be removed
to_remove <- c("TFORF3549-GFP", "TFORF3550-mCherry")
# Subset the data to exclude rows with these TF values
metadata4days <- metadata4days[!metadata4days$TF %in% to_remove, ]
metadata7days <- metadata7days[!metadata7days$TF %in% to_remove, ]
day_4_rownames <- rownames(metadata4days)
OE_PB_day_4 <- OE[, day_4_rownames]
day_7_rownames <- rownames(metadata7days)
OE_PB_day_7 <- OE[, day_7_rownames]

# cluster them together
DefaultAssay(OE_PB_day_4) <- "RNA"
DefaultAssay(OE_PB_day_7) <- "RNA"
OE_PB_day_4 = AggregateExpression(OE_PB_day_4, group.by = c("seurat_clusters", "batch"),
                                  return.seurat = TRUE)
OE_PB_day_7 = AggregateExpression(OE_PB_day_7, group.by = c("seurat_clusters", "batch"),
                                  return.seurat = TRUE)
# update column names
colnames_updated <- make.unique(gsub("_.*", "", colnames(OE_PB_day_4)))
colnames(OE_PB_day_4) <- colnames_updated
# check Idents
Idents(OE_PB_day_4)
# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(OE_PB_day_4)) # Removes ".x" suffix
Idents(OE_PB_day_4) <- group_labels
DefaultAssay(OE_PB_day_4) <- "RNA"
# find differentially expressed genes for each cluster
day4_markers <- FindAllMarkers(OE_PB_day_4, only.pos = TRUE, test.use = "wilcox", assay = "RNA")

# do the same for day 7
colnames_updated <- make.unique(gsub("_.*", "", colnames(OE_PB_day_7)))
colnames(OE_PB_day_7) <- colnames_updated
Idents(OE_PB_day_7)
# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(OE_PB_day_7)) # Removes ".x" suffix
Idents(OE_PB_day_7) <- group_labels
DefaultAssay(OE_PB_day_7) <- "RNA"
# find differentially expressed genes for each cluster
day7_markers <- FindAllMarkers(OE_PB_day_7, only.pos = TRUE, test.use = "wilcox", assay = "RNA")


# Load necessary libraries
library(ggplot2)

# Get list of unique clusters
clusters <- unique(day4_markers$cluster)

# Initialize dataframe to store correlation results
correlation_results <- data.frame(Cluster = character(), Correlation = numeric(), stringsAsFactors = FALSE)

# Loop through each cluster
for (cluster in clusters) {
  # Extract DEGs for current cluster
  day4_genes <- day4_markers$gene[day4_markers$cluster == cluster]
  day7_genes <- day7_markers$gene[day7_markers$cluster == cluster]
  
  # Find common genes
  common_genes <- intersect(day4_genes, day7_genes)
  
  # Skip if no common genes
  if (length(common_genes) == 0) next
  
  # Extract normalized expression data
  day4_data <- OE_PB_day_4@assays$RNA$data[common_genes, cluster]
  day7_data <- OE_PB_day_7@assays$RNA$data[common_genes, cluster]
  
  # Compute correlation
  correlation_value <- cor(day4_data, day7_data, method = "pearson")
  
  # Store results
  correlation_results <- rbind(correlation_results, data.frame(Cluster = cluster, Correlation = correlation_value))
}

# Convert cluster names to factors for proper ordering
correlation_results$Cluster <- factor(correlation_results$Cluster, levels = clusters)
# Remove the "g" prefix and convert clusters to numeric labels
correlation_results$Cluster <- as.numeric(gsub("g", "", correlation_results$Cluster))

# Create bar plot
correlation_bars <- ggplot(correlation_results, aes(x = as.factor(Cluster), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low =  "#FFA500", high = "#E15759") + 
  labs(title = "Pearson Correlation Between Day 4 and Day 7 DEGs\nPer Seurat Cluster",
       x = "Cluster", y = "Correlation") +
  theme_minimal() +
  theme_bw(base_size = 16) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black")
  )


#### give one example scatterplot
day4_g1_genes <- day4_markers$gene[day4_markers$cluster == "g1"]

# Extract genes correctly from OE_markers
day7_g1_genes <- day7_markers$gene[day7_markers$cluster == "g1"]

# find the genes that they share in common
common_genes <- intersect(day4_g1_genes, day7_g1_genes)

# subset data for only cluster 1
OE_PB_day_4_data <- OE_PB_day_4@assays$RNA$data[common_genes, "g1"]
OE_PB_day_7_data <- OE_PB_day_7@assays$RNA$data[common_genes, "g1"]
# Convert to dataframe
df <- data.frame(Day4_Expression = OE_PB_day_4_data, Day7_Expression = OE_PB_day_7_data)

# compute Pearson correlation and R² value
r <- cor(df$Day4_Expression, df$Day7_Expression)
r_squared <- r^2

example_scatterplot <- ggplot(df, aes(x = Day4_Expression, y = Day7_Expression)) +
  geom_point(color = "#1F77B4", alpha = 0.8, size = 2) +  # Remove density coloring
  geom_smooth(method = "lm", color = "#E15759", se = TRUE) +
  labs(
    title = "DEG Correlation Between Seurat Cluster 1 in Day 4 and Day 7",
    x = "Day 4 Expression (log-normalized)",
    y = "Day 7 Expression (log-normalized)"
  ) +
  annotate("text", x = max(df$Day4_Expression) * 0.8, y = max(df$Day7_Expression) * 0.9, 
           label = paste0("R² = ", round(r_squared, 3)), size = 6, color = "black", fontface = "bold") +
  theme_minimal() +
  theme_bw(base_size = 16) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text = element_text(face = "bold", colour = "black")
  )
# use patchwork to collect the figures in one plot
library(patchwork)
cells_per_cluster_wrapped <- wrap_elements(full = cells_per_cluster)
UMAP_days_wrapped <- wrap_elements(full = UMAP_days)
correlation_bars_wrapped <- wrap_elements(full = correlation_bars)
example_scatterplot_wrapped <- wrap_elements(full = example_scatterplot)

final_plot <- (UMAP_days_wrapped | cells_per_cluster_wrapped) / 
  (example_scatterplot_wrapped | correlation_bars_wrapped) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18, colour = "black")
  )

final_plot
#######
# Compare GFP and mCherry with eachother
################# 
# Load Seurat
library(Seurat)

# Read the data
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")

# Create a new column combining GFP, mCherry, and OE status
SHARE_seq_subsample_data@meta.data$control_status <- "OE cells"  # Default to OE cells

# Assign GFP control label
SHARE_seq_subsample_data@meta.data$control_status[SHARE_seq_subsample_data@meta.data$TF == "TFORF3549-GFP"] <- "GFP control cells"

# Assign mCherry control label
SHARE_seq_subsample_data@meta.data$control_status[SHARE_seq_subsample_data@meta.data$TF == "TFORF3550-mCherry"] <- "mCherry control cells"

# Convert to factor for ordered plotting
SHARE_seq_subsample_data@meta.data$control_status <- factor(
  SHARE_seq_subsample_data@meta.data$control_status, 
  levels = c("GFP control cells", "mCherry control cells", "OE cells")  # Order in plot
)

# Generate a single UMAP with all three conditions split
UMAP_of_cell_groups <-  DimPlot(SHARE_seq_subsample_data, reduction = "wnn.umap", split.by = "control_status", group.by = "seurat_clusters") +
  ggtitle("WNN UMAP of GFP, mCherry, and OE Cells in Seurat Clusters") +
  labs(color = "Clusters") +
  theme_bw(base_size = 17) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(t = 10, b = 10)),  # Add top and bottom space
    axis.title = element_text(face = "bold"),  # Bold axis labels
    axis.text = element_text(face = "bold", colour = "black"), 
    strip.text = element_text(face = "bold"),  # Bold facet labels
    plot.margin = margin(t = 10, r = 10, b = 10, l = 30),  # Increase LEFT margin to fix cropping
    legend.position = "right",  legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"))


# subset the data into Overexpressed cells, GFP cells and mCherry cells
GFP_control = subset(SHARE_seq_subsample_data, TF == "TFORF3549-GFP")
mCherry_control = subset(SHARE_seq_subsample_data, TF == "TFORF3550-mCherry")
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")
mCherry_control = SetIdent(mCherry_control, value = mCherry_control$TF)
GFP_control = SetIdent(GFP_control, value = GFP_control$TF)

# Aggregate expression into a single column
GFP_PB_all = AggregateExpression(GFP_control)$RNA
mCherry_PB_all = AggregateExpression(mCherry_control)$RNA
GFP_PB_all <- as.matrix(GFP_PB_all)
mCherry_PB_all <- as.matrix(mCherry_PB_all)
# LogNormalize the counts like seurat does it
mCherry_PB_lognorm <- log1p((mCherry_PB_all / sum(mCherry_PB_all)) * 10000)
GFP_PB_lognorm <- log1p((GFP_PB_all / sum(GFP_PB_all)) * 10000)

# plot the correlation for lognorm
library(ggplot2)

# Create the plot data frame
plot_data <- data.frame(
  GFP = as.data.frame(GFP_PB_lognorm)$all,
  mCherry = as.data.frame(mCherry_PB_lognorm)$all
)

# calculate the correlation and R²
r <- cor(plot_data$mCherry, plot_data$GFP)
r2 <- r^2

# Generate the scatter plot
Control_correlation <- ggplot(plot_data, aes(x = GFP, y = mCherry)) +
  geom_point(color = "#1F77B4", size = 2, alpha = 0.8, shape = 21, fill = "#1F77B4", stroke = 1) +
  geom_smooth(method = "lm", color = "#E15759", fill = "#F28E2B", alpha = 1, size = 1.2, linetype = "solid") +  
  labs(
    x = "GFP Expression\n(log-normalized)",
    y = "mCherry Expression\n(log-normalized)",
    title = "Comparison of GFP and mCherry Expression"
  ) +
  # R² annotation at a custom position
  annotate("text", x = max(plot_data$GFP) * 0.5, 
           y = max(plot_data$mCherry) * 0.9, 
           label = paste0("R² = ", round(r2, 3)),
           color = "black", size = 6, hjust = -0.8, fontface = "bold") +  
  theme_bw(base_size = 17) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title = element_text(face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),
  )

##### enrichment of control cells
# subset data
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
control <- subset(SHARE_seq_subsample_data, TF %in% c("TFORF3550-mCherry", "TFORF3549-GFP"))
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")

# Convert seurat_clusters to numeric
control@meta.data$seurat_clusters <- as.numeric(as.character(control@meta.data$seurat_clusters))
OE@meta.data$seurat_clusters <- as.numeric(as.character(OE@meta.data$seurat_clusters))

# Initialize an empty data frame
enrichment_df <- data.frame(
  cluster = integer(),
  control_enrichment = numeric()
)
# sort the clusters, so it becomes cluster 1, 2, 3 ects. 
clusters <- sort(unique(control@meta.data$seurat_clusters))

for (i in clusters) {
  # Count cells in the control and non-control groups
  n_control_cells_cluster_i <- sum(control@meta.data$seurat_clusters == i)
  n_non_control_cells_cluster_i <- sum(OE@meta.data$seurat_clusters == i)
  n_control_cells_all_other_clusters <- sum(control@meta.data$seurat_clusters)
  n_non_control_cells_all_other_clusters <- sum(OE@meta.data$seurat_clusters)
  
  # Compute control enrichment
  control_enrichment <- log2(
    (n_control_cells_cluster_i / n_non_control_cells_cluster_i) /
      (n_control_cells_all_other_clusters / n_non_control_cells_all_other_clusters)
    +1)
  
  # Store the result in the data frame
  enrichment_df <- rbind(enrichment_df, data.frame(cluster = i, control_enrichment = control_enrichment))
}

library(ggplot2)
library(forcats)  # For factor reordering

Control_enrichment <- ggplot(enrichment_df, aes(x = factor(cluster, levels = sort(unique(cluster))), 
                          y = control_enrichment, fill = control_enrichment)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +  
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  
  geom_text(aes(label = round(control_enrichment, 2)), vjust = -0.5, size = 5) +  
  expand_limits(y = max(enrichment_df$control_enrichment) * 1.2) +
  theme_classic() +  
  labs(
    title = "Control Enrichment in Each Seurat Cluster",
    x = "Cluster",
    y = "log2(Control Enrichment Ratio)"
  ) +
  theme_minimal() + 
  theme_bw(base_size = 17) +
  theme(
    text = element_text(face = "bold", colour = "black"), 
    axis.text = element_text(face = "bold", colour = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  )

library(patchwork)
library(ggplot2)

pC <- Control_enrichment +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(r = 20))

# Create the top row with custom relative widths
top_row <- plot_spacer() + 
  UMAP_of_cell_groups + 
  plot_spacer() +
  plot_layout(widths = c(1, 3, 1))  # Makes UMAP_of_cell_groups wider than the spacers

# Combine with the bottom row
final_plot <- top_row / (Control_correlation + pC) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 18, colour = "black"))

final_plot

########
#### Does any of the OE seurat clusters resemble the controls?
#########
# Find differentially expressed genes for each of them and find their correlation. Does the OE correlate with the control?
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
# subset into overexpressed cells and control
control <- subset(SHARE_seq_subsample_data, TF %in% c("TFORF3550-mCherry", "TFORF3549-GFP"))
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")
# aggregate the control into seurat clusters, and make duplicates based on their batch
DefaultAssay(control) <- "RNA"
control = AggregateExpression(control, group.by = c("seurat_clusters", "batch"),
                              return.seurat = TRUE, normalization.method = "LogNormalize")
colnames_updated <- make.unique(gsub("_.*", "", colnames(control)))
colnames(control) <- colnames_updated
Idents(control)
# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(control)) # Removes ".x" suffix
Idents(control) <- group_labels
DefaultAssay(control) <- "RNA"
# find differentially expressed genes for each cluster of control cells
control_markers <- FindAllMarkers(control, only.pos = TRUE, test.use = "wilcox", assay = "RNA")

# aggregate the overexpressed cells into seurat clusters, and make duplicates based on their batch
OE_PB = AggregateExpression(OE, group.by = c("seurat_clusters", "batch"),
                            normalization.method = "LogNormalize",
                            scale.factor = 10000,
                            return.seurat = TRUE)
colnames_updated <- make.unique(gsub("_.*", "", colnames(OE_PB)))
colnames(OE_PB) <- colnames_updated
# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(OE_PB)) # Removes ".x" suffix
Idents(OE_PB) <- group_labels
DefaultAssay(OE_PB) <- "RNA"
library(devtools)
#devtools::install_github('immunogenomics/presto', force = TRUE)
#install.packages("Rcpp")
library("Rcpp")
library("presto")
# find differentially expressed genes for overexpressed cells for each cluster
OE_markers <- FindAllMarkers(OE_PB, only.pos = TRUE, test.use = "wilcox", assay = "RNA")

# Initialize a data frame to store correlation results
cluster_names <- unique(OE_markers$cluster)
correlation_results <- data.frame(Cluster = character(), Correlation = numeric(), stringsAsFactors = FALSE)

# Loop through each cluster
for (cluster in cluster_names) {
  # Extract genes from this cluster in control_markers
  control_genes <- control_markers$gene[control_markers$cluster == cluster]
  
  # Extract genes from this cluster in OE_markers
  OE_genes <- OE_markers$gene[OE_markers$cluster == cluster]
  
  
  # Find intersection (common genes)
  common_genes <- intersect(control_genes, OE_genes)
  
  # Skip cluster if there are no common genes
  if (length(common_genes) == 0) next
  
  # Extract normalized expression data for common genes
  normdata_OE <- OE_PB@assays$RNA$data[common_genes, cluster]
  normdata_control <- control@assays$RNA$data[common_genes, cluster]
  
  # Compute correlation
  correlation_value <- cor(normdata_OE, normdata_control, method = "pearson")
  
  # Store results
  correlation_results <- rbind(correlation_results, data.frame(Cluster = cluster, Correlation = correlation_value))
}

# Convert cluster names to factors for plotting
correlation_results$Cluster <- factor(correlation_results$Cluster, levels = cluster_names)
correlation_results$Cluster <- as.numeric(gsub("g", "", correlation_results$Cluster))

library(ggplot2)

bar_cor <- ggplot(correlation_results, aes(x = as.factor(Cluster), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#FFA500", high = "#E15759") + 
  labs(title = "Clusterwise Pearson Correlation of DEGs between Control and OE Cells",
       x = "Cluster", y = "Correlation", fill = "Correlation") +  # Add fill label
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 16),  # Bold axis labels
    axis.text = element_text(face = "bold", size = 16, colour = "black"),  # Bold tick labels
    legend.title = element_text(face = "bold", size = 16),  # Bold legend title
    legend.text = element_text(size = 16, face = "bold"),  # legend text
    plot.margin = margin(t = 4, r = 20, b = 20, l = 20)  # change margins
  )


#### find all number of marker genes
# Initialize data frame to store common DEG counts
intersect_DEG_counts <- data.frame(Cluster = character(), Common_DEGs = numeric(), stringsAsFactors = FALSE)

# Loop through each cluster
for (cluster in cluster_names) {
  # Extract differentially expressed genes (DEGs) from control_markers for this cluster
  control_DEGs <- control_markers$gene[control_markers$cluster == cluster]
  
  # Extract DEGs from OE_markers for this cluster
  OE_DEGs <- OE_markers$gene[OE_markers$cluster == cluster]
  
  # Find intersection (common DEGs)
  common_DEGs <- intersect(control_DEGs, OE_DEGs)
  
  # Count the number of common DEGs
  common_DEG_count <- length(common_DEGs)
  
  # Store results
  intersect_DEG_counts <- rbind(intersect_DEG_counts, 
                                data.frame(Cluster = cluster, Common_DEGs = common_DEG_count))
}

# Convert cluster names to numeric factors for ordering
intersect_DEG_counts$Cluster <- as.numeric(gsub("g", "", intersect_DEG_counts$Cluster))

# Load ggplot2
library(ggplot2)

# Calculate percentages
intersect_DEG_counts$Percentage <- round(intersect_DEG_counts$Common_DEGs / sum(intersect_DEG_counts$Common_DEGs) * 100, 1)

# Create the bar plot for common DEG counts with labels
common_DEG <- ggplot(intersect_DEG_counts, aes(x = as.factor(Cluster), y = Common_DEGs, fill = Common_DEGs)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 4, fontface = "bold") +  # Use percentage
  scale_fill_gradient(low = "#1F77B4", high = "#E15759", name = "Common DEGs") +
  labs(title = "Number of Common Differentially Expressed Genes per Seurat Cluster",
       x = "Cluster", y = "Number of Common DEGs") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(t = 20, b = 10)),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 16, colour = "black"),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16, face = "bold")
  )
common_DEG <- common_DEG +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(r = 20))
################# find shared markers in cluster 1
# Load necessary libraries
library(ggplot2)

# Extract genes from cluster g1 in control_markers
control_g1_genes <- control_markers$gene[control_markers$cluster == "g1"]

# Extract genes correctly from OE_markers
OE_g1_genes <- OE_markers$gene[OE_markers$cluster == "g1"]

# Find common DEGs in cluster g1
common_genes <- intersect(control_g1_genes, OE_g1_genes)

# Extract normalized expression data
normdata_OE <- OE_PB@assays$RNA$data[common_genes, "g1"]
normdata_control <- control@assays$RNA$data[common_genes, "g1"]

# Convert to dataframe for ggplot
df <- data.frame(OE_Expression = normdata_OE, Control_Expression = normdata_control)

# Compute correlation and R² value
r <- cor(df$Control_Expression, df$OE_Expression)
r_squared <- r^2

# Scatter plot with R²
example_scatter <- ggplot(df, aes(x = OE_Expression, y = Control_Expression)) +
  geom_point(color = "#1F77B4", alpha = 0.8, size = 2) +  
  geom_smooth(method = "lm", color = "#E15759", linetype = "solid") + 
  labs(
    title = "DEG Correlation Between Cluster 1 in Control and OE cells",
    x = "OE (log-normalized)",
    y = "Control (log-normalized)"
  ) +
  annotate("text", x = max(df$OE_Expression) * 0.7, y = max(df$Control_Expression) * 0.9, 
           label = paste0("R² = ", round(r_squared, 3)), size = 5, color = "black", fontface = "bold") +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 16, colour = "black", face = "bold")
  )

### Selected human embryonic stem cell markers
library(Seurat)
library(ggplot2)
library(patchwork)
library(grid)
# human embryonic stem cell markers (pluripotency)
genes <- c("POU5F1", "SOX2", "MYC", "NANOG")
# Remove "g" from identities
Idents(OE_PB) <- gsub("^g", "", Idents(OE_PB))

# Get list of plots, one per gene
plots <- VlnPlot(OE_PB, features = genes, ncol = 2, combine = FALSE)

# Modify each plot: min, mid (max/2), max
plots <- lapply(plots, function(p) {
  y_vals <- ggplot_build(p)$data[[1]]$y
  min_val <- round(min(y_vals, na.rm = TRUE), 2)
  max_val <- round(max(y_vals, na.rm = TRUE), 2)
  mid_val <- round(max_val / 2, 2)
  
  p + scale_y_continuous(breaks = c(min_val, mid_val, max_val)) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
      strip.text   = element_text(face = "bold"),
      legend.position = "none"
    ) +
    xlab("Cluster")
})

# Combine plots
combined_plot <- wrap_plots(plots, ncol = 2)

# Add shared y-axis
y_label <- grid::textGrob("Expression Level", rot = 90,
                          gp = grid::gpar(fontsize = 14, fontface = "bold"))

# Final layout
stem_cell_markers <- wrap_elements(full = y_label) + combined_plot +
  plot_layout(widths = c(0.05, 1)) +
  plot_annotation(
    title = "Expression of Selected hESC Marker Genes Across Seurat Clusters",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )


# sources for the marker genes
# https://www.nature.com/articles/cr2008309
# https://pubmed.ncbi.nlm.nih.gov/16865673/
# https://pubmed.ncbi.nlm.nih.gov/30031967/
# https://pmc.ncbi.nlm.nih.gov/articles/PMC3629778/
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6268870/
gene_modules <- list(
  Module1 = c("POU5F1", "NANOG", "SOX2", "MYC", "ZFP42", "ZFX", "FOXD3", "HMGA2", "STAT3", "LEF1", 
              "TCF3", "SALL4", "DPPA4", "DPPA5", "DPPA2", "DPPA3", "AK3", "DUSP6", "FRAT2", "GAL", 
              "MYO10", "PLP1", "VRK1", "VSNL1", "FOXH1", "GABRB3", "GAP43", "GRPR", "PHC1", "PRDM14", 
              "PTPRZ1", "SALL1", "SALL2", "THY1", "ZIC2", "ZIC3")


)

DefaultAssay(OE_PB) <- "RNA"
OE_PB2 <- AddModuleScore(
  object = OE_PB,
  features = gene_modules,
  ctrl = 100, # default
  name = "PseudobulkModule"
)

# Generate the violin plot with a title and corrected x/y-axis labels
Module1 <- VlnPlot(OE_PB2, features = "PseudobulkModule1") +
  ggtitle("Module Score (Seurat) of 36 hESC Markers Across Seurat Clusters") +  # Add a title
  xlab("Cluster") +  # Change x-axis label from "Identity" to "Cluster"
  ylab("Relative Expression") +  # Add Y-axis label
  labs(fill = "Clusters") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),  # Rotate and bold x-axis text
    axis.title.x = element_text(face = "bold", size = 16),  # Bold x-axis label
    axis.title.y = element_text(face = "bold", size = 16),  # Bold y-axis label
    axis.text.y = element_text(size = 16, face = "bold"),  # Bold y-axis tick labels
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and bold title
    legend.position = "none",  legend.title = element_text(face = "bold"))

# calculate module score with UCELL
library(UCell)
hESC_markers <- c("POU5F1", "NANOG", "SOX2", "MYC", "ZFP42", "ZFX", "FOXD3", "HMGA2", "STAT3", "LEF1", 
                  "TCF3", "SALL4", "DPPA4", "DPPA5", "DPPA2", "DPPA3", "AK3", "DUSP6", "FRAT2", "GAL", 
                  "MYO10", "PLP1", "VRK1", "VSNL1", "FOXH1", "GABRB3", "GAP43", "GRPR", "PHC1", "PRDM14", 
                  "PTPRZ1", "SALL1", "SALL2", "THY1", "ZIC2", "ZIC3")

data.seurat <- AddModuleScore_UCell(OE_PB, features = hESC_markers,
                                    ncores = 1)

# Identify UCell score columns
UCell_columns <- grep("UCell", colnames(data.seurat@meta.data), value = TRUE)

# Filter out columns where all values are zero
nonzero_UCell <- UCell_columns[colSums(data.seurat@meta.data[, UCell_columns]) > 0]

# Compute mean UCell score only on non-zero scores
data.seurat@meta.data$hESC_UCell <- rowMeans(data.seurat@meta.data[, nonzero_UCell], na.rm = TRUE)

Module2 <- VlnPlot(data.seurat, features = "hESC_UCell") +
  ggtitle("Module Score (UCell) of 36 hESC Markers Across Seurat Clusters") +  # Add a title
  xlab("Cluster") +  # Change x-axis label from "Identity" to "Cluster"
  ylab("Relative Expression") +  # Add Y-axis label
  labs(fill = "Clusters") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Rotate x-axis text
    axis.title.x = element_text(face = "bold", size = 16),  # Bold x-axis label
    axis.title.y = element_text(face = "bold", size = 16),  # Bold y-axis label
    axis.text.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and bold title
    legend.position = "none",  legend.title = element_text(face = "bold"))

library(patchwork)

stem_cell_markers_single <- wrap_elements(full = stem_cell_markers)

final_plot <- (
  (common_DEG + example_scatter) /
    (bar_cor + stem_cell_markers_single) /
    ((Module1 + Module2) + plot_layout(ncol = 2))
) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18)
  )

final_plot

#########
# Determine the cell types
#####
# load data and subset it into only the cells that have been overexpressed
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
# subset the data into only cells that have been overexpressed
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")
# Aggregate the cells
# https://satijalab.org/seurat/articles/de_vignette
OE_PB2 = AggregateExpression(OE, group.by = c("seurat_clusters", "batch"),
                             normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             return.seurat = TRUE)
colnames_updated <- make.unique(gsub("_.*", "", colnames(OE_PB2)))
colnames(OE_PB2) <- colnames_updated
#Idents(OE_PB2)
# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(OE_PB2)) # Removes ".x" suffix
Idents(OE_PB2) <- group_labels
DefaultAssay(OE_PB2) <- "RNA"

# find differentially expressed genes between clusters
library(devtools)
#devtools::install_github('immunogenomics/presto', force = TRUE)
#install.packages("Rcpp")
library("Rcpp")
library("presto")
markers <- FindAllMarkers(OE_PB2, only.pos = TRUE, test.use = "wilcox", assay = "RNA")
# find top 20 markers
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
library(patchwork)

# Convert seurat_clusters to a character, remove the "g" prefix, then convert to numeric
OE_PB2$seurat_clusters <- as.character(OE_PB2$seurat_clusters)
OE_PB2$seurat_clusters <- sub("^g", "", OE_PB2$seurat_clusters)  # remove the 'g' prefix
OE_PB2$seurat_clusters <- as.numeric(OE_PB2$seurat_clusters)

# Make it a factor, sorted in ascending order
OE_PB2$seurat_clusters <- factor(OE_PB2$seurat_clusters,
                                 levels = sort(unique(OE_PB2$seurat_clusters)))
####################################### cluster 0
genes0 <- top20$gene[1:20]
# Use cowplot to get it into one single plot
library(cowplot)
library(remotes)
# Use remotes to install enrichR from GitHub
#remotes::install_github("wjawaid/enrichR", force = TRUE)
library(enrichR)
markers0vln <- VlnPlot(
  OE_PB2, 
  features = genes0, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers0vln_labeled <- ggdraw(markers0vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers0atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g0",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 0")

markers0GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g0",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 0")

library(patchwork)

left_panel <- wrap_plots(markers0vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers0atlas + labs(tag = "B")) /
  (markers0GO    + labs(tag = "C"))

combined_plot0 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot0)

####################### cluster 1
genes1 <- top20$gene[21:40]
markers1vln <- VlnPlot(
  OE_PB2, 
  features = genes1, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers1vln_labeled <- ggdraw(markers1vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers1atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g1",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 1")

markers1GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g1",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 1")

library(patchwork)

left_panel <- wrap_plots(markers1vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers1atlas + labs(tag = "B")) /
  (markers1GO    + labs(tag = "C"))

combined_plot1 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot1)

########################## cluster 2
genes2 <- top20$gene[41:60]
markers2vln <- VlnPlot(
  OE_PB2, 
  features = genes2, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers2vln_labeled <- ggdraw(markers2vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers2atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g2",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 2")

markers2GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g2",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 2")

library(patchwork)

left_panel <- wrap_plots(markers2vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers2atlas + labs(tag = "B")) /
  (markers2GO    + labs(tag = "C"))

combined_plot2 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot2)

###################################### cluster 3
genes3 <- top20$gene[61:80]
markers3vln <- VlnPlot(
  OE_PB2, 
  features = genes3, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers3vln_labeled <- ggdraw(markers3vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers3atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g3",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 3")

markers3GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g3",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 3")

library(patchwork)

left_panel <- wrap_plots(markers3vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers3atlas + labs(tag = "B")) /
  (markers3GO    + labs(tag = "C"))

combined_plot3 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot3)


#################### cluster 4
genes4 <- top20$gene[81:100]
markers4vln <- VlnPlot(
  OE_PB2, 
  features = genes4, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers4vln_labeled <- ggdraw(markers4vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers4atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g4",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 4")

markers4GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g4",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 4")

library(patchwork)

left_panel <- wrap_plots(markers4vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers4atlas + labs(tag = "B")) /
  (markers4GO    + labs(tag = "C"))

combined_plot4 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot4)


#################### cluster 5
# Very Important to write about cells that have been dificult to find the cell type of. Such as this cluster
# the ZG16 suggests goblet cells, while GABRB3 suggests Gabaergic neurons
# according to https://panglaodb.se/ it is goblet cells when all top 20 genes are put in there
genes5 <- top20$gene[101:120]
markers5vln <- VlnPlot(
  OE_PB2, 
  features = genes5, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers5vln_labeled <- ggdraw(markers5vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers5atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g5",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 5")

markers5GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g5",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 5")

library(patchwork)

left_panel <- wrap_plots(markers5vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers5atlas + labs(tag = "B")) /
  (markers5GO    + labs(tag = "C"))

combined_plot5 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot5)

####################### cluster 6
genes6 <- top20$gene[121:140]
markers6vln <- VlnPlot(
  OE_PB2, 
  features = genes6, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers6vln_labeled <- ggdraw(markers6vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers6atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g6",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 6")

markers6GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g6",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 6")

library(patchwork)

left_panel <- wrap_plots(markers6vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers6atlas + labs(tag = "B")) /
  (markers6GO    + labs(tag = "C"))

combined_plot6 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot6)

################ cluster 7
genes7 <- top20$gene[141:160]
markers7vln <- VlnPlot(
  OE_PB2, 
  features = genes7, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers7vln_labeled <- ggdraw(markers7vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers7atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g7",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 7")

markers7GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g7",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 7")

library(patchwork)

left_panel <- wrap_plots(markers7vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers7atlas + labs(tag = "B")) /
  (markers7GO    + labs(tag = "C"))

combined_plot7 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot7)

genes7.2 <- c("CD34", "STAB2", "ESAM", "PECAM1", "CDH5")
markers7.2vln <- VlnPlot(
  OE_PB2, 
  features = genes7.2, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme(
    axis.title.y = element_text(face = "bold"),  
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold", size = 14)
  ) &
  xlab("Cluster")
markers7.2vln

################ cluster 8
# https://panglaodb.se/
# according to that it is fibroblasts. from the gene set: DNAH5, CAMK2D, SHISA2, SLIT2, CHRD, TSPAN12, PKD1L1, NFASC, SULF2, ASCL4, COL4A6, TFF3, SHH
# But ASCL4 and PKD1L1 are not included, so could be them, IT IS FETAL LIVER
# https://www.nature.com/articles/s41392-020-00222-7
# it mentions FGF21
genes8 <- top20$gene[161:180]
markers8vln <- VlnPlot(
  OE_PB2, 
  features = genes8, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers8vln_labeled <- ggdraw(markers8vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers8atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g8",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 8")

markers8GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g8",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 8")

library(patchwork)

left_panel <- wrap_plots(markers8vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers8atlas + labs(tag = "B")) /
  (markers8GO    + labs(tag = "C"))

combined_plot8 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot8)

genes8.2 <- c("FGF21", "FOXA1", "FOXA2", "PROX1", "FGFR2", "GATA4")
markers8.2vln <- VlnPlot(
  OE_PB2, 
  features = genes8.2, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme(
    axis.title.y = element_text(face = "bold"),  
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold", size = 14)
  ) &
  xlab("Cluster")
markers8.2vln

############################ cluster 9
genes9 <- top20$gene[181:200]
markers9vln <- VlnPlot(
  OE_PB2, 
  features = genes9, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers9vln_labeled <- ggdraw(markers9vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers9atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g9",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 9")

markers9GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g9",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 9")

library(patchwork)

left_panel <- wrap_plots(markers9vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers9atlas + labs(tag = "B")) /
  (markers9GO    + labs(tag = "C"))

combined_plot9 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot9)

################################# cluster 10
genes10 <- top20$gene[201:220]
markers10vln <- VlnPlot(
  OE_PB2, 
  features = genes10, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers10vln_labeled <- ggdraw(markers10vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers10atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g10",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 10")

markers10GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g10",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 10")

library(patchwork)

left_panel <- wrap_plots(markers10vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers10atlas + labs(tag = "B")) /
  (markers10GO    + labs(tag = "C"))

combined_plot10 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot10)

####################### cluster 11
genes11 <- top20$gene[221:240]
markers11vln <- VlnPlot(
  OE_PB2, 
  features = genes11, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers11vln_labeled <- ggdraw(markers11vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers11atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g11",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 11")

markers11GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g11",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 11")

library(patchwork)

left_panel <- wrap_plots(markers11vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers11atlas + labs(tag = "B")) /
  (markers11GO    + labs(tag = "C"))

combined_plot11 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot11)

############################## cluster 12
genes12 <- top20$gene[241:260]
markers12vln <- VlnPlot(
  OE_PB2, 
  features = genes12, 
  group.by = "seurat_clusters", 
  assay = "RNA", 
  pt.size = 0, 
  combine = TRUE, 
  ncol = 4
) &
  theme_bw(base_size = 15) &
  theme(
    axis.text = element_text(face = "bold", colour = "black"),
    title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position   = "none"
  ) &
  xlab("Cluster")

markers12vln_labeled <- ggdraw(markers12vln) +               
  draw_label(
    "Expression level", x = 0, y = 0.5, angle = 90,
    vjust = 0.5, hjust = 0.5, fontface = "bold", size = 15
  ) +
  labs(tag = "A")

markers12atlas <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g12",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "Human_Gene_Atlas",
  num.pathway = 10,
  return.gene.list = FALSE  
)&
  theme_bw(base_size = 15) +
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Human Gene Atlas Positive Markers\nCluster 12")

markers12GO <- DEenrichRPlot(
  OE_PB2,
  ident.1 = "g12",
  balanced = FALSE,
  logfc.threshold = 0.1,
  assay = "RNA",
  max.genes = 20,
  test.use = "wilcox",
  p.val.cutoff = 0.01,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2023",
  num.pathway = 10,
  return.gene.list = FALSE
) &
  theme_bw(base_size = 15) &
  theme(
    title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold"),
    axis.text.y = element_blank()
  ) &
  ggtitle("Go Biological Process 2023 Positive Markers\nCluster 12")

library(patchwork)

left_panel <- wrap_plots(markers12vln_labeled) 

left_panel <- left_panel + labs(tag = "A")

right_panel <- (markers12atlas + labs(tag = "B")) /
  (markers12GO    + labs(tag = "C"))

combined_plot12 <- left_panel + right_panel + 
  plot_layout(widths = c(1.85, 1)) & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    title = element_text(face = "bold")  
  )

print(combined_plot12)

######## Marker gene expression comparison for all clusters
library(patchwork)
DefaultAssay(OE_PB2) <- "RNA"
OE_PB2 <- ScaleData(OE_PB2, features = rownames(OE_PB2))
identities <- Idents(OE_PB2)
new_identities <- gsub("^g", "", identities)
Idents(OE_PB2) <- new_identities
genes_of_interest <- c("ZFP42", "SOX2", "DPPA4", 
                       "LYVE1", "CD33", 
                       "COLCA2", "FPR1", 
                       "FOXE1", 
                       "PHOX2B", "NEUROG3", "OLIG2", 
                       "TEX14", "RMRP", "FOLH1", 
                       "ZG16", "GABRB3", "PIP4K2A",
                       "HOXC4", "HOXD4", "HOXD3",
                       "CD34", "STAB2", "ESAM", "COL5A2", "COL5A1",
                       "PKD1L1", "FOXA1", "FOXA2", "PROX1", "GATA4",
                       "NPHS1", "PODXL", "KIRREL2", "SYNPO",
                       "ACTN2", "PAX7", "SOX6",
                       "KRT19", "TP63", "CD44",
                       "PRDM16", "ZNF423", "COL6A1")


the_heatmap <- DoHeatmap(object = OE_PB2,
                         features = genes_of_interest,
                         assay = "RNA",
                         size = 5) + # 3 before
  ggtitle("Expression of Marker Genes in Seurat Clusters") +
  labs(y = "Marker Genes") +   # This adds a y-axis label
  scale_color_discrete(
    name = "Clusters",
    guide = guide_legend(
      override.aes = list(shape = 16, size = 4)
    )
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x  = element_blank(),   # remove text
    panel.border     = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),   # remove tick marks
    axis.title.x = element_blank(),   # remove title
    text = element_text(face = "bold", colour = "black"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(face = "bold", colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold")
  )



# Plot for where we have the 3 hESC samples
SHARE_seq_subsample_data = SetIdent(SHARE_seq_subsample_data, value = SHARE_seq_subsample_data$seurat_clusters)
new.cluster.ids <- c("hESC 1", "hESC 2", "hESC 3",
                     "Neuron like cells", "T/P like cells", "Goblet like cells", 
                     "OP like cells", "VE like cells",
                     "FL like cell", "Podocyte like cells", "SMP like cells", 
                     "BE like cells", "Adipocyte like cells")
names(new.cluster.ids) <- levels(SHARE_seq_subsample_data)
SHARE_seq_subsample_data <- RenameIdents(SHARE_seq_subsample_data, new.cluster.ids)

the_UMAP <- DimPlot(SHARE_seq_subsample_data, reduction = "wnn.umap", label = FALSE, pt.size = 0.5, group.by = "ident", label.color = "black") +
  ggtitle("WNN UMAP of Cell-Like Types in Seurat Clusters") +
  labs(color = "Clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title = element_text(face = "bold"),  # Bold axis titles
        axis.text = element_text(face = "bold", size = 16),   
        legend.position = "right",  legend.title = element_text(face = "bold", size = 17),
        legend.text = element_text(face = "bold", size = 16))


umap_with_labels <- LabelClusters(
  plot     = the_UMAP,
  id       = "ident",          # which metadata column to label
  repel    = TRUE,             # keep labels from overlapping
  size     = 5,                # text size
  fontface = "bold"            # makes the text bold
)

umap_with_labels

the_UMAP_adjusted <- umap_with_labels +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.05))) +  # 5% padding on both sides
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

UMAP_wrapped <- wrap_elements(full = the_UMAP_adjusted)

heatmap_wrapped <- wrap_elements(full = the_heatmap)

final_plot <- (
  heatmap_wrapped | UMAP_wrapped
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag   = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

final_plot
############
## TF enrichment analysis
############
# load data
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
metadata <- SHARE_seq_subsample_data@meta.data
# make it so that only the TF name is there
metadata$TF <- sub("^[^-]*-", "", metadata$TF)
metadata$TF <- as.factor(metadata$TF)
################ TF enrichment relative to GFP AND mCherry control
# Filter the metadata into control and non-control
metadata_control <- metadata[metadata$TF %in% c("GFP", "mCherry"), ]
metadata_no_control <- metadata[!(metadata$TF %in% c("GFP", "mCherry")), ]
# Initialize a data.frame to store results
enrichment_results <- data.frame(
  TF = character(),
  Cluster = character(),
  Enrichment = numeric()
)
# Loop through each unique TF and cluster
unique_TFs <- unique(metadata_no_control$TF)
unique_clusters <- unique(metadata$seurat_clusters)

for (tf in unique_TFs) {
  for (cluster in unique_clusters) {
    # Cells expressing the TF in the cluster
    tf_cluster_count <- sum(metadata_no_control$TF == tf & metadata_no_control$seurat_clusters == cluster)
    
    # Cells expressing the control in the cluster
    control_cluster_count <- sum(metadata_control$seurat_clusters == cluster)
    
    # Cells expressing the TF across all clusters
    tf_total_count <- sum(metadata_no_control$TF == tf)
    
    # Cells expressing the control across all clusters
    control_total_count <- nrow(metadata_control)
    
    # Proportions
    cluster_proportion_tf <- tf_cluster_count / control_cluster_count
    overall_proportion_tf <- tf_total_count / control_total_count
    
    # Log2 enrichment
    enrichment <- log2((cluster_proportion_tf / overall_proportion_tf)+1) 
    
    # Store the result
    enrichment_results <- rbind(enrichment_results, data.frame(
      TF = tf,
      Cluster = cluster,
      Enrichment = enrichment
    ))
  }
}


# Convert results_df to a matrix
enrichment_results <- reshape2::acast(enrichment_results, TF ~ Cluster, value.var = "Enrichment", fill = 0)
enrichment_results <- enrichment_results[,order(as.numeric(colnames(enrichment_results)))]
colnames(enrichment_results) <- c("hESC 1", "hESC 2", "hESC 3",
                                  "Neuron like cells", "T/P like cells", "Goblet like cells", 
                                  "OP like cells", "VE like cells",
                                  "FL like cell", "Podocyte like cells", "SMP like cells", 
                                  "BE like cells", "Adipocyte like cells")
library(ComplexHeatmap)
library(circlize)

# Plot the heatmap
enrichment_heatmap <- ComplexHeatmap::Heatmap(
  enrichment_results,
  name = "Log2 Enrichment",                       # Legend title
  col = c("white", "orange", "red"),                            # Define the color palette
  na_col = "grey",                                # Color for NA values
  cluster_rows = TRUE,                            # Cluster rows
  cluster_columns = TRUE,                         # Cluster columns
  row_names_side = "left",                        # Row labels on the left
  column_names_side = "bottom",                      # Column labels on the top
  column_dend_side = "top",
  show_row_names = FALSE,                         # Hide row names
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),          # Font size for column names
  row_names_gp = gpar(fontsize = 10),             # Font size for row names
  column_title = "TF Enrichment relative to control (log2 ratio) per Seurat Cluster",               # Title for columns
  column_title_side = "top",
  row_title = "Transcription Factors",            # Title for rows
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  row_title_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    title      = "Log2 Ratio",
    title_gp   = gpar(fontsize = 12, fontface = "bold"),     # legend title
    labels_gp  = gpar(fontsize = 11, fontface = "bold")      # legend tick labels
  )
)

################ compare percent of cells with TF ORF with TF enrichment
# Create an empty list to store the results
results <- list()

# Outer loop: Iterate over each unique cluster
for (cluster in unique(metadata$seurat_clusters)) {
  
  # Filter rows for the current cluster
  cluster_rows <- metadata[metadata$seurat_clusters == cluster, ]
  
  # Inner loop: Iterate over each unique TF
  for (tf in unique(metadata$TF)) {
    
    # Filter rows for the current TF within the cluster
    tf_rows <- cluster_rows[cluster_rows$TF == tf, ]
    
    # Calculate the percentage of cells with this TF in the current cluster
    percentage <- (nrow(tf_rows) / nrow(cluster_rows)) * 100
    
    # Store the result as a list entry
    results[[paste0("Cluster_", cluster, "_TF_", tf)]] <- list(
      Cluster = cluster,
      TF = tf,
      Percentage = percentage
    )
  }
}

# Combine the results into a data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))

# Load necessary library
library(pheatmap)

# Convert results_df to a matrix
heatmap_data <- reshape2::acast(results_df, TF ~ Cluster, value.var = "Percentage", fill = 0)
# Reorder columns based on numeric order
heatmap_data <- heatmap_data[, order(as.numeric(colnames(heatmap_data)))]

########### for the enrichment relative to control correlation comparison
# Find common row names
common_rows <- intersect(rownames(heatmap_data), rownames(enrichment_results))

# Subset both matrices to include only common rows
heatmap_data <- heatmap_data[common_rows, ]
enrichment_results <- enrichment_results[common_rows, ]

# calculate pearson correlation
# Compute correlations
cor_mat_log2_percent <- matrix(nrow = nrow(enrichment_results), ncol = 1)
for (i in 1:nrow(enrichment_results)){
  cor_mat_log2_percent[i,] <- cor(enrichment_results[i,], heatmap_data[i,], method = "pearson")
}
rownames(cor_mat_log2_percent) <- rownames(enrichment_results)

# Prepare data for bar plot with sorted values
correlation_data <- data.frame(
  Correlation = cor_mat_log2_percent
)

# Add a Transcription Factor column for sorting
correlation_data$TranscriptionFactor <- seq_along(correlation_data$Correlation)

# Sort the data by correlation values
correlation_data <- correlation_data[order(correlation_data$Correlation), ]

# Calculate the average Pearson correlation
avg_correlation <- mean(correlation_data$Correlation, na.rm = TRUE)

# Create sorted bar plot
library(ggplot2)
# Adjusted bar plot with no white lines
cor_plot <- ggplot(correlation_data, aes(x = factor(TranscriptionFactor, levels = TranscriptionFactor), y = Correlation)) +
  geom_bar(stat = "identity", aes(fill = Correlation > 0), width = 1, show.legend = FALSE) +  # Adjust bar width
  coord_flip() +
  scale_fill_manual(values = c("red", "steelblue"), guide = "none") +  # Suppress guide for `fill`
  geom_hline(aes(yintercept = avg_correlation, color = "Average Correlation"), 
             linetype = "dashed", size = 1.2, show.legend = TRUE) +  # Add horizontal average line with a legend
  scale_color_manual(
    name = "Average Correlation",
    values = c("Average Correlation" = "black"),
    labels = paste0("", round(avg_correlation, 2))  # Format dynamically  Average Correlation:
  ) +  
  theme_minimal() +
  labs(title = "Pearson Correlation between percent cells with TF and TF enrichment relative to control",
       x = "Transcription Factors (Sorted)",  # Update x-axis label
       y = "Pearson Correlation") +
  theme(
    axis.text.y = element_blank(),  # Remove individual TF labels
    axis.ticks.y = element_blank(),  # Remove ticks on the y-axis
    axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 14, face = "bold", colour = "black")
  )

library(ggplot2)

# Create a data frame for ggplot
plot_data <- data.frame(
  HeatmapValue = heatmap_data["HNF4A", ],
  EnrichmentValue = enrichment_results["HNF4A", ]
)
# compute pearson correlation and R² value
r <- cor(plot_data$EnrichmentValue, plot_data$HeatmapValue)
r2 <- r^2

# Create the plot with R²
HNF4A_plot <- ggplot(plot_data, aes(x = HeatmapValue, y = EnrichmentValue)) +
  geom_point(color = "#1F77B4", size = 4, alpha = 0.9, shape = 21, fill = "#1F77B4", stroke = 1) +  
  geom_smooth(method = "lm", color = "#E15759", alpha = 0.15, size = 1.2, se = FALSE) +  
  labs(
    title = "HNF4A Enrichment vs Percent",
    x = "Percent cells with TF ORF",
    y = "Log2 ratio enrichment"
  ) +
  annotate("text", 
           x = max(plot_data$HeatmapValue) * 0.65, 
           y = min(plot_data$EnrichmentValue) + 
             0.8*(max(plot_data$EnrichmentValue)-min(plot_data$EnrichmentValue)), 
           label = paste0("R² = ", round(r2, 3)),
           color = "black", size = 5, hjust = 0, fontface = "bold") +  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    axis.title = element_text(face = "bold", size = 16),  
    axis.text = element_text(color = "black", size = 16, face = "bold"),  
    panel.grid.major = element_line(color = "gray80", linetype = "dotted"),  
    panel.grid.minor = element_blank()
  )

######### Compare Tf enrichment with likelihood of expression
# load data
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
# make it so that only the TF name is there and not the stuff before it
metadata <- SHARE_seq_subsample_data@meta.data
metadata$TF <- sub("^[^-]*-", "", metadata$TF)
metadata$TF <- as.factor(metadata$TF)
# subset metadata into control and non control
metadata_control <- metadata[metadata$TF %in% c("GFP", "mCherry"), ]
metadata_no_control <- metadata[!(metadata$TF %in% c("GFP", "mCherry")), ]
# subset into only unique TFs to loop through
all_TFs <- unique(metadata_no_control$TF)

library(ggplot2)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# Prepare storage for results
results <- data.frame(TF = character(), OE = numeric(), Control = numeric())
# loop through each TF
for (TF in all_TFs){
  # subset for cells with the given TF
  cells_with_TF <- rownames(metadata_no_control[metadata_no_control$TF == TF, ])
  cells_with_control <- rownames(metadata_control)
  
  # Expression data
  tf_expression_in_own_cells <- SHARE_seq_subsample_data@assays[["RNA"]]@data[TF, cells_with_TF]
  tf_expression_in_control <- SHARE_seq_subsample_data@assays[["RNA"]]@data[TF, cells_with_control]
  
  # Compute normalized log2 likelihood
  num_cells_above_0_own <- sum(tf_expression_in_own_cells > 0)
  total_cells_own <- length(tf_expression_in_own_cells)
  likelihood_own <- if (total_cells_own > 0) log2((num_cells_above_0_own / total_cells_own) + 1) else 0
  
  num_cells_above_0_control <- sum(tf_expression_in_control > 0)
  total_cells_control <- length(tf_expression_in_control)
  likelihood_control <- if (total_cells_control > 0) log2((num_cells_above_0_control / total_cells_control) + 1) else 0
  
  # Store results
  results <- rbind(results, data.frame(TF = TF, OE = likelihood_own, Control = likelihood_control))
}

# Plot density distribution of likelihoods
results_long <- melt(results, id.vars = "TF")

self_expression_distribution <- ggplot(results_long, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Self-Expression Likelihoods in OE vs. Control",
       x = "Log2 Likelihood",
       y = "Density",
       fill = "Condition") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Bold and center title
    axis.title = element_text(face = "bold", size = 16),  # Bold axis labels
    axis.text = element_text(face = "bold", size = 16, colour = "black"),  # Bold axis tick labels
    legend.title = element_text(face = "bold", size = 16),  # Bold legend title
    legend.text = element_text(face = "bold", size = 16)  # Bold legend labels
  )

library(ComplexHeatmap)
library(patchwork)
library(grid)

# Convert ComplexHeatmap to grob explicitly
heatmap_grob <- grid.grabExpr(draw(enrichment_heatmap))

# Wrap the grob for patchwork compatibility
heatmap_wrapped <- wrap_elements(full = heatmap_grob)
cor_wrapped <- wrap_elements(full = cor_plot)
self_expression_distribution_wrapped <- wrap_elements(full = self_expression_distribution)
# Now combine all plots correctly:
final_combined_plot <- (
  (heatmap_wrapped | HNF4A_plot) /
    (cor_wrapped | self_expression_distribution_wrapped)
) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

final_combined_plot

##############
# Run chromVAR on all cells
#############
# load data
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
# set identities to be the seurat clusters
SHARE_seq_subsample_data = SetIdent(SHARE_seq_subsample_data, value = SHARE_seq_subsample_data$seurat_clusters)
# makethe identities the cell types
new.cluster.ids <- c("LYVE like cell", "Monocyte like cells", "Thyroid like cells",
                     "Neuron like cells", "T/P like cells", "Goblet like cells", 
                     "OP like cells", "VE like cells",
                     "FL like cell", "Podocyte like cells", "SMP like cells", 
                     "BE like cells", "Adipocyte like cells")
names(new.cluster.ids) <- levels(SHARE_seq_subsample_data)
SHARE_seq_subsample_data <- RenameIdents(SHARE_seq_subsample_data, new.cluster.ids)
# subset into only OE cells
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")

########## Run chromVAR on single cell on all peaks
DefaultAssay(OE) <- "ATAC"
# load motif library
library(chromVARmotifs)
data("human_pwms_v2")
# create the motif matrix
motif.matrix <- CreateMotifMatrix(
  features = granges(OE),
  pwm = human_pwms_v2,
  genome = "BSgenome.Hsapiens.UCSC.hg38"
)
library(BiocParallel)
register(SerialParam())
set.seed(42)
chromVAR_object <- RunChromVAR(
  OE,
  genome = "BSgenome.Hsapiens.UCSC.hg38",
  motif.matrix = motif.matrix,
  assay = "ATAC",
  new.assay.name = "chromvar"
)
#saveRDS(chromVAR_object,"/work/Home/Data_for_masters_project/chromVAR_results/chromVAR_single_cell_all_peaks.rds")
chromVAR_object <- readRDS("/work/Home/Data_for_masters_project/chromVAR_results/chromVAR_single_cell_all_peaks.rds")
############
# Explore GC content in data
###########
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)
#install.packages("jpeg")
library(jpeg)
#BiocManager::install("ShortRead")
library(ShortRead)
#BiocManager::install("EDASeq")
library(EDASeq)
# load data
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
SHARE_seq_subsample_data = SetIdent(SHARE_seq_subsample_data, value = SHARE_seq_subsample_data$seurat_clusters)
new.cluster.ids <- c("LYVE like cell", "Monocyte like cells", "Thyroid like cells",
                     "Neuron like cells", "T/P like cells", "Goblet like cells", 
                     "OP like cells", "VE like cells",
                     "FL like cell", "Podocyte like cells", "SMP like cells", 
                     "BE like cells", "Adipocyte like cells")
names(new.cluster.ids) <- levels(SHARE_seq_subsample_data)
SHARE_seq_subsample_data <- RenameIdents(SHARE_seq_subsample_data, new.cluster.ids)
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")

# Extract GRanges of peaks
gr_peaks <- granges(OE)

# Get the sequences
peak_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_peaks)

# Calculate GC fraction
gc_fraction <- letterFrequency(peak_seqs, letters="GC", as.prob=TRUE)

# Ensure gc_fraction is a numeric vector
gc_fraction <- as.numeric(gc_fraction)  # Convert from matrix to vector

# Create Data Frame
df <- data.frame(GC_Fraction = gc_fraction)

# Check if it worked by Plot
GC_fraction_histogram <- ggplot(df, aes(x = GC_Fraction)) +
  geom_histogram(bins = 50, fill = "#4682B4", color = "black", alpha = 0.75) +
  labs(title = "GC Fraction Across ATAC Peaks",
       x = "GC Fraction",
       y = "Frequency") +
  theme_minimal() +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        text = element_text(face = "bold", colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", colour = "black"),
        panel.grid.major = element_line(linetype = "dashed"))

################ GC content affects differential expression
# Perform a mock null DA analysis
# Perform a mock null DA analysis by splitting cells within the same cluster into two artificial groups.
# Run DA analysis on these groups using edgeR.
# LogFC should already be centered around 0 if no systematic bias exists.
# Plot logFC vs. GC content

# Pick one cluster
OE_subset <- subset(OE, idents = "Neuron like cells")

# Manually assign 8 batches to each group
groupA_batches <- c("0", "1", "2", "3", "4", "5", "6", "7")
groupB_batches <- c("8", "9", "10", "11", "12", "13", "14", "15")

# Create a dataframe with group labels
sample_conditions <- data.frame(
  sample_id = c(groupA_batches, groupB_batches),
  group = rep(c("GroupA", "GroupB"), each = 8)  # Exactly 8 per group
)
rownames(sample_conditions) <- sample_conditions$sample_id

# Aggregate by batch, 16 batches in total, we get 8 vs 8 mock null DA analysis
OE_PB_null <- AggregateExpression(OE_subset, 
                                  group.by = "batch",  
                                  normalization.method = "LogNormalize", 
                                  scale.factor = 10000, 
                                  return.seurat = TRUE)

# Extract pseudobulk counts
counts_mat_null <- as.matrix(OE_PB_null@assays$ATAC$counts)
colnames(counts_mat_null) <- sample_conditions$group
library(edgeR)

# Create DGEList object
dge_null <- DGEList(counts = counts_mat_null, group = sample_conditions$group)

# Normalize counts
dge_null <- calcNormFactors(dge_null)

# Estimate dispersion
design_null <- model.matrix(~ group, data = sample_conditions)
dge_null <- estimateDisp(dge_null, design_null)

# Fit negative binomial distribution for null comparison
fit_null <- glmQLFit(dge_null, design_null)

# Perform DA test (should be null)
da_results_null <- glmQLFTest(fit_null, coef = 2)

# Extract log fold changes (should be centered around 0)
df_lfc_null <- data.frame(peak = rownames(da_results_null$table), 
                          logFC = da_results_null$table$logFC)


# Extract GC fractions for peaks
gc_fraction <- letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, granges(OE_subset)), 
                               letters="GC", as.prob=TRUE)
gc_fraction <- as.numeric(gc_fraction)

# Assign GC bins
df_gc_null <- data.frame(peak = rownames(counts_mat_null), GC_Fraction = gc_fraction)
df_gc_null$GC_Bin <- cut(df_gc_null$GC_Fraction, 
                         breaks = seq(min(df_gc_null$GC_Fraction), max(df_gc_null$GC_Fraction), length.out = 20), 
                         include.lowest = TRUE)

# Merge DA results with GC bins
df_gc_fc_null <- left_join(df_lfc_null, df_gc_null, by = "peak")

# Remove NAs
df_gc_fc_null <- df_gc_fc_null %>% filter(!is.na(GC_Bin))

library(ggplot2)
library(viridis)
#(8 vs. 8 Mock Null DA Test in Neuron Like Cells)
non_GC_adjusted_bins <- ggplot(df_gc_fc_null, aes(x = GC_Bin, y = logFC, fill = GC_Bin)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  geom_smooth(aes(group = 1), 
              method = "gam", formula = y ~ s(x, bs = "cs"), 
              color = "blue", size = 1.2) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_fill_viridis_d(name = "GC Bin") +
  labs(
    title = "Non-GC Adjusted Content vs. log₂FC",
    x = "GC-content bin",
    y = "log₂FC"
  ) +
  theme_minimal() +
  theme_bw(base_size = 18) +
  theme(
    # Make X-axis text vertical if needed:
    axis.text.x   = element_text(angle = 90, hjust = 1, colour = "black", face = "bold"),
    
    axis.text.y   = element_text(angle = 90, hjust = 1, colour = "black", face = "bold"),
    
    # Bold + center the title:
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    
    # Bold axis titles:
    axis.title    = element_text(face = "bold"),
    
    # Bold legend title:
    legend.title  = element_text(face = "bold"),
    legend.position = "none"
  )

########### Test out EDA normalisation on ATAC seq for edgeR
# extract count data
counts_mat <- as.matrix(OE_PB_null@assays$ATAC$counts)
# Calculate GC fraction
gc_fraction <- letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, granges(OE_subset)), 
                               letters="GC", as.prob=TRUE)
gc_fraction <- as.numeric(gc_fraction)
peak_length <- width(granges(OE_subset))

library(Biobase)
# Create SeqExpressionSet
featureData = data.frame(GC = gc_fraction, Length = peak_length)
rownames(featureData) <- rownames(counts_mat)

# create object of the class SeqExpressionSet.
seq_obj <- newSeqExpressionSet(counts = counts_mat,
                               featureData = featureData)

# Within-sample GC correction
seq_obj_gc <- withinLaneNormalization(seq_obj, "GC", which = "full", offset = TRUE)

# Between-sample library size normalization
seq_obj_gc_size <- betweenLaneNormalization(seq_obj_gc, which = "full", offset = TRUE)

# perform differential accesibility analysis with edgeR
# Create DGEList object
dge_null <- DGEList(counts = counts(seq_obj_gc_size), group = sample_conditions$group)
dge_null$offset <- -offst(seq_obj_gc_size)

# Estimate dispersion
design_null <- model.matrix(~ group, data = sample_conditions)
dge_null <- estimateDisp(dge_null, design_null)

# Fit negative binomial distribution for null comparison
fit_null <- glmQLFit(dge_null, design_null)

# Perform DA test (should be null)
da_results_null <- glmQLFTest(fit_null, coef = 2)

# Extract log fold changes (should be centered around 0)
df_lfc_null <- data.frame(peak = rownames(da_results_null$table), 
                          logFC = da_results_null$table$logFC)


# Extract GC fractions for peaks
gc_fraction <- letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, granges(OE_subset)), 
                               letters="GC", as.prob=TRUE)
gc_fraction <- as.numeric(gc_fraction)

# Assign GC bins
df_gc_null <- data.frame(peak = rownames(featureData), GC_Fraction = gc_fraction)
df_gc_null$GC_Bin <- cut(df_gc_null$GC_Fraction, 
                         breaks = seq(min(df_gc_null$GC_Fraction), max(df_gc_null$GC_Fraction), length.out = 20), 
                         include.lowest = TRUE)

# Merge DA results with GC bins
df_gc_fc_null <- left_join(df_lfc_null, df_gc_null, by = "peak")

# Remove NAs
df_gc_fc_null <- df_gc_fc_null %>% filter(!is.na(GC_Bin))

# plot it
GC_adjusted_bins <- ggplot(df_gc_fc_null, aes(x = GC_Bin, y = logFC, fill = GC_Bin)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  geom_smooth(aes(group = 1), 
              method = "gam", formula = y ~ s(x, bs = "cs"), 
              color = "blue", size = 1.2) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_fill_viridis_d(name = "GC Bin") +
  labs(
    title = "GC Adjusted Content vs. log₂FC",
    x = "GC-content bin",
    y = "log₂FC"
  ) +
  theme_minimal() +
  theme_bw(base_size = 18) +
  theme(
  # Make X-axis text vertical if needed:
  axis.text.x   = element_text(angle = 90, hjust = 1, colour = "black", face = "bold"),
  
  axis.text.y   = element_text(angle = 90, hjust = 1, colour = "black", face = "bold"),
  
  # Bold + center the title:
  plot.title    = element_text(hjust = 0.5, face = "bold"),
  
  # Bold axis titles:
  axis.title    = element_text(face = "bold"),
  
  # Bold legend title:
  legend.title  = element_text(face = "bold"),
  legend.position = "right",
  legend.text = element_text(face = "bold")
)
library(patchwork)
library(ggplot2)

# Wrap for patchwork compatibility
histogram_wrapped <- wrap_elements(full = GC_fraction_histogram)
non_GC_wrapped <- wrap_elements(full = non_GC_adjusted_bins)
GC_wrapped <- wrap_elements(full = GC_adjusted_bins)

# Create the top row with custom relative widths
top_row <- plot_spacer() + 
  histogram_wrapped + 
  plot_spacer() +
  plot_layout(widths = c(1, 3, 1)) 

# Combine with the bottom row
final_plot <- top_row / (non_GC_wrapped + GC_wrapped) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 18))

final_plot
##########
# IMAGE analysis and feature selection of peaks, and comparison and validation of IMAGE and chromVAR
#########
# load libraries
library(BSgenome.Hsapiens.UCSC.hg38) 
library(GenomicRanges)
library(Biostrings)
#install.packages("jpeg")
library(jpeg)
#BiocManager::install("ShortRead")
library(ShortRead)
#BiocManager::install("EDASeq")
library(EDASeq)

########## Compute fraction of cells where each peak is detected, 4 replicate
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")

# Compute fraction of cells where each peak is detected, 1%
peaks_to_keep <- c()
for (cl in 1:13) {
  peaks_to_keep <- c(peaks_to_keep, names(which((rowSums(SHARE_seq_subsample_data@assays$ATAC@counts[,SHARE_seq_subsample_data$seurat_clusters == cl] > 0) / sum(SHARE_seq_subsample_data$seurat_clusters == cl)) >= 0.01)))
}
peaks_to_keep <- unique(peaks_to_keep)
# make clusters
SHARE_seq_subsample_data = SetIdent(SHARE_seq_subsample_data, value = SHARE_seq_subsample_data$seurat_clusters)
new.cluster.ids <- c("hESC 1", "hESC 2", "hESC 3",
                     "Neuron like cells", "T/P like cells", "Goblet like cells", 
                     "OP like cells", "VE like cells",
                     "FL like cell", "Podocyte like cells", "SMP like cells", 
                     "BE like cells", "Adipocyte like cells")
names(new.cluster.ids) <- levels(SHARE_seq_subsample_data)
SHARE_seq_subsample_data <- RenameIdents(SHARE_seq_subsample_data, new.cluster.ids)
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")
metadata <- OE@meta.data

# Divide batches into 4 replicates in total
library(dplyr)
# Add a new column based on batch ranges
metadata <- metadata %>%
  mutate(batch_group = case_when(
    batch %in% c(0, 1, 2, 3) ~ 1,
    batch %in% c(4, 5, 6, 7) ~ 2,
    batch %in% c(8, 9, 10, 11) ~ 3,
    batch %in% c(12, 13, 14, 15) ~ 4
  ))
# Add the new batch group column to the Seurat object
OE$batch_group <- metadata$batch_group

# change back to working with peaks instead of gene activities
OE_PB2 = AggregateExpression(OE, group.by = c("ident", "batch_group"),
                             normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             return.seurat = TRUE)

colnames_updated <- make.unique(gsub("_.*", "", colnames(OE_PB2)))
colnames(OE_PB2) <- colnames_updated

# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(OE_PB2)) # Removes ".x" suffix
Idents(OE_PB2) <- group_labels

DefaultAssay(OE_PB2) <- "ATAC"
OE_PB2 <- subset(OE_PB2, features = peaks_to_keep)
bulk_ATAC_data <- OE_PB2@assays$ATAC$counts

# Make bulk data into summarized experiment object for ATAC seq
cd = data.frame(names = colnames(bulk_ATAC_data))
# Create a vector of unique conditions based on unique cell types
unique_cell_types <- unique(sub("\\.\\d+$", "", cd$names))
conditions <- paste0("Cond", seq_along(unique_cell_types))

# Create a mapping between cell types and conditions
cell_type_to_condition <- setNames(conditions, unique_cell_types)

# Assign conditions based on cell type
cd$condition <- cell_type_to_condition[sub("\\.\\d+$", "", cd$names)]

rownames(cd) = cd$names
cd = DataFrame(cd)

library(Matrix)
# Extract row names from ATAC
row_names <- rownames(bulk_ATAC_data)
# Create vectors to store the parsed chromosome, start, and end positions from rownames
chromosomes <- vector("character", length(row_names))
starts <- numeric(length(row_names))
ends <- numeric(length(row_names))
for (i in seq_along(row_names)) {
  parts <- strsplit(row_names[i], "-")[[1]]
  chromosomes[i] <- parts[1]
  starts[i] <- as.numeric(parts[2])
  ends[i] <- as.numeric(parts[3])
}
ATAC = SummarizedExperiment(assays = list(Peaks = bulk_ATAC_data), colData = cd, rowRanges = GRanges(seqnames = chromosomes, IRanges(start = starts, end = ends), peak_id = rownames(bulk_ATAC_data)))
#saveRDS(ATAC, "/work/Home/Data_for_masters_project/ATAC_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_GC_content_4_replicate_fraction_of_cells")
ATAC <- readRDS("/work/Home/Data_for_masters_project/ATAC_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_GC_content_4_replicate_fraction_of_cells")
###### load motifs
library(chromVARmotifs)
data("human_pwms_v2")
###### prepare IMAGE functions for step 1
countMotifs = function(ATAC, genome, pwms, p.cutoff = 5e-05, n_cores) {
  
  ### Import the genome
  genome.ranges = BSgenome::getBSgenome(genome)
  
  ### Check compatibility between ranges and genome
  rangesForMotifs = ATAC@rowRanges
  ATAC.seqnames = as.character(rangesForMotifs@seqnames@values)
  genome.seqnames = seqnames(genome.ranges)
  
  # Adjust for 'chr' prefix
  if (any(grep("chr", ATAC.seqnames))) {
    if (!any(grep("chr", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("chr", "", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("chr", "", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  } else {
    if (any(grep("chr", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(paste0("chr", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = paste0("chr", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  }
  
  # Adjust for chromosome M/MT naming
  if (any(grep("MT", ATAC.seqnames))) {
    if (!any(grep("MT", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("MT", "M", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("MT", "M", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  } else {
    if (any(grep("MT", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("M", "MT", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("M", "MT", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  }
  
  # Setup parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ############################################ IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Set library paths on worker nodes
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load libraries on worker nodes
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
    # Add any other required libraries
  })
  
  # Check ranges against genome bounds
  genome.size <- data.frame(chr = seqnames(genome.ranges), length = seqlengths(genome.ranges))
  keep <- foreach(i = seq_along(rangesForMotifs), .combine = c, .packages = c("GenomicRanges")) %dopar% {
    range_df <- as.data.frame(rangesForMotifs[i])
    chr <- as.character(range_df[1, 1])
    end_position <- range_df[1, 3]
    chr_length <- genome.size[genome.size$chr == chr, "length"]
    if (length(chr_length) > 0 && end_position <= chr_length) i else NULL
  }
  
  # Stop parallel backend
  stopCluster(cl)
  keep <- keep[!sapply(keep, is.null)]
  rangesForMotifs_clean <- rangesForMotifs[keep]
  
  # Match motifs using PWMatrixList
  Countmatrix = motifmatchr::matchMotifs(pwms, rangesForMotifs_clean,
                                         genome = genome,
                                         bg = "subject",
                                         out = "scores",
                                         p.cutoff = p.cutoff)
  
  # Create a count matrix
  counts = motifmatchr::motifCounts(Countmatrix)
  colnames(counts) = colData(Countmatrix)$name
  colnames(counts) <- paste0(colnames(counts), ".motif")
  return(counts)
}
processMotifs = function(motifs) {
  # Normalize motif matrix
  Cols <- Matrix::colMeans(motifs)
  motifs <- t(t(motifs) - Cols)
  motifs <- as.data.frame(as.matrix(motifs))
  
  # Return
  return(motifs)
}
normalizeEnhancers = function(ATAC, normalize = TRUE, standardize = TRUE, genome) {
  # Extract data
  counts_mat <- as.matrix(SummarizedExperiment::assays(ATAC, withDimnames = FALSE)$Peaks)
  
  # calculate GC frquency
  gc_fraction <- letterFrequency(getSeq(genome, granges(ATAC)), 
                                 letters="GC", as.prob=TRUE)
  gc_fraction <- as.numeric(gc_fraction)
  # Extract peak lengths
  peak_length <- width(granges(ATAC))
  
  # Create SeqExpressionSet
  featureData = data.frame(GC = gc_fraction, Length = peak_length)
  rownames(featureData) <- rownames(counts_mat)
  seq_obj <- newSeqExpressionSet(counts = counts_mat,
                                 featureData = featureData)
  
  seq_obj_gc <- withinLaneNormalization(seq_obj, "GC", which = "full")
  
  seq_obj_gc_size <- betweenLaneNormalization(seq_obj_gc, which = "full")
  
  normalized_counts_fqfq <- normCounts(seq_obj_gc_size)
  
  # Normalize
  if (normalize) {
    normalized_counts_fqfq = t(t(normalized_counts_fqfq) / (colSums(normalized_counts_fqfq) / 10000000))
  }
  
  # Standardize
  if (standardize) {
    Cols <- colMeans(normalized_counts_fqfq)
    SD <- apply(normalized_counts_fqfq[,c(1:ncol(normalized_counts_fqfq))],2,FUN="sd")
    normalized_counts_fqfq <- t( ((t(normalized_counts_fqfq) - Cols) / SD) )
  }
  
  # Return
  return(normalized_counts_fqfq)
}
calculateEnhancerMotifActivity <- function(motifs, enhancers, a, seed = 42, ncores) {
  # Define variable with the number of motifs and samples
  nmotifs <- ncol(motifs)
  nsamples <- ncol(enhancers)
  
  common_peaks <- intersect(rownames(enhancers), rownames(motifs))
  
  # subset the objects to have common peaks
  motifs <- motifs[common_peaks,]
  enhancers <- enhancers[common_peaks,]
  
  # Setup parallel computation
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  ############################################ IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Set library paths on worker nodes
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load libraries on worker nodes
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
  })
  # Set seed for reproducibility (this sets the seed for all workers)
  set.seed(seed)

  # Use foreach to loop over samples in parallel
  results <- foreach(samp = 1:nsamples, .combine = cbind, .packages = c("glmnet", "Matrix")) %dopar% {
    # Set seed for each worker to ensure reproducibility
    set.seed(seed + samp)
    
    # Perform 3-fold cross-validation to find the optimal lambda
    cv <- glmnet::cv.glmnet(
      y = as.numeric(enhancers[, samp]),
      x = as.matrix(motifs),
      alpha = a,
      standardize = FALSE,
      intercept = TRUE,
      parallel = FALSE,
      standardize.response = FALSE,
      nfolds = 3
    )
    
    # Extract the optimal lambda
    lambda <- cv$lambda.min
    
    # Capture the regression coefficients (excluding intercept)
    coef_tmp <- as.data.frame(as.matrix(coef(cv, s = "lambda.min")))
    sparse_coef <- Matrix(coef_tmp[c(2:nrow(coef_tmp)), 1], sparse = TRUE)  # Skip intercept
    
    # Return the sparse coefficients for this sample
    sparse_coef
  }
  
  # Insert results into sparse matrix
  coefmat <- results
  
  # Stop the parallel cluster
  stopCluster(cl)
  
  # Assign proper row and column names
  colnames(coefmat) <- colnames(enhancers)
  rownames(coefmat) <- colnames(motifs)
  
  # Return the coefficient matrix
  return(coefmat)
}

### Count Motifs
genome <- "BSgenome.Hsapiens.UCSC.hg38"
start_time1 <- Sys.time()
motif_counts <- countMotifs(ATAC, genome = genome, pwms = human_pwms_v2, p.cutoff = 5e-05, n_cores = 2)
end_time1 <- Sys.time()
time_motif_counts <- end_time1 - start_time1
#saveRDS(motif_counts, "/work/Home/Data_for_masters_project/motif_counts_chromVAR_motifs_4_replicates_GC_fraction_of_cells")
#saveRDS(time_motif_counts, "/work/Home/Data_for_masters_project/time_motif_counts_chromVAR_motifs_4_replicates_GC_fraction_of_cells")

### Normalize Enhancers
genome <- BSgenome.Hsapiens.UCSC.hg38
start_time2 <- Sys.time()
normEnhancers <- normalizeEnhancers(ATAC, normalize = TRUE, standardize = TRUE, genome = genome)
end_time2 <- Sys.time()
time_normEnhancers <- end_time2 - start_time2
#saveRDS(normEnhancers, "/work/Home/Data_for_masters_project/normEnhancers_chromVAR_motifs_4_replicates_GC_fraction_of_cells")
#saveRDS(time_normEnhancers, "/work/Home/Data_for_masters_project/time_normEnhancers_chromVAR_motifs_4_replicates_GC_fraction_of_cells")

### Process Motifs
start_time3 <- Sys.time()
filteredMotifs <- processMotifs(motifs = motif_counts)
end_time3 <- Sys.time()
time_filteredMotifs <- end_time3 - start_time3
#saveRDS(filteredMotifs, "/work/Home/Data_for_masters_project/filteredMotifs_chromVAR_motifs_4_replicates_GC_fraction_of_cells")
#saveRDS(time_filteredMotifs, "/work/Home/Data_for_masters_project/time_filteredMotifs_chromVAR_motifs_4_replicates_GC_fraction_of_cells")

### Calculate Enhancer Motif Activity
alpha <- c(0, 0.25, 0.5, 0.75, 1)
activity_list <- list()

for (a in alpha) {
  start_time4 <- Sys.time()
  activities <- calculateEnhancerMotifActivity(motifs = filteredMotifs, enhancers = normEnhancers, ncores = 2, a = a)
  end_time4 <- Sys.time()
  time_activities <- end_time4 - start_time4

  # Store in a structured list
  activity_list[[as.character(a)]] <- list(
    activities = activities,
    time = time_activities
  )
  #saveRDS(activity_list, "/work/Home/Data_for_masters_project/activities_list_chromVAR_motifs_4_replicates_GC_fraction_of_cells")

}

####################### For all peaks, 4 replicates
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)
#install.packages("jpeg")
library(jpeg)
#BiocManager::install("ShortRead")
library(ShortRead)
#BiocManager::install("EDASeq")
library(EDASeq)
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")

# make clusters, 4 replicate
SHARE_seq_subsample_data = SetIdent(SHARE_seq_subsample_data, value = SHARE_seq_subsample_data$seurat_clusters)
new.cluster.ids <- c("hESC 1", "hESC 2", "hESC 3",
                     "Neuron like cells", "T/P like cells", "Goblet like cells", 
                     "OP like cells", "VE like cells",
                     "FL like cell", "Podocyte like cells", "SMP like cells", 
                     "BE like cells", "Adipocyte like cells")
names(new.cluster.ids) <- levels(SHARE_seq_subsample_data)
SHARE_seq_subsample_data <- RenameIdents(SHARE_seq_subsample_data, new.cluster.ids)
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")
metadata <- OE@meta.data

# Divide batches into 4 replicates
library(dplyr)
# Add a new column based on batch ranges
metadata <- metadata %>%
  mutate(batch_group = case_when(
    batch %in% c(0, 1, 2, 3) ~ 1,
    batch %in% c(4, 5, 6, 7) ~ 2,
    batch %in% c(8, 9, 10, 11) ~ 3,
    batch %in% c(12, 13, 14, 15) ~ 4
  ))
# Add the new batch group column to the Seurat object
OE$batch_group <- metadata$batch_group

# change back to working with peaks instead of gene activities
OE_PB2 = AggregateExpression(OE, group.by = c("ident", "batch_group"),
                             normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             return.seurat = TRUE)

colnames_updated <- make.unique(gsub("_.*", "", colnames(OE_PB2)))
colnames(OE_PB2) <- colnames_updated

# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(OE_PB2)) # Removes ".x" suffix
Idents(OE_PB2) <- group_labels

bulk_ATAC_data <- OE_PB2@assays$ATAC$counts

# Make bulk data into summarized experiment object for ATAC seq
cd = data.frame(names = colnames(bulk_ATAC_data))
# Create a vector of unique conditions based on unique cell types
unique_cell_types <- unique(sub("\\.\\d+$", "", cd$names))
conditions <- paste0("Cond", seq_along(unique_cell_types))

# Create a mapping between cell types and conditions
cell_type_to_condition <- setNames(conditions, unique_cell_types)

# Assign conditions based on cell type
cd$condition <- cell_type_to_condition[sub("\\.\\d+$", "", cd$names)]

rownames(cd) = cd$names
cd = DataFrame(cd)

library(Matrix)
# Extract row names from ATAC
row_names <- rownames(bulk_ATAC_data)
# Create vectors to store the parsed chromosome, start, and end positions from rownames
chromosomes <- vector("character", length(row_names))
starts <- numeric(length(row_names))
ends <- numeric(length(row_names))
for (i in seq_along(row_names)) {
  parts <- strsplit(row_names[i], "-")[[1]]
  chromosomes[i] <- parts[1]
  starts[i] <- as.numeric(parts[2])
  ends[i] <- as.numeric(parts[3])
}
ATAC = SummarizedExperiment(assays = list(Peaks = bulk_ATAC_data), colData = cd, rowRanges = GRanges(seqnames = chromosomes, IRanges(start = starts, end = ends), peak_id = rownames(bulk_ATAC_data)))
#saveRDS(ATAC, "/work/Home/Data_for_masters_project/ATAC_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_GC_content_4_replicate_all_peaks")
ATAC <- readRDS("/work/Home/Data_for_masters_project/ATAC_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_GC_content_4_replicate_all_peaks")

######## Load motifs
library(chromVARmotifs)
data("human_pwms_v2")

########### prepare IMAGE functions for step 1
countMotifs = function(ATAC, genome, pwms, p.cutoff = 5e-05, n_cores) {
  
  ### Import the genome
  genome.ranges = BSgenome::getBSgenome(genome)
  
  ### Check compatibility between ranges and genome
  rangesForMotifs = ATAC@rowRanges
  ATAC.seqnames = as.character(rangesForMotifs@seqnames@values)
  genome.seqnames = seqnames(genome.ranges)
  
  # Adjust for 'chr' prefix
  if (any(grep("chr", ATAC.seqnames))) {
    if (!any(grep("chr", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("chr", "", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("chr", "", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  } else {
    if (any(grep("chr", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(paste0("chr", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = paste0("chr", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  }
  
  # Adjust for chromosome M/MT naming
  if (any(grep("MT", ATAC.seqnames))) {
    if (!any(grep("MT", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("MT", "M", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("MT", "M", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  } else {
    if (any(grep("MT", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("M", "MT", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("M", "MT", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  }
  
  # Setup parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  # Set library paths on worker nodes
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load libraries on worker nodes
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
    # Add any other required libraries
  })
  
  # Check ranges against genome bounds
  genome.size <- data.frame(chr = seqnames(genome.ranges), length = seqlengths(genome.ranges))
  keep <- foreach(i = seq_along(rangesForMotifs), .combine = c, .packages = c("GenomicRanges")) %dopar% {
    range_df <- as.data.frame(rangesForMotifs[i])
    chr <- as.character(range_df[1, 1])
    end_position <- range_df[1, 3]
    chr_length <- genome.size[genome.size$chr == chr, "length"]
    if (length(chr_length) > 0 && end_position <= chr_length) i else NULL
  }
  
  # Stop parallel backend
  stopCluster(cl)
  keep <- keep[!sapply(keep, is.null)]
  rangesForMotifs_clean <- rangesForMotifs[keep]
  
  # Match motifs using PWMatrixList
  Countmatrix = motifmatchr::matchMotifs(pwms, rangesForMotifs_clean,
                                         genome = genome,
                                         bg = "subject",
                                         out = "scores",
                                         p.cutoff = p.cutoff)
  
  # Create a count matrix
  counts = motifmatchr::motifCounts(Countmatrix)
  colnames(counts) = colData(Countmatrix)$name
  colnames(counts) <- paste0(colnames(counts), ".motif")
  return(counts)
}
processMotifs = function(motifs) {
  # Normalize motif matrix
  Cols <- Matrix::colMeans(motifs)
  motifs <- t(t(motifs) - Cols)
  motifs <- as.data.frame(as.matrix(motifs))
  
  # Return
  return(motifs)
}
normalizeEnhancers = function(ATAC, normalize = TRUE, standardize = TRUE, genome) {
  # Extract data
  counts_mat <- as.matrix(SummarizedExperiment::assays(ATAC, withDimnames = FALSE)$Peaks)
  
  # calculate GC frquency
  gc_fraction <- letterFrequency(getSeq(genome, granges(ATAC)), 
                                 letters="GC", as.prob=TRUE)
  gc_fraction <- as.numeric(gc_fraction)
  # Extract peak lengths
  peak_length <- width(granges(ATAC))
  
  # Create SeqExpressionSet
  featureData = data.frame(GC = gc_fraction, Length = peak_length)
  rownames(featureData) <- rownames(counts_mat)
  seq_obj <- newSeqExpressionSet(counts = counts_mat,
                                 featureData = featureData)
  
  seq_obj_gc <- withinLaneNormalization(seq_obj, "GC", which = "full")
  
  seq_obj_gc_size <- betweenLaneNormalization(seq_obj_gc, which = "full")
  
  normalized_counts_fqfq <- normCounts(seq_obj_gc_size)
  
  # # Normalize
  if (normalize) {
    normalized_counts_fqfq = t(t(normalized_counts_fqfq) / (colSums(normalized_counts_fqfq) / 10000000))
  }
  
  # Standardize
  if (standardize) {
    Cols <- colMeans(normalized_counts_fqfq)
    SD <- apply(normalized_counts_fqfq[,c(1:ncol(normalized_counts_fqfq))],2,FUN="sd")
    normalized_counts_fqfq <- t( ((t(normalized_counts_fqfq) - Cols) / SD) )
  }
  
  # Return
  return(normalized_counts_fqfq)
}
calculateEnhancerMotifActivity <- function(motifs, enhancers, a, seed = 42, ncores) {
  # Define variable with the number of motifs and samples
  nmotifs <- ncol(motifs)
  nsamples <- ncol(enhancers)
  
  common_peaks <- intersect(rownames(enhancers), rownames(motifs))
  
  # subset the objects to have common peaks
  motifs <- motifs[common_peaks,]
  enhancers <- enhancers[common_peaks,]
  
  # Setup parallel computation
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  ############################################ IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Set library paths on worker nodes
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load libraries on worker nodes
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
  })
  # Set seed for reproducibility (this sets the seed for all workers)
  set.seed(seed)
  
  # Use foreach to loop over samples in parallel
  results <- foreach(samp = 1:nsamples, .combine = cbind, .packages = c("glmnet", "Matrix")) %dopar% {
    # Set seed for each worker to ensure reproducibility
    set.seed(seed + samp)
    
    # Perform 3-fold cross-validation to find the optimal lambda
    cv <- glmnet::cv.glmnet(
      y = as.numeric(enhancers[, samp]),
      x = as.matrix(motifs),
      alpha = a,
      standardize = FALSE,
      intercept = TRUE,
      parallel = FALSE,
      standardize.response = FALSE,
      nfolds = 3
    )
    
    # Extract the optimal lambda
    lambda <- cv$lambda.min
    
    # Capture the regression coefficients (excluding intercept)
    coef_tmp <- as.data.frame(as.matrix(coef(cv, s = "lambda.min")))
    sparse_coef <- Matrix(coef_tmp[c(2:nrow(coef_tmp)), 1], sparse = TRUE)  # Skip intercept
    
    # Return the sparse coefficients for this sample
    sparse_coef
  }
  
  # Insert results into sparse matrix
  coefmat <- results
  
  # Stop the parallel cluster
  stopCluster(cl)
  
  # Assign proper row and column names
  colnames(coefmat) <- colnames(enhancers)
  rownames(coefmat) <- colnames(motifs)
  
  # Return the coefficient matrix
  return(coefmat)
}

### Count Motifs
genome <- "BSgenome.Hsapiens.UCSC.hg38"
start_time1 <- Sys.time()
motif_counts <- countMotifs(ATAC, genome = genome, pwms = human_pwms_v2, p.cutoff = 5e-05, n_cores = 2)
end_time1 <- Sys.time()
time_motif_counts <- end_time1 - start_time1

#saveRDS(motif_counts, "/work/Home/Data_for_masters_project/motif_counts_chromVAR_motifs_4_replicates_GC_all_peaks")
#saveRDS(time_motif_counts, "/work/Home/Data_for_masters_project/time_motif_counts_chromVAR_motifs_4_replicates_GC_all_peaks")

### Normalize Enhancers
genome <- BSgenome.Hsapiens.UCSC.hg38
start_time2 <- Sys.time()
normEnhancers <- normalizeEnhancers(ATAC, normalize = TRUE, standardize = TRUE, genome = genome)
end_time2 <- Sys.time()
time_normEnhancers <- end_time2 - start_time2
#saveRDS(normEnhancers, "/work/Home/Data_for_masters_project/normEnhancers_chromVAR_motifs_4_replicates_GC_all_peaks")
#saveRDS(time_normEnhancers, "/work/Home/Data_for_masters_project/time_normEnhancers_chromVAR_motifs_4_replicates_GC_all_peaks")

### Process Motifs
start_time3 <- Sys.time()
filteredMotifs <- processMotifs(motifs = motif_counts)
end_time3 <- Sys.time()
time_filteredMotifs <- end_time3 - start_time3
#saveRDS(filteredMotifs, "/work/Home/Data_for_masters_project/filteredMotifs_chromVAR_motifs_4_replicates_GC_all_peaks")
#saveRDS(time_filteredMotifs, "/work/Home/Data_for_masters_project/time_filteredMotifs_chromVAR_motifs_4_replicates_GC_all_peaks")


### Calculate Enhancer Motif Activity
alpha <- c(0, 0.25, 0.5, 0.75, 1)
activity_list <- list()

for (a in alpha) {
  start_time4 <- Sys.time()
  activities <- calculateEnhancerMotifActivity(motifs = filteredMotifs, enhancers = normEnhancers, ncores = 1, a = a)
  end_time4 <- Sys.time()
  time_activities <- end_time4 - start_time4

  # Store in a structured list
  activity_list[[as.character(a)]] <- list(
    activities = activities,
    time = time_activities
  )
  #saveRDS(activity_list, "/work/Home/Data_for_masters_project/activities_list_chromVAR_motifs_4_replicates_GC_all_peaks")

}

################################ For all peaks and then differential accesible peaks, 4 replicates
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
SHARE_seq_subsample_data = SetIdent(SHARE_seq_subsample_data, value = SHARE_seq_subsample_data$seurat_clusters)
new.cluster.ids <- c("hESC 1", "hESC 2", "hESC 3",
                     "Neuron like cells", "T/P like cells", "Goblet like cells", 
                     "OP like cells", "VE like cells",
                     "FL like cell", "Podocyte like cells", "SMP like cells", 
                     "BE like cells", "Adipocyte like cells")
names(new.cluster.ids) <- levels(SHARE_seq_subsample_data)
SHARE_seq_subsample_data <- RenameIdents(SHARE_seq_subsample_data, new.cluster.ids)
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")
metadata <- OE@meta.data

# Divide batches into 4 replicates
library(dplyr)
# Add a new column based on batch ranges
metadata <- metadata %>%
  mutate(batch_group = case_when(
    batch %in% c(0, 1, 2, 3) ~ 1,
    batch %in% c(4, 5, 6, 7) ~ 2,
    batch %in% c(8, 9, 10, 11) ~ 3,
    batch %in% c(12, 13, 14, 15) ~ 4
  ))
# Add the new batch group column to the Seurat object
OE$batch_group <- metadata$batch_group

# change back to working with peaks instead of gene activities
OE_PB2 = AggregateExpression(OE, group.by = c("ident", "batch_group"),
                             normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             return.seurat = TRUE)
colnames_updated <- make.unique(gsub("_.*", "", colnames(OE_PB2)))
colnames(OE_PB2) <- colnames_updated

# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(OE_PB2)) # Removes ".x" suffix
Idents(OE_PB2) <- group_labels


################################ Run chromvar on bulk 4 replicates
DefaultAssay(OE_PB2) <- 'ATAC'
# Store peak coordinates from the original dataset
peak.ranges <- granges(SHARE_seq_subsample_data)
gc()
library(GenomicRanges)
library(Signac)
# Ensure peaks are stored in a proper ATAC assay for the new bulk object
atac_assay <- CreateChromatinAssay(
  counts = GetAssayData(OE_PB2, slot = "counts"),
  ranges = peak.ranges,
  genome = "hg38"
)
# Create a new Seurat object with peaks
OE_PB2_ATAC <- CreateSeuratObject(
  counts = atac_assay,
  assay = "ATAC"
)
## S3 method for class 'Seurat'
#BiocManager::install("chromVAR")
library(chromVAR)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
data("human_pwms_v2")

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(OE_PB2_ATAC) <- 'ATAC'
motif.matrix <- CreateMotifMatrix(
  features = granges(OE_PB2_ATAC),
  pwm = human_pwms_v2,
  genome = "BSgenome.Hsapiens.UCSC.hg38"
)

set.seed(42)
chromVAR_object <- RunChromVAR(
  OE_PB2_ATAC,
  genome = "BSgenome.Hsapiens.UCSC.hg38",
  motif.matrix = motif.matrix,
  assay = "ATAC",
  new.assay.name = "chromvar"
)
#saveRDS(chromVAR_object,"/work/Home/Data_for_masters_project/chromVAR_bulk_4_replicates_all_peaks")

######### make ATAC summerized object and run it with IMAGE
# extract counts
bulk_ATAC_data <- OE_PB2_ATAC@assays$ATAC$counts
# Make bulk data into summarized experiment object for ATAC seq
cd = data.frame(names = colnames(bulk_ATAC_data))
# Create a vector of unique conditions based on unique cell types
unique_cell_types <- unique(sub("\\.\\d+$", "", cd$names))
conditions <- paste0("Cond", seq_along(unique_cell_types))

# Create a mapping between cell types and conditions
cell_type_to_condition <- setNames(conditions, unique_cell_types)

# Assign conditions based on cell type
cd$condition <- cell_type_to_condition[sub("\\.\\d+$", "", cd$names)]

rownames(cd) = cd$names
cd = DataFrame(cd)

library(Matrix)
# Extract row names from ATAC
row_names <- rownames(bulk_ATAC_data)
# Create vectors to store the parsed chromosome, start, and end positions from rownames
chromosomes <- vector("character", length(row_names))
starts <- numeric(length(row_names))
ends <- numeric(length(row_names))
for (i in seq_along(row_names)) {
  parts <- strsplit(row_names[i], "-")[[1]]
  chromosomes[i] <- parts[1]
  starts[i] <- as.numeric(parts[2])
  ends[i] <- as.numeric(parts[3])
}
ATAC = SummarizedExperiment(assays = list(Peaks = bulk_ATAC_data), colData = cd, rowRanges = GRanges(seqnames = chromosomes, IRanges(start = starts, end = ends), peak_id = rownames(bulk_ATAC_data)))
#saveRDS(ATAC, "/work/Home/Data_for_masters_project/ATAC_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_4_replicates_all peaks")

########### prepare IMAGE functions for step 1
ATAC <- readRDS("/work/Home/Data_for_masters_project/ATAC_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_4_replicates_all peaks")
countMotifs = function(ATAC, genome, pwms, p.cutoff = 5e-05, n_cores) {
  
  ### Import the genome
  genome.ranges = BSgenome::getBSgenome(genome)
  
  ### Check compatibility between ranges and genome
  rangesForMotifs = ATAC@rowRanges
  ATAC.seqnames = as.character(rangesForMotifs@seqnames@values)
  genome.seqnames = seqnames(genome.ranges)
  
  # Adjust for 'chr' prefix
  if (any(grep("chr", ATAC.seqnames))) {
    if (!any(grep("chr", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("chr", "", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("chr", "", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  } else {
    if (any(grep("chr", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(paste0("chr", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = paste0("chr", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  }
  
  # Adjust for chromosome M/MT naming
  if (any(grep("MT", ATAC.seqnames))) {
    if (!any(grep("MT", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("MT", "M", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("MT", "M", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  } else {
    if (any(grep("MT", genome.seqnames))) {
      rangesForMotifs@seqnames@values = factor(gsub("M", "MT", as.character(rangesForMotifs@seqnames@values)))
      rangesForMotifs@seqinfo@seqnames = gsub("M", "MT", as.character(rangesForMotifs@seqinfo@seqnames))
    }
  }
  
  # Setup parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ############################################ IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Set library paths on worker nodes
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load libraries on worker nodes
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
    # Add any other required libraries
  })
  
  # Check ranges against genome bounds
  genome.size <- data.frame(chr = seqnames(genome.ranges), length = seqlengths(genome.ranges))
  keep <- foreach(i = seq_along(rangesForMotifs), .combine = c, .packages = c("GenomicRanges")) %dopar% {
    range_df <- as.data.frame(rangesForMotifs[i])
    chr <- as.character(range_df[1, 1])
    end_position <- range_df[1, 3]
    chr_length <- genome.size[genome.size$chr == chr, "length"]
    if (length(chr_length) > 0 && end_position <= chr_length) i else NULL
  }
  
  # Stop parallel backend
  stopCluster(cl)
  keep <- keep[!sapply(keep, is.null)]
  rangesForMotifs_clean <- rangesForMotifs[keep]
  
  # Match motifs using PWMatrixList
  Countmatrix = motifmatchr::matchMotifs(pwms, rangesForMotifs_clean,
                                         genome = genome,
                                         bg = "subject",
                                         out = "scores",
                                         p.cutoff = p.cutoff)
  
  # Create a count matrix
  counts = motifmatchr::motifCounts(Countmatrix)
  colnames(counts) = colData(Countmatrix)$name
  colnames(counts) <- paste0(colnames(counts), ".motif")
  return(counts)
}
processMotifs = function(motifs) {
  # Normalize motif matrix
  Cols <- Matrix::colMeans(motifs)
  motifs <- t(t(motifs) - Cols)
  motifs <- as.data.frame(as.matrix(motifs))
  
  # Return
  return(motifs)
}
normalizeEnhancers = function(ATAC, normalize = TRUE, standardize = TRUE, genome) {
  #### create design matrix
  # Extract sample names
  sample_names <- colnames(ATAC)
  
  # Extract cell type names dynamically (before the first period)
  cell_types <- gsub("\\..*", "", sample_names)  # Removes numeric suffixes like ".1", ".2"
  
  # Initialize replicate vector
  replicate_numbers <- numeric(length(sample_names))
  
  # Assign replicates by looping through unique cell types
  for (ct in unique(cell_types)) {
    # Find indices where this cell type appears
    indices <- which(cell_types == ct)
    # Assign replicate numbers sequentially (1,2,3,4,...)
    replicate_numbers[indices] <- seq_along(indices)
  }
  
  # Create metadata table
  sample_conditions <- data.frame(
    sample_id = sample_names,
    cell_type = factor(cell_types, levels = unique(cell_types)),  # Convert to factor
    replicate = factor(replicate_numbers)  # Convert to factor
  )
  
  rownames(sample_conditions) <- sample_conditions$sample_id  # Match rownames to sample IDs
  
  
  # Extract data
  counts_mat <- as.matrix(SummarizedExperiment::assays(ATAC, withDimnames = FALSE)$Peaks)
  
  # calculate GC frquency
  gc_fraction <- letterFrequency(getSeq(genome, granges(ATAC)), 
                                 letters="GC", as.prob=TRUE)
  gc_fraction <- as.numeric(gc_fraction)
  # Extract peak lengths
  peak_length <- width(granges(ATAC))
  
  # Create SeqExpressionSet
  featureData = data.frame(GC = gc_fraction, Length = peak_length)
  rownames(featureData) <- rownames(counts_mat)
  seq_obj <- newSeqExpressionSet(counts = counts_mat,
                                 featureData = featureData)
  
  seq_obj_gc <- withinLaneNormalization(seq_obj, "GC", which = "full")
  seq_obj_gc_size <- betweenLaneNormalization(seq_obj_gc, which = "full")
  
  normalized_counts_fqfq <- normCounts(seq_obj_gc_size)
  
  # Create DGEList object with FQ-FQ normalized counts and start calculating DA peaks
  dge_fqfq <- DGEList(counts = normalized_counts_fqfq, group = sample_conditions$cell_type)
  
  # Normalize
  dge_fqfq <- calcNormFactors(dge_fqfq)
  
  # Build design matrix with both cell type and replicate
  design <- model.matrix(~ cell_type + replicate, data = sample_conditions)
  
  # Estimate dispersion
  dge_fqfq <- estimateDisp(dge_fqfq, design)
  
  # Fit the negative binomial model
  fit <- glmQLFit(dge_fqfq, design)
  
  # Run DA test
  da_results_gc <- glmQLFTest(fit, coef = 2:ncol(design))
  
  # Extract DA results
  df_da_gc <- data.frame(
    peak = rownames(da_results_gc$table),
    logFC = rowMeans(da_results_gc$table[, grep("^logFC.cell_type", colnames(da_results_gc$table))], na.rm = TRUE),
    PValue = da_results_gc$table$PValue,
    FDR = p.adjust(da_results_gc$table$PValue, method = "BH")
  )
  
  
  # Keep only significantly DA peaks
  significant_peaks <- df_da_gc %>% filter(FDR < 0.05)
  
  # Keep only significantly DA peaks
  normalized_counts_da <- normalized_counts_fqfq[rownames(normalized_counts_fqfq) %in% significant_peaks$peak, ]

  # # Normalize
  if (normalize) {
    normalized_counts_da = t(t(normalized_counts_da) / (colSums(normalized_counts_da) / 10000000))
  }
  
  # Standardize
  if (standardize) {
    Cols <- colMeans(normalized_counts_da)
    SD <- apply(normalized_counts_da[,c(1:ncol(normalized_counts_da))],2,FUN="sd")
    normalized_counts_da <- t( ((t(normalized_counts_da) - Cols) / SD) )
  }
  
  # Return
  return(normalized_counts_da)
}
calculateEnhancerMotifActivity <- function(motifs, enhancers, a, seed = 42, ncores) {
  # Define variable with the number of motifs and samples
  nmotifs <- ncol(motifs)
  nsamples <- ncol(enhancers)
  
  common_peaks <- intersect(rownames(enhancers), rownames(motifs))
  
  # subset the objects to have common peaks
  motifs <- motifs[common_peaks,]
  enhancers <- enhancers[common_peaks,]
  
  # Setup parallel computation
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  ############################################ IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Set library paths on worker nodes
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load libraries on worker nodes
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
    # Add any other required libraries
  })
  # Set seed for reproducibility (this sets the seed for all workers)
  set.seed(seed)

  # Use foreach to loop over samples in parallel
  results <- foreach(samp = 1:nsamples, .combine = cbind, .packages = c("glmnet", "Matrix")) %dopar% {
    # Set seed for each worker to ensure reproducibility
    set.seed(seed + samp)
    
    # Perform 3-fold cross-validation to find the optimal lambda
    cv <- glmnet::cv.glmnet(
      y = as.numeric(enhancers[, samp]),
      x = as.matrix(motifs),
      alpha = a,
      standardize = FALSE,
      intercept = TRUE,
      parallel = FALSE,
      standardize.response = FALSE,
      nfolds = 3
    )
    
    # Extract the optimal lambda
    lambda <- cv$lambda.min
    
    # Capture the regression coefficients (excluding intercept)
    coef_tmp <- as.data.frame(as.matrix(coef(cv, s = "lambda.min")))
    sparse_coef <- Matrix(coef_tmp[c(2:nrow(coef_tmp)), 1], sparse = TRUE)  # Skip intercept
    
    # Return the sparse coefficients for this sample
    sparse_coef
  }
  
  # Insert results into sparse matrix
  coefmat <- results
  
  # Stop the parallel cluster
  stopCluster(cl)
  
  # Assign proper row and column names
  colnames(coefmat) <- colnames(enhancers)
  rownames(coefmat) <- colnames(motifs)
  
  # Return the coefficient matrix
  return(coefmat)
}

# load genome and motifs
genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(chromVARmotifs)
data("human_pwms_v2")

start_time1 <- Sys.time()
motif_counts = countMotifs(ATAC, genome = genome, pwms = human_pwms_v2, p.cutoff = 5e-05, n_cores = 10)
end_time1 <- Sys.time()
time_motif_counts <- end_time1-start_time1
#saveRDS(motif_counts, "/work/Home/Data_for_masters_project/motif_counts_chromVAR_motifs_4_replicates_all_peaks")
#saveRDS(time_motif_counts, "/work/Home/Data_for_masters_project/time_motif_counts_chromVAR_motifs_4_replicates_all_peaks")
#motif_counts <- readRDS("/work/Home/Data_for_masters_project/motif_counts_chromVAR_motifs_4_replicates_all_peaks")
#time_motif_counts <- readRDS("/work/Home/Data_for_masters_project/time_motif_counts_chromVAR_motifs_4_replicates_all_peaks")

start_time2 <- Sys.time()
filteredMotifs = processMotifs(motifs = motif_counts)
end_time2 <- Sys.time()
time_filteredMotifs <- end_time2-start_time2
#saveRDS(filteredMotifs, "/work/Home/Data_for_masters_project/filteredMotifs_chromVAR_motifs_4_replicates_all_peaks")
#saveRDS(time_filteredMotifs, "/work/Home/Data_for_masters_project/time_filteredMotifs_chromVAR_motifs_4_replicates_all_peaks")
#filteredMotifs <- readRDS("/work/Home/Data_for_masters_project/filteredMotifs_chromVAR_motifs_4_replicates_all_peaks")
#time_filteredMotifs <- readRDS("/work/Home/Data_for_masters_project/time_filteredMotifs_chromVAR_motifs_4_replicates_all_peaks")

genome <- BSgenome.Hsapiens.UCSC.hg38
start_time3 <- Sys.time()
normEnhancers <- normalizeEnhancers(ATAC, normalize = TRUE, standardize = TRUE, genome = genome)
end_time3 <- Sys.time()
time_normEnhancers <- end_time3-start_time3
#saveRDS(time_normEnhancers, "/work/Home/Data_for_masters_project/time_normEnhancers_chromVAR_motifs_4_replicates_all_peaks_GC_and_DA")
#saveRDS(normEnhancers, "/work/Home/Data_for_masters_project/normEnhancers_chromVAR_motifs_4_replicates_all_peaks_GC_and_DA")
#normEnhancers <- readRDS("/work/Home/Data_for_masters_project/normEnhancers_chromVAR_motifs_4_replicates_all_peaks_GC_and_DA")

start_time4 <- Sys.time()
activities = calculateEnhancerMotifActivity(motifs = filteredMotifs, enhancers = normEnhancers, ncores = 2, a = 0)
end_time4 <- Sys.time()
time_activities <- end_time4-start_time4
#saveRDS(activities, "/work/Home/Data_for_masters_project/activities_chromVAR_motifs_GC_and_edgeR")
#saveRDS(time_activities, "/work/Home/Data_for_masters_project/time_activities_chromVAR_motifs_GC_and_edgeR")
#activities <- readRDS("/work/Home/Data_for_masters_project/activities_chromVAR_motifs_GC_and_edgeR")

########################################## plot IMAGE motif activities
########## for all peaks
activities_list <- readRDS("/work/Home/Data_for_masters_project/activities_list_chromVAR_motifs_4_replicates_GC_all_peaks")
activities <- activities_list[["0"]]$activities
activities <- as.matrix(activities)
# Load required library
library(Matrix)
# Convert all entries to numeric
activities <- matrix(
  as.numeric(activities),
  nrow = nrow(activities),
  ncol = ncol(activities),
  dimnames = dimnames(activities)  # preserve row/column names
)
# Extract base column names (removing suffixes like .1, .2, etc.)
base_colnames <- sub("\\.\\d+$", "", colnames(activities))

# Compute mean over each group
activities_mean <- sapply(unique(base_colnames), function(cell_type) {
  col_indices <- which(base_colnames == cell_type)
  rowMeans(activities[, col_indices, drop = FALSE])
})

# Assign new column names
colnames(activities_mean) <- c("hESC 1", "hESC 2", "hESC 3",
                               "Neuron like cells", "T/P like cells", "Goblet like cells", 
                               "OP like cells", "VE like cells",
                               "FL like cell", "Podocyte like cells", "SMP like cells", 
                               "BE like cells", "Adipocyte like cells")

### Map the motifs to TFs
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
# These specific motifs were not handled correctly, due to their naming. They are manuallt being handled
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

# Sort the rownames of activity_object
# Create a mapping vector with old rownames as keys and new rownames as values
rename_mapping <- c(
  "LINE11277.motif" = "DUX1.motif",
  "LINE11282.motif" = "DUX3.motif",
  "LINE4118.motif" = "ZNF75C.motif"
)

# Replace matching rownames
rownames(activities_mean) <- ifelse(
  rownames(activities_mean) %in% names(rename_mapping),
  rename_mapping[rownames(activities_mean)],
  rownames(activities_mean)  # Keep other rownames unchanged
)

activities_mean <- activities_mean[order(rownames(activities_mean)), ]
# Sort the rownames of mapping based on motif
mapping <- mapping[order(mapping$V2), ]
# give rownames to activities
rownames(activities_mean) <- mapping$V1

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
# convert it to matrix
heatmap_data <- as.matrix(activities_mean)
# Z-score the data
heatmap_data <- t(scale(t(heatmap_data)))
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]
#col_fun <- colorRamp2(c(min(heatmap_data), 0, max(heatmap_data)), c("#4575B4", "#91BFDB", "white", "#FCAE61", "#D73027"))
# Define breaks manually using the min, 0, and max values with intermediate points
col_fun <- colorRamp2(
  breaks = c(min(heatmap_data),  # lowest
             quantile(heatmap_data, 0.33),  # lower intermediate
             0,                           # mid point
             quantile(heatmap_data, 0.67),  # upper intermediate
             max(heatmap_data)),          # highest
  colors = c("#4575B4", "#91BFDB", "white", "#FCAE61", "#D73027")
)

BASE_FONTSIZE <- 14
All_peaks <- Heatmap(
  heatmap_data,
  col            = col_fun,
  na_col         = "grey",
  ## clustering & dendrogram placement
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  column_dend_side = "top",
  ## row / column labels & titles
  row_names_side    = "left",
  column_names_side = "bottom",
  show_row_names    = FALSE,
  column_names_gp   = gpar(fontsize = BASE_FONTSIZE,   fontface = "bold"),
  row_names_gp      = gpar(fontsize = BASE_FONTSIZE-2),
  column_title      = "IMAGE Enhancer Motif Activity Heatmap\n(4 Replicates, All Peaks, Alpha = 0)",
  row_title         = "Transcription Factors",
  column_title_gp   = gpar(fontsize = BASE_FONTSIZE+4, fontface = "bold"),
  row_title_gp      = gpar(fontsize = BASE_FONTSIZE+4, fontface = "bold"),
  heatmap_legend_param = list(
    title      = "Z-Score Motif Activity",
    title_gp   = gpar(fontsize = BASE_FONTSIZE,   fontface = "bold"),
    labels_gp  = gpar(fontsize = BASE_FONTSIZE-1, fontface = "bold")
  )
)
All_peaks <- grid.grabExpr(draw(All_peaks))

########## for fraction of cells
activities_list <- readRDS("/work/Home/Data_for_masters_project/activities_list_chromVAR_motifs_4_replicates_GC_fraction_of_cells")
activities <- activities_list[["0"]]$activities
activities <- as.matrix(activities)
# Load required library
library(Matrix)
# Convert all entries to numeric
activities <- matrix(
  as.numeric(activities),
  nrow = nrow(activities),
  ncol = ncol(activities),
  dimnames = dimnames(activities)  # preserve row/column names
)
# Extract base column names (removing suffixes like .1, .2, etc.)
base_colnames <- sub("\\.\\d+$", "", colnames(activities))

# Compute mean over each group
activities_mean <- sapply(unique(base_colnames), function(cell_type) {
  col_indices <- which(base_colnames == cell_type)
  rowMeans(activities[, col_indices, drop = FALSE])
})

# Assign new column names
colnames(activities_mean) <- c("hESC 1", "hESC 2", "hESC 3",
                               "Neuron like cells", "T/P like cells", "Goblet like cells", 
                               "OP like cells", "VE like cells",
                               "FL like cell", "Podocyte like cells", "SMP like cells", 
                               "BE like cells", "Adipocyte like cells")

### Map the motifs to TFs
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

# Sort the rownames of activity_object
# Create a mapping vector with old rownames as keys and new rownames as values
rename_mapping <- c(
  "LINE11277.motif" = "DUX1.motif",
  "LINE11282.motif" = "DUX3.motif",
  "LINE4118.motif" = "ZNF75C.motif"
)

# Replace matching rownames
rownames(activities_mean) <- ifelse(
  rownames(activities_mean) %in% names(rename_mapping),
  rename_mapping[rownames(activities_mean)],
  rownames(activities_mean)  # Keep other rownames unchanged
)

activities_mean <- activities_mean[order(rownames(activities_mean)), ]
# Sort the rownames of mapping based on motif
mapping <- mapping[order(mapping$V2), ]
# give rownames to activities
rownames(activities_mean) <- mapping$V1

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
# convert to matrix
heatmap_data <- as.matrix(activities_mean)
# Z-score the data
heatmap_data <- t(scale(t(heatmap_data)))
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]
#col_fun <- colorRamp2(c(min(heatmap_data), 0, max(heatmap_data)), c("#4575B4", "#91BFDB", "white", "#FCAE61", "#D73027"))
# Define breaks manually using the min, 0, and max values with intermediate points
col_fun <- colorRamp2(
  breaks = c(min(heatmap_data),  # lowest
             quantile(heatmap_data, 0.33),  # lower intermediate
             0,                           # mid point
             quantile(heatmap_data, 0.67),  # upper intermediate
             max(heatmap_data)),          # highest
  colors = c("#4575B4", "#91BFDB", "white", "#FCAE61", "#D73027")
)

BASE_FONTSIZE <- 14     
Fraction_peaks <- Heatmap(
  heatmap_data,
  col            = col_fun,
  na_col         = "grey",
  ## clustering & dendrogram placement
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  column_dend_side = "top",
  ## row / column labels & titles
  row_names_side    = "left",
  column_names_side = "bottom",
  show_row_names    = FALSE,
  column_names_gp   = gpar(fontsize = BASE_FONTSIZE,   fontface = "bold"),
  row_names_gp      = gpar(fontsize = BASE_FONTSIZE-2),
  column_title      = "IMAGE Enhancer Motif Activity Heatmap\n(4 Replicates, ≥1% Cells with Peaks, Alpha = 0)",
  row_title         = "Transcription Factors",
  column_title_gp   = gpar(fontsize = BASE_FONTSIZE+4, fontface = "bold"),
  row_title_gp      = gpar(fontsize = BASE_FONTSIZE+4, fontface = "bold"),
  heatmap_legend_param = list(
    title      = "Z-Score Motif Activity",
    title_gp   = gpar(fontsize = BASE_FONTSIZE,   fontface = "bold"),
    labels_gp  = gpar(fontsize = BASE_FONTSIZE-1, fontface = "bold")
  )
)
Fraction_peaks <- grid.grabExpr(draw(Fraction_peaks))

########## for GC differentially accesible peaks
activities <- readRDS("/work/Home/Data_for_masters_project/activities_chromVAR_motifs_GC_and_edgeR")
activities <- as.matrix(activities)
# Load required library
library(Matrix)
# Convert all entries to numeric
activities <- matrix(
  as.numeric(activities),
  nrow = nrow(activities),
  ncol = ncol(activities),
  dimnames = dimnames(activities)  # preserve row/column names
)
# Extract base column names (removing suffixes like .1, .2, etc.)
base_colnames <- sub("\\.\\d+$", "", colnames(activities))

# Compute mean over each group
activities_mean <- sapply(unique(base_colnames), function(cell_type) {
  col_indices <- which(base_colnames == cell_type)
  rowMeans(activities[, col_indices, drop = FALSE])
})

# Assign new column names
colnames(activities_mean) <- c("hESC 1", "hESC 2", "hESC 3",
                               "Neuron like cells", "T/P like cells", "Goblet like cells", 
                               "OP like cells", "VE like cells",
                               "FL like cell", "Podocyte like cells", "SMP like cells", 
                               "BE like cells", "Adipocyte like cells")

### Map the motifs to TFs
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

# Sort the rownames of activity_object
# Create a mapping vector with old rownames as keys and new rownames as values
rename_mapping <- c(
  "LINE11277.motif" = "DUX1.motif",
  "LINE11282.motif" = "DUX3.motif",
  "LINE4118.motif" = "ZNF75C.motif"
)

# Replace matching rownames
rownames(activities_mean) <- ifelse(
  rownames(activities_mean) %in% names(rename_mapping),
  rename_mapping[rownames(activities_mean)],
  rownames(activities_mean)  # Keep other rownames unchanged
)

activities_mean <- activities_mean[order(rownames(activities_mean)), ]
# Sort the rownames of mapping based on motif
mapping <- mapping[order(mapping$V2), ]
# give rownames to activities
rownames(activities_mean) <- mapping$V1

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
# convert to matrix
heatmap_data <- as.matrix(activities_mean)
# Z-score the data
heatmap_data <- t(scale(t(heatmap_data)))
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]

col_fun <- colorRamp2(
  breaks = c(min(heatmap_data),  # lowest
             quantile(heatmap_data, 0.33),  # lower intermediate
             0,                           # mid point
             quantile(heatmap_data, 0.67),  # upper intermediate
             max(heatmap_data)),          # highest
  colors = c("#4575B4", "#91BFDB", "white", "#FCAE61", "#D73027")
)

BASE_FONTSIZE <- 14     
DA_peaks <- Heatmap(
  heatmap_data,
  col            = col_fun,
  na_col         = "grey",
  ## clustering & dendrogram placement
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  column_dend_side = "top",
  ## row / column labels & titles
  row_names_side    = "left",
  column_names_side = "bottom",
  show_row_names    = FALSE,
  column_names_gp   = gpar(fontsize = BASE_FONTSIZE,   fontface = "bold"),
  row_names_gp      = gpar(fontsize = BASE_FONTSIZE-2),
  column_title      = "IMAGE Enhancer Motif Activity Heatmap\n(4 Replicates, Differential Accessible Regions, Alpha = 0)",
  row_title         = "Transcription Factors",
  column_title_gp   = gpar(fontsize = BASE_FONTSIZE+4, fontface = "bold"),
  row_title_gp      = gpar(fontsize = BASE_FONTSIZE+4, fontface = "bold"),
  heatmap_legend_param = list(
    title      = "Z-Score Motif Activity",
    title_gp   = gpar(fontsize = BASE_FONTSIZE,   fontface = "bold"),
    labels_gp  = gpar(fontsize = BASE_FONTSIZE-1, fontface = "bold")
  )
)
DA_peaks <- grid.grabExpr(draw(DA_peaks))

####### process chromVAR
chromVAR_object <- readRDS("/work/Home/Data_for_masters_project/chromVAR_results/chromVAR_single_cell_all_peaks.rds")
chromVAR_object <- chromVAR_object@assays$chromvar$data
# Load mapping of names
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

# Sort the rownames of chromVAR_object
chromVAR_object <- chromVAR_object[order(rownames(chromVAR_object)), ]

# Sort the rownames of mapping
mapping <- mapping[order(rownames(mapping)), ]

rownames(chromVAR_object) <- mapping$V1
# identify clusters based on barcodes
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
metadata <- SHARE_seq_subsample_data@meta.data
metadata$TF <- sub("^[^-]*-", "", metadata$TF)
metadata$TF <- as.factor(metadata$TF)
metadata_no_control <- metadata[!(metadata$TF %in% c("GFP", "mCherry")), ]
# define clusters
clusters <- sort(unique(metadata_no_control$seurat_clusters))

# Store results
mean_list <- list()

for (cluster in clusters){
  # Get barcodes (cells) for the cluster
  cluster_barcodes <- rownames(metadata_no_control[metadata_no_control$seurat_clusters == cluster, ])
  
  # Get only the cells that exist in chromVAR_object
  cells <- intersect(colnames(chromVAR_object), cluster_barcodes)
  
  if (length(cells) > 0) {  # Ensure there are cells
    # Compute mean deviations across cells for the cluster
    mean_list[[as.character(cluster)]] <- rowMeans(chromVAR_object[, cells, drop = FALSE], na.rm = TRUE)
  } else {
    warning(paste("No matching cells found for cluster", cluster))
  }
}

# Convert the list to a data frame
chromvar_avg_motifactivity <- do.call(cbind, mean_list)
colnames(chromvar_avg_motifactivity) <- c("hESC 1", "hESC 2", "hESC 3",
                                          "Neuron like cells", "T/P like cells", "Goblet like cells", 
                                          "OP like cells", "VE like cells",
                                          "FL like cell", "Podocyte like cells", "SMP like cells", 
                                          "BE like cells", "Adipocyte like cells")
# z-score it
chromvar_avg_motifactivity <- t(scale(t(chromvar_avg_motifactivity)))
chromvar_avg_motifactivity <- as.matrix(chromvar_avg_motifactivity)
# Load necessary library
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Create the heatmap
#col_fun <- colorRamp2(c(min(chromvar_avg_motifactivity), 0, max(chromvar_avg_motifactivity)), c("#4575b4", "white", "#d73027"))
col_fun <- colorRamp2(
  breaks = c(min(heatmap_data),  # lowest
             quantile(heatmap_data, 0.33),  # lower intermediate
             0,                           # mid point
             quantile(heatmap_data, 0.67),  # upper intermediate
             max(heatmap_data)),          # highest
  colors = c("#4575B4", "#91BFDB", "white", "#FCAE61", "#D73027")
)
BASE_FONTSIZE <- 14      
chromVAR_peaks <- Heatmap(
  chromvar_avg_motifactivity,
  col            = col_fun,
  na_col         = "grey",
  ## clustering & dendrogram placement
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  column_dend_side = "top",
  ## row / column labels & titles
  row_names_side    = "left",
  column_names_side = "bottom",
  show_row_names    = FALSE,
  column_names_gp   = gpar(fontsize = BASE_FONTSIZE,   fontface = "bold"),
  row_names_gp      = gpar(fontsize = BASE_FONTSIZE-2),
  column_title      = "ChromVAR Motif Activity Heatmap (All Peaks)",
  row_title         = "Transcription Factors",
  column_title_gp   = gpar(fontsize = BASE_FONTSIZE+4, fontface = "bold"),
  row_title_gp      = gpar(fontsize = BASE_FONTSIZE+4, fontface = "bold"),
  heatmap_legend_param = list(
    title      = "Z-Score Motif Activity",
    title_gp   = gpar(fontsize = BASE_FONTSIZE,   fontface = "bold"),
    labels_gp  = gpar(fontsize = BASE_FONTSIZE-1, fontface = "bold")
  )
)
chromVAR_peaks <- grid.grabExpr(draw(chromVAR_peaks))

## plot them all 4
# Wrap for patchwork compatibility
chromVAR_wrapped <- wrap_elements(full = chromVAR_peaks)
DA_wrapped <- wrap_elements(full = DA_peaks)
Fraction_wrapped <- wrap_elements(full = Fraction_peaks)
All_wrapped <- wrap_elements(full = All_peaks)

# Now combine all plots correctly:
final_combined_plot <- (
  (chromVAR_wrapped | All_wrapped) /
    (Fraction_wrapped | DA_wrapped)
) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

final_combined_plot

########### compare IMAGE and chromVAR
# find comment TFs and subset them based on the commen TFs
common_TFs <- intersect(rownames(chromvar_avg_motifactivity), rownames(heatmap_data))
chromvar_avg_motifactivity_subset <- chromvar_avg_motifactivity[common_TFs,]
heatmap_data_subset <- heatmap_data[common_TFs,]
# find correlation
cor_matrix <- cor(as.matrix(chromvar_avg_motifactivity_subset), as.matrix(heatmap_data_subset), method = "pearson")
library(ComplexHeatmap)
library(circlize)
library(grid)
col_fun <- colorRamp2(
  breaks = c(min(cor_matrix),  # lowest
             quantile(cor_matrix, 0.33),  # lower intermediate
             0,                           # mid point
             quantile(cor_matrix, 0.67),  # upper intermediate
             max(cor_matrix)),          # highest
  colors = c("#4575B4", "#91BFDB", "white", "#FCAE61", "#D73027")
)

BASE_FS <- 16
cor_ht <- Heatmap(
  cor_matrix,
  col              = col_fun,
  name             = "Pearson\nCorrelation",   # short legend title
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  # titles
  column_title     = "IMAGE Motif Activity",
  row_title        = "chromVAR Motif Activity",
  column_title_gp  = gpar(fontsize = BASE_FS + 1, fontface = "bold"),
  row_title_gp     = gpar(fontsize = BASE_FS + 1, fontface = "bold"),
  # axis labels
  column_names_gp  = gpar(fontsize = BASE_FS,     fontface = "bold"),
  row_names_gp     = gpar(fontsize = BASE_FS,     fontface = "bold"),
  # add the correlation numbers inside each cell
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cor_matrix[i, j]),
              x, y,
              gp = gpar(fontsize = BASE_FS, fontface = "bold"))
  },
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = BASE_FS,     fontface = "bold"),
    labels_gp = gpar(fontsize = BASE_FS - 1, fontface = "bold")
  )
)

draw(
  cor_ht,
  column_title = "Correlation Between chromVAR and IMAGE (Differential Accessible Regions) Motif Activities",
  column_title_gp = gpar(fontsize = BASE_FS + 2, fontface = "bold")
)
################ validation
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
metadata <- SHARE_seq_subsample_data@meta.data
metadata$TF <- sub("^[^-]*-", "", metadata$TF)
# as.factor to ensure compatability between 
metadata$TF <- as.factor(metadata$TF)
################ TF enrichment relative to GFP AND mCherry control
# Initialize an empty matrix to store the log2 ratios
# Filter the metadata into control and non-control
metadata_control <- metadata[metadata$TF %in% c("GFP", "mCherry"), ]
metadata_no_control <- metadata[!(metadata$TF %in% c("GFP", "mCherry")), ]
# Initialize a data.frame to store results
enrichment_results <- data.frame(
  TF = character(),
  Cluster = character(),
  Enrichment = numeric()
)
# Loop through each unique TF and cluster
unique_TFs <- unique(metadata_no_control$TF)
unique_clusters <- unique(metadata$seurat_clusters)

for (tf in unique_TFs) {
  for (cluster in unique_clusters) {
    # Cells expressing the TF in the cluster
    tf_cluster_count <- sum(metadata_no_control$TF == tf & metadata_no_control$seurat_clusters == cluster)
    
    # Cells expressing the control in the cluster
    control_cluster_count <- sum(metadata_control$seurat_clusters == cluster)
    
    # Cells expressing the TF across all clusters
    tf_total_count <- sum(metadata_no_control$TF == tf)
    
    # Cells expressing the control across all clusters
    control_total_count <- nrow(metadata_control)
    
    # Proportions
    cluster_proportion_tf <- tf_cluster_count / control_cluster_count
    overall_proportion_tf <- tf_total_count / control_total_count
    
    # Log2 enrichment
    enrichment <- log2((cluster_proportion_tf / overall_proportion_tf)+1) 
    
    # Store the result
    enrichment_results <- rbind(enrichment_results, data.frame(
      TF = tf,
      Cluster = cluster,
      Enrichment = enrichment
    ))
  }
}


# Convert results_df to a matrix
enrichment_results <- reshape2::acast(enrichment_results, TF ~ Cluster, value.var = "Enrichment", fill = 0)
enrichment_results <- enrichment_results[,order(as.numeric(colnames(enrichment_results)))]
colnames(enrichment_results) <- c("hESC 1", "hESC 2", "hESC 3",
                                  "Neuron like cells", "T/P like cells", "Goblet like cells", 
                                  "OP like cells", "VE like cells",
                                  "FL like cell", "Podocyte like cells", "SMP like cells", 
                                  "BE like cells", "Adipocyte like cells")

### Expression
mean_list <- list()
clusters <- unique(sort(metadata_no_control$seurat_clusters))
for (cluster in clusters){
  # Get barcodes (cells) for the cluster
  cluster_barcodes <- rownames(metadata_no_control[metadata_no_control$seurat_clusters == cluster, ])
  # Get only the cells that exist in chromVAR_object
  cells <- intersect(colnames(SHARE_seq_subsample_data), cluster_barcodes)
  
  if (length(cells) > 0) {  # Ensure there are cells
    # Compute mean deviations across cells for the cluster
    mean_list[[as.character(cluster)]] <- rowMeans(SHARE_seq_subsample_data@assays$RNA$data[,cells, drop = FALSE], na.rm = TRUE)
  } else {
    warning(paste("No matching cells found for cluster", cluster))
  }
}

Pseudobulk_exprs <- do.call(cbind, mean_list)
colnames(Pseudobulk_exprs) <- c("hESC 1", "hESC 2", "hESC 3",
                                "Neuron like cells", "T/P like cells", "Goblet like cells", 
                                "OP like cells", "VE like cells",
                                "FL like cell", "Podocyte like cells", "SMP like cells", 
                                "BE like cells", "Adipocyte like cells")

############# Compare for selected clusters
clusters = colnames(Pseudobulk_exprs[, 7:13])
enrichment_threshold = log2(2)
no_enrichment_treshold = log2(1)
expression_threshold = 0.2

for (cluster in clusters) {
  # Define enriched and not enriched TFs
  enriched_TFs = names(which(enrichment_results[, cluster] >= enrichment_threshold))
  not_enriched_TFs = names(which(enrichment_results[, cluster] <= no_enrichment_treshold))
  
  # Remove expressed TFs from not enriched TFs
  not_enriched_TFs = not_enriched_TFs[not_enriched_TFs %in% names(which(Pseudobulk_exprs[, cluster] <= expression_threshold))]
  
  # Convert to numeric vector instead of data frame
  IMAGE_activity = abs(as.numeric(heatmap_data[, cluster]))
  CV_activity = abs(as.numeric(chromvar_avg_motifactivity[, cluster]))
  
  # Save results
  res_tmp = data.frame(
    cluster = cluster,
    type = "Enriched",
    method = "IMAGE",
    activity_mean = mean(IMAGE_activity[rownames(heatmap_data) %in% enriched_TFs]),
    activity_median = median(IMAGE_activity[rownames(heatmap_data) %in% enriched_TFs])
  )
  
  res_tmp = rbind(res_tmp, data.frame(
    cluster = cluster,
    type = "Control",
    method = "IMAGE",
    activity_mean = mean(IMAGE_activity[rownames(heatmap_data) %in% not_enriched_TFs]),
    activity_median = median(IMAGE_activity[rownames(heatmap_data) %in% not_enriched_TFs])
  ))
  
  res_tmp = rbind(res_tmp, data.frame(
    cluster = cluster,
    type = "Enriched",
    method = "chromVAR",
    activity_mean = mean(CV_activity[rownames(chromvar_avg_motifactivity) %in% enriched_TFs]),
    activity_median = median(CV_activity[rownames(chromvar_avg_motifactivity) %in% enriched_TFs])
  ))
  
  res_tmp = rbind(res_tmp, data.frame(
    cluster = cluster,
    type = "Control",
    method = "chromVAR",
    activity_mean = mean(CV_activity[rownames(chromvar_avg_motifactivity) %in% not_enriched_TFs]),
    activity_median = median(CV_activity[rownames(chromvar_avg_motifactivity) %in% not_enriched_TFs])
  ))
  
  if (cluster == clusters[1]) {
    res = res_tmp
  } else {
    res = rbind(res, res_tmp)
  }
}

library(dplyr)
library(ggplot2)

# Summarize across TFs or replicate entries, etc.
df_by_cluster <- res %>%
  group_by(cluster, type, method) %>%
  summarize(
    median = activity_median,
    .groups = "drop"
  )

df_by_cluster$method <- factor(df_by_cluster$method, levels = c("IMAGE", "chromVAR"))
df_by_cluster$type <- factor(df_by_cluster$type, levels = c("Enriched", "Control"))

desired_cluster_order <- c("hESC 1", "hESC 2", "hESC 3",
                           "Neuron like cells", "T/P like cells", "Goblet like cells", 
                           "OP like cells", "VE like cells",
                           "FL like cell", "Podocyte like cells", "SMP like cells", 
                           "BE like cells", "Adipocyte like cells")

df_by_cluster$cluster <- factor(df_by_cluster$cluster, levels = desired_cluster_order)

gg1 <- ggplot(df_by_cluster, aes(x = method, y = median, fill = type)) +
  geom_col(position = position_dodge(width = 0.9)) +
  facet_wrap(~ cluster, ncol = 3) +
  labs(
    title = "Comparison of Enriched and Control TF Activities in Cell Types by Method",
    subtitle = paste0(
      "Thresholds:\n",
      "  • Enriched  – log2 ratio ≥ ", enrichment_threshold, " (≥2-fold)\n",
      "  • Control   – log2 ratio ≤ ", round(no_enrichment_treshold, 2),
      " (≤1-fold) and Expression ≤ ", expression_threshold
    ),
    x = NULL,
    y = "Median Activities"
  ) +
  scale_fill_manual(values = c("Enriched" = "#4575B4", "Control" = "#D73027")) +
  theme_bw(base_size = 14) +                
  theme(
    text           = element_text(face = "bold", colour = "black"),  # make *all* text bold
    plot.title     = element_text(size = 18, hjust = 0.5, face = "bold"),      # bump a bit more
    plot.subtitle  = element_text(size = 15, hjust = 0.5),
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 14),
    legend.text    = element_text(size = 14),
    legend.title   = element_text(size = 14),
    strip.text     = element_text(size = 14),      # facet labels
    legend.position = "top",
    axis.text.x  = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 14, face = "bold", colour = "black")
    
  )


############ median of median
clusters = colnames(Pseudobulk_exprs[, 7:13])
enrichment_threshold = log2(2)
no_enrichment_treshold = log2(1)
expression_threshold = 0.2

for (cluster in clusters) {
  # Define enriched and not enriched TFs
  enriched_TFs = names(which(enrichment_results[, cluster] >= enrichment_threshold))
  not_enriched_TFs = names(which(enrichment_results[, cluster] <= no_enrichment_treshold))
  
  # Remove expressed TFs from not enriched TFs
  not_enriched_TFs = not_enriched_TFs[not_enriched_TFs %in% names(which(Pseudobulk_exprs[, cluster] <= expression_threshold))]
  
  # Convert to numeric vector instead of data frame
  IMAGE_activity = abs(as.numeric(heatmap_data[, cluster]))
  CV_activity = abs(as.numeric(chromvar_avg_motifactivity[, cluster]))
  
  # Save results
  res_tmp = data.frame(
    cluster = cluster,
    type = "Enriched",
    method = "IMAGE",
    activity_mean = mean(IMAGE_activity[rownames(heatmap_data) %in% enriched_TFs]),
    activity_median = median(IMAGE_activity[rownames(heatmap_data) %in% enriched_TFs])
  )
  
  res_tmp = rbind(res_tmp, data.frame(
    cluster = cluster,
    type = "Control",
    method = "IMAGE",
    activity_mean = mean(IMAGE_activity[rownames(heatmap_data) %in% not_enriched_TFs]),
    activity_median = median(IMAGE_activity[rownames(heatmap_data) %in% not_enriched_TFs])
  ))
  
  res_tmp = rbind(res_tmp, data.frame(
    cluster = cluster,
    type = "Enriched",
    method = "chromVAR",
    activity_mean = mean(CV_activity[rownames(chromvar_avg_motifactivity) %in% enriched_TFs]),
    activity_median = median(CV_activity[rownames(chromvar_avg_motifactivity) %in% enriched_TFs])
  ))
  
  res_tmp = rbind(res_tmp, data.frame(
    cluster = cluster,
    type = "Control",
    method = "chromVAR",
    activity_mean = mean(CV_activity[rownames(chromvar_avg_motifactivity) %in% not_enriched_TFs]),
    activity_median = median(CV_activity[rownames(chromvar_avg_motifactivity) %in% not_enriched_TFs])
  ))
  
  if (cluster == clusters[1]) {
    res = res_tmp
  } else {
    res = rbind(res, res_tmp)
  }
}

library(dplyr)

df_summary <- res %>%
  group_by(type, method) %>%
  summarize(
    mean_of_median = median(activity_median, na.rm = TRUE),
    .groups = "drop"
  )
# Ensuring method is a factor: IMAGE first, then chromVAR
df_summary$method <- factor(df_summary$method, levels = c("IMAGE", "chromVAR"))

# Make type a factor: Enriched first, then Depleted
df_summary$type <- factor(df_summary$type, levels = c("Enriched", "Control"))

library(ggplot2)
# Final ggplot
gg2 <- ggplot(df_summary, aes(x = method, y = mean_of_median, fill = type)) +
  geom_col(position = position_dodge(width = 0.9)) +
  labs(
    title = "Comparison of Enriched and Control TF Activities in Cell Types by Method",
    subtitle = paste0(
      "Thresholds:\n",
      "  • Enriched  – log2 ratio ≥ ", enrichment_threshold, " (≥2-fold)\n",
      "  • Control   – log2 ratio ≤ ", round(no_enrichment_treshold, 2),
      " (≤1-fold) and Expression ≤ ", expression_threshold
    ),
    x = NULL,
    y = "Median of Median Activities"
  ) +
  scale_fill_manual(values = c("Enriched" = "#4575B4", "Control" = "#D73027")) +
  theme_bw(base_size = 14) +                # base font size ↑
  theme(
    text           = element_text(face = "bold", colour = "black"),  # make *all* text bold
    plot.title     = element_text(size = 18, hjust = 0.5, face = "bold"),      # bump a bit more
    plot.subtitle  = element_text(size = 15, hjust = 0.5),
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 14),
    legend.text    = element_text(size = 14),
    legend.title   = element_text(size = 14),
    strip.text     = element_text(size = 14),      # facet labels
    legend.position = "top",
    axis.text.x  = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 14, face = "bold", colour = "black")
  )

############### line plot
# test enriched thresholds = 2:9 fold
all_folds <- 2:9
expr_threshold <- 0.2
depl_threshold <- 1  # 1-fold => log2(1) = 0   , 2-fold => log2(2) = 1

df_all <- list()
k <- 1
# loop through the fold changes
for (fold in all_folds) {
  # calculate the log2 of the fold
  enrT   <- log2(fold)  # log2 of e.g. 2..9
  noEnrT <- log2(depl_threshold)   # always 1 => 0
  
  for (cluster in clusters) {
    
    # Identify enriched TFs
    enriched_TFs <- names(which(enrichment_results[, cluster] >= enrT))
    
    # Identify not-enriched TFs
    not_enriched_TFs <- names(which(enrichment_results[, cluster] <= noEnrT))
    # Also require that they are low-expressed (<= expr_threshold)
    not_enriched_TFs <- intersect(
      not_enriched_TFs,
      names(which(Pseudobulk_exprs[, cluster] <= expr_threshold))
    )
    
    # Activities for this cluster
    IMAGE_activity <- abs(as.numeric(heatmap_data[, cluster]))
    CV_activity    <- abs(as.numeric(chromvar_avg_motifactivity[, cluster]))
    
    # Identify only the motifs that actually appear in each matrix:
    valid_enriched_IMAGE    <- intersect(enriched_TFs, rownames(heatmap_data))
    valid_control_IMAGE     <- intersect(not_enriched_TFs, rownames(heatmap_data))
    valid_enriched_chromVAR <- intersect(enriched_TFs, rownames(chromvar_avg_motifactivity))
    valid_control_chromVAR  <- intersect(not_enriched_TFs, rownames(chromvar_avg_motifactivity))
    
    # Build rows for (IMAGE.Enriched, IMAGE.Control, chromVAR.Enriched, chromVAR.Control)
    res_tmp <- data.frame(
      cluster         = cluster,
      type            = "Enriched",
      method          = "IMAGE",
      activity_median = median(IMAGE_activity[rownames(heatmap_data) %in% enriched_TFs]),
      n_TFs           = length(valid_enriched_IMAGE)
    )
    
    res_tmp <- rbind(res_tmp, data.frame(
      cluster         = cluster,
      type            = "Control",
      method          = "IMAGE",
      activity_median = median(IMAGE_activity[rownames(heatmap_data) %in% not_enriched_TFs]),
      n_TFs           = length(valid_control_IMAGE)
    ))
    
    res_tmp <- rbind(res_tmp, data.frame(
      cluster         = cluster,
      type            = "Enriched",
      method          = "chromVAR",
      activity_median = median(CV_activity[rownames(chromvar_avg_motifactivity) %in% enriched_TFs]),
      n_TFs           = length(valid_enriched_chromVAR)
    ))
    
    res_tmp <- rbind(res_tmp, data.frame(
      cluster         = cluster,
      type            = "Control",
      method          = "chromVAR",
      activity_median = median(CV_activity[rownames(chromvar_avg_motifactivity) %in% not_enriched_TFs]),
      n_TFs           = length(valid_control_chromVAR)
    ))
    
    if (cluster == clusters[1]) {
      res <- res_tmp
    } else {
      res <- rbind(res, res_tmp)
    }
  }
  
  # Summarize across all clusters for this fold
  df_summary <- res %>%
    dplyr::group_by(type, method) %>%
    dplyr::summarize(
      median_of_median = median(activity_median, na.rm = TRUE),
      total_TFs        = sum(n_TFs),
      .groups          = "drop"
    )
  
  # Record which fold was used
  df_summary$enriched_fold <- fold
  
  # Store in our master list
  df_all[[k]] <- df_summary
  k <- k + 1
}

# Combine into a single data frame
df_all <- do.call(rbind, df_all)

# factorising
df_all$method <- factor(df_all$method, levels = c("IMAGE", "chromVAR"))
df_all$type   <- factor(df_all$type,   levels = c("Enriched", "Control"))

# Now plot
library(ggplot2)

p <- ggplot(df_all, aes(
  x     = enriched_fold, 
  y     = median_of_median,
  color = interaction(method, type),
  group = interaction(method, type)
)) +
  # LOESS smoothing
  geom_smooth(method = "loess", se = FALSE, span = 0.5) + 
  # Keep text labels for the points 2..9
  geom_text(
    data = subset(df_all, enriched_fold %in% 2:9),
    aes(label = paste0("n=", total_TFs)),
    vjust = -1,
    size  = 5,
    show.legend = FALSE,
    colour = "black",
    fontface = "bold") +
  scale_x_continuous(breaks = 2:9) +  # ensure ticks from 2 to 9
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Median of Median TF Enhancer Activities vs. Enrichment Fold",
    subtitle = "Control: ≤ 1-fold & expr ≤ 0.2; Enriched: ≥ 2 to 9-fold",
    x        = "Enrichment Fold (≥)",
    y        = "Median of Median Activities",
    color    = "Method & Type"
  ) +
  theme_bw(base_size = 14) +              
  theme(
    text           = element_text(face = "bold", colour = "black"),  # make *all* text bold
    plot.title     = element_text(size = 18, hjust = 0.5, face = "bold"),      # bump a bit more
    plot.subtitle  = element_text(size = 15, hjust = 0.5),
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 14),
    legend.text    = element_text(size = 14),
    legend.title   = element_text(size = 14),
    strip.text     = element_text(size = 14),      # facet labels
    legend.position = "right",
    axis.text.x  = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 14, face = "bold", colour = "black")
  )



############## boxplots of it
library(dplyr)
library(ggplot2)

# test enriched thresholds = 2 to 9
all_folds <- 2:9
expr_threshold <- 0.2
depl_threshold <- 1

# keep cluster-level results in df_box
df_box <- data.frame()

for (fold in all_folds) {
  
  enrT   <- log2(fold)               
  noEnrT <- log2(depl_threshold)     
  
  # Loop over clusters
  for (cluster in clusters) {
    
    # Define sets of TFs
    enriched_TFs     <- names(which(enrichment_results[, cluster] >= enrT))
    not_enriched_TFs <- names(which(enrichment_results[, cluster] <= noEnrT))
    
    # Constrain the "control" TFs to be low-expressed (≤ expr_threshold)
    not_enriched_TFs <- intersect(
      not_enriched_TFs,
      names(which(Pseudobulk_exprs[, cluster] <= expr_threshold))
    )
    
    # Activities
    IMAGE_activity <- abs(as.numeric(heatmap_data[, cluster]))
    CV_activity    <- abs(as.numeric(chromvar_avg_motifactivity[, cluster]))
    
    # For each cluster, create a small data.frame with 4 rows:
    # (IMAGE, Enriched), (IMAGE, Control), (chromVAR, Enriched), (chromVAR, Control)
    res_tmp <- data.frame(
      cluster = cluster,
      method  = c("IMAGE","IMAGE","chromVAR","chromVAR"),
      type    = c("Enriched","Control","Enriched","Control"),
      activity_median = c(
        median(IMAGE_activity[rownames(heatmap_data) %in% enriched_TFs]),
        median(IMAGE_activity[rownames(heatmap_data) %in% not_enriched_TFs]),
        median(CV_activity[rownames(chromvar_avg_motifactivity) %in% enriched_TFs]),
        median(CV_activity[rownames(chromvar_avg_motifactivity) %in% not_enriched_TFs])
      )
    )
    
    # Add the fold value
    res_tmp$enriched_fold <- fold
    
    # Append to our big data frame
    df_box <- rbind(df_box, res_tmp)
  }
}

# Ensuring method and type are factored 
df_box$method <- factor(df_box$method, levels = c("IMAGE", "chromVAR"))
df_box$type   <- factor(df_box$type,   levels = c("Enriched", "Control"))

# Create a combined factor that explicitly fixes the order
df_box$combo <- interaction(df_box$method, df_box$type)

# Re-level the combined factor so it appears in desired order.
df_box$combo <- factor(
  df_box$combo,
  levels = c("IMAGE.Enriched", "IMAGE.Control",
             "chromVAR.Enriched", "chromVAR.Control")
)

# plot
p_box <- ggplot(df_box, 
                aes(
                  x = factor(enriched_fold), 
                  y = activity_median,
                  fill = combo
                )) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  labs(
    title    = "Distribution of Median TF Enhancer Motif Activities vs. Enrichment Fold",
    subtitle = "Control: ≤ 1-fold and Expression ≤ 0.2, Enriched: ≥ 2 to 9-fold",
    x        = "Enrichment Fold (≥)",
    y        = "Median Activity",
    fill     = "Method & Type"
  ) +
  theme_bw(base_size = 14) +                # base font size ↑
  theme(
    text           = element_text(face = "bold", colour = "black"),  # make *all* text bold
    plot.title     = element_text(size = 18, hjust = 0.5, face = "bold"),      # bump a bit more
    plot.subtitle  = element_text(size = 15, hjust = 0.5),
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 14),
    legend.text    = element_text(size = 14),
    legend.title   = element_text(size = 14),
    strip.text     = element_text(size = 14),      # facet labels
    legend.position = "right",
    axis.text.x  = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 14, face = "bold", colour = "black")
  )

library(patchwork)

combined_plot <- (gg1 + gg2) / (p + p_box) +
  plot_annotation(tag_levels = "A") &  
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

combined_plot

############
# Step 2 of IMAGE analysis
###########
# prepare RNA summerized experiment data
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
SHARE_seq_subsample_data = SetIdent(SHARE_seq_subsample_data, value = SHARE_seq_subsample_data$seurat_clusters)
new.cluster.ids <- c("hESC 1", "hESC 2", "hESC 3",
                     "Neuron like cells", "T/P like cells", "Goblet like cells", 
                     "OP like cells", "VE like cells",
                     "FL like cell", "Podocyte like cells", "SMP like cells", 
                     "BE like cells", "Adipocyte like cells")
names(new.cluster.ids) <- levels(SHARE_seq_subsample_data)
SHARE_seq_subsample_data <- RenameIdents(SHARE_seq_subsample_data, new.cluster.ids)
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")
metadata <- OE@meta.data

# Divide batches into 4 pseudobulks
library(dplyr)
# Add a new column based on batch ranges
metadata <- metadata %>%
  mutate(batch_group = case_when(
    batch %in% c(0, 1, 2, 3) ~ 1,
    batch %in% c(4, 5, 6, 7) ~ 2,
    batch %in% c(8, 9, 10, 11) ~ 3,
    batch %in% c(12, 13, 14, 15) ~ 4
  ))
# Add the new batch group column to the Seurat object
OE$batch_group <- metadata$batch_group

# aggregate
OE_PB2 = AggregateExpression(OE, group.by = c("ident", "batch_group"),
                             normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             return.seurat = TRUE)
colnames_updated <- make.unique(gsub("_.*", "", colnames(OE_PB2)))
colnames(OE_PB2) <- colnames_updated

# Assign Group Labels Based on Identifiers:
group_labels <- sub("\\..*", "", colnames(OE_PB2)) # Removes ".x" suffix
Idents(OE_PB2) <- group_labels

# extract RNA counts
bulk_RNA_data <- OE_PB2@assays$RNA$counts

library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(dplyr)
# establish connection to Ensembl BioMart server and get genes for humans
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
att <- listAttributes(mart)
grep("transcript", att$name, value=TRUE)

# Retrieve Gene Annotation Information
gene_info <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "transcript_count", "chromosome_name", "external_gene_name", "transcript_is_canonical", 
                 "transcript_start", "transcription_start_site", "transcript_end", "transcript_length", "strand"),
  filters = "hgnc_symbol",
  values = rownames(SHARE_seq_subsample_data@assays$RNA),
  mart = mart
)
# remove na or incomplete rows
gene_cleaned <- gene_info[complete.cases(gene_info), ]

# Define standard chromosomes
standard_chromosomes <- c(as.character(1:22), "X", "Y")

# Keep only rows where chromosome_name is in the standard set
gene_standard <- gene_cleaned %>%
  filter(chromosome_name %in% standard_chromosomes)

# Within each gene, keep the row with the highest transcript_count
gene_filtered <- gene_standard %>%
  group_by(external_gene_name) %>%
  arrange(desc(transcript_count)) %>%
  slice(1) %>%  # Keep only the best row per gene
  ungroup()
# convert 1 to + and -1 to -
gene_filtered <- gene_filtered %>%
  mutate(strand = ifelse(strand == 1, "+", "-"))
# add chr
gene_filtered$chromosome_name <- paste0("chr", gene_filtered$chromosome_name)


# Extract the matching rows from annotation_unique based on Name column matching rownames of RNA_summed_sparse
matches <- gene_filtered[gene_filtered$external_gene_name %in% rownames(bulk_RNA_data), ]

# Select only columns of interests
extracted_info <- matches[, c("external_gene_name","chromosome_name", "strand", "transcript_start", "transcript_end", "transcript_length", "ensembl_gene_id")]
colnames(extracted_info) <- c("Name","chr", "strand", "start", "end", "Length", "Ensembl")
# subset so I get the genes names of my data that is found in this database of the extracted genes
bulk_RNA_data_filtered <- bulk_RNA_data[rownames(bulk_RNA_data) %in% extracted_info$Name, ]

# Reorder RNA_summed_subset based on the order of the Name column in extracted_info
bulk_RNA_data_filtered <- bulk_RNA_data_filtered[match(extracted_info$Name, rownames(bulk_RNA_data_filtered)), ]

# Check again if the order now matches
order_matches <- all(rownames(bulk_RNA_data_filtered) == extracted_info$Name)

# Print the result
if (order_matches) {
  print("The order now matches between RNA_summed_subset rownames and the Name column of extracted_info.")
} else {
  print("The order still does not match.")
}
# Define the limits, each range must have a start that is < 2^31 and > - 2^31 in the summerized experiment object
lower_limit <- -2^31
upper_limit <- 2^31-1
valid_ranges <- extracted_info$start < upper_limit & extracted_info$start > lower_limit &
  extracted_info$end < upper_limit & extracted_info$end > lower_limit

extracted_info_filtered <- extracted_info[valid_ranges, ]
bulk_RNA_data_filtered <- bulk_RNA_data_filtered[valid_ranges, ]


library(SummarizedExperiment)
cond = data.frame(names = colnames(bulk_RNA_data_filtered))
cond$condition <- paste0("Cond", seq_len(nrow(cond)))
# Create a vector of unique conditions based on unique cell types
unique_cell_types <- unique(sub("\\.\\d+$", "", cond$names))
conditions <- paste0("Cond", seq_along(unique_cell_types))
# Create a mapping between cell types and conditions
cell_type_to_condition <- setNames(conditions, unique_cell_types)
# Assign conditions based on cell type
cond$condition <- cell_type_to_condition[sub("\\.\\d+$", "", cond$names)]
rownames(cond) = cond$names
cond = DataFrame(cond)

# Create a matrix with count lengths
lengths = matrix(rep(extracted_info_filtered$Length, 52),ncol=52)
colnames(lengths) = colnames(bulk_RNA_data_filtered)
bulk_RNA_data_filtered <- as.matrix(bulk_RNA_data_filtered)
RNA = SummarizedExperiment(assays = list(counts = bulk_RNA_data_filtered, length = lengths),
                           colData = cond,
                           rowRanges = GRanges(seqnames = extracted_info_filtered$chr, strand = extracted_info_filtered$strand, IRanges(start = extracted_info_filtered$start, end = extracted_info_filtered$end), gene_id = extracted_info_filtered$Ensembl, symbol = extracted_info_filtered$Name))

#saveRDS(RNA, "/work/Home/Data_for_masters_project/RNA_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_4_replicates_all peaks")

###################### Run the rest of IMAGE
RNA <- readRDS("/work/Home/Data_for_masters_project/RNA_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_4_replicates_all peaks")
ATAC <- readRDS("/work/Home/Data_for_masters_project/ATAC_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_4_replicates_all peaks")

# Assign the genome variable
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

## Make design vector
# Extract unique cell types from the column names (before the period)
unique_cell_types <- unique(sub("\\..*", "", colnames(RNA)))

# Create a mapping of cell types to numeric values
cell_type_mapping <- setNames(seq_along(unique_cell_types), unique_cell_types)

# Create the numeric vector by mapping the column names to the numeric values
design <- as.numeric(as.character(cell_type_mapping[sub("\\..*", "", colnames(RNA))]))
# load IMAGE functions
analyzeRNA = function(RNA, design, paired = FALSE) {
  ## Get dimensions of samples and conditions
  Samples <- length(design)
  Conditions <- length(unique(design))
  
  # Extract data
  counts_mat <- as.matrix(SummarizedExperiment::assays(RNA, withDimnames = FALSE)$counts)
  
  # calculate GC frquency
  gc_fraction <- letterFrequency(getSeq(genome, granges(RNA)), 
                                 letters="GC", as.prob=TRUE)
  gc_fraction <- as.numeric(gc_fraction)
  # Extract peak lengths
  peak_length <- width(granges(RNA))
  
  # Create SeqExpressionSet
  featureData = data.frame(GC = gc_fraction, Length = peak_length)
  rownames(featureData) <- rownames(counts_mat)
  seq_obj <- newSeqExpressionSet(counts = counts_mat,
                                 featureData = featureData)
  
  seq_obj_gc <- withinLaneNormalization(seq_obj, "GC", which = "full")
  
  seq_obj_gc_size <- betweenLaneNormalization(seq_obj_gc, which = "full")
  
  normalized_counts_fqfq <- normCounts(seq_obj_gc_size)
  
  
  ########### Calculate CPM, filter genes and calculate differential expression
  # Create a DGEList object using the raw counts
  GeneDGE <- edgeR::DGEList(counts = normalized_counts_fqfq)
  
  # Calculate normalization factors for the DGEList object
  # This step adjusts for differences in library sizes among samples
  GeneDGE <- edgeR::calcNormFactors(GeneDGE)
  
  # Compute counts per million (CPM) for each gene, then log-transform the CPM values
  GeneCPM <- as.data.frame(edgeR::cpm(GeneDGE, log = TRUE, normalized.lib.sizes = TRUE))
  
  # Rename the columns of the GeneCPM data frame to include both the sample number and its condition
  for (samp in 1:Samples) {
    colnames(GeneCPM)[samp] <- paste("Expression-Sample_", samp, "-Condition_", design[samp], sep = "")
  }
  
  # For each gene, calculate the maximum logCPM value across all samples
  GeneCPM$Max <- apply(GeneCPM[, 1:ncol(GeneCPM)], 1, FUN = "max")
  
  # Add column to the GeneCPM data frame containing the gene names
  GeneCPM$Factor <- rownames(GeneCPM)
  
  ####### Calculate a logCPM threshold corresponding to a raw count of 10. According to edgeR
  # Identify the smallest library size among all samples
  min_lib_size = min(GeneDGE$samples$lib.size)
  
  # Convert the raw count threshold (10 counts) to CPM by dividing by the library size (in millions)
  # applying a log2 transformation with a pseudo-count of 1.
  threshold = log2(10 / (min_lib_size * 1e-6) + 1)
  
  # filter genes based on the following criteria
  if (!is.null(threshold)) {
    # At least two samples have logCPM values that meet or exceed the threshold.
    # The maximum logCPM value for the gene (across all samples) is at or above the threshold.
    keep <- (rowSums(GeneCPM[, 1:Samples] >= threshold) >= 2) & (GeneCPM$Max >= threshold)
    
    # Subset the GeneCPM data frame to keep only the genes meeting the above criteria.
    GeneCPM <- GeneCPM[keep, ]
  }
  
  # calculate differential expression for different number of condition
  # If there are exactly 2 conditions:
  if (Conditions == 2) {
    # Estimate the common dispersion for the DGEList object.
    GeneDGE <- edgeR::estimateCommonDisp(GeneDGE)
    
    # Estimate the trended dispersion (captures the trend in dispersion with gene abundance).
    GeneDGE <- edgeR::estimateTrendedDisp(GeneDGE)
    
    # Estimate tagwise (gene-specific) dispersion.
    GeneDGE <- edgeR::estimateTagwiseDisp(GeneDGE)
    
    # Perform the exact test for differential expression between the two conditions.
    Test <- as.data.frame(edgeR::topTags(edgeR::exactTest(GeneDGE), nrow(RNA)))
    
    # Save gene identifiers into a new column "Factor" based on the row names.
    Test$Factor <- rownames(Test)
    
    # Rename columns for clarity:
    # Column 1 becomes the log-fold change between condition 2 and condition 1.
    colnames(Test)[1] <- "logFC_Cond2-vs-Cond1"
    # Column 3 is the p-value for the test.
    colnames(Test)[3] <- "Pval_Cond2-vs-Cond1"
    # Column 4 is the false discovery rate (FDR).
    colnames(Test)[4] <- "FDR_Cond2-vs-Cond1"
    
    # Merge the CPM table (GeneCPM) with the differential expression results on "Factor".
    GeneCPM <- merge(GeneCPM, Test[, c(1, 3, 4, 5)], by = "Factor")
  }
  
  
  # If there are three or more conditions:
  if (Conditions >= 3) {
    # The workflow diverges based on whether the experimental design is paired or unpaired.
    
    # For unpaired samples:
    if (!paired) {
      # Construct the design matrix using the 'design' factor.
      DesignMatrix <- model.matrix(~factor(design))
      
      # Estimate dispersion
      GeneDGE <- edgeR::estimateGLMCommonDisp(GeneDGE, DesignMatrix)
      GeneDGE <- edgeR::estimateGLMTrendedDisp(GeneDGE, DesignMatrix)
      GeneDGE <- edgeR::estimateGLMTagwiseDisp(GeneDGE, DesignMatrix)
      
      # Fit negative binomial model to the data.
      GeneFit <- edgeR::glmQLFit(GeneDGE, DesignMatrix)
      
      # Generate a grid of all possible pairwise comparisons using the 'design' factor.
      Tests <- expand.grid(design, design)
      # Remove comparisons where a condition is compared with itself.
      Tests <- Tests[Tests[, 1] != Tests[, 2], ]
      # Create a combined label to identify unique comparison pairs.
      Tests$Combine <- paste(Tests[, 1], Tests[, 2], sep = "-")
      # Remove duplicate comparisons.
      Tests <- Tests[duplicated(Tests$Combine) == F, c(1, 2)]
      # Keep only comparisons where the first condition is less than the second,
      # to establish a consistent contrast direction.
      Tests <- Tests[Tests[, 1] < Tests[, 2], ]
      # Order the comparisons by condition values.
      Tests <- Tests[order(Tests[, 1], Tests[, 2]), ]
      
      # Identify comparisons of condition vs. condition 1.
      Coefficients <- Tests[Tests[, 1] == 1, 2]
      # Identify remaining pairwise contrasts.
      Contrasts <- Tests[Tests[, 1] > 1, ]
      
      # Loop over pairwise contrasts that are not comparisons against condition 1.
      for (i in 1:nrow(Contrasts)) {
        # Create a contrast vector that assigns -1 to the first condition and 1 to the second, zeros otherwise.
        Contrast <- as.vector(rep(0, Conditions))
        Contrast[Contrasts[i, 1]] <- -1
        Contrast[Contrasts[i, 2]] <- 1
        
        # Perform a quasi-likelihood F-test using the defined contrast.
        Test <- as.data.frame(edgeR::topTags(edgeR::glmQLFTest(GeneFit, contrast = Contrast), nrow(RNA)))
        # Store gene identifiers as a new column "Factor".
        Test$Factor <- rownames(Test)
        
        # Rename columns to reflect the contrast between the two conditions.
        colnames(Test)[1] <- paste("logFC_Cond", Contrasts[i, 2], "-vs-Cond", Contrasts[i, 1], sep = "")
        colnames(Test)[4] <- paste("Pval_Cond", Contrasts[i, 2], "-vs-Cond", Contrasts[i, 1], sep = "")
        colnames(Test)[5] <- paste("FDR_Cond", Contrasts[i, 2], "-vs-Cond", Contrasts[i, 1], sep = "")
        
        # Merge the results for this contrast into the GeneCPM table.
        GeneCPM <- merge(GeneCPM, Test[, c(1, 4, 5, 6)], by = "Factor")
      }
      
      # Loop over coefficients comparing conditions against condition 1.
      for (i in Coefficients) {
        # Perform a quasi-likelihood F-test on the specified coefficient.
        Test <- as.data.frame(edgeR::topTags(edgeR::glmQLFTest(GeneFit, coef = i), nrow(RNA)))
        Test$Factor <- rownames(Test)
        
        # Rename columns for the contrast of condition i versus condition 1.
        colnames(Test)[1] <- paste("logFC_Cond", i, "-vs-Cond1", sep = "")
        colnames(Test)[4] <- paste("Pval_Cond", i, "-vs-Cond1", sep = "")
        colnames(Test)[5] <- paste("FDR_Cond", i, "-vs-Cond1", sep = "")
        
        # Merge these results into the GeneCPM table.
        GeneCPM <- merge(GeneCPM, Test[, c(1, 4, 5, 6)], by = "Factor")
      }
    }
    
    # For paired samples:
    if (paired) {
      # Create a frequency table of the 'design' factor.
      Table <- table(design)
      # Initialize an empty vector to store replicate identifiers.
      Paired <- vector()
      # For each condition, create a sequence of replicate indices.
      for (q in 1:length(Table)) {
        Paired <- c(Paired, seq(1, Table[q], by = 1))
      }
      
      # Construct the design matrix with factors for pairing and condition effects.
      DesignMatrix <- model.matrix(~factor(Paired) + factor(design))
      
      # Estimate dispersion
      GeneDGE <- edgeR::estimateGLMCommonDisp(GeneDGE, DesignMatrix)
      GeneDGE <- edgeR::estimateGLMTrendedDisp(GeneDGE, DesignMatrix)
      GeneDGE <- edgeR::estimateGLMTagwiseDisp(GeneDGE, DesignMatrix)
      
      # Fit the model using quasi-likelihood methods.
      GeneFit <- edgeR::glmQLFit(GeneDGE, DesignMatrix)
      
      # Generate pairwise comparisons between conditions.
      Tests <- expand.grid(design, design)
      Tests <- Tests[Tests[, 1] != Tests[, 2], ]
      Tests$Combine <- paste(Tests[, 1], Tests[, 2], sep = "-")
      Tests <- Tests[duplicated(Tests$Combine) == F, c(1, 2)]
      Tests <- Tests[Tests[, 1] < Tests[, 2], ]
      Tests <- Tests[order(Tests[, 1], Tests[, 2]), ]
      
      # Determine the maximum number of replicates across conditions.
      MaxNumberOfReplicates <- max(table(design))
      
      # Adjust coefficients for comparisons involving condition 1 by shifting the index.
      Coefficients <- Tests[Tests[, 1] == 1, 2] + (MaxNumberOfReplicates - 1)
      # For the remaining comparisons, adjust the indices similarly.
      Contrasts <- Tests[Tests[, 1] > 1, ]
      Contrasts[, 3] <- Contrasts[, 1] + (MaxNumberOfReplicates - 1)
      Contrasts[, 4] <- Contrasts[, 2] + (MaxNumberOfReplicates - 1)
      
      # Loop over pairwise contrasts for the paired design.
      for (i in 1:nrow(Contrasts)) {
        # Create a contrast vector of the appropriate length, considering replicate adjustments.
        Contrast <- as.vector(rep(0, (Conditions + (MaxNumberOfReplicates - 1))))
        Contrast[Contrasts[i, 3]] <- -1
        Contrast[Contrasts[i, 4]] <- 1
        
        # Test the specified contrast using a quasi-likelihood F-test.
        Test <- as.data.frame(edgeR::topTags(edgeR::glmQLFTest(GeneFit, contrast = Contrast), nrow(RNA)))
        Test$Factor <- rownames(Test)
        
        # Rename columns to detail the comparison.
        colnames(Test)[1] <- paste("logFC_Cond", Contrasts[i, 2], "-vs-Cond", Contrasts[i, 1], sep = "")
        colnames(Test)[4] <- paste("Pval_Cond", Contrasts[i, 2], "-vs-Cond", Contrasts[i, 1], sep = "")
        colnames(Test)[5] <- paste("FDR_Cond", Contrasts[i, 2], "-vs-Cond", Contrasts[i, 1], sep = "")
        
        # Merge these differential expression results with the CPM table.
        GeneCPM <- merge(GeneCPM, Test[, c(1, 4, 5, 6)], by = "Factor")
      }
      
      # Loop over coefficients comparing each condition versus condition 1 for paired data.
      for (i in Coefficients) {
        Test <- as.data.frame(edgeR::topTags(edgeR::glmQLFTest(GeneFit, coef = i), nrow(RNA)))
        Test$Factor <- rownames(Test)
        
        # Rename columns for clarity.
        colnames(Test)[1] <- paste("logFC_Cond", i, "-vs-Cond1", sep = "")
        colnames(Test)[4] <- paste("Pval_Cond", i, "-vs-Cond1", sep = "")
        colnames(Test)[5] <- paste("FDR_Cond", i, "-vs-Cond1", sep = "")
        
        # Merge the results into the main CPM table.
        GeneCPM <- merge(GeneCPM, Test[, c(1, 4, 5, 6)], by = "Factor")
      }
    }
  }
  
  ## Insert TSS for each gene
  ranges.df = data.frame(
    Factor = rownames(RNA),
    seqnames = GenomicRanges::seqnames(RNA),
    start = GenomicRanges::start(RNA),
    end = GenomicRanges:: end(RNA),
    strand = GenomicRanges::strand(RNA))
  ranges.df$TSS = ifelse(ranges.df$strand == "+", ranges.df$start, ranges.df$end)
  colnames(ranges.df)[2] = "Chr"
  GeneCPM = merge(GeneCPM, ranges.df[,c("Factor","Chr","TSS")], by="Factor")
  # Reorder columns to ensure "Max", "Chr", and "TSS" are not at the back
  desired_order = c("Factor", "Max", "Chr", "TSS", setdiff(colnames(GeneCPM), c("Factor", "Max", "Chr", "TSS")))
  GeneCPM = GeneCPM[, desired_order]
  
  # Return
  return(GeneCPM)
}
findTargetEnhancers = function(motifs, enhancers, activity, raw_motifs, ncores) {
  # find common enhancers
  common_peaks <- intersect(rownames(enhancers), rownames(motifs))
  
  # subset the objects to have common peaks
  motifs <- motifs[common_peaks,]
  enhancers <- enhancers[common_peaks,]
  raw_motifs <- raw_motifs[common_peaks,]
  
  # setup parallel
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  # Set library paths on worker nodes
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load libraries on worker nodes
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
  })
  # Calculate prediction error
  # Setup motifs
  motifs.error = t(motifs)
  motifs.error = motifs.error[ match(rownames(activity), rownames(motifs.error)),]
  # Parallel loop for samples
  FullModel <- foreach(i = 1:(ncol(activity)), .combine = cbind, .packages = c("Matrix")) %dopar% {
    # Calculate predicted signal
    AN = motifs.error * activity[, i]
    AN = colSums(AN)
    
    # Extract observed signal
    response = enhancers[, colnames(enhancers) == colnames(activity)[i]]
    names(response) = rownames(enhancers)
    response = response[match(names(AN), names(response))]
    
    # Calculate error for this sample
    (response - AN)^2
  }
  
  # Set column and row names
  colnames(FullModel) = colnames(activity)
  rownames(FullModel) = rownames(enhancers)
  
  ## Calculate average standard deviation
  Average <- sum(rowSums(FullModel)) / (nrow(FullModel) * ncol(FullModel))
  
  ## Parallelized loop for motifs
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  # Set library paths on worker nodes
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load libraries on worker nodes
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
  })
  score = foreach(motif = colnames(motifs), .combine = cbind, .packages = c("Matrix")) %dopar% {
    # Define object to hold errors
    RedModel <- matrix(nrow=nrow(motifs), ncol=(ncol(activity)))
    
    # Loop across samples
    for (i in 1:(ncol(activity))) {
      # Calculate predicted signal without the target motif
      AN = motifs.error[!(rownames(motifs.error) %in% motif),] * activity[!(rownames(activity) %in% motif),][,i]
      AN = colSums(AN)
      
      # Extract observed signal
      response = enhancers[,colnames(enhancers) == colnames(activity)[i]]
      names(response) = rownames(enhancers)
      response = response[ match(names(AN), names(response))]
      
      # Calculate error
      RedModel[,i] = (response - AN)^2
    }
    
    # Calculate score
    score_result <-rowSums(RedModel - FullModel) / Average
    # Save intermediate result. Now removed, to ensure I am not overwriting the file
    #saveRDS(score_result, file = paste0("/work/Home/Data_for_masters_project/Save_target_enhancers/Score_result_", motif, ".rds"))
    score_result
  }
  
  stopCluster(cl)

  # Assign column names to scores
  colnames(score) = colnames(motifs)
  rownames(score) = rownames(motifs)
  
  # Process scores and define targets
  target.list = list()
  for (motif in colnames(score)) {
    tmp = data.frame(contribution = score[, motif])
    tmp$PeakID = rownames(score)
    motif_count = data.frame(nmotifs = raw_motifs[, motif])
    motif_count$PeakID = rownames(motif_count)
    tmp = merge(tmp, motif_count, by = "PeakID")
    tmp$target = 0
    tmp[tmp$contribution > 0 & tmp$nmotifs > 0, "target"] = 1
    target.list[[motif]] = tmp
  }
  
  # Create a target matrix
  target.matrix = do.call("cbind", lapply(target.list, function(x) x$target))
  rownames(target.matrix) = rownames(motifs)
  # Return
  return(list(matrix = target.matrix, data = target.list))
}
MotifWeightMatrix = function(CPM, ATAC, targets, ncores){
  # Make a variable TSS that contains relevant information from the CPM data.
  TSS <- rbind(CPM[, c("Chr", "TSS", "Factor")])
  # Order the TSS data frame by the chromosome.
  TSS <- TSS[order(TSS$Chr), ]
  # Set row names of TSS to the corresponding gene identifier ("Factor").
  rownames(TSS) <- TSS$Factor
  # Define a region of interest around the TSS: 100 kb upstream.
  TSS$Start <- TSS$TSS - 100000
  # Define 100 kb downstream of TSS.
  TSS$End <- TSS$TSS + 100000
  # Create an identifier for each gene/row.
  TSS$subjectHits <- seq(1, nrow(TSS), by = 1)
  
  # Create a GRanges object for genes using the TSS data.
  GenesRanges <- GRanges(seqnames = TSS[, 1],
                         IRanges(start = TSS$Start, end = TSS$End),
                         strand = rep("+", nrow(TSS)),
                         names = TSS[, 3])
  
  ## Make a GRanges object from the enhancers/ATAC-seq peaks.
  RegionRanges <- ATAC@rowRanges
  Regions = as.data.frame(ATAC@rowRanges)
  # Assign an identifier ("queryHits") for each ATAC-seq peak.
  Regions$queryHits = seq(1, nrow(Regions), 1)
  
  ## Find overlaps between gene regions and ATAC-seq peak regions.
  Overlap <- findOverlaps(query = RegionRanges, subject = GenesRanges, type = "any")
  # Convert the overlap result to a data frame.
  Overlap <- as.data.frame(Overlap)
  
  ## Merge additional information from genes and enhancers into the overlap object.
  Overlap <- merge(Overlap, Regions[, c(1:3, 6, ncol(Regions))], by = "queryHits")
  # Merge the result with information from TSS (gene side), selecting Factor, TSS, and subjectHits.
  Overlap <- merge(Overlap, TSS[, c("Factor", "TSS", "subjectHits")], by = "subjectHits")
  
  ## Calculate scaled regulatory potential.
  Overlap$Center <- (Overlap[, 5] + Overlap[, 4]) / 2
  # Calculate the distance between this center and the gene's TSS.
  Overlap$Distance <- Overlap$Center - Overlap$TSS
  # Keep overlaps only if the absolute distance is within 100 kb.
  Overlap <- Overlap[abs(Overlap$Distance) <= 100000, ]
  # Compute scaled linear regulatory potential
  Overlap$Potential <- (exp(-(0.5 + 4 * (abs(Overlap$Distance) / 100000))) - exp(-1 * (0.5 + 4))) /
    (max(exp(-1 * (0.5 + (4 * (seq(0, 100000, by = 10) / 100000)))) - exp(-1 * (0.5 + 4))))
  
  # Find motifs that have targets in enhancer regions overlapping genes.
  motifs_with_hits = c()
  # Loop through each element (motif) in the targets data list.
  for (i in 1:length(targets$data)) {
    # Extract targets for the current motif.
    current_motif_targets = targets$data[[i]]
    # Filter out targets with no motif count or target count.
    current_motif_targets = current_motif_targets[current_motif_targets$nmotifs > 0 & 
                                                    current_motif_targets$target > 0, ]
    # Check if any of the motif's peaks appear in the Overlap table.
    if (any(current_motif_targets$PeakID %in% Overlap$peak_id) == TRUE) {
      # If so, append the name of the motif to the list.
      motifs_with_hits = c(motifs_with_hits, names(targets$data)[i])
    }
  }
  
  # Create an empty data frame to store weighted motif counts.
  WeightMatrix <- as.data.frame(matrix(0, nrow = length(motifs_with_hits),
                                       ncol = length(unique(Overlap$Factor))))
  colnames(WeightMatrix) <- unique(Overlap$Factor)
  rownames(WeightMatrix) <- motifs_with_hits
  
  # Filter targets for each motif into a list.
  filtered_motif_list = list()
  for (i in 1:nrow(WeightMatrix)) {
    current_motif = rownames(WeightMatrix)[i]
    current_motif_targets = targets$data[[which(names(targets$data) == current_motif)]]
    # Filter out targets with no counts.
    current_motif_targets = current_motif_targets[current_motif_targets$nmotifs > 0 &
                                                    current_motif_targets$target > 0, ]
    filtered_motif_list[[i]] = current_motif_targets
    names(filtered_motif_list)[i] = current_motif
  }
  
  # Create a list of genes with their associated peaks and potential values.
  filtered_gene_list = list()
  counter = 1
  # Loop through each gene in WeightMatrix.
  for (m in 1:ncol(WeightMatrix)) {
    current_gene = colnames(WeightMatrix)[m]
    # Get the ATAC peak IDs associated with the current gene.
    gene_peaks = Overlap[Overlap$Factor == current_gene, "peak_id"]
    if (length(gene_peaks) > 0) {
      # Store the peak IDs and their Potential scores for this gene.
      filtered_gene_list[[counter]] = Overlap[Overlap$Factor == current_gene &
                                                Overlap$peak_id %in% gene_peaks, c("peak_id", "Potential")]
      names(filtered_gene_list)[counter] = current_gene
      counter = counter + 1
    }
  }
  
  # Set up a parallel computing cluster using the available number of cores (ncores).
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  # Set library paths on the worker nodes.
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))
  
  # Load required libraries on each worker node.
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
    # Add any other required libraries here.
  })
  
  # Use parallel processing to compute the weight for each motif across genes.
  WeightMatrix <- foreach(motif_id = seq_along(filtered_motif_list), .combine = 'rbind', .packages = "Matrix") %dopar% {
    # Get the current motif data.
    motif <- filtered_motif_list[[motif_id]]
    # Initialize a numeric vector to store the result for each gene.
    result_list <- vector("numeric", length(filtered_gene_list))
    
    # Loop through each gene in the filtered gene list.
    for (gene_id in seq_along(filtered_gene_list)) {
      gene_matrix <- filtered_gene_list[[gene_id]]
      # Merge motif targets with gene peak information based on peak IDs.
      merged_data <- merge(motif, gene_matrix, by.x = "PeakID", by.y = "peak_id")
      # Sum the product of the motif count and the Potential value for that gene.
      result_list[gene_id] <- sum(merged_data$nmotifs * merged_data$Potential)
    }
    
    # Convert the result list into a sparse matrix row.
    result_matrix <- Matrix(result_list, nrow = 1, sparse = TRUE)
    
    # Set proper row and column names for this row in the result matrix.
    rownames(result_matrix) <- names(filtered_motif_list)[motif_id]
    colnames(result_matrix) <- names(filtered_gene_list)
    
    # Return the result matrix for this motif.
    result_matrix
  }
  
  # Stop the parallel cluster after computations are complete.
  parallel::stopCluster(cl)
  
  # Clean up any lingering parallel objects.
  unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
  }
  unregister_dopar()
  
  # Transpose the weight matrix so that genes become rows (keeping the matrix sparse).
  WeightMatrix <- t(WeightMatrix)
  
  # Compute column means for the sparse matrix using Matrix::colMeans.
  Cols <- Matrix::colMeans(WeightMatrix)
  
  # Normalize the WeightMatrix by subtracting the corresponding column mean from each element.
  WeightMatrix <- WeightMatrix - Matrix(Cols,
                                        nrow = nrow(WeightMatrix),
                                        ncol = ncol(WeightMatrix),
                                        byrow = TRUE,
                                        sparse = TRUE)
  
  # Return the final WeightMatrix.
  return(WeightMatrix)
  
}
calculateMotifActivity = function(design, CPM, WeightMatrix, Pairing = FALSE, mapping, a, ncores, seed = 42){
  # rename the design vector and make it numeric
  RNADesign <- design
  RNADesign <- as.numeric(as.character(RNADesign)) 
  
  # get the number of samples and conditions
  Samples <- length(RNADesign)                     
  Conditions <- length(unique(RNADesign))          
  
  ## extract information
  metadata <- CPM[, 1:4]                           
  DE_statistics <- CPM[, (Samples + 5):ncol(CPM)]    
  CPM <- CPM[, 5:(Samples + 4)]                     
  
  # Compute column means and standard deviations for the gene expression data.
  Cols <- colMeans(CPM)                            
  SD <- apply(CPM, 2, sd)                          
  
  # Standardize the data.
  CPM <- t((t(CPM) - Cols) / SD)                    
  
  # Combine the metadata, standardized expression data, and the DE statistics.
  CPM <- cbind(metadata, CPM, DE_statistics)       
  rownames(CPM) <- CPM$Factor                       
  
  # Find shared row names between WeightMatrix and CPM (common genes).
  shared_rownames <- intersect(rownames(WeightMatrix), rownames(CPM))
  
  # Subset both matrices to include only those shared genes.
  WeightMatrix <- WeightMatrix[shared_rownames, , drop = FALSE]
  CPM <- CPM[shared_rownames, , drop = FALSE]
  
  ### Perform stage 2 ridge regression of all samples
  cat("\n\t\t Calculating factor activity - Full model stage")
  
  # Setup parallel computation
  cl <- parallel::makeCluster(ncores)                        
  doParallel::registerDoParallel(cl)                        
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages"))  
  
  # Load necessary libraries on each worker node.
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
  })
  
  # Set seed for reproducibility (this seed is used for all workers)
  set.seed(seed)
  
  # Parallel loop: Perform ridge regression (with cross validation) across the samples.
  Coef <- foreach(samp = 5:(Samples + 4), 
                  .combine = cbind, 
                  .packages = c("glmnet", "Matrix")) %dopar% {
                    # Set a worker-specific seed for reproducibility.
                    set.seed(seed + samp)
                    # Use 3-fold cross validation to choose lambda (regularization parameter)
                    cv <- cv.glmnet(y = as.matrix(CPM[, samp]), 
                                    x = WeightMatrix, 
                                    alpha = a,                          
                                    standardize = FALSE,                
                                    intercept = TRUE, 
                                    parallel = FALSE,                   
                                    standardize.response = FALSE, 
                                    nfolds = 3)
                    lambda <- cv$lambda.min                   
                    # Extract regression coefficients at lambda.min as a sparse vector (skip the intercept).
                    Tmp <- as.data.frame(as.matrix(coef(cv, s = "lambda.min")))
                    sparse_coef <- Matrix(Tmp[2:nrow(Tmp), 1], sparse = TRUE)
                    sparse_coef                              
                  }
  
  ## Setup objects
  Coef <- data.frame(Coef)                   
  colnames(Coef) <- RNADesign                
  Coef$Motif <- colnames(WeightMatrix)       
  
  ## Convert to standard scores for significance testing
  Significant <- Coef                        
  for (i in 1:Samples) {
    SD <- sd(Significant[, i])               
    Mean <- mean(Significant[, i])           
    Significant[, i] <- (Significant[, i] - Mean) / SD  
  }
  
  ## Set up for all combinations of tests (pairwise comparisons between conditions)
  Tests <- expand.grid(RNADesign, RNADesign)   # Create all combinations of conditions.
  Tests <- Tests[Tests[, 1] != Tests[, 2], ]    # Remove self comparisons.
  Tests$Combine <- paste(Tests[, 1], Tests[, 2], sep = "-")  # Combine conditions into a unique string.
  Tests <- Tests[duplicated(Tests$Combine) == F, c(1, 2)]     # Remove duplicate comparisons.
  Tests <- Tests[Tests[, 1] < Tests[, 2], ]      # Keep only one direction for each pair.
  Tests <- Tests[order(Tests[, 1], Tests[, 2]), ] # Order the tests.
  SampleOutline <- cbind(RNADesign, seq(1, Samples, by = 1))  # Create a mapping between conditions and column indices.
  
  ## Perform T tests for pairwise comparisons of motif activities
  for (q in 1:nrow(Tests)) {
    # Add a new column to Significant for storing the p-value for the current comparison.
    Significant[,(ncol(Significant) + 1)] <- 1
    colnames(Significant)[ncol(Significant)] <- paste("Pvalue_", Tests[q, 1], "_vs_", Tests[q, 2], sep = "")
    if (Pairing == TRUE) {  # If the design is paired:
      for (i in 1:nrow(Significant)) {
        Significant[i, ncol(Significant)] <- t.test(
          as.numeric(Significant[i, SampleOutline[SampleOutline[, 1] == Tests[q, 1], 2]]),
          as.numeric(Significant[i, SampleOutline[SampleOutline[, 1] == Tests[q, 2], 2]]),
          paired = TRUE)$p.value
      }
    } else {  # Otherwise (unpaired design):
      for (i in 1:nrow(Significant)) {
        Significant[i, ncol(Significant)] <- t.test(
          as.numeric(Significant[i, SampleOutline[SampleOutline[, 1] == Tests[q, 1], 2]]),
          as.numeric(Significant[i, SampleOutline[SampleOutline[, 1] == Tests[q, 2], 2]]))$p.value
      }
    }
  }
  
  ## Get the minimal P-values for each motif across comparisons
  if (length(grep("Pvalue", colnames(Significant))) == 1) {
    Significant$Min <- Significant[, grep("Pvalue", colnames(Significant))]
  }
  if (length(grep("Pvalue", colnames(Significant))) >= 2) {
    Significant$Min <- apply(Significant[, grep("Pvalue", colnames(Significant))], 1, FUN = "min")
  }
  
  ## Get the minimal FDR for each gene from the CPM table
  if (length(grep("FDR", colnames(CPM))) == 1) {
    CPM$MinimalFDR <- CPM[, grep("FDR", colnames(CPM))]
  }
  if (length(grep("FDR", colnames(CPM))) >= 2) {
    CPM$MinimalFDR <- apply(CPM[, grep("FDR", colnames(CPM))], 1, FUN = "min")
  }
  
  ## Identify putative regulators by combining motif activity and gene expression data
  colnames(mapping) <- c("Factor", "Motif", "Evidence")  
  Significant <- merge(Significant, mapping, by = "Motif")
  # Merge in minimal FDR and 'Max' (maximum expression?) information from CPM by the gene Factor.
  Significant <- merge(Significant, CPM[, c("MinimalFDR", "Factor", "Max")], all.x = TRUE, by = "Factor")
  # Replace any missing FDR or Max values with 1 (non-significant/default).
  Significant[is.na(Significant$MinimalFDR), "MinimalFDR"] <- 1
  Significant[is.na(Significant$Max), "Max"] <- 1
  
  ## Calculate average motif activity for each condition
  # Create a data frame linking each condition to its corresponding column in the Significant matrix.
  ColIndices <- data.frame(Design = as.numeric(RNADesign), Col = seq(3, (2 + length(RNADesign)), by = 1))
  for (i in unique(RNADesign)) {
    Significant$Tmp <- 1  # Temporary column to hold the activity calculations.
    # Calculate row means for the current condition's columns and store as activity.
    Significant[, ncol(Significant)] <- rowMeans(Significant[, ColIndices[ColIndices$Design == i, 2]])
    # Rename the temporary column to reflect the condition.
    colnames(Significant)[ncol(Significant)] <- paste("Activity_Condition_", i, sep = "")
  }
  
  ## Get maximal activity across all conditions for each motif
  Significant$MaximalActivity <- apply(Significant[, grep("Activity_Condition", colnames(Significant))], 1, FUN = "max")
  
  ## Calculate the weighted P-value to rank regulators
  # The weighted P-value is a product of the ranks of several factors:
  # - Inverse rank of gene expression (Max)
  # - Inverse rank of maximal activity (Activity)
  # - Direct rank of minimal p-value from t-tests
  # - Direct rank of minimal FDR from gene expression
  Significant$WeightedPvalue <- (rank(-Significant$Max) / nrow(Significant)) * 
    (rank(-Significant$MaximalActivity) / nrow(Significant)) * 
    (rank(Significant$Min) / nrow(Significant)) * 
    (rank(Significant$MinimalFDR) / nrow(Significant))
  
  # Save a full copy of the results before filtering.
  Full <- Significant
  # Filter for significance: weighted P-value below threshold and minimal p-value ≤ 0.05.
  Significant <- Full[Full$WeightedPvalue <= 0.001, ]
  Significant <- Significant[Significant$Min <= 0.05, ]
  # Remove duplicate motifs.
  Significant <- Significant[duplicated(Significant$Motif) == FALSE, ]
  
  # Bundle the results into a list.
  Results <- list(Coef = Coef, 
                  Z_norm_CPM = CPM, 
                  Significant = Significant, 
                  Full = Full)
  # Return the final Results list.
  return(Results)
}
findTargetGenes = function(WeightMatrix, MotifActivity, ATAC, CPM, design, TargetsOption = "no", ncores){
  # Initialization of parameters
  # Get the Z-normalized gene expression (CPM) matrix.
  Input <- MotifActivity$Z_norm_CPM       
  # Get the table of significant regulators (from previous analysis).
  Significant <- MotifActivity$Significant  
  # Full table with all regulators before filtering.
  Full <- MotifActivity$Full 
  # Regression coefficients (motif activities) from the previous stage.
  Coef <- MotifActivity$Coef 
  # Use the design vector for the RNA samples.
  RNADesign <- design         
  # Use the weight matrix.
  WeightMatrix <- WeightMatrix   
  # Total number of samples.
  Samples <- length(RNADesign)  
  # Number of unique experimental conditions.
  Conditions <- length(unique(RNADesign))
  # Get ATAC-seq peak ranges (genomic regions) from the ATAC object.
  RegionRanges <- ATAC@rowRanges            
  # Convert ATAC peak ranges into a data frame.
  Regions = as.data.frame(ATAC@rowRanges)   
  # Create an identifier for each ATAC peak.
  Regions$queryHits = seq(1, nrow(Regions), 1)  
  
  # Process target genes - Reduced Model Stage
  cat("\n\t\t Calculating factor activity - Reduced model stage")
  Sys.sleep(5)         
  flush.console()      
  
  ## Calculate Eps - Wpm*Ams (error) and summarize per sample.
  # Subset Input to genes present in WeightMatrix.
  Input <- Input[ Input$Factor %in% rownames(WeightMatrix), ] 
  # Order Input by gene (Factor) name.
  Input <- Input[ order(Input$Factor), ]    
  # Copy the WeightMatrix.
  Tmp3 <- WeightMatrix        
  # Order Tmp3 by the first column.
  Tmp3 <- Tmp3[ order(Tmp3[,1]), ]   
  # Transpose Tmp3 (motifs become rows).
  Tmp3 <- t(as.matrix(Tmp3))       
  # Convert to a data frame.
  Tmp3 <- as.data.frame(Tmp3)  
  # Add a column with motif names.
  Tmp3$Motif <- rownames(Tmp3)  
  # Order Tmp3 by motif name.
  Tmp3 <- Tmp3[ order(Tmp3$Motif), ] 
  # Initialize a matrix (FullModel) to store error differences.
  FullModel <- matrix(nrow = (ncol(Tmp3) - 1), ncol = Samples)   
  # Copy of regression coefficients.
  Tmp4 <- Coef    
  # Remove duplicated motif entries from Tmp4.
  Tmp4 <- Tmp4[ duplicated(Tmp4$Motif) == F, ]                  
  
  # For each sample, calculate the squared error between observed gene expression and model prediction.
  for (samp in 1:Samples) {
    # Temporary matrix to capture per-sample weighted motif contributions.
    Tmp2 <- matrix(ncol = (ncol(Tmp3) - 1), nrow = nrow(Tmp3))  
    # Extract the coefficient for the current sample along with Motif name.
    Tmp <- Tmp4[, c((samp), ncol(Tmp4))]  
    # Order Tmp by motif name.
    Tmp <- Tmp[ order(Tmp$Motif), ]                              
    for (i in 1:(ncol(Tmp3) - 1)) {
      # Multiply each motif weight by its coefficient.
      Tmp2[, i] <- as.numeric(as.character(Tmp3[, i])) * Tmp[, 1]
    }
    # Compute squared error: (observed gene expression - sum of weighted motif contributions)^2.
    FullModel[, samp] <- ( as.numeric(as.character(Input[, (4 + samp)])) - colSums(Tmp2) )^2
  }
  
  ## Calculate average standard deviation (used for later scaling).
  Average <- sum(rowSums(FullModel)) / (nrow(FullModel) * ncol(FullModel))
  
  ## Define the motifs to test.
  TestMotifs <- Coef
  if (TargetsOption == "no") {
    # If TargetsOption is "no", only test motifs that are already deemed significant.
    TestMotifs <- TestMotifs[TestMotifs$Motif %in% Significant$Motif, ]
  }
  
  ## Setup parallel processing to compute target gene scores for each motif
  cl <- parallel::makeCluster(ncores)                      
  doParallel::registerDoParallel(cl)                      
  clusterEvalQ(cl, .libPaths("/work/Home/Data_for_masters_project/R packages")) 
  clusterEvalQ(cl, {
    library(Matrix)
    library(doParallel)
  })
  
  ## Loop through all motifs in TestMotifs and determine target genes in parallel.
  Targets <- foreach(motf = 1:nrow(TestMotifs), .inorder = TRUE, .packages = c("Matrix")) %dopar% {
    # Ensure strings are not converted to factors.
    options(stringsAsFactors = FALSE)    
    # Define the current motif.
    Motif <- TestMotifs[motf, "Motif"]  
    # Initialize matrix to capture reduced model errors.
    RedModel <- matrix(nrow = (nrow(WeightMatrix)), ncol = Samples)  
    # Remove the current motif column from the WeightMatrix.
    Tmp3 <- WeightMatrix
    Tmp3 <- Tmp3[, colnames(Tmp3) != Motif]
    Tmp3 <- Tmp3[ order(Tmp3[, 1]), ]
    Tmp3 <- t(as.matrix(Tmp3))
    Tmp3 <- as.data.frame(Tmp3)
    Tmp3$Motif <- rownames(Tmp3)
    Tmp3 <- Tmp3[ order(Tmp3$Motif), ]
    # For each sample, compute reduced model error without the current motif.
    for (samp in 1:Samples) {
      Tmp2 <- matrix(ncol = (ncol(Tmp3) - 1), nrow = nrow(Tmp3))
      Tmp <- Tmp4[, c((samp), ncol(Tmp4))]
      Tmp <- Tmp[Tmp$Motif != Motif, ]
      Tmp <- Tmp[ order(Tmp$Motif), ]
      for (i in 1:(ncol(Tmp3) - 1)) {
        Tmp2[, i] <- as.numeric(as.character(Tmp3[, i])) * Tmp[, 1]
      }
      # Compute squared error for the reduced model.
      RedModel[, samp] <- ( as.numeric(as.character(Input[, (4 + samp)])) - colSums(Tmp2) )^2
    }
    ## Calculate a motif gene contribution score: sum the difference between the reduced model and the full model, normalized by Average.
    Score <- matrix(nrow = nrow(RedModel), ncol = 1)
    for (site in 1:nrow(RedModel)) {
      Score[site, 1] <- sum(RedModel[site, ] - FullModel[site, ]) / Average
    }
    ## Prepare a data frame to combine weight information for the current motif:
    MatrixData <- data.frame(as.matrix(WeightMatrix))
    # Subset the data to the current motif (and the "Factor" column).
    MatrixData <- MatrixData[, colnames(MatrixData) == Motif | colnames(MatrixData) == "Factor"]
    # Combine the motif's weight data with the calculated Score.
    Score <- as.data.frame(cbind(MatrixData, Score))
    ## Calculate a p-value based on ranking: combine inverse rank of weight and score.
    Score$Pval <- (rank(-Score[, 1]) / nrow(Score)) * (rank(-Score[, 2]) / nrow(Score))
    rownames(Score) <- rownames(WeightMatrix)
    ## Save the Score result for this motif to disk.
    saveRDS(Score, file = paste0("/work/Home/Data_for_masters_project/save_target_genes/Score_result_", Motif, ".rds"))
    Score   # Return the Score for this motif.
  }
  
  stopCluster(cl)   
  
  # Process the results and finalize output
  cat("\n\t\t Processing results")
  
  ## Setup objects for final output
  # Gene (motif) activity table.
  Gene_Activity <- Coef                      
  Enhancers <- Regions[, c(1:(ncol(Regions) - 1))]  
  
  ## Process target genes: assign 'Target' flag based on p-value threshold.
  for (motf in 1:nrow(TestMotifs)) {
    Motif <- TestMotifs[motf, "Motif"]
    Tmp <- Targets[[motf]]
    Tmp$Target <- 0                        # Initialize target flag to 0.
    Tmp[Tmp$Pval <= 0.005, "Target"] <- 1   # Mark as target if p-value is <= 0.005.
    Targets[[motf]] <- Tmp
    names(Targets)[[motf]] <- Motif          # Rename list element using motif name.
  }
  Target_Genes <- Targets                  # Save the list of target genes.
  
  ## Calculate average motif activity for each sample
  ColIndices <- data.frame(Design = as.numeric(RNADesign), Col = seq(1, (length(RNADesign)), by = 1))
  for (i in unique(RNADesign)) {
    Coef$Tmp <- 1
    # Compute row means for coefficients corresponding to each condition.
    Coef[, ncol(Coef)] <- rowMeans(Coef[, ColIndices[ColIndices$Design == i, 2]])
    colnames(Coef)[ncol(Coef)] <- paste("Average_Motif_Activity_Condition_", i, sep = "")
  }
  
  ## Merge additional information: combine Coef with selected columns from Full.
  Coef <- suppressWarnings(merge(Coef, Full[, c(1, 2, grep("Min", colnames(Full)), grep("WeightedPvalue", colnames(Full)))], by = "Motif"))
  Coef <- Coef[, grep("MinimalFDR", colnames(Coef), invert = TRUE)] 
  
  ## Calculate average gene expression strength per condition:
  ColIndices <- data.frame(Design = as.numeric(RNADesign), 
                           Col = seq(5, (4 + length(RNADesign)), by = 1))
  for (i in unique(RNADesign)) {
    # Append a new column to CPM with average gene expression for condition i.
    CPM[, ncol(CPM) + 1] <- rowMeans(CPM[, ColIndices[ColIndices$Design == i, "Col"]])
    colnames(CPM)[ncol(CPM)] <- paste("Average_Gene_Expression_Condition_", i, sep = "")
  }
  
  ## Perform correlations between motif activity and gene expression (per condition)
  Coef <- merge(Coef, CPM[, c(1, grep("Expression-Sample", colnames(CPM)))], by = "Factor")
  Coef$Pearsons <- 0
  # For each motif, compute the Pearson correlation between its activity and gene expression.
  for (i in 1:nrow(Coef)) { 
    Coef[i, "Pearsons"] <- cor(as.numeric(Coef[i, c(3:(2 + Samples))]), 
                               as.numeric(Coef[i, grep("Expression-Sample", colnames(Coef))])) 
  }
  # Remove rows with missing Pearson correlation.
  Coef <- na.omit(Coef)
  
  ## Clean up columns by removing those related to "Expression-Sample"
  Coef <- Coef[, grep("Expression-Sample", colnames(Coef), invert = TRUE)]
  Coef <- Coef[, c(1, 2, (Samples + 3):ncol(Coef))]
  colnames(Coef)[grep("Min", colnames(Coef))] <- "MinimalActivityPvalue"
  colnames(Coef)[grep("Combined", colnames(Coef))] <- "WeightedPvalue"
  
  ## Identify putative regulatory factors:
  # Create a flag "CausalTF" based on thresholds for WeightedPvalue and MinimalActivityPvalue.
  Coef$CausalTF <- 0
  Coef[Coef$WeightedPvalue <= 0.005 & Coef$MinimalActivityPvalue <= 0.05, "CausalTF"] <- 2
  Coef[Coef$WeightedPvalue <= 0.001 & Coef$MinimalActivityPvalue <= 0.05, "CausalTF"] <- 1
  
  ## Prepare final results excluding columns related to "Evidence"
  Results <- list(MotifActivities = Coef, target_Genes = Target_Genes)
  cat("\n\t\t Outputting results\n\n")
  
  # Return the final Results list.
  return(Results)
}

start_time2 <- Sys.time()
CPM <- analyzeRNA(RNA, paired = FALSE, design = design)
end_time2 <- Sys.time()
time_CPM <- end_time2-start_time2
#saveRDS(CPM, "/work/Home/Data_for_masters_project/CPM_chromVAR_motifs_GC_and_edgeR")
#saveRDS(time_CPM, "/work/Home/Data_for_masters_project/time_CPM_chromVAR_motifs_GC_and_edgeR")

start_time6 <- Sys.time()
target_enhancers = findTargetEnhancers(motifs = filteredMotifs, enhancers = normEnhancers, activity = activities, raw_motifs = motif_counts, ncores = 2)
end_time6 <- Sys.time()
time_target_enhancers <- end_time6-start_time6 
#saveRDS(target_enhancers, "/work/Home/Data_for_masters_project/target_enhancers_chromVAR_motifs_GC_and_edgeR")
#saveRDS(time_target_enhancers, "/work/Home/Data_for_masters_project/time_target_enhancers_chromVAR_motifs_GC_and_edgeR")
target_enhancers <- readRDS("/work/Home/Data_for_masters_project/target_enhancers_chromVAR_motifs_GC_and_edgeR")
start_time7 <- Sys.time()
WeightMatrix = MotifWeightMatrix(CPM = CPM, ATAC = ATAC, targets = target_enhancers, ncores = 2)
end_time7 <- Sys.time()
time_WeightMatrix <- end_time7-start_time7
#saveRDS(WeightMatrix, "/work/Home/Data_for_masters_project/WeightMatrix_chromVAR_motifs_GC_and_edgeR")
#saveRDS(time_WeightMatrix, "/work/Home/Data_for_masters_project/time_WeightMatrix_chromVAR_motifs_GC_and_edgeR")

start_time8 <- Sys.time()
motifactivity = calculateMotifActivity(design = design, CPM = CPM, WeightMatrix = WeightMatrix, mapping = mapping, a = 0, ncores = 2)
end_time8 <- Sys.time()
time_motifactivity <- end_time8-start_time8
#saveRDS(motifactivity, "/work/Home/Data_for_masters_project/motifactivity_chromVAR_motifs_GC_and_edgeR")
#saveRDS(time_motifactivity, "/work/Home/Data_for_masters_project/time_motifactivity_chromVAR_motifs_GC_and_edgeR")

start_time9 <- Sys.time()
target_genes = findTargetGenes(WeightMatrix = WeightMatrix, MotifActivity = motifactivity, ATAC = ATAC, CPM = CPM, design = design, TargetsOption = "no", ncores = 2)
end_time9 <- Sys.time()
time_target_genes <- end_time9-start_time9
#saveRDS(target_genes, "/work/Home/Data_for_masters_project/target_genes_chromVAR_motifs_GC_and_edgeR")
#saveRDS(time_target_genes, "/work/Home/Data_for_masters_project/time_target_genes_chromVAR_motifs_GC_and_edgeR")

########### plot the Gene motif activity
motifactivity <- readRDS("/work/Home/Data_for_masters_project/motifactivity_chromVAR_motifs_GC_and_edgeR")
activities <- as.matrix(motifactivity$Coef)
rownames(activities) <- activities[, ncol(activities)]
activities <- activities[, 1:ncol(activities)-1]
activities <- as.matrix(activities)
# Load required library
library(Matrix)
# Convert all entries to numeric
activities <- matrix(
  as.numeric(activities),
  nrow = nrow(activities),
  ncol = ncol(activities),
  dimnames = dimnames(activities)  # preserve row/column names
)
# Extract base column names (removing suffixes like .1, .2, etc.)
base_colnames <- sub("\\.\\d+$", "", colnames(activities))

# Compute mean over each group
activities_mean <- sapply(unique(base_colnames), function(cell_type) {
  col_indices <- which(base_colnames == cell_type)
  rowMeans(activities[, col_indices, drop = FALSE])
})

# Assign new column names
colnames(activities_mean) <- c("hESC 1", "hESC 2", "hESC 3",
                               "Neuron like cells", "T/P like cells", "Goblet like cells", 
                               "OP like cells", "VE like cells",
                               "FL like cell", "Podocyte like cells", "SMP like cells", 
                               "BE like cells", "Adipocyte like cells")

### Map the motifs to TFs
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

# Sort the rownames of activity_object
# Create a mapping vector with old rownames as keys and new rownames as values
rename_mapping <- c(
  "LINE11277.motif" = "DUX1.motif",
  "LINE11282.motif" = "DUX3.motif",
  "LINE4118.motif" = "ZNF75C.motif"
)

# Replace matching rownames
rownames(activities_mean) <- ifelse(
  rownames(activities_mean) %in% names(rename_mapping),
  rename_mapping[rownames(activities_mean)],
  rownames(activities_mean)  # Keep other rownames unchanged
)

activities_mean <- activities_mean[order(rownames(activities_mean)), ]
# Sort the rownames of mapping based on motif
mapping <- mapping[order(mapping$V2), ]
# give rownames to activities
rownames(activities_mean) <- mapping$V1

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
# convert to matrix
heatmap_data <- as.matrix(activities_mean)
# Z-score the data
heatmap_data <- t(scale(t(heatmap_data)))
IMAGE_motif_activity <- heatmap_data[complete.cases(heatmap_data), ]
col_fun <- colorRamp2(
  breaks = c(min(IMAGE_motif_activity),  # lowest
             quantile(IMAGE_motif_activity, 0.33),  # lower intermediate
             0,                           # mid point
             quantile(IMAGE_motif_activity, 0.67),  # upper intermediate
             max(IMAGE_motif_activity)),          # highest
  colors = c("#4575B4", "#91BFDB", "white", "#FCAE61", "#D73027")
)

BASE_FS <- 18   
ht <- Heatmap(
  IMAGE_motif_activity,
  col              = col_fun,
  na_col           = "grey",
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  column_dend_side = "top",
  row_names_side   = "left",
  column_names_side = "bottom",
  show_row_names   = FALSE,
  
  # axis labels
  column_names_gp  = gpar(fontsize = BASE_FS,     fontface = "bold"),
  row_names_gp     = gpar(fontsize = BASE_FS - 4),
  
  # titles next to the matrix
  row_title        = "Transcription Factors",
  column_title_side = "top",
  row_title_gp     = gpar(fontsize = BASE_FS + 2, fontface = "bold"),
  column_title_gp  = gpar(fontsize = BASE_FS + 2, fontface = "bold"),
  
  # legend
  heatmap_legend_param = list(
    title      = "Z-Score Motif Activity",
    title_gp   = gpar(fontsize = BASE_FS,     fontface = "bold"),
    labels_gp  = gpar(fontsize = BASE_FS - 1, fontface = "bold")
  )
)

draw(
  ht,
  column_title    = "IMAGE Gene Motif Activity Heatmap (4 Replicates, Differential Accessible Regions, Alpha = 0)",
  column_title_gp = gpar(fontsize = BASE_FS + 4, fontface = "bold")
)
### validate and compare enhancer motif activity and gene motif activity with eachother
activities <- readRDS("/work/Home/Data_for_masters_project/activities_chromVAR_motifs_GC_and_edgeR")
activities <- as.matrix(activities)
# Load required library
library(Matrix)
# Convert all entries to numeric
activities <- matrix(
  as.numeric(activities),
  nrow = nrow(activities),
  ncol = ncol(activities),
  dimnames = dimnames(activities)  
)
# Extract base column names (removing suffixes like .1, .2, etc.)
base_colnames <- sub("\\.\\d+$", "", colnames(activities))

# Compute mean over each group
activities_mean <- sapply(unique(base_colnames), function(cell_type) {
  col_indices <- which(base_colnames == cell_type)
  rowMeans(activities[, col_indices, drop = FALSE])
})

# Assign new column names
colnames(activities_mean) <- c("hESC 1", "hESC 2", "hESC 3",
                               "Neuron like cells", "T/P like cells", "Goblet like cells", 
                               "OP like cells", "VE like cells",
                               "FL like cell", "Podocyte like cells", "SMP like cells", 
                               "BE like cells", "Adipocyte like cells")

### Map the motifs to TFs
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

# Sort the rownames of activity_object
# Create a mapping vector with old rownames as keys and new rownames as values
rename_mapping <- c(
  "LINE11277.motif" = "DUX1.motif",
  "LINE11282.motif" = "DUX3.motif",
  "LINE4118.motif" = "ZNF75C.motif"
)

# Replace matching rownames
rownames(activities_mean) <- ifelse(
  rownames(activities_mean) %in% names(rename_mapping),
  rename_mapping[rownames(activities_mean)],
  rownames(activities_mean)  # Keep other rownames unchanged
)

activities_mean <- activities_mean[order(rownames(activities_mean)), ]
# Sort the rownames of mapping based on motif
mapping <- mapping[order(mapping$V2), ]
# give rownames to activities
rownames(activities_mean) <- mapping$V1

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
# convert to matrix
heatmap_data <- as.matrix(activities_mean)
# Z-score the data
heatmap_data <- t(scale(t(heatmap_data)))
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]

######## enrichment
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
metadata <- SHARE_seq_subsample_data@meta.data
metadata$TF <- sub("^[^-]*-", "", metadata$TF)
# as.factor to ensure compatability between 
metadata$TF <- as.factor(metadata$TF)
################ TF enrichment relative to GFP AND mCherry control
# Initialize an empty matrix to store the log2 ratios
# Filter the metadata into control and non-control
metadata_control <- metadata[metadata$TF %in% c("GFP", "mCherry"), ]
metadata_no_control <- metadata[!(metadata$TF %in% c("GFP", "mCherry")), ]
# Initialize a data.frame to store results
enrichment_results <- data.frame(
  TF = character(),
  Cluster = character(),
  Enrichment = numeric()
)
# Loop through each unique TF and cluster
unique_TFs <- unique(metadata_no_control$TF)
unique_clusters <- unique(metadata$seurat_clusters)

for (tf in unique_TFs) {
  for (cluster in unique_clusters) {
    # Cells expressing the TF in the cluster
    tf_cluster_count <- sum(metadata_no_control$TF == tf & metadata_no_control$seurat_clusters == cluster)
    
    # Cells expressing the control in the cluster
    control_cluster_count <- sum(metadata_control$seurat_clusters == cluster)
    
    # Cells expressing the TF across all clusters
    tf_total_count <- sum(metadata_no_control$TF == tf)
    
    # Cells expressing the control across all clusters
    control_total_count <- nrow(metadata_control)
    
    # Proportions
    cluster_proportion_tf <- tf_cluster_count / control_cluster_count
    overall_proportion_tf <- tf_total_count / control_total_count
    
    # Log2 enrichment
    enrichment <- log2((cluster_proportion_tf / overall_proportion_tf)+1) 
    
    # Store the result
    enrichment_results <- rbind(enrichment_results, data.frame(
      TF = tf,
      Cluster = cluster,
      Enrichment = enrichment
    ))
  }
}


# Convert results_df to a matrix
enrichment_results <- reshape2::acast(enrichment_results, TF ~ Cluster, value.var = "Enrichment", fill = 0)
enrichment_results <- enrichment_results[,order(as.numeric(colnames(enrichment_results)))]
colnames(enrichment_results) <- c("hESC 1", "hESC 2", "hESC 3",
                                  "Neuron like cells", "T/P like cells", "Goblet like cells", 
                                  "OP like cells", "VE like cells",
                                  "FL like cell", "Podocyte like cells", "SMP like cells", 
                                  "BE like cells", "Adipocyte like cells")

### Expression
mean_list <- list()
clusters <- unique(sort(metadata_no_control$seurat_clusters))
for (cluster in clusters){
  # Get barcodes (cells) for the cluster
  cluster_barcodes <- rownames(metadata_no_control[metadata_no_control$seurat_clusters == cluster, ])
  # Get only the cells that exist in chromVAR_object
  cells <- intersect(colnames(SHARE_seq_subsample_data), cluster_barcodes)
  # Ensure there are cells
  if (length(cells) > 0) {  
    # Compute mean deviations across cells for the cluster
    mean_list[[as.character(cluster)]] <- rowMeans(SHARE_seq_subsample_data@assays$RNA$data[,cells, drop = FALSE], na.rm = TRUE)
  } else {
    warning(paste("No matching cells found for cluster", cluster))
  }
}

Pseudobulk_exprs <- do.call(cbind, mean_list)
colnames(Pseudobulk_exprs) <- c("hESC 1", "hESC 2", "hESC 3",
                                "Neuron like cells", "T/P like cells", "Goblet like cells", 
                                "OP like cells", "VE like cells",
                                "FL like cell", "Podocyte like cells", "SMP like cells", 
                                "BE like cells", "Adipocyte like cells")

##### compare motif activity between IMAGE enhancer motif activity and IMAGE gene expression motif activity
#### line plot
clusters = colnames(Pseudobulk_exprs[, 7:13])
all_folds <- 2:9
expr_threshold <- 0.2
depl_threshold <- 1  

df_all <- list()
k <- 1

for (fold in all_folds) {
  
  enrT   <- log2(fold)
  noEnrT <- log2(depl_threshold)
  
  # collect cluster-level data in 'res'
  for (cluster in clusters) {
    
    # Identify enriched TFs
    enriched_TFs <- names(which(enrichment_results[, cluster] >= enrT))
    
    # Identify not-enriched TFs
    not_enriched_TFs <- names(which(enrichment_results[, cluster] <= noEnrT))
    # Also require that they are low-expressed (<= expr_threshold)
    not_enriched_TFs <- intersect(
      not_enriched_TFs,
      names(which(Pseudobulk_exprs[, cluster] <= expr_threshold))
    )
    
    # Activities for this cluster
    IMAGE_activity <- abs(as.numeric(heatmap_data[, cluster]))
    IM_activity    <- abs(as.numeric(IMAGE_motif_activity[, cluster]))
    
    # Identify only the motifs that actually appear in each matrix:
    valid_enriched_IMAGE_enhancer    <- intersect(enriched_TFs, rownames(heatmap_data))
    valid_control_IMAGE_enhancer     <- intersect(not_enriched_TFs, rownames(heatmap_data))
    valid_enriched_IMAGE_gene <- intersect(enriched_TFs, rownames(IMAGE_motif_activity))
    valid_control_IMAGE_gene  <- intersect(not_enriched_TFs, rownames(IMAGE_motif_activity))
    
    # Build rows for (IMAGE.Enriched, IMAGE.Control, chromVAR.Enriched, chromVAR.Control)
    res_tmp <- data.frame(
      cluster         = cluster,
      type            = "Enriched",
      method          = "IMAGE_enhancer",
      activity_median = median(IMAGE_activity[rownames(heatmap_data) %in% enriched_TFs]),
      n_TFs           = length(valid_enriched_IMAGE_enhancer)
    )
    
    res_tmp <- rbind(res_tmp, data.frame(
      cluster         = cluster,
      type            = "Control",
      method          = "IMAGE_enhancer",
      activity_median = median(IMAGE_activity[rownames(heatmap_data) %in% not_enriched_TFs]),
      n_TFs           = length(valid_control_IMAGE_enhancer)
    ))
    
    res_tmp <- rbind(res_tmp, data.frame(
      cluster         = cluster,
      type            = "Enriched",
      method          = "IMAGE_gene",
      activity_median = median(IM_activity[rownames(IMAGE_motif_activity) %in% enriched_TFs]),
      n_TFs           = length(valid_enriched_IMAGE_gene)
    ))
    
    res_tmp <- rbind(res_tmp, data.frame(
      cluster         = cluster,
      type            = "Control",
      method          = "IMAGE_gene",
      activity_median = median(IM_activity[rownames(IMAGE_motif_activity) %in% not_enriched_TFs]),
      n_TFs           = length(valid_control_IMAGE_gene)
    ))
    
    if (cluster == clusters[1]) {
      res <- res_tmp
    } else {
      res <- rbind(res, res_tmp)
    }
  } # end for(cluster in clusters)
  
  # Summarize across all clusters for this fold
  df_summary <- res %>%
    dplyr::group_by(type, method) %>%
    dplyr::summarize(
      median_of_median = median(activity_median, na.rm = TRUE),
      total_TFs        = sum(n_TFs),
      .groups          = "drop"
    )
  
  # Record which fold was used
  df_summary$enriched_fold <- fold
  
  # Store in our master list
  df_all[[k]] <- df_summary
  k <- k + 1
}

# Combine into a single data frame
df_all <- do.call(rbind, df_all)

# factorise
df_all$method <- factor(df_all$method, levels = c("IMAGE_enhancer", "IMAGE_gene"))
df_all$type   <- factor(df_all$type,   levels = c("Enriched", "Control"))

# Now plot
library(ggplot2)

p2 <- ggplot(df_all, aes(
  x     = enriched_fold, 
  y     = median_of_median,
  color = interaction(method, type),
  group = interaction(method, type)
)) +
  # LOESS smoothing
  geom_smooth(method = "loess", se = FALSE, span = 0.5) + 
  # Keep text labels for the points 2..9
  geom_text(
    data = subset(df_all, enriched_fold %in% 2:9),
    aes(label = paste0("n=", total_TFs)),
    vjust = -1,
    size  = 5,
    show.legend = FALSE,
    colour = "black",
    fontface = "bold"
  ) +
  scale_x_continuous(breaks = 2:9) +  # ensure ticks from 2 to 9
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Median of Median TF Enhancer Activities vs. Enrichment Fold",
    subtitle = "Control: ≤ 1-fold & expr ≤ 0.2; Enriched: ≥ 2 to 9-fold",
    x        = "Enrichment Fold (≥)",
    y        = "Median of Median Activities",
    color    = "Method & Type"
  ) +
  theme_bw(base_size = 14) +              
  theme(
    text           = element_text(face = "bold", colour = "black"),  # make *all* text bold
    plot.title     = element_text(size = 18, hjust = 0.5, face = "bold"),      # bump a bit more
    plot.subtitle  = element_text(size = 15, hjust = 0.5),
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 14),
    legend.text    = element_text(size = 14),
    legend.title   = element_text(size = 14),
    strip.text     = element_text(size = 14),      # facet labels
    legend.position = "right",
    axis.text.x  = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 14, face = "bold", colour = "black")
  )


########## boxplot
library(dplyr)
library(ggplot2)

# test enriched thresholds = 2 to 9
all_folds <- 2:9
expr_threshold <- 0.2
depl_threshold <- 1  

# keep cluster-level results in df_box
df_box <- data.frame()

for (fold in all_folds) {
  enrT   <- log2(fold)              
  noEnrT <- log2(depl_threshold)
  
  # Loop over clusters
  for (cluster in clusters) {
    
    # Define sets of TFs
    enriched_TFs     <- names(which(enrichment_results[, cluster] >= enrT))
    not_enriched_TFs <- names(which(enrichment_results[, cluster] <= noEnrT))
    
    # Constrain the "control" TFs to be low-expressed (≤ expr_threshold)
    not_enriched_TFs <- intersect(
      not_enriched_TFs,
      names(which(Pseudobulk_exprs[, cluster] <= expr_threshold))
    )
    
    # Activities
    IMAGE_activity <- abs(as.numeric(heatmap_data[, cluster]))
    IM_activity    <- abs(as.numeric(IMAGE_motif_activity[, cluster]))
    
    # For each cluster, create a small data.frame with 4 rows:
    # (IMAGE, Enriched), (IMAGE, Control), (chromVAR, Enriched), (chromVAR, Control)
    res_tmp <- data.frame(
      cluster = cluster,
      method  = c("IMAGE_enhancer","IMAGE_enhancer","IMAGE_gene","IMAGE_gene"),
      type    = c("Enriched","Control","Enriched","Control"),
      activity_median = c(
        median(IMAGE_activity[rownames(heatmap_data) %in% enriched_TFs]),
        median(IMAGE_activity[rownames(heatmap_data) %in% not_enriched_TFs]),
        median(IM_activity[rownames(IMAGE_motif_activity) %in% enriched_TFs]),
        median(IM_activity[rownames(IMAGE_motif_activity) %in% not_enriched_TFs])
      )
    )
    
    # Add the fold value
    res_tmp$enriched_fold <- fold
    
    # Append to our big data frame
    df_box <- rbind(df_box, res_tmp)
  }
}

# factorize the method
df_box$method <- factor(df_box$method, levels = c("IMAGE_enhancer", "IMAGE_gene"))
df_box$type   <- factor(df_box$type,   levels = c("Enriched", "Control"))

# Create a "combined" factor that explicitly fixes the order
df_box$combo <- interaction(df_box$method, df_box$type)

# factorise
df_box$combo <- factor(
  df_box$combo,
  levels = c("IMAGE_enhancer.Enriched", "IMAGE_enhancer.Control",
             "IMAGE_gene.Enriched", "IMAGE_gene.Control")
)

p_box2 <- ggplot(df_box, 
                 aes(
                   x = factor(enriched_fold), 
                   y = activity_median,
                   fill = combo
                 )) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  labs(
    title    = "Distribution of Median TF Activities vs. Enrichment Fold",
    subtitle = "Control: ≤ 1-fold and Expression ≤ 0.2, Enriched: ≥ 2 to 9-fold",
    x        = "Enrichment Fold (≥)",
    y        = "Median Activity",
    fill     = "Method & Type"
  ) +
  theme_bw(base_size = 14) +              
  theme(
    text           = element_text(face = "bold", colour = "black"),  # make *all* text bold
    plot.title     = element_text(size = 18, hjust = 0.5, face = "bold"),      # bump a bit more
    plot.subtitle  = element_text(size = 15, hjust = 0.5),
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 14),
    legend.text    = element_text(size = 14),
    legend.title   = element_text(size = 14),
    strip.text     = element_text(size = 14),      # facet labels
    legend.position = "right",
    axis.text.x  = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 14, face = "bold", colour = "black")
  )

library(patchwork)
# plot them
combined_plot <- p2/p_box2 +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

combined_plot

########## 
# subset based on motif activities 5th and 95th percentiles. Characterize the TFs
##########
# load data and prepare it
motifactivity <- readRDS("/work/Home/Data_for_masters_project/motifactivity_chromVAR_motifs_GC_and_edgeR")
activities <- as.matrix(motifactivity$Coef)
rownames(activities) <- activities[, ncol(activities)]
activities <- activities[, 1:ncol(activities)-1]
activities <- as.matrix(activities)
# Load required library
library(Matrix)
# Convert all entries to numeric
activities <- matrix(
  as.numeric(activities),
  nrow = nrow(activities),
  ncol = ncol(activities),
  dimnames = dimnames(activities)  # preserve row/column names
)
# Extract base column names (removing suffixes like .1, .2, etc.)
base_colnames <- sub("\\.\\d+$", "", colnames(activities))

# Compute mean over each group
activities_mean <- sapply(unique(base_colnames), function(cell_type) {
  col_indices <- which(base_colnames == cell_type)
  rowMeans(activities[, col_indices, drop = FALSE])
})

# Assign new column names
colnames(activities_mean) <- c("hESC 1", "hESC 2", "hESC 3",
                               "Neuron like cells", "T/P like cells", "Goblet like cells", 
                               "OP like cells", "VE like cells",
                               "FL like cell", "Podocyte like cells", "SMP like cells", 
                               "BE like cells", "Adipocyte like cells")

### Map the motifs to TFs
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

# Sort the rownames of activity_object
# Create a mapping vector with old rownames as keys and new rownames as values
rename_mapping <- c(
  "LINE11277.motif" = "DUX1.motif",
  "LINE11282.motif" = "DUX3.motif",
  "LINE4118.motif" = "ZNF75C.motif"
)

# Replace matching rownames
rownames(activities_mean) <- ifelse(
  rownames(activities_mean) %in% names(rename_mapping),
  rename_mapping[rownames(activities_mean)],
  rownames(activities_mean)  # Keep other rownames unchanged
)

activities_mean <- activities_mean[order(rownames(activities_mean)), ]
# Sort the rownames of mapping based on motif
mapping <- mapping[order(mapping$V2), ]
# give rownames to activities
rownames(activities_mean) <- mapping$V1

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
# convert to matrix
IMAGE_motif_activity <- as.matrix(activities_mean)

#### enhancer activity
activities <- readRDS("/work/Home/Data_for_masters_project/activities_chromVAR_motifs_GC_and_edgeR")
activities <- as.matrix(activities)
# Load required library
library(Matrix)
# Convert all entries to numeric
activities <- matrix(
  as.numeric(activities),
  nrow = nrow(activities),
  ncol = ncol(activities),
  dimnames = dimnames(activities)  # preserve row/column names
)
# Extract base column names (removing suffixes like .1, .2, etc.)
base_colnames <- sub("\\.\\d+$", "", colnames(activities))

# Compute mean over each group
activities_mean <- sapply(unique(base_colnames), function(cell_type) {
  col_indices <- which(base_colnames == cell_type)
  rowMeans(activities[, col_indices, drop = FALSE])
})

# Assign new column names
colnames(activities_mean) <- c("hESC 1", "hESC 2", "hESC 3",
                               "Neuron like cells", "T/P like cells", "Goblet like cells", 
                               "OP like cells", "VE like cells",
                               "FL like cell", "Podocyte like cells", "SMP like cells", 
                               "BE like cells", "Adipocyte like cells")

### Map the motifs to TFs
library(chromVARmotifs)
data("human_pwms_v2")
# Initialize an empty list to store results
results <- list()

# Loop through each element in the list
for (name in names(human_pwms_v2@listData)) {
  # Extract the TF Name
  tf_name <- as.character(human_pwms_v2@listData[[name]]@name)
  
  # Determine if the motif is Direct or Inferred based on the fourth character
  motif_type <- ifelse(strsplit(name, "_")[[1]][4] == "D", "Direct", "Inferred")
  
  # Construct the second column value
  second_column <- paste0(tf_name, ".motif")
  
  # Append the row to results
  results[[name]] <- c(V1 = tf_name, V2 = second_column, V3 = motif_type)
}

# Convert results to a data.frame
mapping <- do.call(rbind, results)
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
mapping[851,1] <- "ZNF75C"
mapping[851,2] <- "ZNF75C.motif"
mapping[859,1] <- "DUX1"
mapping[859,2] <- "DUX1.motif"
mapping[860,1] <- "DUX3"
mapping[860,2] <- "DUX3.motif"

# Sort the rownames of activity_object
# Create a mapping vector with old rownames as keys and new rownames as values
rename_mapping <- c(
  "LINE11277.motif" = "DUX1.motif",
  "LINE11282.motif" = "DUX3.motif",
  "LINE4118.motif" = "ZNF75C.motif"
)

# Replace matching rownames
rownames(activities_mean) <- ifelse(
  rownames(activities_mean) %in% names(rename_mapping),
  rename_mapping[rownames(activities_mean)],
  rownames(activities_mean)  # Keep other rownames unchanged
)

activities_mean <- activities_mean[order(rownames(activities_mean)), ]
# Sort the rownames of mapping based on motif
mapping <- mapping[order(mapping$V2), ]
# give rownames to activities
rownames(activities_mean) <- mapping$V1

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
# convert to matrix
heatmap_data <- as.matrix(activities_mean)



### plot it
library(ggplot2)
library(tidyr)
library(patchwork)
heatmap_data_long <- pivot_longer(as.data.frame(heatmap_data), 
                                  cols = everything(),
                                  names_to = "Condition",
                                  values_to = "Value")
desired_order <- c("hESC 1", "hESC 2", "hESC 3", "Neuron like cells", 
                   "T/P like cells", "Goblet like cells", "OP like cells",
                   "VE like cells", "FL like cell", "Podocyte like cells", 
                   "SMP like cells", "BE like cells", "Adipocyte like cells")

heatmap_data_long$Condition <- factor(heatmap_data_long$Condition, levels = desired_order)
IMAGE_motif_activity_long <- pivot_longer(as.data.frame(IMAGE_motif_activity), 
                                          cols = everything(),
                                          names_to = "Condition",
                                          values_to = "Value")
IMAGE_motif_activity_long$Condition <- factor(IMAGE_motif_activity_long$Condition, levels = desired_order)

library(dplyr)

# calculate quantiles for the violin plots
df_percentiles <- heatmap_data_long %>%
  group_by(Condition) %>%
  summarize(
    p5 = quantile(Value, 0.05),
    p95 = quantile(Value, 0.95)
  )
df_percentiles_IMAGE <- IMAGE_motif_activity_long %>%
  group_by(Condition) %>%
  summarize(
    p5 = quantile(Value, 0.05),
    p95 = quantile(Value, 0.95)
  )

violin1 <- ggplot(heatmap_data_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  # Draw small horizontal segments at p5 and p95
  geom_segment(
    data = df_percentiles,
    aes(
      x = as.numeric(Condition) - 0.3,
      xend = as.numeric(Condition) + 0.3,
      y = p5,
      yend = p5
    ),
    color = "#D73027",
    size = 0.8
  ) +
  geom_segment(
    data = df_percentiles,
    aes(
      x = as.numeric(Condition) - 0.3,
      xend = as.numeric(Condition) + 0.3,
      y = p95,
      yend = p95
    ),
    color = "#4575B4",
    size = 0.8
  ) +
  theme_bw(base_size = 16) +
  coord_flip() +
  ggtitle("Raw Enhancer Motif Activity for each Cell Type") +
  theme(
    text = element_text(face = "bold", colour = "black"),
    # Center & bold title
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text = element_text(face = "bold", colour = "black")
  )

# For the second plot (Gene Motif Activity)
violin2 <- ggplot(IMAGE_motif_activity_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  # Draw small horizontal segment at the 5th percentile with a mapped color
  geom_segment(
    data = df_percentiles_IMAGE,
    aes(
      x = as.numeric(Condition) - 0.3,
      xend = as.numeric(Condition) + 0.3,
      y = p5,
      yend = p5,
      color = "5th Percentile"   # Legend label for 5th percentile
    ),
    size = 0.8
  ) +
  # Draw small horizontal segment at the 95th percentile with a mapped color
  geom_segment(
    data = df_percentiles_IMAGE,
    aes(
      x = as.numeric(Condition) - 0.3,
      xend = as.numeric(Condition) + 0.3,
      y = p95,
      yend = p95,
      color = "95th Percentile"  # Legend label for 95th percentile
    ),
    size = 0.8
  ) +
  # Define the same custom colors
  scale_color_manual(
    name = "Percentile Lines",
    values = c("5th Percentile" = "#D73027", 
               "95th Percentile" = "#4575B4")
  ) +
  theme_bw(base_size = 16) +
  coord_flip() +
  ggtitle("Raw Gene Motif Activity for each Cell Type") + 
  theme(
    text = element_text(face = "bold", colour = "black"),
    # Center & bold title
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(),
    axis.text = element_text(face = "bold", colour = "black")
  )



library(patchwork)
# plot them
combined_plot <- violin1 + violin2 +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18)
  )

combined_plot

######## For IMAGE_motif_activity
extreme_image_sets <- list()

for (colname in colnames(IMAGE_motif_activity)) {
  # Compute the 5th and 95th percentiles
  p05 <- quantile(IMAGE_motif_activity[, colname], probs = 0.05, na.rm = TRUE)
  p95 <- quantile(IMAGE_motif_activity[, colname], probs = 0.95, na.rm = TRUE)
  
  # Identify TFs in the bottom 5% or top 5%
  bottom_tf_set <- rownames(IMAGE_motif_activity)[
    which(IMAGE_motif_activity[, colname] < p05)
  ]
  top_tf_set <- rownames(IMAGE_motif_activity)[
    which(IMAGE_motif_activity[, colname] > p95)
  ]
  
  # Combine top and bottom
  tf_set <- c(bottom_tf_set, top_tf_set)
  
  # Store it in the list
  extreme_image_sets[[colname]] <- tf_set
}

# Combine into a single vector and remove duplicates
extreme_image_sets_unlisted <- unlist(extreme_image_sets)
extreme_image_sets_unlisted_unique <- unique(extreme_image_sets_unlisted)

# For heatmap_data
extreme_heatmap_sets <- list()

for (colname in colnames(heatmap_data)) {
  # Compute the 5th and 95th percentile for the column
  p05 <- quantile(heatmap_data[, colname], probs = 0.05, na.rm = TRUE)
  p95 <- quantile(heatmap_data[, colname], probs = 0.95, na.rm = TRUE)
  
  # Find rownames (TFs) for bottom 5% or top 5%
  bottom_tf_set <- rownames(heatmap_data)[
    which(heatmap_data[, colname] < p05)
  ]
  top_tf_set <- rownames(heatmap_data)[
    which(heatmap_data[, colname] > p95)
  ]
  
  # Combine TFs from both sets
  tf_set <- c(bottom_tf_set, top_tf_set)
  
  # Store them in the list
  extreme_heatmap_sets[[colname]] <- tf_set
}

# Combine into a single vector and remove duplicates
extreme_heatmap_sets_unlisted <- unlist(extreme_heatmap_sets)
extreme_heatmap_sets_unlisted_unique <- unique(extreme_heatmap_sets_unlisted)

#######  Classify sign usage in IMAGE enhancer motif activity and IMAGE motif activity data
# Row-wise min & max
row_mins_heat_enh <- apply(heatmap_data[extreme_heatmap_sets_unlisted_unique,], 1, min)
row_maxs_heat_enh <- apply(heatmap_data[extreme_heatmap_sets_unlisted_unique,], 1, max)

# Classify each TF
cat_levels <- c("Negative effectors", "Dual function", "Positive effectors")
cat_heat_enh <- factor(
  ifelse(
    row_mins_heat_enh > 0, 
    "Positive effectors",
    ifelse(
      row_maxs_heat_enh < 0,
      "Negative effectors",
      "Dual function"
    )
  ),
  levels = cat_levels
)

# Count how many TFs in each category
tab_heat_enh <- table(cat_heat_enh)
df_heat_enh <- data.frame(
  Category = names(tab_heat_enh),
  Count    = as.vector(tab_heat_enh)
)
df_heat_enh$Category <- factor(df_heat_enh$Category, levels = cat_levels)
library(ggplot2)

gg1 <- ggplot(df_heat_enh, aes(x = Category, y = Count, fill = Category)) +
  geom_col() +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  expand_limits(y = max(df_heat_enh$Count) * 1.2) +
  scale_fill_manual(values = c(
    "Negative effectors" = "#D73027",
    "Dual function" = "#00BA3E",
    "Positive effectors" = "#4575B4"
  )) +
  labs(
    title = "Number of TFs in the 5th and 95th Percentiles of Raw Motif Activity\nGlobally (IMAGE Enhancer)",
    x = NULL,
    y = "Number of TFs"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # Center & bold title
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.text = element_text(face = "bold", colour = "black"),
    text = element_text(face = "bold", colour = "black")
  )


# Row-wise min & max
row_mins_heat_img <- apply(IMAGE_motif_activity[extreme_image_sets_unlisted_unique,], 1, min)
row_maxs_heat_img <- apply(IMAGE_motif_activity[extreme_image_sets_unlisted_unique,], 1, max)

cat_heat_img <- factor(
  ifelse(
    row_mins_heat_img > 0, 
    "Positive effectors",
    ifelse(
      row_maxs_heat_img < 0,
      "Negative effectors",
      "Dual function"
    )
  ),
  levels = cat_levels
)

tab_heat_img <- table(cat_heat_img)
df_heat_img <- data.frame(
  Category = names(tab_heat_img),
  Count    = as.vector(tab_heat_img)
)
df_heat_img$Category <- factor(df_heat_img$Category, levels = cat_levels)
gg2 <- ggplot(df_heat_img, aes(x = Category, y = Count, fill = Category)) +
  geom_col() +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  expand_limits(y = max(df_heat_img$Count) * 1.2) +
  scale_fill_manual(values = c(
    "Negative effectors" = "#D73027",
    "Dual function" = "#00BA3E",
    "Positive effectors" = "#4575B4"
  )) +
  labs(
    title = "Number of TFs in the 5th and 95th Percentiles of Raw Motif Activity\nGlobally (IMAGE Gene)",
    x = NULL,
    y = "Number of TFs"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.text = element_text(face = "bold", colour = "black"),
    text = element_text(face = "bold", colour = "black")
  )
############ number of repressors, activators and dual function for cell types incrementally

cat_levels <- c("Negative effectors", "Dual function", "Positive effectors")
clusters   <- colnames(heatmap_data)
n_clusters <- length(clusters)

df_increments <- data.frame()

for (i in seq_len(n_clusters)) {
  
  # Gather all TFs from columns 1..i in 'extreme_heatmap_sets'
  union_tfs <- unique(unlist(extreme_heatmap_sets[1:i]))
  
  # Subset 'heatmap_data' by these TFs (rows) and columns 1..i
  sub_mat <- heatmap_data[union_tfs, seq_len(i), drop = FALSE]
  
  # Compute row-wise min & max
  row_mins <- apply(sub_mat, 1, min)
  row_maxs <- apply(sub_mat, 1, max)
  
  # Classify each TF
  cat_img <- factor(
    ifelse(
      row_mins > 0, 
      "Positive effectors",
      ifelse(
        row_maxs < 0,
        "Negative effectors",
        "Dual function"
      )
    ),
    levels = cat_levels
  )
  
  # Count how many TFs in each category
  tab_img <- table(cat_img)
  
  # Convert to data frame, store the step i
  df_tmp <- data.frame(
    step     = i,
    Category = names(tab_img),
    Count    = as.vector(tab_img)
  )
  
  # Accumulate in df_increments
  df_increments <- rbind(df_increments, df_tmp)
}
# Now df_increments has a row for each step i and each Category.
library(ggplot2)
df_increments$Category <- factor(df_increments$Category, levels = cat_levels)
gg5 <- ggplot(df_increments, aes(x = factor(step), y = Count, fill = Category)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    title = "Number of TFs in the 5th and 95th Percentiles of Raw Motif Activity\nAcross Incrementally Added Cell Types (IMAGE Enhancer)",
    x = "Number of Cell Types (Incrementally Included)",
    y = "Number of TFs",
    fill = "Activity type"
  ) +
  scale_fill_manual(values = c(
    "Negative effectors" = "#D73027",
    "Dual function" = "#00BA3E",
    "Positive effectors" = "#4575B4"
  )) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    text = element_text(face = "bold", colour = "black")
  )

cat_levels <- c("Negative effectors", "Dual function", "Positive effectors")
clusters   <- colnames(IMAGE_motif_activity)
n_clusters <- length(clusters)

df_increments <- data.frame()

for (i in seq_len(n_clusters)) {
  
  # Union of TFs from columns 1..i
  union_tfs <- unique(unlist(extreme_image_sets[1:i]))
  
  # Subset IMAGE_motif_activity by these TFs and columns 1..i
  sub_mat <- IMAGE_motif_activity[union_tfs, seq_len(i), drop = FALSE]
  
  # Compute row-wise min & max
  row_mins <- apply(sub_mat, 1, min)
  row_maxs <- apply(sub_mat, 1, max)
  
  # Classify each TF
  cat_img <- factor(
    ifelse(
      row_mins > 0, 
      "Positive effectors",
      ifelse(
        row_maxs < 0,
        "Negative effectors",
        "Dual function"
      )
    ),
    levels = cat_levels
  )
  
  # Count how many TFs in each category
  tab_img <- table(cat_img)
  
  # Convert to data frame, store the step i
  df_tmp <- data.frame(
    step     = i,
    Category = names(tab_img),
    Count    = as.vector(tab_img)
  )
  
  # Accumulate
  df_increments <- rbind(df_increments, df_tmp)
}

# Finally, plot
library(ggplot2)
df_increments$Category <- factor(df_increments$Category, levels = cat_levels)
gg6 <- ggplot(df_increments, aes(x = factor(step), y = Count, fill = Category)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    title = "Number of TFs in the 5th and 95th Percentiles of Raw Motif Activity\nAcross Incrementally Added Cell Types (IMAGE Gene)",
    x     = "Number of Cell Types Incrementally Included",
    y     = "Number of TFs",
    fill  = "Activity type"
  ) +
  scale_fill_manual(values = c(
    "Negative effectors" = "#D73027",
    "Dual function" = "#00BA3E",
    "Positive effectors" = "#4575B4"
  )) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    # Bold axis titles
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    text = element_text(face = "bold", colour = "black")
  )


##### for each cell type
library(ggplot2)
library(dplyr)

# Define categories
cat_levels <- c("Negative effectors", "Dual function", "Positive effectors")

# Get cluster names in the original order
clusters <- colnames(heatmap_data)
n_clusters <- length(clusters)

# Create an empty dataframe to store results
df_per_cluster <- data.frame()

# Iterate over each cluster independently
for (i in clusters) {
  sub_mat <- heatmap_data[extreme_heatmap_sets[[i]], i, drop = FALSE]  # Take only 1 column at a time
  
  # Compute min/max for each row
  row_mins <- apply(sub_mat, 1, min)
  row_maxs <- apply(sub_mat, 1, max)
  
  # Classify TFs
  cat_img <- factor(
    ifelse(
      row_mins > 0, "Positive effectors",
      ifelse(row_maxs < 0, "Negative effectors", "Dual function")
    ),
    levels = cat_levels
  )
  
  # Count TFs per category
  tab_img <- table(cat_img)
  
  # Convert to dataframe
  df_tmp <- data.frame(
    Cluster  = i,
    Category = names(tab_img),
    Count    = as.vector(tab_img)
  )
  
  # Append results
  df_per_cluster <- rbind(df_per_cluster, df_tmp)
}

# Ensure clusters appear in their original order
df_per_cluster$Cluster <- factor(df_per_cluster$Cluster, levels = clusters)
df_per_cluster$Category <- factor(df_per_cluster$Category, levels = cat_levels)
# Plot per cluster
gg3 <- ggplot(df_per_cluster, aes(x = Cluster, y = Count, fill = Category)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    title = "Number of TFs in the 5th and 95th Percentiles of Raw Motif Activity\nPer Cell Type (IMAGE Enhancer)",
    x = "Cell Type",
    y = "Number of TFs",
    fill = "Activity type"
  ) +
  scale_fill_manual(values = c(
    "Negative effectors" = "#D73027",
    "Dual function" = "#00BA3E",
    "Positive effectors" = "#4575B4"
  )) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    # Bold axis titles
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    text = element_text(face = "bold", colour = "black")
  )

library(ggplot2)
library(dplyr)

# Define categories
cat_levels <- c("Negative effectors", "Dual function", "Positive effectors")

# Get cluster names in the original order
clusters <- colnames(IMAGE_motif_activity)
n_clusters <- length(clusters)

df_per_cluster <- data.frame()

for (i in clusters) {
  sub_mat <- IMAGE_motif_activity[extreme_image_sets[[i]], i, drop = FALSE]  # Take only 1 column at a time
  
  row_mins <- apply(sub_mat, 1, min)
  row_maxs <- apply(sub_mat, 1, max)
  
  cat_img <- factor(
    ifelse(
      row_mins > 0, "Positive effectors",
      ifelse(row_maxs < 0, "Negative effectors", "Dual function")
    ),
    levels = cat_levels
  )
  
  tab_img <- table(cat_img)
  
  df_tmp <- data.frame(
    Cluster  = i,
    Category = names(tab_img),
    Count    = as.vector(tab_img)
  )
  
  df_per_cluster <- rbind(df_per_cluster, df_tmp)
}

# Ensure clusters appear in their original order
df_per_cluster$Cluster <- factor(df_per_cluster$Cluster, levels = clusters)
df_per_cluster$Category <- factor(df_per_cluster$Category, levels = cat_levels)
# Plot per cluster
gg4 <- ggplot(df_per_cluster, aes(x = Cluster, y = Count, fill = Category)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    title = "Number of TFs in the 5th and 95th Percentiles of Raw Motif Activity\nPer Cell Type (IMAGE Gene)",
    x = "Cell Type",
    y = "Number of TFs",
    fill = "Activity type"
  ) + 
  scale_fill_manual(values = c(
    "Negative effectors" = "#D73027",
    "Dual function" = "#00BA3E",
    "Positive effectors" = "#4575B4"
  )) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    # Bold axis titles
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold", colour = "black"),
    text = element_text(face = "bold", colour = "black")
  )


library(ggplot2)
library(patchwork)

# Arrange the plots in 3 rows, 2 columns
final_plot <- (gg1 + gg2) / (gg3 + gg4) / (gg5 + gg6) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18)
  )
  

# Display the combined plot
print(final_plot)


####################### does any of them change sign from enhancer motif activity to motif activity?
# map numeric values to +1, -1, or 0
clusters <- colnames(IMAGE_motif_activity)
sign2 <- function(x) {
  ifelse(x > 0, 1, ifelse(x < 0, -1, 0))
}

# Desired cluster order
desired_order <- c(
  "hESC 1", "hESC 2", "hESC 3",
  "Neuron like cells", "T/P like cells", "Goblet like cells", 
  "OP like cells", "VE like cells", "FL like cell", 
  "Podocyte like cells", "SMP like cells", 
  "BE like cells", "Adipocyte like cells"
)

library(dplyr)
library(ggplot2)
library(ggrepel)

df_per_cluster_detailed <- data.frame()
df_tf_details <- data.frame()
for (cl in clusters) {
  
  # Intersection of TFs 
  tf_enh <- extreme_heatmap_sets[[cl]]
  tf_img <- extreme_image_sets[[cl]]
  tf_common <- intersect(tf_enh, tf_img)
  
  # If none, store zero rows (remove them later since we don't want 0% slices)
  if (length(tf_common) == 0) {
    df_per_cluster_detailed <- rbind(
      df_per_cluster_detailed,
      data.frame(Cluster=cl, Difference="Same positive effector",         Count=0),
      data.frame(Cluster=cl, Difference="Same negative effector",         Count=0),
      data.frame(Cluster=cl, Difference="Positive effector → Negative effector",  Count=0),
      data.frame(Cluster=cl, Difference="Negative effector → Positive effector",  Count=0)
    )
    next
  }
  
  # Single column values
  enh_vals <- heatmap_data[tf_common, cl]
  img_vals <- IMAGE_motif_activity[tf_common, cl]
  
  # Convert to ±1 or 0
  enh_signs <- sign2(enh_vals)
  img_signs <- sign2(img_vals)
  # Per-TF labeling: did it stay pos, neg, or switched
  library(dplyr)
  difference_labels <- case_when(
    enh_signs ==  1 & img_signs ==  1 ~ "Same positive effector",
    enh_signs == -1 & img_signs == -1 ~ "Same negative effector",
    enh_signs ==  1 & img_signs == -1 ~ "Positive effector → Negative effector",
    enh_signs == -1 & img_signs ==  1 ~ "Negative effector → Positive effector",
    TRUE                              ~ "No change or zero sign"
  )
  # Store the per-TF details in df_tf_details
  df_tf_details <- rbind(
    df_tf_details,
    data.frame(
      TF        = tf_common,
      Cluster   = cl,
      Enh_sign  = enh_signs,
      Gene_sign = img_signs,
      Difference= difference_labels,
      stringsAsFactors = FALSE
    )
  )
  # sign combos
  same_activator       <- sum(enh_signs == 1  & img_signs == 1)
  same_repressor       <- sum(enh_signs == -1 & img_signs == -1)
  switch_to_repressor  <- sum(enh_signs == 1  & img_signs == -1)
  switch_to_activator  <- sum(enh_signs == -1 & img_signs == 1)
  
  # Store results
  df_per_cluster_detailed <- rbind(
    df_per_cluster_detailed,
    data.frame(Cluster=cl, Difference="Same positive effector",         Count=same_activator),
    data.frame(Cluster=cl, Difference="Same negative effector",         Count=same_repressor),
    data.frame(Cluster=cl, Difference="Positive effector → Negative effector",  Count=switch_to_repressor),
    data.frame(Cluster=cl, Difference="Negative effector → Positive effector",  Count=switch_to_activator)
  )
}

library(dplyr)
library(scales)
# compute total, percentage, and a label with Count + Percent
df_plot <- df_per_cluster_detailed %>%
  group_by(Cluster) %>%
  mutate(
    total = sum(Count),
    perc  = Count / total * 100,
    label = paste0(Count, " (", round(perc, 1), "%)")
  ) %>%
  ungroup()

# desired cluster order:
desired_order <- c(
  "hESC 1", "hESC 2", "hESC 3",
  "Neuron like cells", "T/P like cells", "Goblet like cells", 
  "OP like cells", "VE like cells", "FL like cell", 
  "Podocyte like cells", "SMP like cells", 
  "BE like cells", "Adipocyte like cells"
)

library(ggplot2)
library(scales)
df_plot <- df_plot %>%
  filter(Count > 0)
# desired cluster order
df_plot$Cluster <- factor(df_plot$Cluster, levels = desired_order)

bar_change <- ggplot(data = df_plot, aes(x = Cluster, y = Count, fill = Difference)) +
  # Stacked bars from raw counts
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = c("Same positive effector" = "#4575B4",
                               "Same negative effector" = "#D73027",
                               "Positive effector → Negative effector" = "darkgreen",
                               "Negative effector → Positive effector" = "forestgreen")) +
  # Label with precomputed percent column
  geom_text(
    aes(label = label),               
    position = position_fill(vjust=0.5),
    color    = "white",
    size     = 4,
    angle = 90,
    fontface = "bold"
  ) +
  
  # Make the y-axis show 0..100%
  scale_y_continuous(labels = percent_format()) +
  
  # Titles, themes, etc.
  labs(
    title    = "Raw TF Motif Activity Shifts from Enhancers to Genes\n(5th & 95th Percentiles)",
    x        = "Cluster",
    y        = "Percentage of TFs",
    fill     = "Activity Shift"
  ) +
  theme_minimal() +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    axis.text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    # Bold axis titles
    legend.title = element_text(face = "bold")
    
  )

# Which TFs switched from positive to negative in Neuron like cells?
df_tf_details %>%
  filter(Cluster == "Neuron like cells", 
         Difference == "Positive effector → Negative effector")
# Which TFs switch sign in any cluster?
df_tf_details %>%
  filter(Difference %in% c("Positive effector → Negative effector",
                           "Negative effector → Positive effector"))
# Which TFs never changed signs in all clusters
always_same <- df_tf_details %>%
  group_by(TF) %>%
  # Count how many times was each TF was in switched categories.
  summarise(n_switches = sum(Difference %in% c("Positive effector → Negative effector",
                                               "Negative effector → Positive effector"))) %>%
  # Keep TFs with zero sign changes
  filter(n_switches == 0)

print(always_same, n = Inf)

### make it into a heatmap
all_TFs <- unique(df_tf_details$TF)
all_clusters <- desired_order

# Create a dataframe that has one row per TF and Cluster combination
df_all <- expand.grid(TF = all_TFs, Cluster = all_clusters, stringsAsFactors = FALSE) %>%
  left_join(df_tf_details, by = c("TF", "Cluster"))

# Replace any missing difference values with Not present
df_all$Difference[is.na(df_all$Difference)] <- "Not present"

# Order the Cluster factor according to the desired order
df_all$Cluster <- factor(df_all$Cluster, levels = desired_order)

# Define ordering for the Difference category
df_all$Difference <- factor(df_all$Difference,
                            levels = c("Same positive effector", "Same negative effector",
                                       "Positive effector → Negative effector",
                                       "Negative effector → Positive effector",
                                       "Not present")
)
shift_heatmap <- ggplot(df_all, aes(x = Cluster, y = TF, fill = Difference)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("Same positive effector" = "#4575B4",
                               "Same negative effector" = "#D73027",
                               "Positive effector → Negative effector" = "darkgreen",
                               "Negative effector → Positive effector" = "forestgreen",
                               "Not present" = "grey90")) +
  labs(title = "TF Sign Change Patterns Across Clusters\n(5th & 95th Percentiles)",
       x = "Cluster",
       y = "TF",
       fill = "Activity Shift") +
  theme_minimal() +
  theme_bw(base_size = 14) +
  theme(text = element_text(face = "bold", colour = "black"),
        axis.text = element_text(face = "bold", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right",
        # Bold axis titles
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

library(ggplot2)
library(patchwork)

# Wrap for patchwork compatibility
bar_wrapped <- wrap_elements(full = bar_change)
shifts_wrapped <- wrap_elements(full = shift_heatmap)

# Combine with the bottom row
final_plot <- (bar_wrapped | shifts_wrapped) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18))

final_plot

########### motif activities correlation with their own expression
CPM <- readRDS("/work/Home/Data_for_masters_project/CPM_chromVAR_motifs_GC_and_edgeR")
gene_expression <- CPM[, 5:56]

# Extract condition labels from column names
condition_labels <- sub(".*-Condition_", "Condition_", colnames(gene_expression))

# Create an empty matrix to store results
unique_conditions <- unique(condition_labels)
gene_expression_means <- matrix(NA, nrow = nrow(gene_expression), ncol = length(unique_conditions))

# Compute row means for each condition using a loop
for (i in seq_along(unique_conditions)) {
  condition <- unique_conditions[i]
  cols <- which(condition_labels == condition)  # Get columns for the condition
  gene_expression_means[, i] <- rowMeans(gene_expression[, cols, drop = FALSE])
}
rownames(gene_expression_means) <- CPM$Factor

# find the intersect of TFs and subset based on those
common_TFs <- intersect(rownames(IMAGE_motif_activity), rownames(gene_expression_means))
common_TFs <- intersect(common_TFs, extreme_image_sets_unlisted_unique)
IMAGE_motif_activity_subset <- IMAGE_motif_activity[common_TFs,]
gene_expression_means_subset <- gene_expression_means[common_TFs,]

# Compute self-correlations (each TF's own motif activity vs. its own gene expression)
self_correlations <- sapply(rownames(IMAGE_motif_activity_subset), function(tf) {
  cor(IMAGE_motif_activity_subset[tf, ], gene_expression_means_subset[tf, ], method = "pearson")
})

# Convert to a data frame
self_correlation_df <- data.frame(Correlation = self_correlations)

# Histogram of self-correlations
gg_self_correlation <- ggplot(self_correlation_df, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = "#4575B4", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of Self-Correlations (Motif Activity vs Gene Expression)",
       x = "Pearson Correlation",
       y = "Frequency") +
  theme_bw(base_size = 20) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    axis.text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

#### example plot
# Load required libraries
library(ggplot2)
library(dplyr)

# Compute self-correlations (each TF's own motif activity vs. its own gene expression)
self_correlations <- sapply(rownames(IMAGE_motif_activity_subset), function(tf) {
  cor(IMAGE_motif_activity_subset[tf, ], gene_expression_means_subset[tf, ], method = "pearson")
})

# Convert the correlations into a data frame for easier manipulation
cor_df <- data.frame(TF = names(self_correlations), Correlation = self_correlations)

# Identify TFs with correlation closest to -1, 0, and 1, or find specific TFs of interest
# Sort the TFs by correlation to find the best matches
cor_df_sorted <- cor_df[order(cor_df$Correlation), ]

# Select the TF with the **lowest** correlation (closest to -1)
tf_neg1 <- "RBPJ"#cor_df_sorted$TF[which.min(cor_df_sorted$Correlation)]

# Select the TF with the **highest** correlation (closest to 1)
tf_1 <- "POU5F1"#cor_df_sorted$TF[which.max(cor_df_sorted$Correlation)]

# Select the TF with the correlation **closest to 0** (no correlation)
tf_0 <- "TBX3"#cor_df_sorted$TF[which.min(abs(cor_df_sorted$Correlation))]


# Create a new data frame for plotting
plot_data <- list()

# Function to extract and format data for a given TF
extract_tf_data <- function(tf, correlation_value) {
  data.frame(
    Motif_Activity = IMAGE_motif_activity_subset[tf, ],
    Gene_Expression = gene_expression_means_subset[tf, ],
    TF_Label = paste0(tf, " (Pearson Correlation = ", round(correlation_value, 2), ")")  # Label for facet
  )
}

# Get data for the three selected TFs
plot_data[[1]] <- extract_tf_data(tf_neg1, self_correlations[tf_neg1])
plot_data[[2]] <- extract_tf_data(tf_0, self_correlations[tf_0])
plot_data[[3]] <- extract_tf_data(tf_1, self_correlations[tf_1])

# Combine into a single data frame for plotting
plot_data_df <- bind_rows(plot_data)


library(viridis)
# Create the scatter plot 
gg_example_scatter <- ggplot(plot_data_df, aes(x = Motif_Activity, y = Gene_Expression, color = Motif_Activity)) +
  geom_point(alpha = 0.7, size = 2) +  # Scatter plot points with size and transparency
  geom_smooth(method = "lm", color = "#440154FF", se = FALSE, linewidth = 1.2) +  # Dark purple regression line
  scale_color_viridis(option = "magma", direction = 1) +  # Use a modern color gradient
  facet_wrap(~TF_Label, scales = "free") +  # Create separate plots for each TF
  theme_minimal() +
  labs(
    title = "Examples of Motif Activity and Gene Expression Self-Correlations",
    x = "Motif Activity",
    y = "Gene Expression",
    color = "Motif Activity"  # Color legend title
  ) +
  theme_bw(base_size = 20) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    axis.text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"  
  )

################################### Do the correlation for target genes
target_genes <- readRDS("/work/Home/Data_for_masters_project/target_genes_chromVAR_motifs_GC_and_edgeR")
MotifActivities <- target_genes[["MotifActivities"]]
# Filter rows where CausalTF column contains 1 or 2
subset_MotifActivities <- MotifActivities[MotifActivities$CausalTF %in% c(0, 1, 2), ]
# Extract motifs from the MotifActivities data frame
motifs <- subset_MotifActivities$Motif

# Append ".motif" to match with the gene list names
extreme_image_sets_unlisted_unique_motifs <- paste0(extreme_image_sets_unlisted_unique, ".motif")
subset_motifs <- motifs[motifs %in% extreme_image_sets_unlisted_unique_motifs]
# Get the names of the target_Genes list
target_genes_names <- names(target_genes[["target_Genes"]])

# Find matching motifs
matching_motifs <- subset_motifs[subset_motifs %in% target_genes_names]

# Subset target_Genes based on matching motifs (get target genes for the motifs)
filtered_target_Genes <- target_genes[["target_Genes"]][matching_motifs]

# Filter the list to include only elements where the "Target" column has at least one '1'
filtered_target_Genes_subset <- lapply(filtered_target_Genes, function(df) {
  subset_df <- df[df$Target == 1, ]
  if (nrow(subset_df) > 0) {
    return(subset_df)
  } else {
    return(NULL)
  }
})

IMAGE_motif_activity_correlation <- IMAGE_motif_activity[extreme_image_sets_unlisted_unique,]
rownames(IMAGE_motif_activity_correlation) <- extreme_image_sets_unlisted_unique_motifs

# Initialize an empty list to store results
motif_gene_correlations <- list()

# Loop over each motif
for (motif in rownames(IMAGE_motif_activity_correlation)) {
  
  # Extract motif activity for the current motif
  motif_activity <- as.numeric(IMAGE_motif_activity_correlation[motif, ])
  
  # Get target genes for the current motif
  if (!is.null(filtered_target_Genes_subset[[motif]])) {
    target_genes <- rownames(filtered_target_Genes_subset[[motif]])
  } else {
    next  # Skip if no target genes
  }
  
  # Find common genes between target genes and expression data
  common_genes <- intersect(target_genes, rownames(gene_expression_means))
  
  if (length(common_genes) == 0) next  # Skip motifs with no matching genes
  
  # Subset gene expression data for common genes
  target_gene_expression <- gene_expression_means[common_genes, ]
  
  # Compute correlations
  correlations <- apply(target_gene_expression, 1, function(expr_values) {
    cor(expr_values, motif_activity, method = "pearson")
  })
  
  # Store results in a dataframe
  motif_gene_correlations[[motif]] <- data.frame(
    TF = motif,
    Gene = common_genes,
    Correlation = correlations
  )
}

# Combine results into a single dataframe
final_correlation_results <- do.call(rbind, motif_gene_correlations)

# Sort by absolute correlation strength
final_correlation_results_sorted <- final_correlation_results[order(-abs(final_correlation_results$Correlation)), ]

# Inspect the top correlations
head(final_correlation_results_sorted)

# remove .motif
final_correlation_results_sorted$TF <- gsub("\\.motif$", "", final_correlation_results_sorted$TF)
rownames(IMAGE_motif_activity_correlation) <- gsub("\\.motif$", "", rownames(IMAGE_motif_activity_correlation))

#### scatterplots of 3 examples
## select TFs of interest or the highest, lowest or 0 correlation
# Select the TF with the **lowest** correlation (closest to -1)
tf_neg1 <- "RBPJ.motif.SNAI2"#rownames(final_correlation_results_sorted)[which.min(final_correlation_results_sorted$Correlation)]#"NANOG"#

# Select the TF with the **highest** correlation (closest to 1)
tf_1 <- "POU5F1.motif.DPPA4"#rownames(final_correlation_results_sorted)[which.max(final_correlation_results_sorted$Correlation)] #"POU5F1"#

# Select the TF with the correlation **closest to 0** (no correlation)
tf_0 <- "TBX3.motif.CRHR2"#rownames(final_correlation_results_sorted)[which.min(abs(final_correlation_results_sorted$Correlation))]#"TBX3"#

# subset the gene expression to match with the genes that are found to be targets
gene_expression_means_subset <- gene_expression_means[final_correlation_results_sorted$Gene,]
colnames(gene_expression_means_subset) <- colnames(IMAGE_motif_activity_correlation)
tf_list <- c(tf_neg1, tf_1, tf_0)

# Extract TF and Gene names
split_names <- strsplit(tf_list, "\\.motif\\.")
TF_names    <- sapply(split_names, `[`, 1)
Gene_names  <- sapply(split_names, `[`, 2)

# Get corresponding correlation values
cor_values <- final_correlation_results_sorted[tf_list, "Correlation"]

# Create a clear data frame
tf_gene_corr_df <- data.frame(TF = TF_names, Gene = Gene_names, Correlation = cor_values)

library(dplyr)
library(tidyr)

# Extract motif activities for each TF
motif_activity_subset <- IMAGE_motif_activity_correlation[TF_names, ]

# Extract gene expression for each Gene
gene_expression_subset <- gene_expression_means_subset[Gene_names, ]

# Tidy format
motif_df <- as.data.frame(motif_activity_subset) %>%
  mutate(TF = TF_names) %>%
  pivot_longer(-TF, names_to = "Condition", values_to = "Motif_Activity")

gene_expr_df <- as.data.frame(gene_expression_subset) %>%
  mutate(Gene = Gene_names) %>%
  pivot_longer(-Gene, names_to = "Condition", values_to = "Gene_Expression")

# Merge data
plot_data_df <- motif_df %>%
  left_join(tf_gene_corr_df, by = "TF") %>%
  left_join(gene_expr_df, by = c("Gene", "Condition"))

# Add a label with TF, Gene, and correlation
plot_data_df$TF_Label <- paste0(plot_data_df$TF, " - ", plot_data_df$Gene, 
                                "\nPearson Correlation = ", round(plot_data_df$Correlation, 3))


# Load required library for better color palettes
library(viridis)

# Create the scatter plot 
gg_motif_gene_scatter <- ggplot(plot_data_df, aes(x = Motif_Activity, y = Gene_Expression, color = Motif_Activity)) +
  geom_point(alpha = 0.7, size = 2) +  # Scatter plot points with size and transparency
  geom_smooth(method = "lm", color = "#440154FF", se = FALSE, linewidth = 1.2) +  # Dark purple regression line
  scale_color_viridis(option = "magma", direction = 1) +  # Use a modern color gradient
  facet_wrap(~TF_Label, scales = "free") +  # Create separate plots for each TF
  theme_minimal() +
  labs(
    title = "Examples of Motif Activity and TF-Specific Target Gene Correlations",
    x = "Raw Motif Activity",
    y = "Gene Expression",
    color = "Raw Motif Activity"  # Color legend title
  ) +
  theme_bw(base_size = 20) + 
  theme(
    text = element_text(face = "bold", colour = "black"),
    axis.text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"  
  )

#### histogram of the correlations
gg_motif_gene_correlation <- ggplot(final_correlation_results_sorted, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.05, fill = "#4575B4", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of Motifs and Target Gene Correlations (Motif Activity vs Gene Expression)",
       x = "Pearson Correlation",
       y = "Frequency") +
  theme_bw(base_size = 20) + 
  theme(
    text = element_text(face = "bold", colour = "black"),
    axis.text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )
library(patchwork)

# Adjust relative heights and widths explicitly
combined_plot1 <- gg_example_scatter/gg_self_correlation +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

combined_plot2 <- gg_motif_gene_scatter/gg_motif_gene_correlation +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
# Display combined plot
combined_plot1
combined_plot2
################
# IMAGE gene regulatory network processing
###############
target_enhancers <- readRDS("/work/Home/Data_for_masters_project/target_enhancers_chromVAR_motifs_GC_and_edgeR")
target_genes <- readRDS("/work/Home/Data_for_masters_project/target_genes_chromVAR_motifs_GC_and_edgeR")
normEnhancers <- readRDS("/work/Home/Data_for_masters_project/normEnhancers_chromVAR_motifs_4_replicates_all_peaks_GC_and_DA")


MotifActivities <- target_genes[["MotifActivities"]]
subset_MotifActivities <- MotifActivities#MotifActivities[MotifActivities$CausalTF %in% c(1, 2), ]
# Extract motifs from the MotifActivities data frame
motifs <- subset_MotifActivities$Motif
rownames(subset_MotifActivities) <- subset_MotifActivities$Motif

# Get the names of the target_Genes list
target_genes_names <- names(target_genes[["target_Genes"]])

# Find matching motifs
matching_motifs <- motifs[motifs %in% target_genes_names]

# Subset target_Genes based on matching motifs (get target genes for the motifs)
filtered_target_Genes_subset <- target_genes[["target_Genes"]][matching_motifs]
subset_MotifActivities <- subset_MotifActivities[names(filtered_target_Genes_subset),]
filtered_target_enhancers_subset <- target_enhancers[["data"]][matching_motifs]

ATAC <- readRDS("/work/Home/Data_for_masters_project/ATAC_summarizedExperimentObject_TFatlas_bulk_for_seurat_clusters_4_replicates_all peaks")
CPM <- readRDS("/work/Home/Data_for_masters_project/CPM_chromVAR_motifs_GC_and_edgeR")

# Make variable that contains relevent information regarding the transcription start site
TSS <- rbind(CPM[,c("Chr", "TSS", "Factor")])
TSS <- TSS[order(TSS$Chr), ]
rownames(TSS) <- TSS$Factor
TSS$Start <- TSS$TSS - 100000
TSS$End <- TSS$TSS + 100000
# Subject hits is used to give an "identifier" to each gene/row
TSS$subjectHits <- seq(1, nrow(TSS), by=1)
# Make a gene ranges object, that  contains relevent information about the TSS
GenesRanges <- GRanges(seqnames=TSS[,1], IRanges(start=TSS$Start, end=TSS$End), strand=rep("+", nrow(TSS)), names=TSS[,3])

## Make a GRange object from the enhancers
# Quiry hits is used to give an "identifier" to each ATAC seq peak ID
RegionRanges <- ATAC@rowRanges
Regions = as.data.frame(ATAC@rowRanges)
Regions$queryHits = seq(1, nrow(Regions),1)

## find gene regions and enhancers (open) regions that overlap
Overlap <- findOverlaps(query=RegionRanges, subject=GenesRanges, type="any")
Overlap <- as.data.frame(Overlap)

## Merge the appropriate information from gene and enhancers
Overlap <- merge(Overlap, Regions[,c(1:3,6,ncol(Regions))], by="queryHits")
Overlap <- merge(Overlap, TSS[,c("Factor","TSS","subjectHits")], by="subjectHits")

## Calculate scaled regulatory potential
Overlap$Center <- (Overlap[,5]+Overlap[,4])/2
Overlap$Distance <- Overlap$Center - Overlap$TSS
Overlap <- Overlap[ abs(Overlap$Distance) <= 100000,]
## Consider alternative formulas
Overlap$Potential <- (exp(-(0.5+4*(abs(Overlap$Distance)/100000)))-exp(-1 * (0.5 + 4)))/(max((exp(-1 * (0.5 + (4 * (seq(0,100000,by=10)/100000))))-exp(-1 * (0.5 + 4)))))


############################## Get all enhancers to each gene
# find all unique target genes
unique_target_genes <- list()
for (motif in names(filtered_target_Genes_subset)){
  unique_target_genes[[motif]] <- rownames(filtered_target_Genes_subset[[motif]])
}

unique_target_genes <- unlist(unique_target_genes)
unique_target_genes <- unique(unique_target_genes)

# subset overlap based on it
Overlap_subset <- Overlap[Overlap$Factor %in% unique_target_genes,]

library(dplyr)

# Initialize an empty list to store gene-specific data
gene_list <- list()

# Loop through motifs
for (motif in names(filtered_target_Genes_subset)) {
  
  # Find target genes for the motif
  gene_rows <- filtered_target_Genes_subset[[motif]]
  
  # Find target enhancers for the motif
  enhancer_rows <- filtered_target_enhancers_subset[[motif]]$PeakID
  
  # Loop through genes
  for (gene in unique_target_genes) {
    
    # Check if the gene exists in the rownames of gene_rows
    if (gene %in% rownames(gene_rows)) {
      
      # Subset Overlap data for the specific gene
      Overlap_subset_gene_specific <- Overlap_subset[Overlap_subset$Factor %in% gene, ]
      
      # Add motif information
      Overlap_subset_gene_specific$Motif <- motif
      
      # find motif gene contribution score
      motifs_gene_contribution_score <- gene_rows[gene,"V2"] # V2 is contribution score column
      
      # insert gene motif contribution score
      Overlap_subset_gene_specific$Motif_gene_contribution_score <- motifs_gene_contribution_score
      
      # is the gene predicted to be a target gene?
      Target_gene_prediction <- gene_rows[gene,"Target"]
      
      # insert it into object
      Overlap_subset_gene_specific$Target_gene_prediction <- Target_gene_prediction
      
      if (any(Overlap_subset_gene_specific$peak_id %in% enhancer_rows)) {
        
        # Filter based on enhancer rows
        Overlap_subset_gene_specific <- Overlap_subset_gene_specific[Overlap_subset_gene_specific$peak_id %in% enhancer_rows, ]
        
        # Find and get enhancer motif contribution score
        specific_enhancers <- Overlap_subset_gene_specific$peak_id
        motifs_enhancer_contrib <- filtered_target_enhancers_subset[[motif]][
          filtered_target_enhancers_subset[[motif]]$PeakID %in% specific_enhancers,
        ]$contribution
        
        # insert them in that order
        Overlap_subset_gene_specific$Motifs_enhancer_contribution_score <- motifs_enhancer_contrib
        
        # is the enhancer predicted to be a target enhancer?
        Target_enhancer_prediction <- filtered_target_enhancers_subset[[motif]][
          filtered_target_enhancers_subset[[motif]]$PeakID %in% specific_enhancers,
        ]$target
        
        # insert them into object
        Overlap_subset_gene_specific$Target_enhancer_prediction <- Target_enhancer_prediction
        
        # If the gene already exists in the list, bind new rows
        if (gene %in% names(gene_list)) {
          gene_list[[gene]] <- bind_rows(gene_list[[gene]], Overlap_subset_gene_specific)
        } else {
          gene_list[[gene]] <- Overlap_subset_gene_specific
        }
      }
    }
  }
}

#saveRDS(gene_list,"/work/Home/Data_for_masters_project/gene_list_2_for_locus_plot_IMAGE")

###################
# setup for scE2G
###################
library(Signac)
library(Seurat)
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")

# Paths
fragment_file <- "/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_fragments.tsv.gz"  # Path to the input fragment file
output_dir <- "/work/Home/Data_for_masters_project/scE2G_cluster_files/"  # Directory to save pseudobulk files

# Read the fragment file
fragment_data <- read.table(gzfile(fragment_file), header = FALSE, stringsAsFactors = FALSE)
colnames(fragment_data) <- c("chr", "start", "end", "cell_name", "read_count")

# Read metadata from Seurat object
metadata <- SHARE_seq_subsample_data@meta.data

# Get unique clusters from metadata
clusters <- unique(metadata$seurat_clusters)

# Process fragments for each cluster
for (cluster in clusters) {
  # Subset metadata for the current cluster
  cluster_metadata <- metadata[metadata$seurat_clusters == cluster, ]
  
  # Get the cell barcodes for the current cluster
  cluster_cells <- rownames(cluster_metadata)
  
  # Subset fragment_data for these cell barcodes
  cluster_fragments <- fragment_data[fragment_data$cell_name %in% cluster_cells, ]
  gc()
  # Save the pseudobulk fragment file
  output_file <- paste0(output_dir, "atac_fragments_cluster_", cluster, ".tsv")
  write.table(cluster_fragments, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Compress and index the file
  # Compress and index the file
  system(paste0("sort -k1,1 -k2,2n \"", output_file, "\" | bgzip > \"", output_file, ".gz\""))
  system(paste("tabix -p bed", paste0("\"", output_file, ".gz\"")))
  
}

# Process RNA data for each cluster
for (cluster in clusters) {
  # Subset metadata for the current cluster
  cluster_metadata <- metadata[metadata$seurat_clusters == cluster, ]
  
  # Get the cell barcodes for the current cluster
  cluster_cells <- rownames(cluster_metadata)
  
  # Subset RNA counts for these cell barcodes
  cluster_RNA <- SHARE_seq_subsample_data@assays$RNA@counts[, cluster_cells]
  
  # Convert to data frame for better readability in CSV
  cluster_RNA_df <- as.data.frame(as.matrix(cluster_RNA))
  
  # Save the file as .csv.gz
  output_file <- paste0(output_dir, "RNA_cluster_", cluster, ".csv.gz")
  write.csv(cluster_RNA_df, gzfile(output_file), row.names = TRUE)  # Save with row names (gene names)
  
  # Free up memory
  gc()
}
# What I did in BASH:
# git clone --recurse-submodules https://github.com/EngreitzLab/scE2G.git
# cd scE2G
# git config --global submodule.recurse true
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# bash miniconda.sh
# source ~/miniconda3/bin/activate
# conda config --set channel_priority flexible
# conda create -n mamba -c conda-forge mamba=1.5.11 "python<=3.11"
# conda activate mamba
# mamba install -c conda-forge -c bioconda snakemake=7

# To configure the pipeline:
# nano /work/scE2G/config/config.yaml
# nano /work/scE2G/config/config_cell_clusters.tsv
# snakemake -s workflow/Snakefile -j1 --use-conda


# verify before deletion:  ls -lh /work
# remove scE2G:  rm -rf /work/scE2G/
# verify after deletion:  ls -lh /work
# Reinstall scE2G
#  cd /work
# cp /work/Home/config_cell_clusters.tsv /work/Home/config_models.tsv /work/Home/config_training.yaml /work/Home/config.yaml /work/scE2G/config/

##########################
# locus plots and confusion matrices and recall and precision plot
#######################
# Memory heavy to combine scE2G for all clusters
library(readr)
library(tidyr)
library(reshape2)
# Read the gzipped TSV file and make a combined enhancer-gene link dataframe for scE2G
sce2g_0 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_0/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_1 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_1/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_combined <- rbind(sce2g_0, sce2g_1)
rm(sce2g_0)
rm(sce2g_1)
gc()
sce2g_2 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_2/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_3 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_3/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_4 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_4/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_combined <- rbind(sce2g_combined, sce2g_2, sce2g_3, sce2g_4)
rm(sce2g_2)
rm(sce2g_3)
rm(sce2g_4)
gc()
sce2g_5 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_5/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_6 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_6/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_7 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_7/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_combined <- rbind(sce2g_combined, sce2g_5, sce2g_6, sce2g_7)
rm(sce2g_5)
rm(sce2g_6)
rm(sce2g_7)
gc()
sce2g_8 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_8/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_combined <- rbind(sce2g_combined, sce2g_8)
rm(sce2g_8)
gc()
sce2g_9 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_9/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_combined <- rbind(sce2g_combined, sce2g_9)
rm(sce2g_9)
gc()
sce2g_10 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_10/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
sce2g_combined <- rbind(sce2g_combined, sce2g_10)
rm(sce2g_10)
gc()
#saveRDS(sce2g_combined,"/work/Home/Data_for_masters_project/sce2g_combined")
sce2g_combined <- readRDS("/work/Home/Data_for_masters_project/sce2g_combined")

sce2g_11 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_11/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
gc()
sce2g_combined <- rbind(sce2g_combined, sce2g_11)
rm(sce2g_11)
gc()
#saveRDS(sce2g_combined,"/work/Home/Data_for_masters_project/sce2g_combined")
sce2g_combined <- readRDS("/work/Home/Data_for_masters_project/sce2g_combined")
sce2g_12 <- read_tsv("/work/Home/Data_for_masters_project/scE2G_results/Cluster_12/multiome_powerlaw_v2/encode_e2g_predictions.tsv.gz")
library(data.table)
# Convert data frames to data.table for efficient binding. I had memory issues with this step.
# so I tried data.table and it worked
setDT(sce2g_combined)
setDT(sce2g_12)
# Efficiently combine large data
sce2g_combined <- rbindlist(list(sce2g_combined, sce2g_12), use.names = TRUE, fill = TRUE)
# Remove unused object and free memory
rm(sce2g_12)
gc()
sce2g_combined <- as.data.frame(sce2g_combined)
gc()
#saveRDS(sce2g_combined,"/work/Home/Data_for_masters_project/sce2g_combined")

library(GenomicRanges)
library(dplyr)

# load enhancer-gene links for IMAGE
gene_list <- readRDS("/work/Home/Data_for_masters_project/gene_list_2_for_locus_plot_IMAGE")
# load enhancer-gene links for scE2G
sce2g_combined <- readRDS("/work/Home/Data_for_masters_project/sce2g_combined")
sce2g_combined <- sce2g_combined %>%
  mutate(Coordinate = paste(chr, start, end, sep = "-"))
# subset both based on common genes
sce2g_combined_subset <- sce2g_combined[sce2g_combined$TargetGene %in% names(gene_list), ]
gene_list_subset <- gene_list[names(gene_list) %in% sce2g_combined_subset$TargetGene]
#saveRDS(sce2g_combined_subset,"/work/Home/Data_for_masters_project/sce2g_combined_subset_2_based_on_genes")
#saveRDS(gene_list_subset,"/work/Home/Data_for_masters_project/gene_list_subset_2_based_on_genes_for_locus_plot_IMAGE")

# Loaddata
sce2g <- readRDS("/work/Home/Data_for_masters_project/sce2g_combined_subset_2_based_on_genes")
gene_list <- readRDS("/work/Home/Data_for_masters_project/gene_list_subset_2_based_on_genes_for_locus_plot_IMAGE")

# Binarize IMAGE prediction
df_image <- bind_rows(gene_list, .id = "Gene")
df_image$is_IMAGE <- ifelse(df_image$Target_gene_prediction == 1 & 
                              df_image$Target_enhancer_prediction == 1, 1, 0)
rm(gene_list)
gc()
# Binarize scE2G prediction using threshold on E2G.Score.qnorm. The threshold is according to the paper
sce2g$is_sce2g <- ifelse(sce2g$E2G.Score.qnorm >= 0.164, 1, 0)
# subset sce2g to only containing rows that are withing 100 kb of gene
sce2g <- sce2g[sce2g$distance <= 1e5, ]
gc()
# Keep only the shared genes 
shared_genes <- intersect(unique(df_image$Gene), unique(sce2g$TargetGene))
df_image <- df_image[df_image$Gene %in% shared_genes, ]
sce2g   <- sce2g[sce2g$TargetGene %in% shared_genes, ]

# setup for collect row indices that pass the 60% overlap.
df_image$row_idx <- as.numeric(rownames(df_image))
sce2g$row_idx    <- as.numeric(rownames(sce2g))

# Split into lists keyed by gene, for faster comparisons
image_by_gene <- split(df_image, df_image$Gene)
sce2g_by_gene <- split(sce2g,   sce2g$TargetGene)

rows_to_keep_image <- integer(0)
rows_to_keep_sce2g <- integer(0)
# counter to keep track of the process.
counter <- 0
# Loop over shared genes, find intervals with ≥40% overlap
for (gene in shared_genes) {
  image_gene <- image_by_gene[[gene]]
  sce_gene   <- sce2g_by_gene[[gene]]
  
  # Build GRanges for IMAGE
  gr_image <- GRanges(
    seqnames = image_gene$seqnames,
    ranges   = IRanges(start = image_gene$start, end = image_gene$end)
  )
  mcols(gr_image)$image_row_idx <- image_gene$row_idx
  
  # Build GRanges for scE2G
  gr_sce <- GRanges(
    seqnames = sce_gene$chr,
    ranges   = IRanges(start = sce_gene$start, end = sce_gene$end)
  )
  mcols(gr_sce)$sce2g_row_idx <- sce_gene$row_idx
  
  # Find all overlaps
  ov <- findOverlaps(gr_image, gr_sce)
  if (length(ov) > 0) {
    # For each overlapping pair, compute fraction of scE2G range that is overlapped
    q <- gr_image[queryHits(ov)]
    s <- gr_sce[subjectHits(ov)]
    overlap_fraction <- width(pintersect(q, s)) / width(s)
    
    # Indices of the pairs passing >= 40% overlap
    keep_idx <- which(overlap_fraction >= 0.4)
    
    if (length(keep_idx) > 0) {
      # Which row indices from each data set do these correspond to?
      keep_rows_image <- mcols(q)$image_row_idx[keep_idx]
      keep_rows_sce2g <- mcols(s)$sce2g_row_idx[keep_idx]
      
      rows_to_keep_image <- c(rows_to_keep_image, keep_rows_image)
      rows_to_keep_sce2g <- c(rows_to_keep_sce2g, keep_rows_sce2g)
    }
  }
  print(counter)
  counter <- counter + 1
}


# Subset original data to the intervals that pass 40% overlap
rows_to_keep_image <- unique(rows_to_keep_image)
rows_to_keep_sce2g <- unique(rows_to_keep_sce2g)

df_image_40pct <- df_image[df_image$row_idx %in% rows_to_keep_image, ]
sce2g_40pct    <- sce2g[sce2g$row_idx %in% rows_to_keep_sce2g, ]

#saveRDS(df_image_40pct, "/work/Home/Data_for_masters_project/df_image_40pct_overlap")
#saveRDS(sce2g_40pct, "/work/Home/Data_for_masters_project/sce2g_40pct_overlap")

####################################### Make locus plots. Visualize IMAGE and sce2g enhancer-gene links
df_image_40pct <- readRDS("/work/Home/Data_for_masters_project/df_image_40pct_overlap")
sce2g_40pct <- readRDS("/work/Home/Data_for_masters_project/sce2g_40pct_overlap")
# subset based on only those that are predicted to be positive links
sce2g_40pct <- sce2g_40pct[sce2g_40pct$is_sce2g == 1,]

library(Signac)
library(GenomicRanges)
library(dplyr)
SHARE_seq_subsample_data <- readRDS("/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_subsample.rds")
SHARE_seq_subsample_data = SetIdent(SHARE_seq_subsample_data, value = SHARE_seq_subsample_data$seurat_clusters)
new.cluster.ids <- c("hESC 1", "hESC 2", "hESC 3",
                     "Neuron like cells", "T/P like cells", "Goblet like cells", 
                     "OP like cells", "VE like cells",
                     "FL like cell", "Podocyte like cells", "SMP like cells", 
                     "BE like cells", "Adipocyte like cells")
names(new.cluster.ids) <- levels(SHARE_seq_subsample_data)
SHARE_seq_subsample_data <- RenameIdents(SHARE_seq_subsample_data, new.cluster.ids)
OE = subset(SHARE_seq_subsample_data, TF != "TFORF3550-mCherry")
OE = subset(OE, TF != "TFORF3549-GFP")
rm(SHARE_seq_subsample_data)
gc()
# Build a Seurat object with an "ATAC" assay
chrom_assay <- CreateChromatinAssay(
  counts    = OE@assays[["ATAC"]]@data,
  sep       = c("-", "-"),  
  genome    = "hg38",       
  fragments = "/work/Home/Data_for_masters_project/GSE217215_201218_ATAC_fragments.tsv.gz",
)
atac_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay  = "ATAC",
  meta.data = OE@meta.data
)
DefaultAssay(atac_obj) <- "ATAC"
gc()

###  Attach a gene annotation for hg38 (EnsDb.Hsapiens.v86 supports hg38)
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# Convert to UCSC style ("chr1" naming)
seqlevelsStyle(annotation) <- "UCSC"

# Attach to the ChromatinAssay
Annotation(atac_obj[["ATAC"]]) <- annotation

# Define the gene of interest
gene_of_interest <- "DPPA4"
# Get GRanges object
gene_coords <- LookupGeneCoords(atac_obj, gene_of_interest)

# Extract seqnames, start, and end
gene_region <- paste0(seqnames(gene_coords), "-", start(gene_coords), "-", end(gene_coords))

# subset for gene of interest
sce2g_combined_subset_gene <- sce2g_40pct[sce2g_40pct$TargetGene == gene_of_interest,]
# create links between enhancer and gene. Using gene TSS as anchor
my_links <- sce2g_combined_subset_gene %>%
  rename(
    enhancer_start = start,
    enhancer_end   = end
  ) %>%
  mutate(
    seqnames = chr,
    score    = E2G.Score.qnorm,
    start    = pmin(TargetGeneTSS, enhancer_start),
    end      = pmax(TargetGeneTSS, enhancer_end)
  )

# convert to a GRanges
gr_links <- makeGRangesFromDataFrame(
  df = my_links,
  seqnames.field = "seqnames",
  start.field    = "start",
  end.field      = "end",
  keep.extra.columns = TRUE
)

Links(atac_obj) <- gr_links
DefaultAssay(atac_obj) <- "ATAC"

library(Seurat)
library(Signac)
library(patchwork)

# Generate Coverage Plot
cov_plot <- CoveragePlot(
  object     = atac_obj,
  region     = gene_region,
  annotation = TRUE,
  peaks      = TRUE,
  links      = TRUE,
  expression.assay = NULL,
  features   = NULL,
  extend.upstream = 100000,
  extend.downstream = 100000
) &
  theme(
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  )
cov_plot[[1]][[4]] <- cov_plot[[1]][[4]] &
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  )


# Add a title
Locus_plot_sce2g <- cov_plot + plot_annotation(
  title = bquote(bold("Chromatin Accessibility and Enhancer–Gene Links for " ~ .(gene_of_interest) ~
                        " from scE2G (Region: " ~ .(gene_region) ~ " ± 100 kb)")),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 18))
)

### locus plot for IMAGE
# Define the gene of interest
gene_of_interest <- "DPPA4"
# Get GRanges object
gene_coords <- LookupGeneCoords(atac_obj, gene_of_interest)

# Extract seqnames, start, and end
gene_region <- paste0(seqnames(gene_coords), "-", start(gene_coords), "-", end(gene_coords))

# subset for gene of interest
IMAGE_subset_gene <- df_image_40pct[df_image_40pct$Gene == gene_of_interest,]
# create links between enhancer and gene. Using gene TSS as anchor
my_links <- IMAGE_subset_gene %>%
  rename(
    enhancer_start = start,
    enhancer_end   = end
  ) %>%
  mutate(
    seqnames = seqnames,
    score    = Potential,
    start    = pmin(TSS, enhancer_start),
    end      = pmax(TSS, enhancer_end)
  )

# convert to a GRanges
gr_links <- makeGRangesFromDataFrame(
  df = my_links,
  seqnames.field = "seqnames",
  start.field    = "start",
  end.field      = "end",
  keep.extra.columns = TRUE
)

Links(atac_obj) <- gr_links
DefaultAssay(atac_obj) <- "ATAC"

library(Seurat)
library(Signac)
library(patchwork)

# Generate Coverage Plot
cov_plot <- CoveragePlot(
  object     = atac_obj,
  region     = gene_region,
  annotation = TRUE,
  peaks      = TRUE,
  links      = TRUE,
  expression.assay = NULL,
  features   = NULL,
  extend.upstream = 100000,
  extend.downstream = 100000
) &
  theme(
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  )
cov_plot[[1]][[4]] <- cov_plot[[1]][[4]] &
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  )

Locus_plot_IMAGE <- cov_plot + plot_annotation(
  title = bquote(bold("Chromatin Accessibility and Enhancer–Gene Links for " ~ .(gene_of_interest) ~
                        " from IMAGE (Region: " ~ .(gene_region) ~ " ± 100 kb)")),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 18))
)

library(Seurat)
library(Signac)
library(patchwork)

# IMAGE plot with explicit title
IMAGE_plot_final <- wrap_elements(full = Locus_plot_IMAGE, clip = FALSE) +
  plot_annotation(
    title = paste0("Chromatin Accessibility and Enhancer–Gene Links for ", gene_of_interest,
                   " from IMAGE (Region: ", gene_region, " ± 100 kb)")
  )

# sce2g plot with explicit title
sce2g_plot_final <- wrap_elements(full = Locus_plot_sce2g) +
  plot_annotation(
    title = paste0("Chromatin Accessibility and Enhancer–Gene Links for ", gene_of_interest,
                   " from scE2G (Region: ", gene_region, " ± 100 kb)")
  )

# Combine, label (A/B), and center titles clearly
final_locus_plot <- (IMAGE_plot_final / sce2g_plot_final) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.tag = element_text(size = 20, face = "bold")
  )

# ensure titles appear clearly by setting clip = FALSE
wrap_elements(final_locus_plot, clip = FALSE)

####################### make confusion matrix and recall and precision curve
library(GenomicRanges)
library(dplyr)
library(future.apply)

# Parallel backend
plan(multisession, workers = 5)

df_image_40pct <- readRDS("/work/Home/Data_for_masters_project/df_image_40pct_overlap")
sce2g_40pct <- readRDS("/work/Home/Data_for_masters_project/sce2g_40pct_overlap")
length(unique(df_image_40pct$Gene))
# returns 12771
length(unique(sce2g_40pct$TargetGene))
# returns 12771
identical(sort(unique(sce2g_40pct$TargetGene)), sort(unique(df_image_40pct$Gene)))
# returns TRUE
shared_genes <- intersect(df_image_40pct$Gene, sce2g_40pct$TargetGene)
image_by_gene <- split(df_image_40pct, df_image_40pct$Gene)
sce2g_by_gene <- split(sce2g_40pct, sce2g_40pct$TargetGene)
# Run in parallel to compute confusion matrix
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
confusion_list <- future_lapply(shared_genes, function(gene) {
  image_gene <- image_by_gene[[gene]]
  sce_gene <- sce2g_by_gene[[gene]]
  
  if (is.null(image_gene) || is.null(sce_gene)) return(NULL)
  
  # GRanges from IMAGE
  gr_image <- GRanges(
    seqnames = image_gene$seqnames,
    ranges = IRanges(start = image_gene$start, end = image_gene$end)
  )
  mcols(gr_image)$is_IMAGE <- image_gene$is_IMAGE
  mcols(gr_image)$is_sce2g <- 0  # default to 0 (not predicted by scE2G)
  
  # GRanges from scE2G
  gr_sce <- GRanges(
    seqnames = sce_gene$chr,
    ranges = IRanges(start = sce_gene$start, end = sce_gene$end),
    is_sce2g = sce_gene$is_sce2g
  )
  
  # Find any overlaps. we already only have peaks that overlap by 40%
  ov <- findOverlaps(gr_image, gr_sce)
  
  # if there are overlaps, get the hits for IMAGE
  # get the binerisation from sce2g where there was a hit from sce2g
  # give the sce2g binerisation to IMAGE. We map sce2g binerisation to IMAGE
  if (length(ov) > 0) {
    gr_image_hits <- queryHits(ov)
    sce2g_flags <- mcols(gr_sce)$is_sce2g[subjectHits(ov)]
    mcols(gr_image)$is_sce2g[gr_image_hits] <- sce2g_flags
  }
  
  # create confusion matrix of counts
  tab_actual <- table(
    factor(mcols(gr_image)$is_IMAGE, levels = c(0, 1)),
    factor(mcols(gr_image)$is_sce2g, levels = c(0, 1))
  )
  
  # extract confusion matrix metrics
  TP <- tab_actual["1", "1"]
  FP <- tab_actual["1", "0"]
  FN <- tab_actual["0", "1"]
  TN <- tab_actual["0", "0"]
  
  return(data.frame(Gene = gene, TP = TP, FP = FP, FN = FN, TN = TN))
})

# Combine results
confusion_matrix_per_gene <- bind_rows(confusion_list)
colSums(confusion_matrix_per_gene[, c("TP", "FP", "FN", "TN")])
conf_matrix_total <- colSums(confusion_matrix_per_gene[, c("TP", "FP", "FN", "TN")])
library(reshape2)
# Fix the order for the reduced confusion matrix (Ensuring FN, TP, FP order)
conf_matrix_reduced <- matrix(c(
  conf_matrix_total["TP"], conf_matrix_total["FN"],
  conf_matrix_total["FP"], conf_matrix_total["TN"]
), nrow = 2, byrow = TRUE)

# Assign row and column names
rownames(conf_matrix_reduced) <- c("scE2G Predicted Positive", "scE2G Predicted Negative")
colnames(conf_matrix_reduced) <- c("IMAGE Predicted Positive", "IMAGE Predicted Negative")
# Convert matrices to long format for ggplot
conf_matrix_long2 <- melt(conf_matrix_reduced)

# Compute percentages correctly
conf_matrix_long2$Percentage <- round((conf_matrix_long2$value / sum(conf_matrix_total[c("TP", "FP", "FN", "TN")])) * 100, 2)

# Create the second confusion matrix plot
gg_confusion <- ggplot(conf_matrix_long2, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = ifelse(value > 0, paste0(value, "\n", Percentage, "%"), "")), 
            color = "white", size = 6) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  labs(title = "Confusion Matrix for Enhancer-Gene Links Between IMAGE and scE2G",
       x = "IMAGE",
       y = "scE2G",
       fill = "Number of Links") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold", colour = "black"),
        axis.text = element_text(face  = "bold", colour = "black"),
        legend.title    = element_text(face = "bold"),
        legend.text = element_text(face = "bold", colour = "black"),
        legend.position = "right")


############## Precision recall
# Requires 128 gb memory to run
library(precrec)
library(GenomicRanges)
df_image_40pct <- readRDS("/work/Home/Data_for_masters_project/df_image_40pct_overlap")
sce2g_40pct <- readRDS("/work/Home/Data_for_masters_project/sce2g_40pct_overlap")


gr_image <- GRanges(
  seqnames = df_image_40pct$seqnames,
  ranges   = IRanges(df_image_40pct$start, df_image_40pct$end),
  is_IMAGE  = df_image_40pct$is_IMAGE,
  Potential = df_image_40pct$Potential,
  MotifScore = df_image_40pct$Motifs_enhancer_contribution_score
)

gr_sce2g <- GRanges(
  seqnames = sce2g_40pct$chr,
  ranges   = IRanges(sce2g_40pct$start, sce2g_40pct$end),
  is_sce2g = sce2g_40pct$is_sce2g
)

# findOverlaps
ov <- findOverlaps(gr_image, gr_sce2g)
# make a new column that only contains 0
df_image_40pct$label_sce2g <- 0
# find indices where sce2g have a predicted link
pos_idx <- which(gr_sce2g$is_sce2g == 1)
# get the hits for where we have an overlap for IMAGE
df_hits <- queryHits(ov)
# get the hits for where we have an overlap for sce2g
sce2g_hits <- subjectHits(ov)
# fill in 1 where there is an overlap and where sce2g have predicted it to be a link
df_image_40pct$label_sce2g[df_hits[sce2g_hits %in% pos_idx]] <- 1

# Compare two scores as predictors of label_sce2g
label_vec <- df_image_40pct$label_sce2g
score1 <- df_image_40pct$Potential
score2 <- df_image_40pct$Motifs_enhancer_contribution_score
score3 <- df_image_40pct$Motif_gene_contribution_score
score4 <- df_image_40pct$Motifs_enhancer_contribution_score + df_image_40pct$Motif_gene_contribution_score
score5 <- df_image_40pct$Motifs_enhancer_contribution_score * df_image_40pct$Motif_gene_contribution_score
score6 <- df_image_40pct$Motifs_enhancer_contribution_score * df_image_40pct$Motif_gene_contribution_score * df_image_40pct$Potential
score7 <- df_image_40pct$Motifs_enhancer_contribution_score + df_image_40pct$Motif_gene_contribution_score + df_image_40pct$Potential

mm_obj <- mmdata(
  scores   = list(score1, score2, score3, score4, score5, score6, score7),
  labels   = list(label_vec, label_vec, label_vec, label_vec, label_vec, label_vec, label_vec),
  modnames = c("Regulatory Potential Score", "Motif Enhancer Contribution Score", "Motif Gene Contribution Score", "Motif Contribution Scores Added", "Motif Contribution Scores Multiplied",
               "Motif Contribution Scores and Potential Multiplied", "Motif Contribution Scores and Potential Added")
)
mm_curves <- evalmod(mm_obj)

myplot <- autoplot(mm_curves, "prc")

PRcurve <- myplot + 
  labs(
    title    = "Precision-Recall Plot",
    x        = "Recall", 
    y        = "Precision",
    color    = "Score Metric",
    linetype = "Score Metric" 
  ) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    axis.text = element_text(face = "bold", colour = "black"),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    legend.text = element_text(face = "bold", colour = "black")
  )

# Get AUC
auc_list <- auc(mm_curves)
# subset it to AUPRC
auc_list <- auc_list[auc_list$curvetypes == "PRC",]
# plot it
library(ggplot2)
auc <- ggplot(auc_list, aes(x = modnames, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    x = "Score Metric",
    y = "AUPRC",
    title = "AUPRC for Each Score Metric"
  ) +
  theme_minimal() +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", colour = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold", colour = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)
  )

auc_wrapped <- wrap_elements(full = auc)
PRcurve_wrapped <- wrap_elements(full = PRcurve)
gg_confusion_wrapped <- wrap_elements(full = gg_confusion)

# Create the top row with custom relative widths
bottom_row <- plot_spacer() + 
  auc_wrapped + 
  plot_spacer() +
  plot_layout(widths = c(1, 4, 1))

combined_plot <- (gg_confusion_wrapped + PRcurve_wrapped)/(bottom_row) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )
# Display combined plot
combined_plot
