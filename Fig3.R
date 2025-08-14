# Updated Fig 3.
library(ArchR)
library(grid)
library(gridExtra)
library(Seurat)
library(Signac)

plot.dir <- "/gstore/data/project/fibro_multiome/Fig3plots"
archr_proj <- loadArchRProject("/gstore/project/fibro_multiome/Save-ProjMulti2")
addArchRGenome("hg38")
group <- factor(archr_proj$Sample)
treatment <- factor(archr_proj$Sample)
levels(group) <- c(rep("Normal",4),rep("IPF",4))
levels(treatment) <- c("Normal-NoTreatment","Normal-OSM","Normal-TGF",
                       "Normal-IL13","IPF-NoTreatment","IPF-OSM",
                       "IPF-TGF","IPF-IL13")

archr_proj$Group <- as.character(group)
archr_proj$Treatment <- as.character(treatment)


# Panel 1: TF motif enrichment in IPF vs. normal HLFs ---------
NoTreatment <- archr_proj[archr_proj$Treatment %in%  c("Normal-NoTreatment","IPF-NoTreatment"), ]

## Load motifs ---------
library(TFBSTools)
library(universalmotif)

# Read the MEME file
motifs <- read_meme("/gstore/scratch/u/lucast3/fibroMultiome/jaspar_2024_clusters.meme")

# Convert the motifs to a PWMMatrixList
pwms <- list()
for (i in c(1:length(motifs))) {
  motif <- motifs[[i]]
  pwm <- convert_motifs(motif, "TFBSTools-PWMatrix")
  pwm_name <- gsub(":", "|", pwm@name)
  pwms[[pwm_name]] <- pwm
}

# Convert the list to a PWMatrixList
pwm_list <- do.call(PWMatrixList, pwms)

# Add motif annotations
NoTreatment <- addMotifAnnotations(
  NoTreatment,
  motifPWMs = pwm_list,
  annoName = "jaspar_2024_clusters",
  species = "Homo sapiens"
)

##Markers ---------
markersDA <- getMarkerFeatures(ArchRProj = NoTreatment, useMatrix = "PeakMatrix",
                               groupBy = "Group", testMethod = "wilcoxon",
                               bias = c("TSSEnrichment", "log10(nFrags)")) #A positive log2FC means the peak is more accessible in "IPF-NoTreatment"
enrichedMotifs <- peakAnnoEnrichment(seMarker = markersDA, ArchRProj = NoTreatment,
                                     peakAnnotation = "jaspar_2024_clusters",
                                     cutOff = "FDR <= 0.1 & Log2FC >= 0")
### plot by Pvalue ---------------
heatmap_plot <- plotEnrichHeatmap(enrichedMotifs,
                                  n = 10,
                                  pal =RColorBrewer::brewer.pal(n=7,"YlOrRd")
)

draw(heatmap_plot, heatmap_legend_side = "bottom", )# Adjusts x-axis (column) labels)

### plot by enrichment ----------
library(ComplexHeatmap)
library(circlize)
library(SummarizedExperiment)
enrichment_matrix <- assay(enrichedMotifs, "Enrichment")  # or "mlog10p", depending on what you used

top_n <- 10
top_motifs <- unique(c(
  rownames(enrichment_matrix)[order(enrichment_matrix[, "IPF"], decreasing = TRUE)[1:top_n]],
  rownames(enrichment_matrix)[order(enrichment_matrix[, "Normal"], decreasing = TRUE)[1:top_n]]
))

filtered_matrix <- enrichment_matrix[top_motifs, ]
z <- Heatmap(
  scale(filtered_matrix),
  name = "Norm. Enrichment",
  col = RColorBrewer::brewer.pal(n=7,"YlOrRd"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 12),  # Adjusts y-axis (row) labels
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(direction = "horizontal")
)
draw(z, heatmap_legend_side = "bottom")


# Enriched TFs show increased activity ? ----------------------
NoTreatment <- addDeviationsMatrix(
  ArchRProj = NoTreatment,
  peakAnnotation = "jaspar_2024_clusters",  # must match your custom annotation name
  force = TRUE  # optional, but ensures it re-runs if needed
)

deviation_matrix <- getMatrixFromProject(NoTreatment, useMatrix = "jaspar_2024_clustersMatrix")
deviation_scores <- assay(deviation_matrix,i=2)  # rows = motifs, cols = cells
motif_names <- rownames(deviation_scores)

# Step 2: Get condition metadata
meta <- getCellColData(NoTreatment, select = "Group") |> as.data.frame()
meta$Cell <- rownames(meta)

# Ensure metadata matches order of deviation matrix
meta <- meta[colnames(deviation_scores), , drop = FALSE]

# Step 3: Define groups
ipf_cells <- rownames(meta)[meta$Group == "IPF"]
normal_cells <- rownames(meta)[meta$Group == "Normal"]

# Step 4: Calculate average deviation per group
mean_ipf <- rowMeans(deviation_scores[, ipf_cells], na.rm = TRUE)
mean_normal <- rowMeans(deviation_scores[, normal_cells], na.rm = TRUE)

# Step 5: Calculate delta and assemble into data frame
delta_dev <- mean_ipf - mean_normal

motif_diff_df <- data.frame(
  Motif = motif_names,
  Deviation_IPF = mean_ipf,
  Deviation_Normal = mean_normal,
  Delta = delta_dev
)

## Plot activity ----------
enrichment_matrix <- assay(enrichedMotifs, "mlog10Padj")  # or "mlog10p", depending on what you used

top_n <- 10
top_motifs <- unique(c(
  rownames(enrichment_matrix)[order(enrichment_matrix[, "IPF"], decreasing = TRUE)[1:top_n]],
  rownames(enrichment_matrix)[order(enrichment_matrix[, "Normal"], decreasing = TRUE)[1:top_n]]
))

filtered_diff <- motif_diff_df[top_motifs, ]

ggplot(filtered_diff, aes(x = reorder(Motif, Delta), y = Delta)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top Differential Motif Activities (IPF vs Normal)",
    x = "Motif",
    y = "Deviation Score"
  )

# Increasesd Expression? --------------
tf_families <- sapply(strsplit(top_motifs, "\\|"), function(x) x[2])

# Split multi-TF families into individual gene names
tf_genes <- unique(unlist(strsplit(tf_families, "_")))

# Get gene names from Expression matrix
se <- getMatrixFromProject(NoTreatment, useMatrix = "GeneExpressionMatrix")
se <- logNormCounts(se,assay.type = "GeneExpressionMatrix")

archr_genes <-rowData(se)$name
# Match by startsWith (e.g., all "FOXO", "ATF" genes)
matched_tfs <- unique(unlist(sapply(tf_genes, function(tf) grep(paste0("^", tf), archr_genes, value = TRUE))))

expr_data <- assay(se,i = "GeneExpressionMatrix")
rownames(expr_data) <- rowData(se)$name
expr_data <- expr_data[matched_tfs, ]

# Get condition metadata
meta <- getCellColData(NoTreatment, select = "Group") |> as.data.frame()
meta$Cell <- rownames(meta)

# Match metadata to expression matrix
meta <- meta[colnames(expr_data), , drop = FALSE]

ipf_cells <- rownames(meta)[meta$Group == "IPF"]
normal_cells <- rownames(meta)[meta$Group == "Normal"]

mean_expr_ipf <- rowMeans(expr_data[, ipf_cells], na.rm = TRUE)
mean_expr_normal <- rowMeans(expr_data[, normal_cells], na.rm = TRUE)

expr_comparison <- data.frame(
  Gene = matched_tfs,
  Mean_IPF = mean_expr_ipf,
  Mean_Normal = mean_expr_normal,
  Log2FC = log2((mean_expr_ipf + 1e-3) / (mean_expr_normal + 1e-3))
)

# Rank by differential expression
expr_comparison <- expr_comparison[order(-abs(expr_comparison$Log2FC)), ]
top_n <- 15
top_motifs <- unique(c(
  rownames(expr_comparison[order(expr_comparison$Log2FC, decreasing = TRUE)[1:top_n],]),
  rownames(expr_comparison[order(expr_comparison$Log2FC, decreasing = F)[1:top_n],])))

ggplot(expr_comparison[top_motifs,], aes(x = reorder(Gene, Log2FC), y = Log2FC)) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  theme_minimal() +
  labs(title = "TF Expression Difference (IPF vs Normal)",
       x = "Gene", y = "Log2 Fold Change (IPF / Normal)")


# Correaltion -------------------
library(Hmisc)
library(dplyr)
library(ggplot2)

# Get motif names like "C012|FOXO_FOXP|Forkhead"
motif_names <- rownames(deviation_scores)

# Extract TF family portion: second field after splitting on "|"
tf_families <- sapply(strsplit(motif_names, "\\|"), function(x) x[2])

# Expand multi-TF families (e.g., "FOXO_FOXP") to individual TFs
motif_to_tfs <- strsplit(tf_families, "_")
names(motif_to_tfs) <- motif_names

# For each motif, find the first TF that exists in the expression matrix
expr_genes <- rownames(expr_data)

motif_to_gene <- sapply(motif_to_tfs, function(tfs) {
  match <- intersect(tfs, expr_genes)
  if (length(match) > 0) match[1] else NA
})

# Filter to only motifs with an expressed TF
motif_to_gene <- motif_to_gene[!is.na(motif_to_gene)]
common_motifs <- names(motif_to_gene)
common_genes <- unname(motif_to_gene)

### Subset -------
# Subset chromVAR and expression matrices
dev_filtered <- deviation_scores[common_motifs, ]
expr_filtered <- expr_data[common_genes, ]

# Rename motif rows by gene symbol to match
rownames(dev_filtered) <- common_genes

# Transpose for per-cell comparison (rows = cells)
dev_df <- as.data.frame(t(dev_filtered))
expr_df <- as.data.frame(t(expr_filtered))

# Align by cell (row)
shared_cells <- intersect(rownames(expr_df), rownames(dev_df))
expr_df <- expr_df[shared_cells, ]
dev_df <- dev_df[shared_cells, ]

# Ensure same TF columns
common_tfs <- intersect(colnames(expr_df), colnames(dev_df))
expr_df <- expr_df[, common_tfs]
dev_df <- dev_df[, common_tfs]

### Correlation ------------

cor_results <- data.frame()

for (tf in common_tfs) {
  r <- rcorr(expr_df[[tf]], dev_df[[tf]], type = "pearson")
  cor_results <- rbind(cor_results, data.frame(
    TF = tf,
    Correlation = r$r[1, 2],
    P_value = r$P[1, 2]
  ))
}

# Rank results
cor_results <- cor_results %>% arrange(-abs(Correlation))



### plot -------
ggplot(cor_results, aes(x = reorder(TF, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(title = "TF Expression vs Motif Activity (chromVAR)",
       x = "TF", y = "Pearson Correlation")


# Regulon Activity ------------------------

motif_families <- c("FOXO_FOXP", "CEBP_ATF", "TGIF_MYF", "JUN_FOS", "NFI_TLX", "HIC")
# Assume collectri_regulons is a named list like:
# collectri_regulons$JUN = c("IL6", "STAT3", ...)
# Map motif families to TFs used in regulon list
motif_to_tfs <- list(
  FOXO_FOXP = c("FOXO", "FOXP"),
  CEBP_ATF = c("CEBPB", "ATF"),
  TGIF_MYF = c("TGIF", "MYF"),
  JUN_FOS = c("JUN", "FOS"),
  NFI_TLX = c("NFIX", "TLX"),
  HIC = c("HIC")
)
expanded_tfs <- lapply(motif_to_tfs, function(roots) {
  unique(unlist(sapply(roots, function(root) grep(paste0("^", root), expr_genes, value = TRUE))))
})
expanded_tfs[["JUN_FOS"]] <- c("JUN", "FOS")
expanded_tfs[["NFI_TLX"]] <- c("NFIX", "TLX1","TLX2","TLX3")

# Flatten list into regulon targets
selected_regulons <- list()

for (fam in motif_families) {
  tf_candidates <- expanded_tfs[[fam]]
  valid_tfs <- intersect(tf_candidates, names(collections))
  regulon_genes <- unique(unlist(collections[valid_tfs]))
  selected_regulons[[fam]] <- intersect(regulon_genes, rownames(mat))
}

## Score ---------
for (fam in names(selected_regulons)) {
  genes <- selected_regulons[[fam]]
  if (length(genes) >= 5) {
    NoTreatment <- addModuleScore(
      ArchRProj = NoTreatment,
      useMatrix = "GeneExpressionMatrix",
      name = paste0("Regulon_", fam),
      features = list(genes),
    )
  } else {
    message("Skipping ", fam, ": not enough genes in regulon")
  }
}


group_col <- "Group"
meta <- getCellColData(NoTreatment, select = group_col) |> as.data.frame()
meta$Cell <- rownames(meta)
library(tidyverse)
# Extract regulon scores from cellColData
regulon_cols <- grep("^Regulon_", colnames(getCellColData(NoTreatment)), value = TRUE)
regulon_scores <- getCellColData(NoTreatment, select = regulon_cols[52:length(regulon)]) |> as.data.frame()
regulon_scores$Group <- meta[rownames(regulon_scores), group_col]

# Average per group
regulon_summary <- regulon_scores %>%
  pivot_longer(cols = starts_with("Regulon_"), names_to = "Regulon", values_to = "Score") %>%
  group_by(Regulon, Group) %>%
  summarise(MeanScore = mean(Score, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = MeanScore) %>%
  mutate(Log2FC = log2((`IPF` + 1e-3) / (`Normal` + 1e-3))) %>%
  arrange(-Log2FC)


ggplot(regulon_summary, aes(x = reorder(Regulon, Log2FC), y = Log2FC)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Regulon Activity (IPF vs Normal)",
       x = "TF Motif Family", y = "Log2 Fold Change")
