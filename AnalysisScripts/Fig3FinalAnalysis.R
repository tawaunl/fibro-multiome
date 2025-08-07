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

NoTreatment <- addMotifAnnotations(ArchRProj = NoTreatment, motifSet = "cisbp", name = "Motif")
markersDA <- getMarkerFeatures(ArchRProj = NoTreatment, useMatrix = "PeakMatrix",
                               groupBy = "Group", testMethod = "wilcoxon",
                               bias = c("TSSEnrichment", "log10(nFrags)")) #A positive log2FC means the peak is more accessible in "IPF-NoTreatment"
enrichedMotifs <- peakAnnoEnrichment(seMarker = markersDA, ArchRProj = NoTreatment,
                                     peakAnnotation = "Motif",
                                     cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
library(RColorBrewer)
## Heatmap of top motifs ---------
heatmap_plot <- plotEnrichHeatmap(enrichedMotifs,
                                  n = 15,
                                  pal =RColorBrewer::brewer.pal(n=7,"YlOrRd"))

draw(heatmap_plot, heatmap_legend_side = "bottom")

## Priming Signature -------------
priming_genes <- list(
  ## Transcription factors (priming regulators)
  TFs=c("FOXC2", "FOXP1", "FOXP2", "FOXA1", "FOXA2", "CEBPB", "CEBPD", "STAT3",
        "IRF1", "ATF3", "JUN", "FOS", "SMAD3", "SMAD2", "RELA", "KLF6", "NR3C1"),

  ## Canonical profibrotic effectors
  PFE=c("TGFB1", "TGFBR1", "TGFBR2", "IL6", "IL6ST", "IL13RA1", "PDGFRA", "PDGFRB",
        "FN1", "COL1A1", "COL3A1", "COL5A1", "COL6A1", "POSTN", "SPARC", "THBS1"),

  ## Fibroblast contractility and survival
  FCS=c("ACTA2", "TAGLN", "ITGA5", "ITGB1", "CD44", "CDH11", "MMP2", "MMP9", "TIMP1"),

  ## Senescence / stress response
  SSE=c("CDKN1A", "GDF15", "SERPINE1", "SOD2", "HMOX1", "PLK2", "IGFBP3"),

  ## EMT and matrix remodeling
  Matrix=c("ZEB1", "ZEB2", "TWIST1", "VIM", "SNAI1", "SNAI2", "LOX", "LOXL2"),

  ## Proinflammatory cytokines & chemokines
  CYT=c("CXCL12", "CCL2", "CCL7", "CXCL1", "IL1B", "TNFAIP3"),

  ## Chromatin remodelers & pioneer activity
  Chrom=c("HMGA1", "HMGA2", "EZH2", "BRD4", "KDM6B")
)

priming_genes[["all"]] <- as.vector(unlist(priming_genes))
NoTreatment <- addModuleScore(
  ArchRProj = NoTreatment,
  useMatrix = "GeneExpressionMatrix",
  name = "PrimingSignature",
  features = priming_genes
)

### Plot Modules -----------------
module_names <- c("PrimingSignature.TFs", "PrimingSignature.PFE", "PrimingSignature.FCS",
                  "PrimingSignature.SSE", "PrimingSignature.Matrix", "PrimingSignature.CYT",
                  "PrimingSignature.Chrom","PrimingSignature.all")
df <- getCellColData(NoTreatment, select = c("Sample","Group", module_names)) |>
  as.data.frame()
avg <- df %>%
  group_by(Sample, Group) %>%
  summarise(across(all_of(module_names), mean, na.rm = TRUE), .groups = "drop")

df_long <- df %>% pivot_longer(cols = module_names, names_to = "Module", values_to = "Score")

module_labels <- c(
  PrimingSignature.TFs = "Priming TFs",
  PrimingSignature.PFE = "Profibrotic Effectors",
  PrimingSignature.FCS = "Fibroblast Contractility",
  PrimingSignature.SSE = "Senescence/Stress",
  PrimingSignature.Matrix = "EMT/Matrix Remodeling",
  PrimingSignature.CYT = "Cytokines & Chemokines",
  PrimingSignature.Chrom = "Chromatin Remodelers",
  PrimingSignature.all = "All Priming Genes"
)

df_long$Module <- factor(df_long$Module, levels = module_names, labels = module_labels[module_names])

ggplot(df_long, aes(x = Group, y = Score, fill = Group)) +
  ggrastr::geom_jitter_rast(width = 0.2, size = 0.6, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~ Module, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 12) +
  labs(title = "Priming Module Scores by Condition",
       x = NULL, y = "Module Score") +
  scale_fill_brewer(palette = "Set2") +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


## adding Deviation Scores
NoTreatment <- addDeviationsMatrix(
  ArchRProj = NoTreatment,
  peakAnnotation = "Motif",
  name = "MotifMatrix"
)

TFs_of_interest <- c("CEBPA_155","CEBPB_140", "FOXP2_356", "FOXN3_829", "FOXN2_854",
                     "FOXL2_855","FOXD4_830",
                     "FOXJ1_853","FOXP3_348","FOXD4L3_832")
TF_devs <- getMatrixFromProject(NoTreatment, useMatrix = "MotifMatrix")

TF_scores <- assay(TF_devs)[TFs_of_interest, ]
colnames(TF_scores) <- colnames(NoTreatment)  # Match cell IDs
TF_scores <- t(TF_scores)

TF_df <- as.data.frame(TF_scores)
TF_df$Group <- getCellColData(NoTreatment, select = "Group")[,1]

mod_scores <- getCellColData(NoTreatment, select = module_names) |> as.data.frame()

activity_df <- cbind(TF_df, mod_scores)
activity_df$Group <- getCellColData(NoTreatment, select = "Group")[,1]

activity_df$PrimingComposite <- rowMeans(activity_df[, c(TFs_of_interest,module_names)], na.rm = TRUE)

NoTreatment$PrimingComposite <- activity_df$PrimingComposite

plotEmbedding(
  ArchRProj = NoTreatment,
  colorBy = "cellColData",
  name = "PrimingComposite",
  embedding = "UMAP_Combined"
)


# Panel 3: Treatment-Induced Transcriptional Activation in IPF Fibroblasts --------

ipf <- archr_proj[archr_proj$Group == "IPF", ]
markersTFs <- getMarkerFeatures(ArchRProj = ipf, useMatrix = "GeneExpressionMatrix",
                                groupBy = "Treatment", testMethod = "wilcoxon",
                                bias = c("TSSEnrichment", "log10(nFrags)"))
# Read in TFs from Lambert et. al
tfs <- openxlsx::readWorkbook("/gstore/data/project/fibro_multiome/Lambert_TFs.xlsx",
                              sheet = 2)
colnames(tfs) <- tfs[1,]
colnames(tfs)[4] <- "status"
tfs <- tfs %>% filter(status=="Yes")

tf_names <- tfs$Name
rownames(markersTFs) <- markersTFs@NAMES
tf_markers <- markersTFs[rownames(markersTFs) %in% tf_names,]
marker_list <- getMarkers(tf_markers)  # no filter yet
## Profibrotic TFs--------
profibrotic_TFs <- c(

  ## TGFβ/SMAD axis
  "SMAD2", "SMAD3", "SMAD4", "RUNX1", "RUNX2",

  ## JAK/STAT and cytokine-responsive TFs
  "STAT1", "STAT3", "STAT5A", "STAT6",
  "IRF1", "IRF5", "IRF7",

  ## AP-1 complex and stress response
  "JUN", "JUNB", "FOS", "FOSB", "ATF3", "BATF", "MAF",

  ## EMT, ECM remodeling, mesenchymal transition
  "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "GLI1",

  ## Senescence & inflammation
  "CEBPB", "CEBPD", "RELA", "RELB", "NFKB1", "NFKB2", "ETS2",

  ## Fibroblast activation & matrix production
  "KLF4", "KLF6", "FOXC2", "FOXP1", "NFATC2", "MEF2A", "MEF2C",

  ## Pro-survival and proliferation
  "MYC", "EGR1", "ELK1", "SP1", "YY1", "NR3C1", "BCL6"
)

marker_fibrotic <- lapply(marker_list, function(df) {
  data.frame(df) %>%
    filter(name %in% profibrotic_TFs) %>%
    arrange(desc(Log2FC))
})
top_markers <- lapply(marker_list, function(df) {
  df <- data.frame(df) %>%
    arrange(desc(Log2FC)) %>%
    filter(!is.na(Log2FC)) %>%
    slice_head(n = 7)
  df$group <- df$name  # keep gene symbol
  return(df)
})
top_markers_fibrotic <- lapply(marker_fibrotic, function(df) {
  df %>% slice_head(n = 7)
})
top_TFs_fibrotic <- unique(unlist(lapply(top_markers_fibrotic, \(df) df$name)))

top_markers_df <- bind_rows(top_markers, .id = "Treatment")
top_TFs <- unique(top_markers_df$name)

## Plot Heatmap of DE TFs -------------
plotMarkerHeatmap(seMarker = tf_markers,
                  labelMarkers = top_TFs,
                  transpose = TRUE)
plotMarkerHeatmap(
  seMarker = tf_markers,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
  labelMarkers = top_TFs_fibrotic,
  transpose = TRUE
)

### Confirm TF Activity via chromVAR Deviation Scores -------
deviationMat <- getMatrixFromProject(ipf, useMatrix = "MotifMatrix")
motif_names <- rownames(getMatrixFromProject(ipf, useMatrix = "MotifMatrix"))
matched_motifs <- motif_names[sapply(motif_names, function(motif) {
  any(startsWith(motif, top_TFs_fibrotic))
})]

dev_scores <- assay(deviationMat)[matched_motifs,]
marker_list <- getMarkers(markersTFs, cutOff = "FDR <= 1")

# Get one combined table with treatment labels
logfc_df <- bind_rows(lapply(names(marker_list), function(treat) {
  data.frame(marker_list[[treat]]) %>%
    filter(name %in% profibrotic_TFs) %>%
    mutate(Treatment = treat)
}), .id = "TreatmentGroup")

cell_metadata <- getCellColData(ipf, select = c("Treatment")) |> as.data.frame()

chromvar_df <- as.data.frame(t(dev_scores))
chromvar_df$Treatment <- cell_metadata[, "Treatment"]

chromvar_avg <- chromvar_df %>%
  group_by(Treatment) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(-Treatment, names_to = "TF", values_to = "DeviationScore")
chromvar_avg$TF_clean <- sapply(chromvar_avg$TF, function(x) strsplit(x, "_")[[1]][1])
logfc_df$TF <- logfc_df$name

merged_df <- left_join(logfc_df, chromvar_avg,
                       by = c("TF" = "TF_clean", "Treatment" = "Treatment"))

library(ggrepel)
merged_df$Log2FC <- merged_df$Log2FC*-1
ggplot(merged_df, aes(x = Log2FC, y = DeviationScore, color = Treatment, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(max.overlaps = 20, size = 3.5) +
  facet_wrap(~ Treatment, ncol = 4) +
  theme_minimal() +
  labs(title = "Profibrotic TF Expression vs. Motif Activity",
       x = "Log2 Fold Change (Expression)",
       y = "chromVAR Deviation Score (Motif Activity)")

# Panel 4: Differential Response of Normal Fibroblasts to Profibrotic Cytokines--------


normal <- archr_proj[archr_proj$Group == "Normal", ]
markersTFs <- getMarkerFeatures(ArchRProj = normal, useMatrix = "GeneExpressionMatrix",
                                groupBy = "Treatment", testMethod = "wilcoxon",
                                bias = c("TSSEnrichment", "log10(nFrags)"))
rownames(markersTFs) <- rowData(markersTFs)$name
tf_markers <- markersTFs[rownames(markersTFs) %in% tf_names,]

marker_list <- getMarkers(tf_markers,cutOff = "FDR <= 1")  # no filter yet

marker_fibrotic <- lapply(marker_list, function(df) {
  data.frame(df) %>%
    filter(name %in% profibrotic_TFs) %>%
    arrange(desc(Log2FC))
})
top_markers <- lapply(marker_list, function(df) {
  df <- data.frame(df) %>%
    arrange(desc(Log2FC)) %>%
    filter(!is.na(Log2FC)) %>%
    slice_head(n = 7)
  df$group <- df$name  # keep gene symbol
  return(df)
})
top_markers_fibrotic <- lapply(marker_fibrotic, function(df) {
  df %>% slice_head(n = 7)
})
top_TFs_fibrotic <- unique(unlist(lapply(top_markers_fibrotic, \(df) df$name)))

top_markers_df <- bind_rows(top_markers, .id = "Treatment")
top_TFs <- unique(top_markers_df$name)

## Plot Heatmap of DE TFs -------------
plotMarkerHeatmap(seMarker = tf_markers,
                  labelMarkers = top_TFs,
                  transpose = TRUE)
plotMarkerHeatmap(
  seMarker = tf_markers,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
  labelMarkers = top_TFs_fibrotic,
  transpose = TRUE
)

### Confirm TF Activity via chromVAR Deviation Scores -------
deviationMat <- getMatrixFromProject(normal, useMatrix = "MotifMatrix")
motif_names <- rownames(getMatrixFromProject(normal, useMatrix = "MotifMatrix"))
matched_motifs <- motif_names[sapply(motif_names, function(motif) {
  any(startsWith(motif, top_TFs_fibrotic))
})]

dev_scores <- assay(deviationMat)[matched_motifs,]
marker_list <- getMarkers(markersTFs, cutOff = "FDR <= 1")

# Get one combined table with treatment labels
logfc_df <- bind_rows(lapply(names(marker_list), function(treat) {
  data.frame(marker_list[[treat]]) %>%
    filter(name %in% profibrotic_TFs) %>%
    mutate(Treatment = treat)
}), .id = "TreatmentGroup")

cell_metadata <- getCellColData(normal, select = c("Treatment")) |> as.data.frame()

chromvar_df <- as.data.frame(t(dev_scores))
chromvar_df$Treatment <- cell_metadata[, "Treatment"]

chromvar_avg <- chromvar_df %>%
  group_by(Treatment) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(-Treatment, names_to = "TF", values_to = "DeviationScore")
chromvar_avg$TF_clean <- sapply(chromvar_avg$TF, function(x) strsplit(x, "_")[[1]][1])
logfc_df$TF <- logfc_df$name

merged_df <- left_join(logfc_df, chromvar_avg,
                       by = c("TF" = "TF_clean", "Treatment" = "Treatment"))

library(ggrepel)
merged_df$Log2FC <- merged_df$Log2FC*-1
ggplot(merged_df, aes(x = Log2FC, y = DeviationScore, color = Treatment, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(max.overlaps = 20, size = 3.5) +
  facet_wrap(~ Treatment, ncol = 4) +
  theme_minimal() +
  labs(title = "Profibrotic TF Expression vs. Motif Activity",
       x = "Log2 Fold Change (Expression)",
       y = "chromVAR Deviation Score (Motif Activity)")


# Wrapping it all together ------------
expr_mat <- getMatrixFromProject(archr_proj, useMatrix = "GeneExpressionMatrix")
expr_df <- as.data.frame(assay(expr_mat))
expr_df$TF <- rowData(expr_mat)$name

# Keep only profibrotic TFs
expr_df <- expr_df[expr_df$TF %in% profibrotic_TFs, ]

expr_long <- expr_df %>%
  pivot_longer(cols = -TF, names_to = "Cell", values_to = "Expression") %>%
  mutate(TF = factor(TF)) %>%
  left_join(getCellColData(archr_proj)[, c("Sample", "Treatment", "Group")] |>
              as.data.frame() |>
              tibble::rownames_to_column("Cell"))

dev_mat <- getMatrixFromProject(archr_proj, useMatrix = "MotifMatrix")
motif_names <- rownames(dev_mat)
tf_motifs <- motif_names[sapply(motif_names, function(m) any(startsWith(m, profibrotic_TFs)))]

# Create lookup
motif_to_tf <- sapply(tf_motifs, function(m) strsplit(m, "_")[[1]][1])
dev_df <- as.data.frame(assay(dev_mat)[tf_motifs, ])
dev_df$Motif <- tf_motifs
dev_df$TF <- motif_to_tf

dev_long <- dev_df %>%
  pivot_longer(cols = -c(TF, Motif),
               names_to = "Cell", values_to = "DeviationScore") %>%
  left_join(getCellColData(archr_proj)[, c("Sample", "Treatment", "Group")] |>
              as.data.frame() |>
              tibble::rownames_to_column("Cell"))

tf_combined <- left_join(expr_long, dev_long,
                         by = c("TF", "Cell", "Treatment", "Sample", "Group"))
mean_expr <- tf_combined %>%
  group_by(TF, Treatment) %>%
  summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
ref_expr <- mean_expr %>% filter(Treatment == "Normal-NoTreatment") %>%
  rename(RefExpr = MeanExpression)

# Join to compute log2FC for all treatments relative to baseline
expr_logfc <- mean_expr %>%
  left_join(ref_expr, by = "TF") %>%
  mutate(Log2FC = log2((MeanExpression + 1e-3) / (RefExpr + 1e-3)))  # small offset to avoid divide-by-zero
colnames(expr_logfc)[2] <- "Treatment"
mean_dev <- tf_combined %>%
  group_by(TF, Treatment) %>%
  summarise(MeanDeviation = mean(DeviationScore, na.rm = TRUE), .groups = "drop")
tf_summary <- left_join(expr_logfc, mean_dev, by = c("TF", "Treatment"))
tf_summary_clean <- tf_summary %>%
  filter(!is.na(Log2FC), !is.na(MeanDeviation))
tf_summary_scaled <- tf_summary_clean %>%
  mutate(
    Log2FC_scaled = scale(Log2FC)[,1],
    Deviation_scaled = scale(MeanDeviation)[,1],
    TFPrimingScore = Log2FC_scaled + Deviation_scaled
  )
summary_stats <- tf_summary_scaled %>%
  group_by(Treatment) %>%
  summarise(
    MeanLog2FC = mean(Log2FC, na.rm = TRUE),
    MeanDeviation = mean(MeanDeviation, na.rm = TRUE),
    TFPrimingIndex = mean(TFPrimingScore, na.rm = TRUE)
  )


ggplot(summary_stats, aes(x = Treatment, y = TFPrimingIndex, fill = Treatment)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_minimal() +
  labs(title = "TF Priming Index by Treatment",
       y = "TF Priming Index (scaled expression + motif activity)",
       x = NULL)

long_stats <- pivot_longer(summary_stats, cols = colnames(summary_stats)[c(2:4)],
                           values_to = "value",names_to ="variable")

ggplot(long_stats, aes(x = Treatment, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean TF Activation Metrics by Treatment",
       y = "Mean Value", fill = "Metric")


## Stats ------------
# Assuming 'SampleType' column exists: IPF or Normal
tf_combined$SampleType <- ifelse(grepl("IPF", tf_combined$Treatment), "IPF", "Normal")

t.test(log2FC ~ SampleType, data = tf_combined)
t.test(DeviationScore ~ SampleType, data = tf_combined)

