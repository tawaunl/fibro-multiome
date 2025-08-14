# Download CollecTRI
#
# Author:  Oriol Fornes
# Date:    Oct 1, 2024

## Get CollecTRI
library(decoupleR)
library(ArchR)
library(grid)
library(gridExtra)
library(Seurat)
library(Signac)
if(file.exists("/gstore/project/fibro_multiome/CollecTRI.rds")){
  collections <- readRDS("/gstore/project/fibro_multiome/CollecTRI.rds")
} else{
  human_collectri <- get_collectri(organism="human", split_complexes=FALSE)

  ## Group by `source_genesymbol` and create a list
  library(dplyr)
  library(purrr)
  library(tibble)
  human_collectri_regulons <- human_collectri %>%
    group_by(source) %>%
    summarise(target = list(target)) %>%
    deframe()

  ## Convert the list to a `CharacterList`
  library(IRanges)
  collections <- CharacterList(human_collectri_regulons)

  ## Save
  saveRDS(collections, "/gstore/project/fibro_multiome/CollecTRI.rds")

}

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

se <- getMatrixFromProject(NoTreatment,useMatrix = "GeneExpressionMatrix")
library(scater)
se <- logNormCounts(se,assay.type = "GeneExpressionMatrix")
mat <- se@assays@data@listData[["logcounts"]]

rownames(mat) <- rowData(se)$name
mat <- as.matrix(mat)

net <- decoupleR::get_collectri(organism = 'human',
                                split_complexes = FALSE)

# Run ulm
acts <- decoupleR::run_ulm(mat =mat ,
                           net = net,
                           .source = 'source',
                           .target = 'target',
                           .mor='mor',
                           minsize = 5)
saveRDS(acts,"/gstore/scratch/u/lucast3/fibroMultiome/acts.rds")
saveRDS(se,"/gstore/scratch/u/lucast3/fibroMultiome/se_Expression.rds")

acts <- readRDS("/gstore/scratch/u/lucast3/fibroMultiome/acts.rds")
se <- readRDS("/gstore/scratch/u/lucast3/fibroMultiome/se_Expression.rds")

library(SingleCellExperiment)
se <- as(se, "SingleCellExperiment")
assayNames(se)[1] <- "counts"
rownames(se) <- rowData(se)$name

data <- as.Seurat(se)
data[['tfsulm']] <- acts %>%
  tidyr::pivot_wider(id_cols = 'source',
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = data) <- "tfsulm"

# Scale the data
data <- Seurat::ScaleData(data)
data@assays$tfsulm@data <- data@assays$tfsulm@scale.data
Embeddings(data)


n_tfs <- 25

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = data$Group) %>%
  tidyr::pivot_longer(cols = -cluster,
                      names_to = "source",
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))


# Get top tfs with more variable means across clusters
tfs <- df %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)


# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  dplyr::filter(source %in% tfs) %>%
  tidyr::pivot_wider(id_cols = 'cluster',
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   treeheight_row = 20,
                   treeheight_col = 20)

## get TFs of interest -----------
tfs <- c("JUN","FOS","CEBPA","CEBPB","CEBPD","FOXO1","FOXO3","FOXO4","FOXO6",
         "ATF1","ATF2","ATF3","ATF4","ATF5","ATF6","MYF5","MYF6")
# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  dplyr::filter(source %in% tfs) %>%
  tidyr::pivot_wider(id_cols = 'cluster',
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames('cluster') %>%
  as.matrix()

pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   treeheight_row = 20,
                   treeheight_col = 20)
