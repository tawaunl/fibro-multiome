## ----setup, include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,message=FALSE, warning=FALSE,
                      out.width = '150%',out.height = '150%')


## ----Enviroment----------------------------------------------------------------
#load the ArchR library
library(ArchR)
library(grid)
library(gridExtra)
library(dsassembly)
#set a seed to facilitate replication of operations requiring randomization
set.seed(824)

#default number of Parallel threads is 16
# working on a local computer, 1 thread works best
addArchRThreads(14)

#Before we begin, we need add a reference genome annotation for ArchR 
#to have access to chromosome and gene information. ArchR supports hg19, hg38, mm9, and mm10.
addArchRGenome("hg38")
data.dir <- "/gstore/scratch/u/lucast3/fibroticmemory" # setdirectory for files 



## ----Load Project--------------------------------------------------------------
projMulti2 <- loadArchRProject(file.path(data.dir,"Save-ProjMulti2"),
                        showLogo = FALSE)



## ----collapse=TRUE-------------------------------------------------------------
projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)


## ----collapse=TRUE-------------------------------------------------------------
projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)


## ------------------------------------------------------------------------------
projMulti2 <- addCombinedDims(projMulti2, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined",corCutOff = 0.75)


## ------------------------------------------------------------------------------
# Need to order colnames so that findclusters will run
colnames(projMulti2@reducedDims$LSI_Combined$matRD)[1:29]<-1:29
colnames(projMulti2@reducedDims$LSI_Combined$matRD)[30:59]<-1:30
projMulti2@reducedDims$LSI_Combined$matRD <- projMulti2@reducedDims$LSI_Combined$matRD[,order(as.numeric(colnames(projMulti2@reducedDims$LSI_Combined$matRD)))]


## ----collapse=TRUE-------------------------------------------------------------
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)


## ----collapse=TRUE-------------------------------------------------------------
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.4, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.4, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.4, force = TRUE)


## ----collapse=TRUE-------------------------------------------------------------
p1 <- plotEmbedding(projMulti2, name = "Clusters_Combined", embedding = "UMAP_ATAC", size = 1, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(projMulti2, name = "Clusters_Combined", embedding = "UMAP_RNA", size = 1, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(projMulti2, name = "Clusters_Combined", embedding = "UMAP_Combined", size = 1, labelAsFactors=F, labelMeans=F)

p <- lapply(list(p1,p2,p3), function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})


## ----collapse=TRUE-------------------------------------------------------------
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)


## ----collapse=TRUE-------------------------------------------------------------
cM_atac_rna <- confusionMatrix(paste0(projMulti2$Clusters_ATAC), paste0(projMulti2$Clusters_RNA))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)

library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
plotPDF(p_atac_rna,name = "ConfusionMatrix-scATAC-scRNA", addDOC = FALSE)
p_atac_rna


## ----collapse=TRUE-------------------------------------------------------------
require(BSgenome.Hsapiens.UCSC.hg38)
projMulti2 <- addGroupCoverages(ArchRProj = projMulti2, groupBy = "Clusters_Combined",force = TRUE)

projMulti2 <- addReproduciblePeakSet(
  ArchRProj = projMulti2, groupBy = "Clusters_Combined",
  pathToMacs2 = "/gstore/home/lucast3/.conda/envs/MACS/bin/macs2",
  force = TRUE)

projMulti2 <- addPeakMatrix(ArchRProj = projMulti2,force = TRUE)
projMulti2 <- addPeak2GeneLinks(
  ArchRProj = projMulti2, 
  reducedDims = "LSI_Combined",
  useMatrix = "GeneExpressionMatrix")

# examine merged peak set.
getPeakSet(projMulti2)

p2g <- getPeak2GeneLinks(ArchRProj = projMulti2)

#examine links
p2g[[1]]


## ----collapse=TRUE-------------------------------------------------------------
se <- getMarkerFeatures(
  ArchRProj = projMulti2,
  groupBy = "Clusters_Combined",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"))

heatmap_gex <- plotMarkerHeatmap(
  seMarker = se, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.5",
  nLabel = 4,
  transpose = TRUE
)

plotPDF(draw(heatmap_gex, heatmap_legend_side = "bot", annotation_legend_side = "bot"),
        name="Heatmap-ClusterMarkers",addDOC = FALSE)

draw(heatmap_gex, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#We can get a list of DataFrame objects, one for each of our clusters, 
#containing the relevant marker features using the getMarkers() function.
markerList <- getMarkers(se, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")

#save the marker genes object by cluster
saveRDS(se, file = paste0("/gstore/scratch/u/lucast3/fibroticmemory/Save-ProjMulti2/", "markers_Cluster.rds"))


## ------------------------------------------------------------------------------
group <- factor(projMulti2$Sample)
treatment <- factor(projMulti2$Sample)
levels(group) <- c(rep("Normal",4),rep("IPF",4))
levels(treatment) <- c("Normal-NoTreatment","Normal-OSM","Normal-TGF",
                       "Normal-IL13","IPF-NoTreatment","IPF-OSM",
                       "IPF-TGF","IPF-IL13")

projMulti2$Group <- as.character(group)
projMulti2$Treatment <- as.character(treatment)

markersGS <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Group",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

heatmap_gex <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.5",
  nLabel = 4,plotLog2FC = TRUE,
  transpose = TRUE
)

plotPDF(draw(heatmap_gex, heatmap_legend_side = "bot", annotation_legend_side = "bot"),
        name="HeatmapGroupmarkers",addDOC = FALSE)
draw(heatmap_gex, heatmap_legend_side = "bot", annotation_legend_side = "bot")


#save the marker genes object by samples
saveRDS(markersGS, file = paste0("/gstore/scratch/u/lucast3/fibroticmemory/Save-ProjMulti2/", "markers_Group.rds"))


## ------------------------------------------------------------------------------
markersGS <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

heatmap_gex <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  nLabel = 4,plotLog2FC = TRUE,
  transpose = TRUE
)

plotPDF(draw(heatmap_gex, heatmap_legend_side = "bot", annotation_legend_side = "bot"),
        name="HeatmapTreamentMarkers",addDOC=FALSE)

draw(heatmap_gex, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#save the marker genes object by samples
saveRDS(markersGS,
        file = paste0("/gstore/scratch/u/lucast3/fibroticmemory/Save-ProjMulti2/",
                      "markers_Treatment.rds"))


## ------------------------------------------------------------------------------
projMulti2 <- addPeakMatrix(projMulti2,force=TRUE)
getAvailableMatrices(projMulti2)
markersPeaks <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Combined_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)","log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

#DE Peaks by Sample
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
saveRDS(markerList,
        file = paste0("/gstore/scratch/u/lucast3/fibroticmemory/Save-ProjMulti2/",
                      "markerPeaks.rds"))


## ------------------------------------------------------------------------------
#visualizing marker peaks as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.001 & Log2FC >= 3",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projMulti2, addDOC = FALSE)


## ------------------------------------------------------------------------------
##need to install chromVARmotifs
#devtools::install_github("GreenleafLab/chromVARmotifs")
library(chromVARmotifs)
require(BSgenome.Hsapiens.UCSC.hg38)
projMulti2 <- addMotifAnnotations(ArchRProj = projMulti2, motifSet = "cisbp",
                                  name = "Motif",force=TRUE)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projMulti2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

#visualize using heatmap
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projMulti2, addDOC = FALSE)



## ------------------------------------------------------------------------------
projMulti2 <- addBgdPeaks(projMulti2,force=TRUE)

projMulti2 <- addDeviationsMatrix(
  ArchRProj = projMulti2, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projMulti2, name = "MotifMatrix", plot = TRUE)


## ------------------------------------------------------------------------------
markerMotifs <- getFeatures(projMulti2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

projMulti2 <- addImputeWeights(projMulti2,reducedDims = "LSI_Combined")
p <- plotGroups(ArchRProj = projMulti2, 
                groupBy = "Combined_Clusters", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(projMulti2)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})

png("/gstore/scratch/u/lucast3/fibroticmemory/Plots/MotifMarkers.png")
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
dev.off()

## ------------------------------------------------------------------------------
saveArchRProject(ArchRProj = projMulti2,  load = FALSE)


## ------------------------------------------------------------------------------
sessionInfo()

