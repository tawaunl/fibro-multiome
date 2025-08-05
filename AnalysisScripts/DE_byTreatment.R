# DE analysis
# Load Libraries and set up enviroment ======================================

library(ArchR)
library(grid)
library(gridExtra)
library(dsassembly)
library(gg4way)
#set a seed to facilitate replication of operations requiring randomization
set.seed(824)

#default number of Parallel threads is 16
# working on a local computer, 1 thread works best
addArchRThreads(1)

#Before we begin, we need add a reference genome annotation for ArchR 
addArchRGenome("hg38")
data.dir <- "/gstore/scratch/u/lucast3/fibroticmemory" # setdirectory for files ]


# Load Data ============================================================
projMulti2 <- loadArchRProject(file.path(data.dir,"Save-ProjMulti2"),
                                showLogo = FALSE)
library("org.Hs.eg.db")

# IL13 =========================================================================
markersIL13_IPF <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("IPF-IL13"),
  bgdGroups = c("IPF-NoTreatment"),
  useMatrix = "GeneExpressionMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

rowData(markersIL13_IPF)$ID <- mapIds(org.Hs.eg.db,
                             keys = markersIL13_IPF@elementMetadata@listData[["name"]],
                             column = c("ENSEMBL"),
                             keytype = "SYMBOL")

pv <- markerPlot(seMarker = markersIL13_IPF, 
                 name = "IPF-IL13", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "IPF-IL13-vs-IPF-NoTreatment-Markers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)

saveRDS(markersIL13_IPF,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/markers_IPF_IL13.rds")))

## Normal =====================================
markersIL13_Normal <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("Normal-IL13"),
  bgdGroups = c("Normal-NoTreatment"),
  useMatrix = "GeneExpressionMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
rowData(markersIL13_Normal)$ID <- mapIds(org.Hs.eg.db,
                                      keys = markersIL13_Normal@elementMetadata@listData[["name"]],
                                      column = c("ENSEMBL"),
                                      keytype = "SYMBOL")
pv <- markerPlot(seMarker = markersIL13_Normal, 
                 name = "Normal-IL13", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "Normal-IL13-vs-Normal-NoTreatment-Markers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)
saveRDS(markersIL13_Normal,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/markers_Normal_IL13.rds")))

### 4way =====================================================================
IPF <- as.data.frame(markersIL13_IPF@assays@data@listData)
colnames(IPF) <- names(markersIL13_IPF@assays@data@listData)
IPF$symbol <- rowData(markersIL13_IPF)$name
IPF$ID <- rowData(markersIL13_IPF)$ID

Normal <- as.data.frame(markersIL13_Normal@assays@data@listData)
colnames(Normal) <- names(markersIL13_Normal@assays@data@listData)
Normal$symbol <- rowData(markersIL13_Normal)$name
Normal$ID <- rowData(markersIL13_Normal)$ID

x <- list("IL13vsNotreatment_IPF"=IPF,"IL13vsNotreatment_Normal"=Normal)
png(file.path(data.dir,"Save-ProjMulti2/Plots/IL13_4way.png"))
gg4way(DGEdata = x,
       x = "IL13vsNotreatment_IPF",
       y = "IL13vsNotreatment_Normal", sep = "vs",
       logFC = "Log2FC",FDR="FDR", label=F, textSize =10,) 
dev.off()

p <- gg4way(DGEdata = x,
            x = "IL13vsNotreatment_IPF",
            y = "IL13vsNotreatment_Normal", sep = "vs",
            logFC = "Log2FC",FDR="FDR", label=F, textSize =10,) 
data <- p$data
uniqueIL13 <- data[which(data$Significant=="Significant in IL13vsNotreatment_IPF"),]

# OSM =========================================================================
markersOSM_IPF <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("IPF-OSM"),
  bgdGroups = c("IPF-NoTreatment"),
  useMatrix = "GeneExpressionMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
pv <- markerPlot(seMarker = markersOSM_IPF, 
                 name = "IPF-OSM", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "IPF-OSM-vs-IPF-NoTreatment-Markers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)

saveRDS(markersOSM_IPF,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/markers_IPF_OSM.rds")))

markersOSM_Normal <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("Normal-OSM"),
  bgdGroups = c("Normal-NoTreatment"),
  useMatrix = "GeneExpressionMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
pv <- markerPlot(seMarker = markersOSM_Normal, 
                 name = "Normal-OSM", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "Normal-OSM-vs-Normal-NoTreatment-Markers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)
saveRDS(markersOSM_Normal,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/markers_Normal_OSM.rds")))
### 4way =====================================================================
IPF <- as.data.frame(markersOSM_IPF@assays@data@listData)
colnames(IPF) <- names(markersOSM_IPF@assays@data@listData)
IPF$symbol <- rowData(markersOSM_IPF)$name
IPF$ID <- rowData(markersOSM_IPF)$ID

Normal <- as.data.frame(markersOSM_Normal@assays@data@listData)
colnames(Normal) <- names(markersOSM_Normal@assays@data@listData)
Normal$symbol <- rowData(markersOSM_Normal)$name
Normal$ID <- rowData(markersOSM_Normal)$ID

x <- list("OSMvsNotreatment_IPF"=IPF,"OSMvsNotreatment_Normal"=Normal)

png(file.path(data.dir,"Save-ProjMulti2/Plots/OSM_4way.png"))
gg4way(DGEdata = x,
       x = "OSMvsNotreatment_IPF",
       y = "OSMvsNotreatment_Normal", sep = "vs",
       logFC = "Log2FC",FDR="FDR", label=F, textSize =10,) 
dev.off()


p <- gg4way(DGEdata = x,
            x = "OSMvsNotreatment_IPF",
            y = "OSMvsNotreatment_Normal", sep = "vs",
            logFC = "Log2FC",FDR="FDR", label=F, textSize =10,) 
data <- p$data
uniqueOSM <- data[which(data$Significant=="Significant in OSMvsNotreatment_IPF"),]
# TGF =========================================================================
markersTGF_IPF <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("IPF-TGF"),
  bgdGroups = c("IPF-NoTreatment"),
  useMatrix = "GeneExpressionMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
pv <- markerPlot(seMarker = markersTGF_IPF, 
                 name = "IPF-TGF", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "IPF-TGF-vs-IPF-NoTreatment-Markers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)

saveRDS(markersTGF_IPF,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/markers_IPF_TGF.rds")))

markersTGF_Normal <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("Normal-TGF"),
  bgdGroups = c("Normal-NoTreatment"),
  useMatrix = "GeneExpressionMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
pv <- markerPlot(seMarker = markersTGF_Normal, 
                 name = "Normal-TGF", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "Normal-TGF-vs-Normal-NoTreatment-Markers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)
saveRDS(markersTGF_Normal,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/markers_Normal_TGF.rds")))

### 4way =====================================================================
IPF <- as.data.frame(markersTGF_IPF@assays@data@listData)
colnames(IPF) <- names(markersTGF_IPF@assays@data@listData)
IPF$symbol <- rowData(markersTGF_IPF)$name
IPF$ID <- rowData(markersTGF_IPF)$ID

Normal <- as.data.frame(markersTGF_Normal@assays@data@listData)
colnames(Normal) <- names(markersTGF_Normal@assays@data@listData)
Normal$symbol <- rowData(markersTGF_Normal)$name
Normal$ID <- rowData(markersTGF_Normal)$ID

x <- list("TGFvsNotreatment_IPF"=IPF,"TGFvsNotreatment_Normal"=Normal)
png(file.path(data.dir,"Save-ProjMulti2/Plots/TGF_4way.png"))
gg4way(DGEdata = x,
       x = "TGFvsNotreatment_IPF",
       y = "TGFvsNotreatment_Normal", sep = "vs",
       logFC = "Log2FC",FDR="FDR", label=F, textSize =10) 
dev.off()

p <- gg4way(DGEdata = x,
            x = "TGFvsNotreatment_IPF",
            y = "TGFvsNotreatment_Normal", sep = "vs",
            logFC = "Log2FC",FDR="FDR", label=F, textSize =10) 
data <- p$data
uniqueTGF <- data[which(data$Significant=="Significant in TGFvsNotreatment_IPF"),]

# Get Overlapping Genes =======================================================
library(VennDiagram)

pos_list <- list(#IL13=uniqueIL13$symbol[which(uniqueIL13$`IL13vsNotreatment_IPF Direction`=="Up")],
                 OSM=uniqueOSM$symbol[which(uniqueOSM$`OSMvsNotreatment_IPF Direction`=="Up")],
                 TGF=uniqueTGF$symbol[which(uniqueTGF$`TGFvsNotreatment_IPF Direction`=="Up")])

neg_list <- list(#IL13=uniqueIL13$symbol[which(uniqueIL13$`IL13vsNotreatment_IPF Direction`=="Down")],
                 OSM=uniqueOSM$symbol[which(uniqueOSM$`OSMvsNotreatment_IPF Direction`=="Down")],
                 TGF=uniqueTGF$symbol[which(uniqueTGF$`TGFvsNotreatment_IPF Direction`=="Down")])
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

## Chart Positive Enrichment ==============================================
venn.diagram(
  x = pos_list,
  category.names = names(pos_list),
  filename = file.path(data.dir,"Save-ProjMulti2/Plots/Positivegenes.png"),
  output=TRUE,
  disable.logging = TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

## Chart Negative Enrichment ==============================================
venn.diagram(
  x = neg_list,
  category.names = names(neg_list),
  filename = file.path(data.dir,"Save-ProjMulti2/Plots/Negativegenes.png"),
  output=TRUE,
  disable.logging = TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

Reduce(intersect, pos_list)
Reduce(intersect, neg_list)
save(list=c("pos_list","neg_list"),file = file.path(data.dir,"overlaps.Rdata"))


se <- getMarkerFeatures(ArchRProj = projMulti2, useMatrix = "GeneExpressionMatrix",
                        groupBy = "Treatment")
subsetSE <- se[which(rowData(se)$name %in% c(Reduce(intersect, pos_list)
                                             ,Reduce(intersect, neg_list))),]
plotMarkerHeatmap(seMarker = subsetSE)


