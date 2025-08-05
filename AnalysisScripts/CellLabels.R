# Script to add Fibroblast Subtype scoring to dataset

# Load Libraries and set up enviroment ======================================

library(ArchR)
library(grid)
library(gridExtra)
library(dsassembly)
#set a seed to facilitate replication of operations requiring randomization
set.seed(824)

#default number of Parallel threads is 16
# working on a local computer, 1 thread works best
addArchRThreads(1)

#Before we begin, we need add a reference genome annotation for ArchR 
addArchRGenome("hg38")
data.dir <- "/gstore/scratch/u/lucast3/fibroticmemory" # setdirectory for files ]


# Load Data ============================================================
projMulti2 <- loadprojMulti2ect(file.path(data.dir,"Save-ProjMulti2"),
                               showLogo = FALSE)

# Mesenchymal markers
mesenchymal_list <- list(
                         PanFibroblast=c("DCN", "DPT", "FBLN1", "FN1", "HAS2", 
                                         "LUM", "PDGFRA", "PDGFRB", "PDPN"),
                         Adventitial=c("CD34", "PI16", "CADM3", "CLEC3B", "HAS1", "MFAP5", 
                                       "PDGFRL", "SFRP2", "PLA2G2A", "VIT"),
                         Alveolar=c("CD82", "CES1", "FGFR4", "GPC3", "GPM6B", "INMT", 
                                    "ITGA8", "LIMCH1", "NPNT", "SPINT2"),
                         Myofibroblast=c("CDH11", "ASPN", "CTHRC1", "COMP", "FAP", 
                                         "ELN", "POSTN", "TNC"),
                         SMC=c("ACTA2", "ACTG2", "CNN1", "DES",  
                               "MYH11", "MYL9", "NOTCH3", "TAGLN", "TPM2"),
                         Pericyte=c("COX4I2", "CSPG4", "RGS5-ENSG00000143248", "TRPC6"),
                         Lipofibroblast=c("APOE", "ALDH1A3", "ANGPTL4", "FGF10", "FST", 
                                          "LEP", "PLIN2"),
                         Mesothelial=c("CD44", "VCAM1", "WT1", "CALB2", "KRT18", "KRT19"),
                         Collagens=c("COL1A1", "COL3A1", "COL4A1", "COL6A1", "COL8A1", "COL10A1", 
                                     "COL13A1", "COL14A1", "COL15A1", "COL16A1"),
                         TGFb=c("TGFB1", "TGFB2", "TGFB3", "TGFBI", "TGFBR1", "TGFBR2", "TGFBR3"),
                         Fibrosis=c("LRRC15", "LRRC17", "LRRC32", "LOXL1", "MMP11",
                                    "PITX2", "PTX3", "SERPINE1", "SERPINH1"),
                         CSF=c("CSF1", "CSF2", "CSF3"),
                         gp130=c("IL6", "LIF", "OSM", "IL6R", "LIFR", "OSMR", "OSMR-AS1"))
#Heatmap of markers modules
projMulti2 <- addModuleScore(projMulti2, features = mesenchymal_list, useMatrix = "GeneExpressionMatrix")

se <- readRDS(paste0("/gstore/scratch/u/lucast3/fibroticmemory/Save-ProjMulti2/", "markers_Cluster.rds"))

subsetSE <- se[which(rowData(se)$name %in% c(mesenchymal_list$PanFibroblast,
                                             mesenchymal_list$Adventitial,mesenchymal_list$Alveolar,
                                             mesenchymal_list$Myofibroblast,mesenchymal_list$SMC,
                                             mesenchymal_list$Lipofibroblast,mesenchymal_list$Mesothelial,
                                             mesenchymal_list$Collagens,mesenchymal_list$Fibrosis,mesenchymal_list$CSF,
                                             mesenchymal_list$gp130)),]
heatmap_gex <- plotMarkerHeatmap(
  seMarker = subsetSE, 
  cutOff = "FDR <= 1.5 & Log2FC >= 0",
  nLabel = 25,
  transpose = TRUE,
  clusterCols = FALSE
)

draw(heatmap_gex, heatmap_legend_side = "bot", annotation_legend_side = "bot")


# umaps of moudule scores ==============================================
useMatrix = "GeneScoreMatrix"
threads=1
nBgd=100
features <- mesenchymal_list
featureDF <- ArchR:::.getFeatureDF(head(getArrowFiles(projMulti2),2), subGroup="GeneScoreMatrix")

rownames(featureDF) <- paste0(featureDF$seqnames, ":", featureDF$idx)
featureDF$Match <- seq_len(nrow(featureDF))

if(useMatrix %ni% getAvailableMatrices(projMulti2)){
  stop("useMatrix not in available matrices! See getAvailableMatrices!")
}

matrixClass <- h5read(getArrowFiles(projMulti2)[1], paste0("GeneScoreMatrix", "/Info/Class"))

if(matrixClass == "Sparse.Assays.Matrix"){
  if(!all(unlist(lapply(unlist(features), function(x) grepl(":",x))))){
    .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(projMulti2, useMatrix) to list out available formats for input!", logFile = logFile)
    stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(projMulti2, useMatrix) to list out available formats for input!")
  }
}

if(grepl(":",unlist(features)[1])){
  
  sname <- stringr::str_split(unlist(features),pattern=":",simplify=TRUE)[,1]
  name <- stringr::str_split(unlist(features),pattern=":",simplify=TRUE)[,2]
  
  idx <- lapply(seq_along(name), function(x){
    ix <- intersect(which(tolower(name[x]) == tolower(featureDF$name)), BiocGenerics::which(tolower(sname[x]) == tolower(featureDF$seqnames)))
    if(length(ix)==0){
      .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", name[x]), logFile = logFile)
    }
    ix
  }) %>% unlist
  
}else{
  
  idx <- lapply(seq_along(unlist(features)), function(x){
    ix <- which(tolower(unlist(features)[x]) == tolower(featureDF$name))[1]
    if(length(ix)==0){
      .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", unlist(features)[x]), logFile = logFile)
    }
    ix
  }) %>% unlist
  
}

if(is.null(names(features))){
  names(features) <- paste0(name, seq_along(features))
}else{
  names(features) <- paste0(name, ".", names(features))
}

featuresUse <- featureDF[idx[which(is.na(idx)==FALSE)],]
featuresUse$Module <- Rle(stack(features)[,2][-49])

#Get Averages
rS <- ArchR:::.getRowSums(ArrowFiles = getArrowFiles(projMulti2), useMatrix = "GeneScoreMatrix")
rS <- rS[order(rS[,3]), ]
rS$Bins <- Rle(ggplot2::cut_number(x = rS[,3] + rnorm(length(rS[,3]))/1e30, n = 25, labels = FALSE, right = FALSE))
rS$Match <- match(paste0(rS$seqnames, ":", rS$idx), rownames(featureDF))

if(100 > min(rS$Bins@lengths)){
  stop("nBgd must be lower than ", min(rS$Bins@lengths), "!")
}

idxMatch <- match(paste0(featuresUse$seqnames, ":", featuresUse$idx), paste0(rS$seqnames, ":", rS$idx))
featuresUse$Bins <- as.vector(rS$Bins[idxMatch])

#MakeLists
featureList <- split(featuresUse$Match, featuresUse$Module)
moduleList <- split(featuresUse$Bins, featuresUse$Module)
binList <- split(rS$Match, rS$Bins)

dfM <- lapply(seq_along(featureList), function(x){
  message("Computing Module ",x, " of ", length(featureList))
  binx <- binList[moduleList[[x]]]
  idxFgd <- featureList[[x]]
  idxBgd <- unlist(lapply(binx, function(x) sample(x, nBgd)), use.names=FALSE)
  m <- ArchR:::.getPartialMatrix(
    ArrowFiles = getArrowFiles(projMulti2),
    featureDF = featureDF[c(idxFgd, idxBgd), ],
    useMatrix = useMatrix,
    cellNames = projMulti2$cellNames,
    threads = threads,
    verbose = FALSE,
    doSampleCells = FALSE
  )
  Matrix::colMeans(m[seq_along(idxFgd), ]) - Matrix::colMeans(m[-seq_along(idxFgd), ])
}) %>% Reduce("cbind", .)

for(x in seq_len(ncol(dfM))){
  projMulti2 <- addCellColData(projMulti2, data = dfM[,x], name=names(featureList)[x], cells=rownames(dfM), force = TRUE)
}

projMulti2

projMulti2 <- addImputeWeights(projMulti2,reducedDims = "LSI_Combined")


plots <- lapply(names(featureList), FUN=function(x){
  plotEmbedding(projMulti2, name=x,embedding = "UMAP_Combined",
                imputeWeights = getImputeWeights(projMulti2))
})

plots[[14]] <- plotEmbedding(projMulti2, name = "Clusters_Combined",
                             embedding = "UMAP_Combined", size = 1,
                             labelAsFactors=F, labelMeans=T)

plotPDF(ggAlignPlots(plots,draw=F,type="h"))
png("/gstore/scratch/u/lucast3/fibroticmemory/Plots/CellScores.png",
    width=2400,height=1400)
plot_grid(plotlist = plots,align = "h",
          nrow = 3)
dev.off()


map <- setNames(paste0("C",1:13),
                c("Lipofibroblast/Mesothelial","Other","Adventitial/Pericyte","SMC",
                  "Myofibroblast","Lipofibroblast","Mesotheial","Mesotheial",
                  "Mesotheial","Lipofibroblast",
                  "Other","Alveolar","Alveolar"))
celltype <- factor(projMulti2$Clusters_Combined)
levels(celltype)[match(map,levels(celltype))] <- names(map)

projMulti2$celltype <- unfactor(celltype)

png("/gstore/scratch/u/lucast3/fibroticmemory/Plots/Celltypes.png",
    width=500,height=500)
plotEmbedding(projMulti2, name = "celltype",
                             embedding = "UMAP_Combined", size = 1,
                             labelAsFactors=F, labelMeans=T)
dev.off()
saveArchRProject(ArchRProj = projMulti2,  load = FALSE)

markerGenes <- c(
  "PI16", #Adventitial
  "CTHRC1", #Myofibroblast
  "CES1"
  
)

p <- plotEmbedding(
  ArchRProj = projMulti2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP_Combined",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
