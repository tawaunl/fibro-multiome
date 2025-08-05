#PIM1 Analysis

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
data.dir <- "/gstore/scratch/u/lucast3/fibroticmemory" # set directory for files ]


# Load Data ============================================================
projMulti2 <- loadArchRProject(file.path(data.dir,"Save-ProjMulti2"),
                               showLogo = FALSE)

load(file = file.path(data.dir,"overlaps.Rdata"))
library("org.Hs.eg.db")

markerGenes  <- c(Reduce(intersect, pos_list)
                  ,Reduce(intersect, neg_list))
level_order <- c("Normal-NoTreatment","Normal-OSM","Normal-TGF","Normal-IL13",
                 "IPF-NoTreatment","IPF-OSM","IPF-TGF","IPF-IL13")
x <- factor(projMulti2$Treatment,levels = level_order)

new_levels <- c("1_Normal-NoTreatment","2_Normal-OSM","3_Normal-TGF","4_Normal-IL13",
                "5_IPF-NoTreatment","6_IPF-OSM","7_IPF-TGF","8_IPF-IL13")
levels(x) <- new_levels


projMulti2$Treatment <- as.vector(x)
p <- plotBrowserTrack(
  ArchRProj = projMulti2, 
  groupBy = "Treatment", 
  geneSymbol = markerGenes, 
  upstream = 20000,
  downstream = 20000
)

p <- plotBrowserTrack(
  ArchRProj = projMulti2, 
  groupBy = "Treatment", 
  geneSymbol = "PIM1", 
  upstream = 30000,
  downstream = 30000
)
grid::grid.newpage()
grid::grid.draw(p$PIM1)

p <- plotBrowserTrack(
  ArchRProj = projMulti2, 
  groupBy = "Treatment", 
  geneSymbol = "IL17B", 
  upstream = 60000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$IL17B)

plotGroups(
  ArchRProj = projMulti2, 
  groupBy = "Treatment", 
  colorBy = "GeneScoreMatrix", 
  name = "PIM1",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

grid::grid.draw(p$PAQR5)
grid::grid.draw(p$CYP2S1)

plotGroups(
  ArchRProj = projMulti2, 
  groupBy = "Treatment", 
  colorBy = "GeneScoreMatrix", 
  name = "PIM1",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotGroups(
  ArchRProj = projMulti2, 
  groupBy = "Treatment", 
  colorBy = "GeneScoreMatrix", 
  name = "NFATC1",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)



plotTSSEnrichment(projMulti2,
                  groupBy="Treatment")
plotTSSEnrichment(projMulti2,
                  groupBy="Group")