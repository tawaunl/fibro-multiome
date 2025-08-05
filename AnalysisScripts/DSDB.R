# Place dataset into DSDB
library(ArchR)
library(DataSetDB)
library(dsassembly)
library(Matrix)
library(SingleCellExperiment)

#set a seed to facilitate replication of operations requiring randomization
set.seed(824)

#default number of Parallel threads is 16
# working on a local computer, 1 thread works best
addArchRThreads(7)

#Before we begin, we need add a reference genome annotation for ArchR 
#to have access to chromosome and gene information. ArchR supports hg19, hg38, mm9, and mm10.
addArchRGenome("hg38")
data.dir <- "/gstore/scratch/u/lucast3/fibroticmemory" # setdirectory for files 



## ----Load Project--------------------------------------------------------------
projMulti2 <- loadArchRProject(file.path(data.dir,"Save-ProjMulti2"),
                               showLogo = FALSE)


# Get GeneExpressionMatrix =================================
ge <- getMatrixFromProject(projMulti2,
                           useMatrix = "GeneExpressionMatrix",
                           threads = getArchRThreads())

names(assays(ge)) <- "counts"
assay(ge, "log2cpm") <- gp.sa.diff::normalizedCPM(assay(ge, "counts"), log=TRUE) # add log2CPM

ge <- as(ge, "SingleCellExperiment")

ge <- annotateExperiment(ge,
                         title = "GeneExpression Matrix",
                         description = "Gene Expression from 10x Arc workflow",
                         annotation =NULL,
                         organism = "Homo sapiens",
                         namespace = list(list(type="genome",id="GRCh38")),
                         technology = list(name = "scRNA-seq",
                                           terms=c("EFO:0008913", #scRNA
                                                   "EFO:0030059",#10xMultiome
                                                   "EFO:0008563"), #Hiseq 4000
                                           details = "Illumina HiSeq"),
                         sources = list(list(name="FireDB", id="FRS18379"))
)
ge$Sample <- as.vector(ge$Sample)
rowData(ge)$seqnames <- as.vector(rowData(ge)$seqnames)

## Add GeneScoreMatrix as Alternative experiment ---------------
gs <- getMatrixFromProject(projMulti2,useMatrix = "GeneScoreMatrix",
                           threads = getArchRThreads())
gs <- annotateExperiment(gs,
                         title = "ArchR GeneScore Matrix",
                         description = "GeneScoreMatrix from the scATAC-seq ArchR workflow",
                         annotation =NULL,
                         organism = "Homo sapiens",
                         namespace = list(list(type="genome",id="GRCh38")),
                         technology = list(name = "scATAC-seq",
                                           terms=c( "EFO:0030007",#scATAC
                                                    "EFO:0030059",#10xMultiome
                                                    "EFO:0008563"), #Hiseq 4000, 
                                           details = "Illumina HiSeq"),
                         sources = list(list(name="FireDB", id="FRS18379")))

gs$Sample <- as.vector(gs$Sample)
rowData(gs)$seqnames <- as.vector(rowData(gs)$seqnames)


altExp(ge, "ArchR GeneScoreMatrix") <- gs

## Add Tile Matrix==========================================
tilematrix <- getMatrixFromProject(
  ArchRProj = projMulti2,
  useMatrix = "TileMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

tilematrix <- annotateExperiment(tilematrix,
                                 title = "ArchR TileMatrix",
                                 description = "TileMatrix from the scATAC-seq ArchR workflow",
                                 annotation =NULL,
                                 organism = "Homo sapiens",
                                 namespace = list(list(type="genome",id="GRCh38")),
                                 technology = list(name = "scATAC-seq",
                                                   terms=c( "EFO:0030007",#scATAC
                                                            "EFO:0030059",#10xMultiome
                                                            "EFO:0008563"), #Hiseq 4000,),
                                                   details = "Illumina HiSeq"),
                                 sources = list(list(name="FireDB", id="FRS18379"))
)
tilematrix$Sample <- as.vector(tilematrix$Sample)
rowData(tilematrix)$seqnames <- as.vector(rowData(tilematrix)$seqnames)

altExp(ge, "ArchR TileMatrix") <- tilematrix

# Add in peaks and motifs as seperate assays
peaks <- getMatrixFromProject(
  ArchRProj = projMulti2,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

Motifs <- getMatrixFromProject(
  ArchRProj = projMulti2,
  useMatrix = "MotifMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

peaks <- annotateExperiment(peaks,
                            title = "ArchR Peaks Matrix",
                            description = "Peaks counts from the scATAC-seq ArchR workflow",
                            annotation =NULL,
                            organism = "Homo sapiens",
                            namespace = list(list(type="genome",id="GRCh38")),
                            technology = list(name = "scATAC-seq",
                                              terms=c( "EFO:0030007",#scATAC
                                                       "EFO:0030059",#10xMultiome
                                                       "EFO:0008563"), #Hiseq 4000,),
                                              details = "Illumina HiSeq"),
                            sources = list(list(name="FireDB", id="FRS18379")))


Motifs <- annotateExperiment(Motifs,
                             title = "ArchR Motifs Matrix",
                             description = "Motif annotations from the scATAC-seq ArchR workflow",
                             annotation =NULL,
                             organism = "Homo sapiens",
                             namespace = list(list(type="genome",id="GRCh38")),
                             technology = list(name = "scATAC-seq",
                                               terms=c( "EFO:0030007",#scATAC
                                                        "EFO:0030059",#10xMultiome
                                                        "EFO:0008563"), #Hiseq 4000,),
                                               details = "Illumina HiSeq"),
                             sources = list(list(name="FireDB", id="FRS18379")))






peaks$Sample <- as.vector(peaks$Sample)
rowData(peaks)$seqnames <- as.vector(rowData(peaks)$seqnames)

Motifs$Sample <- as.vector(Motifs$Sample)
rowData(Motifs)$seqnames <- as.vector(rowData(Motifs)$seqnames)

altExp(ge, "ArchR PeaksMatrix") <- peaks
altExp(ge, "ArchR MotifsMatrix") <- Motifs


atac.mae <- MultiAssayExperiment(experiments= list(multiome=ge))

desc=c("Prcessed data from DS000015342,
       
       Idiopathic pulmonary fibrosis (IPF) is a progressive interstitial lung disease characterized by scaring of the lungs that eventually leads to respiratory failure. Repetitive lung injury due to chronic exposure to toxic chemicals, pathogens and various environmental factors can lead to fibrotic lung condition. However, the pathogenesis of IPF is still unclear. Preliminary studies in our lab showed that the IPF-donor derived lung fibroblasts had elevated levels of total STAT3 compared to healthy-donor derived lung fibroblasts, which can induce several inflammatory signaling pathways downstream of STAT3. By performing scRNA-sequencing and scATAC-sequencing, we aim to determine whether the IPF-donor derived lung fibroblasts have epigenetic memory that allows these cells to robustly induce inflammatory signals when stimulated with cytokines including OSM, TGF-beta and IL-13.")
atac.mae <- annotateDataset(atac.mae, 
                            title = "NGS4921 - IPF FIBROBLASTS VS NORMAL FIBROBLAST TRANSCRIPTIONAL AND EPIGENETIC PROFILE WITH OR WITHOUT TREATMENT (ASSOCIATED WITH NGS4701) [GARFIELD,SIMOPOUC]",
                            description = desc,
                            authors = c("lucast3","kharwadr"))
options(ArtifactDB.upload.size.limit = 5)
#redirecting from 'GMTY157:GRCh38_tiles500_filtered@REVISION-7' to 'GMTY157:GRCh38_tiles500_filtered@60cc587226e68983d031e6a79532486aff136cd7'
#redirecting from 'GMTY159:GRCh38_genescorefeatures@REVISION-5' to 'GMTY159:GRCh38_genescorefeatures@89ead38f0bd668cd7776b40969f85c5877ad100a'
#redirecting from 'GMTY167:GRCh38-exonic@REVISION-3' to 'GMTY167:GRCh38-exonic@54a8c2c0ebf8e1540acf6f75b479223f6c8a2a1c'

saveDataset(atac.mae)


