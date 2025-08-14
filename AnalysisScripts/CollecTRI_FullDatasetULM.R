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

se <- getMatrixFromProject(archr_proj,useMatrix = "GeneExpressionMatrix")
library(scater)
se <- logNormCounts(se,assay.type = "GeneExpressionMatrix")
mat <- se@assays@data@listData[["logcounts"]]

rownames(mat) <- rowData(se)$name
mat <- as.matrix(mat)

net <- decoupleR::get_collectri(organism = 'human',
                                split_complexes = FALSE)

# Run ulm
acts <- decoupleR::run_ulm(mat = mat ,
                           net = net,
                           .source = 'source',
                           .target = 'target',
                           .mor='mor',
                           minsize = 5)
saveRDS(acts,"/gstore/scratch/u/lucast3/fibroMultiome/actsFull.rds")
saveRDS(se,"/gstore/scratch/u/lucast3/fibroMultiome/se_ExpressionFull.rds")
