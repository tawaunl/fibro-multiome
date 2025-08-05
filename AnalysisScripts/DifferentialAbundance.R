#load the ArchR library
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
#to have access to chromosome and gene information. ArchR supports hg19, hg38, mm9, and mm10.
addArchRGenome("hg38")
data.dir <- "/gstore/scratch/u/lucast3/fibroticmemory" # setdirectory for files 

## ----Load Project--------------------------------------------------------------
projMulti2 <- loadArchRProject(file.path(data.dir,"Save-ProjMulti2"),
                               showLogo = FALSE)

#Refine Clusters post Harmony ===============================
## iterate over a few different resoultions for clustering 

resoultions = c(0.05,0.1,0.15,0.2,0.25,0.3)
umaps <- list()
for(resolution in resoultions){
  projMulti2 <- addClusters(projMulti2, reducedDims = "Harmony",
                            name = paste0("Clusters_Combined_Harmony_res",resolution),
                                         resolution = resolution,
                                         force = TRUE)
  umaps[[paste0("Res_",resolution)]] <- plotEmbedding(projMulti2,
                                                      name = paste0("Clusters_Combined_Harmony_res",resolution),
                                                      embedding = "UMAP_Combined_Harmony",
                                                      size = 1,
                                                      labelAsFactors=F, labelMeans=F)
  
}

do.call(cowplot::plot_grid, c(list(ncol = 3),umaps))

# Stacked Bar plots
## ignoring treatment ===========================
library(ggplot2)
library(wesanderson)
library(grid)
library(gridExtra)
library(tidyverse)
library(dplyr)
coldata <- as.data.frame(getCellColData(ArchRProj = projMulti2))
study_cluster_tmp_by_study <- coldata %>%
  group_by(Group, Clusters_Combined_Harmony_res0.05) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

cluster1 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$Clusters_Combined_Harmony_res0.05=="C1"),]
cluster2 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$Clusters_Combined_Harmony_res0.05=="C2"),]
cluster3 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$Clusters_Combined_Harmony_res0.05=="C3"),]

cluster1$scaled <- cluster1$freq/sum(cluster1$freq)
cluster2$scaled <- cluster2$freq/sum(cluster2$freq)
cluster3$scaled <- cluster3$freq/sum(cluster3$freq)


study_cluster_tmp_by_study_scaled <- rbind(cluster1,cluster2,cluster3)


ggplot(study_cluster_tmp_by_study_scaled, aes(fill=Group, y=scaled, x=Clusters_Combined_Harmony_res0.05)) +
  geom_bar(position="stack", stat="identity") + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within condition") +
  theme(axis.title.y =  element_text(angle = 90))

## Stack by Treatment =============================

study_cluster_tmp_by_study <- coldata %>%
  group_by(Treatment, Clusters_Combined_Harmony_res0.05) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

cluster1 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$Clusters_Combined_Harmony_res0.05=="C1"),]
cluster2 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$Clusters_Combined_Harmony_res0.05=="C2"),]
cluster3 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$Clusters_Combined_Harmony_res0.05=="C3"),]

cluster1$scaled <- cluster1$freq/sum(cluster1$freq)
cluster2$scaled <- cluster2$freq/sum(cluster2$freq)
cluster3$scaled <- cluster3$freq/sum(cluster3$freq)


study_cluster_tmp_by_study_scaled <- rbind(cluster1,cluster2,cluster3)
ggplot(study_cluster_tmp_by_study_scaled, aes(fill=Treatment, y=scaled, x=Clusters_Combined_Harmony_res0.05)) +
  geom_bar(position="stack", stat="identity") + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within condition") +
  theme(axis.title.y =  element_text(angle = 90))


## subset just by no treatment=========
selection <- c("IPF-NoTreatment","Normal-NoTreatment")
subset <- study_cluster_tmp_by_study_scaled[which(study_cluster_tmp_by_study_scaled$Treatment %in% selection),]
ggplot(subset, aes(fill=Treatment, y=scaled, x=Clusters_Combined_Harmony_res0.05)) +
  geom_bar(position="stack", stat="identity") + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within condition") +
  theme(axis.title.y =  element_text(angle = 90))

# C. Abundance Plots ===========================================================
library(edgeR)
library(RColorBrewer)
library(viridis)
abundances <- table(projMulti2$Clusters_Combined_Harmony_res0.05,projMulti2$Sample) 
abundances <- unclass(abundances) 

extra.info <- coldata[match(colnames(abundances), projMulti2$Sample),]
d <- DGEList(abundances, samples=extra.info)
d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)

norm_counts <- as.data.frame(t(d$counts)) 
colnames(norm_counts) <- paste("Cluster ", colnames(norm_counts), sep="")
norm_counts <- norm_counts %>% mutate(Sample = rownames(.))

coldata_short <- coldata %>% dplyr::select(Sample,Treatment,Group) %>% unique


df_long_final <- left_join(norm_counts,coldata_short, by="Sample") 

percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))

sums <- data.frame( clustersums=rowSums(df_long_final[,1:length(levels(factor(projMulti2$Clusters_Combined_Harmony_res0.05)))]),
                    ident=df_long_final$Sample
)

for (clust in 1:length(levels(factor(projMulti2$Clusters_Combined_Harmony_res0.05)))) {
  for (sample in 1:length(df_long_final[,clust])) {
    percent <- df_long_final[sample,clust]/sums[sample,1]*100
    percentages[sample,clust]<- percent
  }
}

percentages[,(length(levels(factor(projMulti2$Clusters_Combined_Harmony_res0.05)))+1):(dim(df_long_final)[2])] <- df_long_final[,(length(levels(factor(projMulti2$Clusters_Combined_Harmony_res0.05)))+1):(dim(df_long_final)[2])]


colnames(percentages) <-  colnames(df_long_final)

df_long <- percentages %>% 
  pivot_longer(c(c(paste0("Cluster C",1:length(levels(factor(projMulti2$Clusters_Combined_Harmony_res0.05)))))))

colnames(df_long) <- c("Sample","Treatment","Group","Cluster","Percent")

level_order <- c("Normal-NoTreatment","IPF-NoTreatment",
                 "Normal-OSM","IPF-OSM",
                 "Normal-TGF","IPF-TGF",
                 "Normal-IL13","IPF-IL13")
df_long$Treatment <- factor(df_long$Treatment,levels = level_order) 
df_long$Cluster <- factor(df_long$Cluster,levels =c(paste0("Cluster C",1:(length(levels(factor(projMulti2$Clusters_Combined_Harmony_res0.05)))))))


ggplot(df_long, aes(x=Treatment, y=Percent,fill=Treatment)) + theme_classic() +
  geom_col() + scale_fill_viridis(discrete = TRUE, alpha=1)   +
  facet_wrap(~factor(Cluster),scales = "free", ncol=4) + xlab("Treatment") +
  theme(axis.text.x = element_text(angle=75,vjust = 0.5),axis.text.y = element_text(size=20)) +
  ylab("Cellularity Percent") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),legend.text = element_text(size=14),
        legend.title = element_text(size=16)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50,75)) # Change font size

