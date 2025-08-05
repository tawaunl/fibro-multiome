# Get marker Peaks

# Load Libraries and set up enviroment ======================================

library(ArchR)
library(grid)
library(gridExtra)
#library(dsassembly)
#library(gg4way)
#set a seed to facilitate replication of operations requiring randomization
set.seed(824)

#default number of Parallel threads is 16
# working on a local computer, 1 thread works best
addArchRThreads(14)

#Before we begin, we need add a reference genome annotation for ArchR 
addArchRGenome("hg38")
data.dir <- "/gstore/scratch/u/lucast3/fibroticmemory" # set directory for files ]


# Load Data ============================================================
projMulti2 <- loadArchRProject(file.path(data.dir,"Save-ProjMulti2"),
                               showLogo = FALSE)

#Treatment ====================================================================
markersPeaks <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Treatment",
  bias = c("TSSEnrichment", "log10(nFrags)","log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
saveRDS(markersPeaks,file.path(data.dir,"Save-ProjMulti2/MarkerPeaks_Treatment.rds"))
saveRDS(markerList,file.path(data.dir,"Save-ProjMulti2/MarkerPeaks_Treatment_list.rds"))

# Group ======================================================================
markersPeaks <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Group",
  bias = c("TSSEnrichment", "log10(nFrags)","log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")

saveRDS(markersPeaks,file.path(data.dir,"Save-ProjMulti2/MarkerPeaks_Group.rds"))

#Plotting ========================================================
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot"),
        name="Heatmap-TreatmentPeaks",addDOC = FALSE)


heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,plotLog2FC = TRUE
)

plotPDF(plot(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot"),
        name="Heatmap-GroupPeaks",addDOC = FALSE)

library("org.Hs.eg.db")

# DA analysis between groups =================================

markersPeaks <- readRDS(file.path(data.dir,"Save-ProjMulti2/MarkerPeaks_Treatment.rds"))

## Add TFBS from encode so we can find overlaps in interesting TFS
projMulti2 <- addArchRAnnotations(ArchRProj = projMulti2, collection = "EncodeTFBS")
# Pair wise testing ==========================
# IL13 =========================================================================
markersIL13_IPF <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("IPF-IL13"),
  bgdGroups = c("IPF-NoTreatment"),
  useMatrix = "PeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
list <- getMarkers(markersIL13_IPF, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.25")
IL13_IPF <- as.data.frame(unlist(list))

pv <- plotMarkers(seMarker = markersIL13_IPF, 
                 name = "IPF-IL13", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "IPF-IL13-vs-IPF-NoTreatment-PeakMarkers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)

saveRDS(markersIL13_IPF,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/Peakmarkers_IPF_IL13.rds")))

motifsUpIL13_IPF <- peakAnnoEnrichment(
  seMarker = markersIL13_IPF,
  ArchRProj = projMulti2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUpIL13_IPF), mlog10Padj = assay(motifsUpIL13_IPF)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

## Normal =====================================
markersIL13_Normal <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("Normal-IL13"),
  bgdGroups = c("Normal-NoTreatment"),
  useMatrix = "PeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
list <- getMarkers(markersIL13_Normal, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.25")
IL13_Normal <- as.data.frame(unlist(list))
pv <- plotMarkers(seMarker = markersIL13_Normal, 
                 name = "Normal-IL13", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "Normal-IL13-vs-Normal-NoTreatment-PeakMarkers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)
saveRDS(markersIL13_Normal,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/Peakmarkers_Normal_IL13.rds")))


### 4way =====================================================================
IL13_IPF$peak <- paste0(IL13_IPF$seqnames,":",
                        IL13_IPF$idx,":",
                        IL13_IPF$start,"-",IL13_IPF$end)
IL13_Normal$peak <- paste0(IL13_Normal$seqnames,":",
                           IL13_Normal$idx,":",
                           IL13_Normal$start,"-",IL13_Normal$end)

commonGenes <- intersect(IL13_IPF$peak,IL13_Normal$peak)
sub1 <- IL13_IPF[which(IL13_IPF$peak %in% commonGenes==TRUE),]
sub2 <- IL13_Normal[which(IL13_Normal$peak %in% commonGenes==TRUE),]
df <- data.frame(sub1,sub2)
df$gene <- df$peak
rsq <- round(cor(x=df$Log2FC,y=df$Log2FC.1,use= "complete.obs" )^2,3)

comparisonName="IL13_IPF vs. IL13_Normal"
xlab= "IPF_IL13  vs. IPF_NoTreatment"
ylab <- "Normal_IL13 vs. Normal_NoTreament"
xsig <- df %>% filter(FDR<=0.1,abs(Log2FC)>=1)
ysig <- df %>% filter(FDR.1<=0.1,abs(Log2FC.1)>=1)

both.sig <- intersect(xsig,ysig)

x.uniq <- setdiff(xsig,ysig)
y.uniq <- setdiff(ysig,xsig)

cols <- c("red","blue","green","darkgrey")
cols <-setNames(cols,
                c("Sig-Both",
                  paste0("Sig-",strsplit(comparisonName," ")[[1]][1]),
                  paste0("Sig-",strsplit(comparisonName," ")[[1]][3]),
                  "NS"))


plot <- ggplot(df,aes(y=Log2FC.1,x=Log2FC,colour="NS")) + geom_point(size=3) 

if(dim(x.uniq)[1]>0){
  plot <- plot + geom_point(data=x.uniq, 
                            aes(x=Log2FC,
                                color=paste0("Sig-",strsplit(comparisonName," ")[[1]][1])),
                            size=4)
}
if(dim(y.uniq)[1]>0){
  plot <- plot + geom_point(data=y.uniq, 
                            aes(y=Log2FC.1,
                                color=paste0("Sig-",strsplit(comparisonName," ")[[1]][3])),
                            size=4)
}
if(dim(both.sig)[1]>0){
  plot <- plot + geom_point(data=both.sig, 
                            aes(y=Log2FC.1,x=Log2FC, color="Sig-Both"),
                            size=5)
}

plot <- plot + scale_colour_manual(name="Significance",
                                   values=cols,
                                   guide = guide_legend(fill = NULL,colour = NULL,
                                                        override.aes = list(size = 6))) +
  geom_label_repel(data = subset(df,subset= c(abs(Log2FC) >0.5,abs(Log2FC.1)>0.5)),
                   aes(label=gene,color="black"),
                   color="black",size =4,
                   box.padding   = 0.35,
                   point.padding = 0.5) +
  labs(title=paste0(comparisonName,' comparison'),
       subtitle= paste0('R2 = ',rsq,"\n",
                        "# of DEGs = ",dim(both.sig)[1],"\n",
                        paste("# of",strsplit(comparisonName," ")[[1]][1],"DEGs =" ,dim(xsig)[1]),
                        "\n",
                        paste("# of",strsplit(comparisonName," ")[[1]][3],"DEGs =" ,dim(ysig)[1])),
       x= paste("log FoldChange",xlab), y= paste("log FoldChange",ylab),
       caption ="Signicance: -1 \u2264 LogFC \u2265 1 , Pvalue \u2264 0.05 ") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.background = element_rect(fill= "white",colour = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        plot.caption = element_text(size=14),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

uniqueIL13 <- x.uniq
png(file.path(data.dir,"Save-ProjMulti2/Plots/IL13_4wayPeaks.png"))
print(plot)
dev.off()

# OSM =========================================================================
markersOSM_IPF <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("IPF-OSM"),
  bgdGroups = c("IPF-NoTreatment"),
  useMatrix = "PeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

list <- getMarkers(markersOSM_IPF, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.25")
OSM_IPF <- as.data.frame(unlist(list))
pv <- plotMarkers(seMarker = markersOSM_IPF, 
                 name = "IPF-OSM", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "IPF-OSM-vs-IPF-NoTreatment-PeakMarkers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)

saveRDS(markersOSM_IPF,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/Peakmarkers_IPF_OSM.rds")))

motifsUpOSM_IPF <- peakAnnoEnrichment(
  seMarker = markersOSM_IPF,
  ArchRProj = projMulti2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)


markersOSM_Normal <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("Normal-OSM"),
  bgdGroups = c("Normal-NoTreatment"),
  useMatrix = "PeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

list <- getMarkers(markersOSM_Normal, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.25")
OSM_Normal <- as.data.frame(unlist(list))
pv <- plotMarkers(seMarker = markersOSM_Normal, 
                 name = "Normal-OSM", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "Normal-OSM-vs-Normal-NoTreatment-PeakMarkers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)
saveRDS(markersOSM_Normal,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/Peakmarkers_Normal_OSM.rds")))
motifsUpOSM_IPF <- peakAnnoEnrichment(
  seMarker = markersOSM_Normal,
  ArchRProj = projMulti2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
### 4way =====================================================================
OSM_IPF$peak <- paste0(OSM_IPF$seqnames,":",
                       OSM_IPF$idx,":",
                       OSM_IPF$start,"-",OSM_IPF$end)
OSM_Normal$peak <- paste0(OSM_Normal$seqnames,":",
                          OSM_Normal$idx,":",
                          OSM_Normal$start,"-",OSM_Normal$end)

commonGenes <- intersect(OSM_IPF$peak,OSM_Normal$peak)
sub1 <- OSM_IPF[which(OSM_IPF$peak %in% commonGenes==TRUE),]
sub2 <- OSM_Normal[which(OSM_Normal$peak %in% commonGenes==TRUE),]
df <- data.frame(sub1,sub2)
df$gene <- df$peak
rsq <- round(cor(x=df$Log2FC,y=df$Log2FC.1,use= "complete.obs" )^2,3)

comparisonName="OSM_IPF vs. OSM_Normal"
xlab= "IPF_OSM  vs. IPF_NoTreatment"
ylab <- "Normal_OSM vs. Normal_NoTreament"
xsig <- df %>% filter(FDR<=0.1,abs(Log2FC)>=1)
ysig <- df %>% filter(FDR.1<=0.1,abs(Log2FC.1)>=1)

both.sig <- intersect(xsig,ysig)

x.uniq <- setdiff(xsig,ysig)
y.uniq <- setdiff(ysig,xsig)

cols <- c("red","blue","green","darkgrey")
cols <-setNames(cols,
                c("Sig-Both",
                  paste0("Sig-",strsplit(comparisonName," ")[[1]][1]),
                  paste0("Sig-",strsplit(comparisonName," ")[[1]][3]),
                  "NS"))


plot <- ggplot(df,aes(y=Log2FC.1,x=Log2FC,colour="NS")) + geom_point(size=3) 

if(dim(x.uniq)[1]>0){
  plot <- plot + geom_point(data=x.uniq, 
                            aes(x=Log2FC,
                                color=paste0("Sig-",strsplit(comparisonName," ")[[1]][1])),
                            size=4)
}
if(dim(y.uniq)[1]>0){
  plot <- plot + geom_point(data=y.uniq, 
                            aes(y=Log2FC.1,
                                color=paste0("Sig-",strsplit(comparisonName," ")[[1]][3])),
                            size=4)
}
if(dim(both.sig)[1]>0){
  plot <- plot + geom_point(data=both.sig, 
                            aes(y=Log2FC.1,x=Log2FC, color="Sig-Both"),
                            size=5)
}

plot <- plot + scale_colour_manual(name="Significance",
                                   values=cols,
                                   guide = guide_legend(fill = NULL,colour = NULL,
                                                        override.aes = list(size = 6))) +
  geom_label_repel(data = subset(df,subset= c(abs(Log2FC) >0.5,abs(Log2FC.1)>0.5)),
                   aes(label=gene,color="black"),
                   color="black",size =4,
                   box.padding   = 0.35,
                   point.padding = 0.5) +
  labs(title=paste0(comparisonName,' comparison'),
       subtitle= paste0('R2 = ',rsq,"\n",
                        "# of DEGs = ",dim(both.sig)[1],"\n",
                        paste("# of",strsplit(comparisonName," ")[[1]][1],"DEGs =" ,dim(xsig)[1]),
                        "\n",
                        paste("# of",strsplit(comparisonName," ")[[1]][3],"DEGs =" ,dim(ysig)[1])),
       x= paste("log FoldChange",xlab), y= paste("log FoldChange",ylab),
       caption ="Signicance: -1 \u2264 LogFC \u2265 1 , Pvalue \u2264 0.05 ") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.background = element_rect(fill= "white",colour = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        plot.caption = element_text(size=14),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

uniqueOSM <- x.uniq
png(file.path(data.dir,"Save-ProjMulti2/Plots/OSM_4wayPeaks.png"))
print(plot)
dev.off()

# TGF =========================================================================
markersTGF_IPF <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("IPF-TGF"),
  bgdGroups = c("IPF-NoTreatment"),
  useMatrix = "PeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
list <- getMarkers(markersTGF_IPF, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.25")
TGF_IPF <- as.data.frame(unlist(list))

pv <- plotMarkers(seMarker = markersTGF_IPF, 
                 name = "IPF-TGF", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "IPF-TGF-vs-IPF-NoTreatment-PeakMarkers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)

saveRDS(markersTGF_IPF,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/Peakmarkers_IPF_TGF.rds")))

markersTGF_Normal <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("Normal-TGF"),
  bgdGroups = c("Normal-NoTreatment"),
  useMatrix = "PeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
list <- getMarkers(markersTGF_Normal, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.25")
TGF_Normal <- as.data.frame(unlist(list))
pv <- plotMarkers(seMarker = markersTGF_Normal, 
                 name = "Normal-TGF", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                 plotAs = "Volcano")
plotPDF(pv, name = "Normal-TGF-vs-Normal-NoTreatment-PeakMarkers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)
saveRDS(markersTGF_Normal,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/Peakmarkers_Normal_TGF.rds")))

### 4way =====================================================================
TGF_IPF$peak <- paste0(TGF_IPF$seqnames,":",
                       TGF_IPF$idx,":",
                       TGF_IPF$start,"-",TGF_IPF$end)
TGF_Normal$peak <- paste0(TGF_Normal$seqnames,":",
                          TGF_Normal$idx,":",
                          TGF_Normal$start,"-",TGF_Normal$end)

commonGenes <- intersect(TGF_IPF$peak,TGF_Normal$peak)
sub1 <- TGF_IPF[which(TGF_IPF$peak %in% commonGenes==TRUE),]
sub2 <- TGF_Normal[which(TGF_Normal$peak %in% commonGenes==TRUE),]
df <- data.frame(sub1,sub2)
df$gene <- df$peak
rsq <- round(cor(x=df$Log2FC,y=df$Log2FC.1,use= "complete.obs" )^2,3)

comparisonName="TGF_IPF vs. TGF_Normal"
xlab= "IPF_TGF  vs. IPF_NoTreatment"
ylab <- "Normal_TGF vs. Normal_NoTreament"
xsig <- df %>% filter(FDR<=0.1,abs(Log2FC)>=1)
ysig <- df %>% filter(FDR.1<=0.1,abs(Log2FC.1)>=1)

both.sig <- intersect(xsig,ysig)

x.uniq <- setdiff(xsig,ysig)
y.uniq <- setdiff(ysig,xsig)

cols <- c("red","blue","green","darkgrey")
cols <-setNames(cols,
                c("Sig-Both",
                  paste0("Sig-",strsplit(comparisonName," ")[[1]][1]),
                  paste0("Sig-",strsplit(comparisonName," ")[[1]][3]),
                  "NS"))


plot <- ggplot(df,aes(y=Log2FC.1,x=Log2FC,colour="NS")) + geom_point(size=3) 

if(dim(x.uniq)[1]>0){
  plot <- plot + geom_point(data=x.uniq, 
                            aes(x=Log2FC,
                                color=paste0("Sig-",strsplit(comparisonName," ")[[1]][1])),
                            size=4)
}
if(dim(y.uniq)[1]>0){
  plot <- plot + geom_point(data=y.uniq, 
                            aes(y=Log2FC.1,
                                color=paste0("Sig-",strsplit(comparisonName," ")[[1]][3])),
                            size=4)
}
if(dim(both.sig)[1]>0){
  plot <- plot + geom_point(data=both.sig, 
                            aes(y=Log2FC.1,x=Log2FC, color="Sig-Both"),
                            size=5)
}

plot <- plot + scale_colour_manual(name="Significance",
                                   values=cols,
                                   guide = guide_legend(fill = NULL,colour = NULL,
                                                        override.aes = list(size = 6))) +
  geom_label_repel(data = subset(df,subset= c(abs(Log2FC) >0.5,abs(Log2FC.1)>0.5)),
                   aes(label=gene,color="black"),
                   color="black",size =4,
                   box.padding   = 0.35,
                   point.padding = 0.5) +
  labs(title=paste0(comparisonName,' comparison'),
       subtitle= paste0('R2 = ',rsq,"\n",
                        "# of DEGs = ",dim(both.sig)[1],"\n",
                        paste("# of",strsplit(comparisonName," ")[[1]][1],"DEGs =" ,dim(xsig)[1]),
                        "\n",
                        paste("# of",strsplit(comparisonName," ")[[1]][3],"DEGs =" ,dim(ysig)[1])),
       x= paste("log FoldChange",xlab), y= paste("log FoldChange",ylab),
       caption ="Signicance: -1 \u2264 LogFC \u2265 1 , Pvalue \u2264 0.05 ") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.background = element_rect(fill= "white",colour = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        plot.caption = element_text(size=14),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

uniqueTGF <- x.uniq

png(file.path(data.dir,"Save-ProjMulti2/Plots/TGF_4wayPeaks.png"))
print(plot)
dev.off()
# Get Overlapping Genes =======================================================
library(VennDiagram)

pos_list <- list(IL13=uniqueIL13$gene[which(uniqueIL13$Log2FC > 0)],
                 OSM=uniqueOSM$gene[which(uniqueOSM$Log2FC > 0 )],
                 TGF=uniqueTGF$gene[which(uniqueTGF$Log2FC > 0)])

neg_list <- list(IL13=uniqueIL13$gene[which(uniqueIL13$Log2FC < 0)],
                 OSM=uniqueOSM$gene[which(uniqueOSM$Log2FC < 0)],
                 TGF=uniqueTGF$gene[which(uniqueTGF$Log2FC < 0)])
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


se <- getMarkerFeatures(ArchRProj = projMulti2, useMatrix = "PeakMatrix",
                        groupBy = "Treatment")
subsetSE <- se[which(rowData(se)$name %in% c(Reduce(intersect, pos_list)
                                             ,Reduce(intersect, neg_list))),]
plotMarkerHeatmap(seMarker = subsetSE)
gr <- GRanges(
  seqnames = Rle(c("chr5"), c(1)),
  ranges = IRanges(start=149347487, end = 149447987))


plotGroups(
  ArchRProj = projMulti2, 
  groupBy = "Treatment", 
  colorBy = "GeneScoreMatrix", 
  name = "IL17B",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
# IPF vs Normal no treatment ==================
markers <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("IPF-NoTreatment"),
  bgdGroups = c("Normal-NoTreatment"),
  useMatrix = "PeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)

pv <- plotMarkers(seMarker = markers, 
                  name = "IPF-NoTreatment", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1",
                  plotAs = "Volcano")
plotPDF(pv, name = "IPF-NoTreatment-vs-Normal-NoTreatment-PeakMarkers-Volcano", 
        width = 5, height = 5, ArchRProj = projMulti2, addDOC = FALSE)
saveRDS(markers,
        file = paste0(file.path(data.dir,"Save-ProjMulti2/Peakmarkers_IPFvsNormal_NoTreatment.rds")))
list <- getMarkers(markers, cutOff = "FDR <= 1 & abs(Log2FC) >= 0.25")
markers.df <- as.data.frame(unlist(list))

markers.df <- markers.df[order(markers.df$Log2FC,decreasing = T),]

topPeaks <- markers.df[1:10,1:4]
set <- getPeakSet(projMulti2)

gr <- GRanges(
  seqnames = Rle(markers.df$seqnames),
  ranges = IRanges(start=markers.df$start, end = as.vector(markers.df$end)))
overlaps <- subsetByOverlaps(set, gr)

overlaps.df <- data.frame(chr=as.vector(overlaps@seqnames),
                          start=overlaps@ranges@start,
                          end=overlaps@ranges@start+500,
                          nearestGene=overlaps$nearestGene)

markers.df$peak <- paste0(markers.df$seqnames,":",
                          markers.df$start,"-",markers.df$end)
overlaps.df$peak <- paste0(overlaps.df$chr,":",
                           overlaps.df$start,"-",overlaps.df$end)
markers.df <- markers.df[order(markers.df$peak,decreasing = T),]

overlaps.df <- overlaps.df[order(overlaps.df$peak,decreasing = T),]
markers.df$nearestGene <- overlaps.df$nearestGene
markers.df$peakConfirm <- overlaps.df$peak

markers.df <- markers.df[order(markers.df$Log2FC,decreasing = T),]

write.csv(markers.df,
          file =file.path(data.dir,"Save-ProjMulti2/Plots/IPFvsNorma_NoTreatment_DA.csv"))
p <- plotBrowserTrack(projMulti2,groupBy = "Treatment",
                      useGroups = c("IPF-NoTreatment","Normal-NoTreatment"),
                      region = gr)

grid::grid.newpage()
grid::grid.draw(p)

## TFBS enrichment ==============================
projMulti2 <- addArchRAnnotations(ArchRProj = projMulti2,
                                  collection = "EncodeTFBS",
                                  db="/gstore/scratch/u/lucast3/fibroticmemory/ArchR-Hg38-v1.Anno")
enrichEncode <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projMulti2,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEncode <- plotEnrichHeatmap(enrichEncode,
                                   n = 14, transpose = TRUE)

ComplexHeatmap::draw(heatmapEncode,
                     heatmap_legend_side = "bot",
                     annotation_legend_side = "bot")

# Volcano Plot of DA ==================

#FDR <= 0.1 & abs(Log2FC) >= 1
markers.df <- markers.df %>% mutate(sig = ifelse((FDR<=0.1 | FDR <= 0.1) & abs(Log2FC) >= 1, "yes", "no"))


markers.df$diffexpressed <- "unchanged"

markers.df$diffexpressed[markers.df$Log2FC >= 1 &
                           markers.df$FDR<=0.1 &
                           markers.df$sig=="yes" ] <- "up"

markers.df$diffexpressed[markers.df$Log2FC <= -1 &
                           markers.df$FDR<=0.1 &
                           markers.df$sig=="yes" ] <- "down"

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("down", "up", "unchanged")

png(file.path(data.dir,"Save-ProjMulti2/Plots/DAR_IPFvsNormal_NoTreatment.png"),
    width=1000,height=1000)
ggplot(data=markers.df, aes(x=Log2FC, y=-log10(FDR), col=diffexpressed, label=nearestGene)) + 
  geom_point() + 
  theme_classic() +
  geom_text_repel(data= subset(markers.df,diffexpressed=="up"|diffexpressed=="down"),
                  show.legend = FALSE)+ 
  scale_colour_manual(values = mycolors) + ggtitle("DAR IPF-NoTreatment vs Normal-NoTreatment ")+
  theme(legend.position = "none",
        axis.title = element_text(size=18,face='bold'),
        plot.title = element_text(size=22,face='bold',hjust = 0.5))
dev.off()


# Heatmap ==================================
heat.up <- markers.df[order(markers.df$Log2FC,decreasing = T),]
heat.up <- heat.up[1:10,]
heat.down <- markers.df[order(markers.df$Log2FC,decreasing = F),]
heat.down <- heat.down[1:10,]

heatmap_data <- rbind(heat.up,heat.down)
h.data <- matrix(data.frame(gene=heatmap_data$nearestGene,FC=heatmap_data$Log2FC))
ComplexHeatmap::draw(heatmap_data,
                     heatmap_legend_side = "bot",
                     annotation_legend_side = "bot")

DE_markers <- getMarkerFeatures(
  ArchRProj = projMulti2, 
  groupBy = "Treatment",
  useGroups = c("IPF-NoTreatment"),
  bgdGroups = c("Normal-NoTreatment"),
  useMatrix = "GeneExpressionMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
