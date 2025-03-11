setwd("~/Desktop/CCM_iEC_bulk_RNAseq/DataDelivery_2024-07-24_14-35-35_snpseq00893/files/XD-3981/results/star_salmon")

library(DESeq2)
library(ggplot2)
library(tximport)
library(readr)
library(openxlsx)
library(ggrepel)
library(cowplot)

count.raw<-read.table("salmon.merged.gene_counts.tsv",sep = "\t",header = T,row.names = 1)

sample.name<-colnames(count.raw)[2:19]
sample.name<-gsub("\\.","-",sample.name)

group<-data.frame("Sample"=sample.name,
                  "group"=c(rep("Brain.slides.CCM01",3),rep("Brain.slides.HD01",3),
                                      rep("iECs.CCM01",3),rep("iECs.CCM02",3),rep("iECs.CCM03",3),
                                      rep("iECs.HD01",3)),
                                      row.names = sample.name)

#Load Data
files <- file.path(".", group$Sample, "quant.sf")
tx2gene <- read.table("/Users/fanya562/Desktop/CCM_iEC_bulk_RNAseq/DataDelivery_2024-07-24_14-35-35_snpseq00893/files/XD-3981/results/genome/salmon_tx2gene.tsv",sep = "\t",header = F)

txi <- tximport(files, type="salmon", tx2gene=tx2gene[,-2])

dds <- DESeqDataSetFromTximport(txi,
                                   colData = group,
                                   design = ~ group)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 20) >= smallestGroupSize
dds <- dds[keep,]

dds

library(RColorBrewer)
library(pheatmap)

#get vst normalized data
vsd <- vst(dds, blind=FALSE)

#Sample distance with only iECs samples
sampleDists <- dist(t(assay(vsd[,c(7:18)])))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


pdf("./Result/Sample_distance_heatmap.pdf",width = 8,height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

tiff("iECs_Sample_distance_heatmap.tiff",width = 2000,height = 1200,compression = "none",res = 400)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

plotPCA(vsd[,c(7:18)], intgroup=c("group"))

pcaData <- plotPCA(vsd[,c(7:18)], intgroup=c("group"), returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
tiff("iECs_Sample_PCA.tiff",width = 2000,height = 2000,compression = "none",res = 400)
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme_bw()+theme(panel.grid.major.y=element_blank(), 
                                 panel.grid.minor.y=element_blank(), 
                                 panel.grid.minor.x=element_blank(),
                                 panel.grid.major.x=element_blank())
dev.off()


pcaData1 <- plotPCA(vsd[,c(1:9,16:18)], intgroup=c("group"), returnData=T)
percentVar1 <- round(100 * attr(pcaData1, "percentVar"))
tiff("CCM01_HD01_Sample_PCA.tiff",width = 2000,height = 2000,compression = "none",res = 400)
ggplot(pcaData1, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) + 
  coord_fixed()+theme_bw()+theme(panel.grid.major.y=element_blank(), 
                                 panel.grid.minor.y=element_blank(), 
                                 panel.grid.minor.x=element_blank(),
                                 panel.grid.major.x=element_blank())
dev.off()

#DEGs analysis
dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
pdf("./Result/Sample_PCA.pdf",width = 6,height = 4)
plotPCA( DESeqTransform(se),intgroup=c("group"))+ 
  coord_fixed()+theme_bw()+theme(panel.grid.major.y=element_blank(), 
                                 panel.grid.minor.y=element_blank(), 
                                 panel.grid.minor.x=element_blank(),
                                 panel.grid.major.x=element_blank())
dev.off()

dds<-DESeq(dds)

#BS.CCM01.vs.BS.HD01 -------------------------------------------------
BS.CCM01.vs.BS.HD01.resLFC <- lfcShrink(dds, contrast = c("group","Brain.slides.CCM01","Brain.slides.HD01"), type="ashr")
BS.CCM01.vs.BS.HD01.resLFC<-data.frame(BS.CCM01.vs.BS.HD01.resLFC)
BS.CCM01.vs.BS.HD01.resLFC.f<-subset(BS.CCM01.vs.BS.HD01.resLFC,abs(BS.CCM01.vs.BS.HD01.resLFC$log2FoldChange)>1.5 & BS.CCM01.vs.BS.HD01.resLFC$padj<0.05)
write.xlsx(BS.CCM01.vs.BS.HD01.resLFC.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.CCM01.vs.BS.HD01/BS.CCM01.vs.BS.HD01.DEG.xlsx",rowNames=T,colNames=T)
library(clusterProfiler)
library(org.Hs.eg.db)
gene.list1<-rownames(subset(BS.CCM01.vs.BS.HD01.resLFC.f,BS.CCM01.vs.BS.HD01.resLFC.f$log2FoldChange>0))
eg1 = bitr(gene.list1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg1)
BS.CCM01.vs.BS.HD01.up.ego <- enrichGO(gene= eg1$ENTREZID,
                                      OrgDb         = org.Hs.eg.db,
                                      ont           = "BP",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)

BS.CCM01.vs.BS.HD01.up.ego.res<-BS.CCM01.vs.BS.HD01.up.ego@result
BS.CCM01.vs.BS.HD01.up.ego.res.f<-subset(BS.CCM01.vs.BS.HD01.up.ego.res,BS.CCM01.vs.BS.HD01.up.ego.res$p.adjust<0.05)
BS.CCM01.vs.BS.HD01.up.ego.res.f$log2foldchange<-apply(BS.CCM01.vs.BS.HD01.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})

write.xlsx(BS.CCM01.vs.BS.HD01.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.CCM01.vs.BS.HD01/BS.CCM01.vs.BS.HD01_up_GO.xlsx",rowNames=T,colNames=T)

#Go plot
BS.CCM01.vs.BS.HD01.up.ego.res.f<-BS.CCM01.vs.BS.HD01.up.ego.res.f[order(BS.CCM01.vs.BS.HD01.up.ego.res.f$Count,decreasing = T),]
GO.plot.data<-BS.CCM01.vs.BS.HD01.up.ego.res.f[c(2,12:14,20,23,26,28,32,38,43,61,92,144,272,302,338),]
GO.plot.data$GeneRatio<-as.numeric(GO.plot.data$Count)/537
GO.plot.data<-GO.plot.data[order(GO.plot.data$GeneRatio,decreasing = F),]
GO.plot.data$Description<-factor(GO.plot.data$Description,levels = GO.plot.data$Description)
color =colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50)

tiff("BS.CCM01.vs.BS.HD01.up.GO.dotplot.tiff",width = 3800,height = 2500,compression = "none",res = 400)
ggplot(GO.plot.data,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

pdf("BS.CCM01.vs.BS.HD01.up.GO.dotplot.pdf",width = 10,height = 6)
ggplot(GO.plot.data,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

gene.list2<-rownames(subset(BS.CCM01.vs.BS.HD01.resLFC.f,BS.CCM01.vs.BS.HD01.resLFC.f$log2FoldChange<0))
eg2 = bitr(gene.list2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg2)
BS.CCM01.vs.BS.HD01.down.ego <- enrichGO(gene= eg2$ENTREZID,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05,
                                    readable      = TRUE)

BS.CCM01.vs.BS.HD01.down.ego.res<-BS.CCM01.vs.BS.HD01.down.ego@result
BS.CCM01.vs.BS.HD01.down.ego.res.f<-subset(BS.CCM01.vs.BS.HD01.down.ego.res,BS.CCM01.vs.BS.HD01.down.ego.res$p.adjust<0.05)
BS.CCM01.vs.BS.HD01.down.ego.res.f$log2foldchange<-apply(BS.CCM01.vs.BS.HD01.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(BS.CCM01.vs.BS.HD01.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.CCM01.vs.BS.HD01/BS.CCM01.vs.BS.HD01_down_GO.xlsx",rowNames=T,colNames=T)

#vocanol plot
vocanol.color.scale<-scale_color_gradientn(
  colours=c("#0C2C84","green","yellow","#CE1256"),
  values=c(0,0.3,0.4,1),
  guide = guide_colorbar(title = "-log10 adj.pval")
)


BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj<-(-log10(BS.CCM01.vs.BS.HD01.resLFC.f$padj))
BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj<-BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj/10

max_size=0.5
min_size=0.1
BS.CCM01.vs.BS.HD01.resLFC.f<-BS.CCM01.vs.BS.HD01.resLFC.f[order(BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj,decreasing = T),]
BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj[1:3]<-28
BS.CCM01.vs.BS.HD01.resLFC.f$size=sqrt(min_size+(max_size-min_size)*(BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj-min(BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj))/(28-min(BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj)))
BS.CCM01.vs.BS.HD01.resLFC.f<-BS.CCM01.vs.BS.HD01.resLFC.f[order(BS.CCM01.vs.BS.HD01.resLFC.f$log10.p.adj,decreasing = T),]

BS.CCM01.vs.BS.HD01.resLFC.f$gene<-rownames(BS.CCM01.vs.BS.HD01.resLFC.f)
BS.CCM01.vs.BS.HD01.resLFC.f <- BS.CCM01.vs.BS.HD01.resLFC.f %>% 
  mutate(label=ifelse(gene %in% c("ADAMTS1","ADAMTS4","CAV1","ANGPT2","SEMA3C",
                                  "SEMA3A","WNT5A","WNT7B","MMP10","CXCL8",
                                  "CLDN5","SEMA6D","OCLN","PALMD","WNT5B",
                                  "SEMA5A","MMP9"),gene,""))


tiff("BS.CCM01.vs.BS.HD01.DEGs.VolcanoPlot.tiff",width = 3000,height = 2200,compression = "none",res = 400)
ggplot(BS.CCM01.vs.BS.HD01.resLFC.f,aes(x=log2FoldChange,y=log10.p.adj))+
  geom_point(aes(color=log10.p.adj))+
  geom_point(aes(color=log10.p.adj),stroke = 1.2,shape=21)+vocanol.color.scale+
  labs(y="-log10(adjust p value)")+
  theme_cowplot()+
  theme(panel.border = element_rect(colour = "black",fill = NA,linewidth = 1))+
  geom_vline(xintercept = c(-1.5,1.5),linetype="dotted",color="black")+
  geom_label_repel(data = BS.CCM01.vs.BS.HD01.resLFC.f,aes(label=label),
                   size=3,box.padding = 0.6,max.overlaps = 5000)
dev.off()




#BS.CCM01.vs.iEC.CCM01 -------------------------------------------------
BS.CCM01.vs.iEC.CCM01<-lfcShrink(dds, contrast = c("group","Brain.slides.CCM01","iECs.CCM01"), type="ashr")
BS.CCM01.vs.iEC.CCM01.DEG<-data.frame(BS.CCM01.vs.iEC.CCM01)
BS.CCM01.vs.iEC.CCM01.DEG<-BS.CCM01.vs.iEC.CCM01.DEG[order(BS.CCM01.vs.iEC.CCM01.DEG$log2FoldChange,decreasing = T),]
BS.CCM01.vs.iEC.CCM01.DEG.f<-subset(BS.CCM01.vs.iEC.CCM01.DEG,BS.CCM01.vs.iEC.CCM01.DEG$padj<0.05&abs(BS.CCM01.vs.iEC.CCM01.DEG$log2FoldChange)>1.5)
write.xlsx(BS.CCM01.vs.iEC.CCM01.DEG.f,file ="~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.CCM01.vs.iEC.CCM01/BS.CCM01.vs.iEC.CCM01.DEG.xlsx",rowNames=T,colNames=T )

gene.list3<-rownames(subset(BS.CCM01.vs.iEC.CCM01.DEG.f,BS.CCM01.vs.iEC.CCM01.DEG.f$log2FoldChange>0))
eg3 = bitr(gene.list3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg3)
BS.CCM01.vs.iEC.CCM01.up.ego <- enrichGO(gene= eg3$ENTREZID,
                                       OrgDb         = org.Hs.eg.db,
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)

BS.CCM01.vs.iEC.CCM01.up.ego.res<-BS.CCM01.vs.iEC.CCM01.up.ego@result
BS.CCM01.vs.iEC.CCM01.up.ego.res.f<-subset(BS.CCM01.vs.iEC.CCM01.up.ego.res,BS.CCM01.vs.iEC.CCM01.up.ego.res$p.adjust<0.05)
BS.CCM01.vs.iEC.CCM01.up.ego.res.f$log2foldchange<-apply(BS.CCM01.vs.iEC.CCM01.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(BS.CCM01.vs.iEC.CCM01.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.CCM01.vs.iEC.CCM01/BS.CCM01.vs.iEC.CCM01_up_GO.xlsx",rowNames=T,colNames=T)


#Go plot
BS.CCM01.vs.iEC.CCM01.up.ego.res.f<-BS.CCM01.vs.iEC.CCM01.up.ego.res.f[order(BS.CCM01.vs.iEC.CCM01.up.ego.res.f$Count,decreasing = T),]
GO.plot.data5<-BS.CCM01.vs.iEC.CCM01.up.ego.res.f[c(1:4,6,13,14,18,21,22,24:28,34),]
GO.plot.data5$GeneRatio<-as.numeric(GO.plot.data5$Count)/102
GO.plot.data5<-GO.plot.data5[order(GO.plot.data5$GeneRatio,decreasing = F),]
GO.plot.data5$Description<-factor(GO.plot.data5$Description,levels = GO.plot.data5$Description)

tiff("BS.CCM01.vs.iEC.CCM01.up.GO.dotplot.tiff",width = 3400,height = 2500,compression = "none",res = 400)
ggplot(GO.plot.data5,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

pdf("BS.CCM01.vs.iEC.CCM01.up.GO.dotplot.pdf",width = 10,height = 6)
ggplot(GO.plot.data5,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

gene.list4<-rownames(subset(BS.CCM01.vs.iEC.CCM01.DEG.f,BS.CCM01.vs.iEC.CCM01.DEG.f$log2FoldChange<0))
eg4 = bitr(gene.list4, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg4)
BS.CCM01.vs.iEC.CCM01.down.ego <- enrichGO(gene= eg4$ENTREZID,
                                         OrgDb         = org.Hs.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE)

BS.CCM01.vs.iEC.CCM01.down.ego.res<-BS.CCM01.vs.iEC.CCM01.down.ego@result
BS.CCM01.vs.iEC.CCM01.down.ego.res.f<-subset(BS.CCM01.vs.iEC.CCM01.down.ego.res,BS.CCM01.vs.iEC.CCM01.down.ego.res$p.adjust<0.05)
BS.CCM01.vs.iEC.CCM01.down.ego.res.f$log2foldchange<-apply(BS.CCM01.vs.iEC.CCM01.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(BS.CCM01.vs.iEC.CCM01.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.CCM01.vs.iEC.CCM01/BS.CCM01.vs.iEC.CCM01_down_GO.xlsx",rowNames=T,colNames=T)




#iEC.CCM01.vs.iEC.HD01 -------------------------------------------------
iEC.CCM01.vs.iEC.HD01<-lfcShrink(dds, contrast = c("group","iECs.CCM01","iECs.HD01"), type="ashr")
iEC.CCM01.vs.iEC.HD01.DEG<-data.frame(iEC.CCM01.vs.iEC.HD01)
iEC.CCM01.vs.iEC.HD01.DEG.f<-subset(iEC.CCM01.vs.iEC.HD01.DEG,iEC.CCM01.vs.iEC.HD01.DEG$padj<0.05&abs(iEC.CCM01.vs.iEC.HD01.DEG$log2FoldChange)>1.5)
write.xlsx(iEC.CCM01.vs.iEC.HD01.DEG.f,file ="~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.HD01/iEC.CCM01.vs.iEC.HD01.DEG.xlsx",rowNames=T,colNames=T )

gene.list5<-rownames(subset(iEC.CCM01.vs.iEC.HD01.DEG.f,iEC.CCM01.vs.iEC.HD01.DEG.f$log2FoldChange>0))
eg5 = bitr(gene.list5, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg5)
iEC.CCM01.vs.iEC.HD01.up.ego <- enrichGO(gene= eg5$ENTREZID,
                                         OrgDb         = org.Hs.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE)

iEC.CCM01.vs.iEC.HD01.up.ego.res<-iEC.CCM01.vs.iEC.HD01.up.ego@result
iEC.CCM01.vs.iEC.HD01.up.ego.res.f<-subset(iEC.CCM01.vs.iEC.HD01.up.ego.res,iEC.CCM01.vs.iEC.HD01.up.ego.res$p.adjust<0.05)
iEC.CCM01.vs.iEC.HD01.up.ego.res.f$log2foldchange<-apply(iEC.CCM01.vs.iEC.HD01.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM01.vs.iEC.HD01.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.HD01/iEC.CCM01.vs.iEC.HD01_up_GO.xlsx",rowNames=T,colNames=T)

gene.list6<-rownames(subset(iEC.CCM01.vs.iEC.HD01.DEG.f,iEC.CCM01.vs.iEC.HD01.DEG.f$log2FoldChange<0))
eg6 = bitr(gene.list6, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg6)
iEC.CCM01.vs.iEC.HD01.down.ego <- enrichGO(gene= eg6$ENTREZID,
                                           OrgDb         = org.Hs.eg.db,
                                           ont           = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)

iEC.CCM01.vs.iEC.HD01.down.ego.res<-iEC.CCM01.vs.iEC.HD01.down.ego@result
iEC.CCM01.vs.iEC.HD01.down.ego.res.f<-subset(iEC.CCM01.vs.iEC.HD01.down.ego.res,iEC.CCM01.vs.iEC.HD01.down.ego.res$p.adjust<0.05)
iEC.CCM01.vs.iEC.HD01.down.ego.res.f$log2foldchange<-apply(iEC.CCM01.vs.iEC.HD01.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})

write.xlsx(iEC.CCM01.vs.iEC.HD01.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.HD01/iEC.CCM01.vs.iEC.HD01_down_GO.xlsx",rowNames=T,colNames=T)

#Go plot
GO.plot.data1<-iEC.CCM01.vs.iEC.HD01.up.ego.res.f[c(11,55,71,74,83,85,93,94,110,123,136,137,141,159,171,188,226,240,388,467,548,559,561),]
GO.plot.data1$GeneRatio<-as.numeric(GO.plot.data1$Count)/450
GO.plot.data1<-GO.plot.data1[order(GO.plot.data1$GeneRatio,decreasing = F),]
GO.plot.data1$Description<-factor(GO.plot.data1$Description,levels = GO.plot.data1$Description)

tiff("iEC.CCM01.vs.iEC.HD01.up.GO.dotplot.tiff",width = 3400,height = 2500,compression = "none",res = 400)
ggplot(GO.plot.data1,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
 scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

pdf("iEC.CCM01.vs.iEC.HD01.up.GO.dotplot.pdf",width = 10,height = 6)
ggplot(GO.plot.data1,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

#iEC.CCM01.vs.iEC.CCM02 -------------------------------------------------
iEC.CCM01.vs.iEC.CCM02 <-lfcShrink(dds, contrast = c("group","iECs.CCM01","iECs.CCM02"), type="ashr")
iEC.CCM01.vs.iEC.CCM02.DEG<-data.frame(iEC.CCM01.vs.iEC.CCM02)
iEC.CCM01.vs.iEC.CCM02.DEG.f<-subset(iEC.CCM01.vs.iEC.CCM02.DEG,iEC.CCM01.vs.iEC.CCM02.DEG$padj<0.05&abs(iEC.CCM01.vs.iEC.CCM02.DEG$log2FoldChange)>1.5)
write.xlsx(iEC.CCM01.vs.iEC.CCM02.DEG.f,file ="~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.CCM02/iEC.CCM01.vs.iEC.CCM02.DEG.xlsx",rowNames=T,colNames=T )

gene.list7<-rownames(subset(iEC.CCM01.vs.iEC.CCM02.DEG.f,iEC.CCM01.vs.iEC.CCM02.DEG.f$log2FoldChange>0))
eg7 = bitr(gene.list7, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg7)
iEC.CCM01.vs.iEC.CCM02.up.ego <- enrichGO(gene= eg7$ENTREZID,
                                         OrgDb         = org.Hs.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE)

iEC.CCM01.vs.iEC.CCM02.up.ego.res<-iEC.CCM01.vs.iEC.CCM02.up.ego@result
iEC.CCM01.vs.iEC.CCM02.up.ego.res.f<-subset(iEC.CCM01.vs.iEC.CCM02.up.ego.res,iEC.CCM01.vs.iEC.CCM02.up.ego.res$p.adjust<0.05)
iEC.CCM01.vs.iEC.CCM02.up.ego.res.f$log2foldchange<-apply(iEC.CCM01.vs.iEC.CCM02.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM01.vs.iEC.CCM02.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.CCM02/iEC.CCM01.vs.iEC.CCM02_up_GO.xlsx",rowNames=T,colNames=T)

gene.list8<-rownames(subset(iEC.CCM01.vs.iEC.CCM02.DEG.f,iEC.CCM01.vs.iEC.CCM02.DEG.f$log2FoldChange<0))
eg8 = bitr(gene.list8, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg8)
iEC.CCM01.vs.iEC.CCM02.down.ego <- enrichGO(gene= eg8$ENTREZID,
                                           OrgDb         = org.Hs.eg.db,
                                           ont           = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)

iEC.CCM01.vs.iEC.CCM02.down.ego.res<-iEC.CCM01.vs.iEC.CCM02.down.ego@result
iEC.CCM01.vs.iEC.CCM02.down.ego.res.f<-subset(iEC.CCM01.vs.iEC.CCM02.down.ego.res,iEC.CCM01.vs.iEC.CCM02.down.ego.res$p.adjust<0.05)
iEC.CCM01.vs.iEC.CCM02.down.ego.res.f$log2foldchange<-apply(iEC.CCM01.vs.iEC.CCM02.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM01.vs.iEC.CCM02.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.CCM02/iEC.CCM01.vs.iEC.CCM02_down_GO.xlsx",rowNames=T,colNames=T)



#iEC.CCM01.vs.iEC.CCM03 -------------------------------------------------
iEC.CCM01.vs.iEC.CCM03 <-lfcShrink(dds, contrast = c("group","iECs.CCM01","iECs.CCM03"), type="ashr")
iEC.CCM01.vs.iEC.CCM03.DEG<-data.frame(iEC.CCM01.vs.iEC.CCM03)
iEC.CCM01.vs.iEC.CCM03.DEG.f<-subset(iEC.CCM01.vs.iEC.CCM03.DEG,iEC.CCM01.vs.iEC.CCM03.DEG$padj<0.05&abs(iEC.CCM01.vs.iEC.CCM03.DEG$log2FoldChange)>1.5)
write.xlsx(iEC.CCM01.vs.iEC.CCM03.DEG.f,file ="~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.CCM03/iEC.CCM01.vs.iEC.CCM03.DEG.xlsx",rowNames=T,colNames=T )

gene.list9<-rownames(subset(iEC.CCM01.vs.iEC.CCM03.DEG.f,iEC.CCM01.vs.iEC.CCM03.DEG.f$log2FoldChange>0))
eg9 = bitr(gene.list9, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg9)
iEC.CCM01.vs.iEC.CCM03.up.ego <- enrichGO(gene= eg9$ENTREZID,
                                          OrgDb         = org.Hs.eg.db,
                                          ont           = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)

iEC.CCM01.vs.iEC.CCM03.up.ego.res<-iEC.CCM01.vs.iEC.CCM03.up.ego@result
iEC.CCM01.vs.iEC.CCM03.up.ego.res.f<-subset(iEC.CCM01.vs.iEC.CCM03.up.ego.res,iEC.CCM01.vs.iEC.CCM03.up.ego.res$p.adjust<0.05)
iEC.CCM01.vs.iEC.CCM03.up.ego.res.f$log2foldchange<-apply(iEC.CCM01.vs.iEC.CCM03.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM01.vs.iEC.CCM03.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.CCM03/iEC.CCM01.vs.iEC.CCM03_up_GO.xlsx",rowNames=T,colNames=T)

gene.list10<-rownames(subset(iEC.CCM01.vs.iEC.CCM03.DEG.f,iEC.CCM01.vs.iEC.CCM03.DEG.f$log2FoldChange<0))
eg10 = bitr(gene.list10, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg10)
iEC.CCM01.vs.iEC.CCM03.down.ego <- enrichGO(gene= eg10$ENTREZID,
                                           OrgDb         = org.Hs.eg.db,
                                           ont           = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)

iEC.CCM01.vs.iEC.CCM03.down.ego.res<-iEC.CCM01.vs.iEC.CCM03.down.ego@result
iEC.CCM01.vs.iEC.CCM03.down.ego.res.f<-subset(iEC.CCM01.vs.iEC.CCM03.down.ego.res,iEC.CCM01.vs.iEC.CCM03.down.ego.res$p.adjust<0.05)
iEC.CCM01.vs.iEC.CCM03.down.ego.res.f$log2foldchange<-apply(iEC.CCM01.vs.iEC.CCM03.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM01.vs.iEC.CCM03.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM01.vs.iEC.CCM03/iEC.CCM01.vs.iEC.CCM03_down_GO.xlsx",rowNames=T,colNames=T)


#iEC.CCM02.vs.iEC.CCM03 -------------------------------------------------
iEC.CCM02.vs.iEC.CCM03 <-lfcShrink(dds, contrast = c("group","iECs.CCM02","iECs.CCM03"), type="ashr")
iEC.CCM02.vs.iEC.CCM03.DEG<-data.frame(iEC.CCM02.vs.iEC.CCM03)
iEC.CCM02.vs.iEC.CCM03.DEG.f<-subset(iEC.CCM02.vs.iEC.CCM03.DEG,iEC.CCM02.vs.iEC.CCM03.DEG$padj<0.05&abs(iEC.CCM02.vs.iEC.CCM03.DEG$log2FoldChange)>1.5)
write.xlsx(iEC.CCM02.vs.iEC.CCM03.DEG.f,file ="~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM02.vs.iEC.CCM03/iEC.CCM02.vs.iEC.CCM03.DEG.xlsx",rowNames=T,colNames=T )

gene.list11<-rownames(subset(iEC.CCM02.vs.iEC.CCM03.DEG.f,iEC.CCM02.vs.iEC.CCM03.DEG.f$log2FoldChange>0))
eg11 = bitr(gene.list11, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg11)
iEC.CCM02.vs.iEC.CCM03.up.ego <- enrichGO(gene= eg11$ENTREZID,
                                          OrgDb         = org.Hs.eg.db,
                                          ont           = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)

iEC.CCM02.vs.iEC.CCM03.up.ego.res<-iEC.CCM02.vs.iEC.CCM03.up.ego@result
iEC.CCM02.vs.iEC.CCM03.up.ego.res.f<-subset(iEC.CCM02.vs.iEC.CCM03.up.ego.res,iEC.CCM02.vs.iEC.CCM03.up.ego.res$p.adjust<0.05)
iEC.CCM02.vs.iEC.CCM03.up.ego.res.f$log2foldchange<-apply(iEC.CCM02.vs.iEC.CCM03.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM02.vs.iEC.CCM03.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM02.vs.iEC.CCM03/iEC.CCM02.vs.iEC.CCM03_up_GO.xlsx",rowNames=T,colNames=T)

gene.list12<-rownames(subset(iEC.CCM02.vs.iEC.CCM03.DEG.f,iEC.CCM02.vs.iEC.CCM03.DEG.f$log2FoldChange<0))
eg12 = bitr(gene.list12, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg12)
iEC.CCM02.vs.iEC.CCM03.down.ego <- enrichGO(gene= eg12$ENTREZID,
                                            OrgDb         = org.Hs.eg.db,
                                            ont           = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            qvalueCutoff  = 0.05,
                                            readable      = TRUE)

iEC.CCM02.vs.iEC.CCM03.down.ego.res<-iEC.CCM02.vs.iEC.CCM03.down.ego@result
iEC.CCM02.vs.iEC.CCM03.down.ego.res.f<-subset(iEC.CCM02.vs.iEC.CCM03.down.ego.res,iEC.CCM02.vs.iEC.CCM03.down.ego.res$p.adjust<0.05)
iEC.CCM02.vs.iEC.CCM03.down.ego.res.f$log2foldchange<-apply(iEC.CCM02.vs.iEC.CCM03.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})

write.xlsx(iEC.CCM02.vs.iEC.CCM03.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM02.vs.iEC.CCM03/iEC.CCM02.vs.iEC.CCM03_down_GO.xlsx",rowNames=T,colNames=T)

#iEC.CCM02.vs.iEC.HD01 -------------------------------------------------

iEC.CCM02.vs.iEC.HD01 <-lfcShrink(dds, contrast = c("group","iECs.CCM02","iECs.HD01"), type="ashr")
iEC.CCM02.vs.iEC.HD01.DEG<-data.frame(iEC.CCM02.vs.iEC.HD01)
iEC.CCM02.vs.iEC.HD01.DEG.f<-subset(iEC.CCM02.vs.iEC.HD01.DEG,iEC.CCM02.vs.iEC.HD01.DEG$padj<0.05&abs(iEC.CCM02.vs.iEC.HD01.DEG$log2FoldChange)>1.5)
write.xlsx(iEC.CCM02.vs.iEC.HD01.DEG.f,file ="~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM02.vs.iEC.HD01/iEC.CCM02.vs.iEC.HD01.DEG.xlsx",rowNames=T,colNames=T )

gene.list13<-rownames(subset(iEC.CCM02.vs.iEC.HD01.DEG.f,iEC.CCM02.vs.iEC.HD01.DEG.f$log2FoldChange>0))
eg13 = bitr(gene.list13, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg13)
iEC.CCM02.vs.iEC.HD01.up.ego <- enrichGO(gene= eg13$ENTREZID,
                                          OrgDb         = org.Hs.eg.db,
                                          ont           = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)

iEC.CCM02.vs.iEC.HD01.up.ego.res<-iEC.CCM02.vs.iEC.HD01.up.ego@result
iEC.CCM02.vs.iEC.HD01.up.ego.res.f<-subset(iEC.CCM02.vs.iEC.HD01.up.ego.res,iEC.CCM02.vs.iEC.HD01.up.ego.res$p.adjust<0.05)
iEC.CCM02.vs.iEC.HD01.up.ego.res.f$log2foldchange<-apply(iEC.CCM02.vs.iEC.HD01.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM02.vs.iEC.HD01.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM02.vs.iEC.HD01/iEC.CCM02.vs.iEC.HD01_up_GO.xlsx",rowNames=T,colNames=T)

gene.list14<-rownames(subset(iEC.CCM02.vs.iEC.HD01.DEG.f,iEC.CCM02.vs.iEC.HD01.DEG.f$log2FoldChange<0))
eg14 = bitr(gene.list14, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg14)
iEC.CCM02.vs.iEC.HD01.down.ego <- enrichGO(gene= eg14$ENTREZID,
                                            OrgDb         = org.Hs.eg.db,
                                            ont           = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            qvalueCutoff  = 0.05,
                                            readable      = TRUE)

iEC.CCM02.vs.iEC.HD01.down.ego.res<-iEC.CCM02.vs.iEC.HD01.down.ego@result
iEC.CCM02.vs.iEC.HD01.down.ego.res.f<-subset(iEC.CCM02.vs.iEC.HD01.down.ego.res,iEC.CCM02.vs.iEC.HD01.down.ego.res$p.adjust<0.05)
iEC.CCM02.vs.iEC.HD01.down.ego.res.f$log2foldchange<-apply(iEC.CCM02.vs.iEC.HD01.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM02.vs.iEC.HD01.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM02.vs.iEC.HD01/iEC.CCM02.vs.iEC.HD01_down_GO.xlsx",rowNames=T,colNames=T)


GO.plot.data2<-iEC.CCM02.vs.iEC.HD01.up.ego.res.f[c(1,2,4,11,13,15,16,19,22,30,31,33,39,42),]
GO.plot.data2$GeneRatio<-as.numeric(GO.plot.data2$Count)/288
GO.plot.data2<-GO.plot.data2[order(GO.plot.data2$GeneRatio,decreasing = F),]
GO.plot.data2$Description<-factor(GO.plot.data2$Description,levels = GO.plot.data2$Description)

tiff("iEC.CCM02.vs.iEC.HD01.up.GO.dotplot.tiff",width = 3800,height = 2000,compression = "none",res = 400)
ggplot(GO.plot.data2,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

pdf("iEC.CCM02.vs.iEC.HD01.up.GO.dotplot.pdf",width = 12,height = 6)
ggplot(GO.plot.data2,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

#iEC.CCM03.vs.iEC.HD01 -------------------------------------------------
iEC.CCM03.vs.iEC.HD01 <-lfcShrink(dds, contrast = c("group","iECs.CCM03","iECs.HD01"), type="ashr")
iEC.CCM03.vs.iEC.HD01.DEG<-data.frame(iEC.CCM03.vs.iEC.HD01)
iEC.CCM03.vs.iEC.HD01.DEG.f<-subset(iEC.CCM03.vs.iEC.HD01.DEG,iEC.CCM03.vs.iEC.HD01.DEG$padj<0.05&abs(iEC.CCM03.vs.iEC.HD01.DEG$log2FoldChange)>1.5)
write.xlsx(iEC.CCM03.vs.iEC.HD01.DEG.f,file ="~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM03.vs.iEC.HD01/iEC.CCM03.vs.iEC.HD01.DEG.xlsx",rowNames=T,colNames=T )

gene.list15<-rownames(subset(iEC.CCM03.vs.iEC.HD01.DEG.f,iEC.CCM03.vs.iEC.HD01.DEG.f$log2FoldChange>0))
eg15 = bitr(gene.list15, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg15)
iEC.CCM03.vs.iEC.HD01.up.ego <- enrichGO(gene= eg15$ENTREZID,
                                         OrgDb         = org.Hs.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE)

iEC.CCM03.vs.iEC.HD01.up.ego.res<-iEC.CCM03.vs.iEC.HD01.up.ego@result
iEC.CCM03.vs.iEC.HD01.up.ego.res.f<-subset(iEC.CCM03.vs.iEC.HD01.up.ego.res,iEC.CCM03.vs.iEC.HD01.up.ego.res$p.adjust<0.05)
iEC.CCM03.vs.iEC.HD01.up.ego.res.f$log2foldchange<-apply(iEC.CCM03.vs.iEC.HD01.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM03.vs.iEC.HD01.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM03.vs.iEC.HD01/iEC.CCM03.vs.iEC.HD01_up_GO.xlsx",rowNames=T,colNames=T)

gene.list16<-rownames(subset(iEC.CCM03.vs.iEC.HD01.DEG.f,iEC.CCM03.vs.iEC.HD01.DEG.f$log2FoldChange<0))
eg16 = bitr(gene.list16, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg16)
iEC.CCM03.vs.iEC.HD01.down.ego <- enrichGO(gene= eg16$ENTREZID,
                                           OrgDb         = org.Hs.eg.db,
                                           ont           = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)

iEC.CCM03.vs.iEC.HD01.down.ego.res<-iEC.CCM03.vs.iEC.HD01.down.ego@result
iEC.CCM03.vs.iEC.HD01.down.ego.res.f<-subset(iEC.CCM03.vs.iEC.HD01.down.ego.res,iEC.CCM03.vs.iEC.HD01.down.ego.res$p.adjust<0.05)
iEC.CCM03.vs.iEC.HD01.down.ego.res.f$log2foldchange<-apply(iEC.CCM03.vs.iEC.HD01.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(iEC.CCM03.vs.iEC.HD01.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/iEC.CCM03.vs.iEC.HD01/iEC.CCM03.vs.iEC.HD01_down_GO.xlsx",rowNames=T,colNames=T)


GO.plot.data3<-iEC.CCM03.vs.iEC.HD01.up.ego.res.f
GO.plot.data3$GeneRatio<-as.numeric(GO.plot.data3$Count)/150
GO.plot.data3<-GO.plot.data3[order(GO.plot.data3$GeneRatio,decreasing = F),]
GO.plot.data3$Description<-factor(GO.plot.data3$Description,levels = GO.plot.data3$Description)

tiff("iEC.CCM03.vs.iEC.HD01.up.GO.dotplot.tiff",width = 3800,height = 1200,compression = "none",res = 400)
ggplot(GO.plot.data3,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

pdf("iEC.CCM03.vs.iEC.HD01.up.GO.dotplot.pdf",width = 12,height = 6)
ggplot(GO.plot.data3,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

#BS.HD01.vs.iEC.HD01 -------------------------------------------------
BS.HD01.vs.iEC.HD01 <-lfcShrink(dds, contrast = c("group","Brain.slides.HD01","iECs.HD01"), type="ashr")
BS.HD01.vs.iEC.HD01.DEG<-data.frame(BS.HD01.vs.iEC.HD01)
BS.HD01.vs.iEC.HD01.DEG.f<-subset(BS.HD01.vs.iEC.HD01.DEG,BS.HD01.vs.iEC.HD01.DEG$padj<0.05&abs(BS.HD01.vs.iEC.HD01.DEG$log2FoldChange)>1.5)
write.xlsx(BS.HD01.vs.iEC.HD01.DEG.f,file ="~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.HD01.vs.iEC.HD01/BS.HD01.vs.iEC.HD01.DEG.xlsx",rowNames=T,colNames=T )

gene.list17<-rownames(subset(BS.HD01.vs.iEC.HD01.DEG.f,BS.HD01.vs.iEC.HD01.DEG.f$log2FoldChange>0))
eg17 = bitr(gene.list17, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg17)
BS.HD01.vs.iEC.HD01.up.ego <- enrichGO(gene= eg17$ENTREZID,
                                         OrgDb         = org.Hs.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE)

BS.HD01.vs.iEC.HD01.up.ego.res<-BS.HD01.vs.iEC.HD01.up.ego@result
BS.HD01.vs.iEC.HD01.up.ego.res.f<-subset(BS.HD01.vs.iEC.HD01.up.ego.res,BS.HD01.vs.iEC.HD01.up.ego.res$p.adjust<0.05)
BS.HD01.vs.iEC.HD01.up.ego.res.f$log2foldchange<-apply(BS.HD01.vs.iEC.HD01.up.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(BS.HD01.vs.iEC.HD01.up.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.HD01.vs.iEC.HD01/BS.HD01.vs.iEC.HD01_up_GO.xlsx",rowNames=T,colNames=T)

gene.list18<-rownames(subset(BS.HD01.vs.iEC.HD01.DEG.f,BS.HD01.vs.iEC.HD01.DEG.f$log2FoldChange<0))
eg18 = bitr(gene.list18, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg18)
BS.HD01.vs.iEC.HD01.down.ego <- enrichGO(gene= eg18$ENTREZID,
                                           OrgDb         = org.Hs.eg.db,
                                           ont           = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)

BS.HD01.vs.iEC.HD01.down.ego.res<-BS.HD01.vs.iEC.HD01.down.ego@result
BS.HD01.vs.iEC.HD01.down.ego.res.f<-subset(BS.HD01.vs.iEC.HD01.down.ego.res,BS.HD01.vs.iEC.HD01.down.ego.res$p.adjust<0.05)
BS.HD01.vs.iEC.HD01.down.ego.res.f$log2foldchange<-apply(BS.HD01.vs.iEC.HD01.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})
write.xlsx(BS.HD01.vs.iEC.HD01.down.ego.res.f,file = "~/Desktop/CCM_iEC_bulk_RNAseq/Result/BS.HD01.vs.iEC.HD01/BS.HD01.vs.iEC.HD01_down_GO.xlsx",rowNames=T,colNames=T)

mytheme=theme(axis.text.x=element_text(hjust=0.5,size=10),
                axis.text.y=element_text(size=10),
                axis.title.x=element_text(size=10),
                axis.title.y=element_text(size=10),
                axis.line=element_line(size=1),
                plot.margin=unit(c(1,1,1,1),"cm"),
                plot.title=element_text(hjust=0.5,size=12),
                legend.title=element_text(size=12),
                legend.text=element_text(size=12),
                legend.position="right",
                panel.background = element_rect(fill="white"), 
                legend.background=element_rect(fill='transparent'))

GO.plot.data4<-BS.HD01.vs.iEC.HD01.up.ego.res.f
GO.plot.data4$GeneRatio<-as.numeric(GO.plot.data4$Count)/87
GO.plot.data4<-GO.plot.data4[order(GO.plot.data4$GeneRatio,decreasing = F),]
GO.plot.data4$Description<-factor(GO.plot.data4$Description,levels = GO.plot.data4$Description)
tiff("BS.HD01.vs.iEC.HD01.up.GO.dotplot.tiff",width = 3200,height = 1400,compression = "none",res = 400)
ggplot(GO.plot.data4,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()

pdf("BS.HD01.vs.iEC.HD01.up.GO.dotplot.pdf",width = 8,height = 4)
ggplot(GO.plot.data4,aes(x=GeneRatio,y=Description,size=log2foldchange,color=-log10(qvalue)))+
  geom_point()+scale_size(range = c(3, 5),breaks = c(1,1.5,2,2.5))+mytheme+
  scale_color_gradientn(colours=color)+ylab("")+
  labs(color = "-Log10(qvalue)",size="Log2 Fold Enrichment")
dev.off()



library(ggvenn)
iECs.common.up.genes<-list("iEC1"=gene.list5,
                           "iEC2"=gene.list13,
                           "iEC3"=gene.list15)

Royal2<-c("#e4c9b2","#f1c2a5","#f49d98","#fcd68f","#629076")
Darjeeling1<-c("#fb0007","#f56f08","#4caecc")
fill_colors = Darjeeling1[1:length(iECs.common.up.genes)]
venn.diagram(iECs.common.up.genes,fill=fill_colors,filename = "venn.tiff",width = 1700,height = 1600)

Overlap <- calculate.overlap(x = iECs.common.up.genes)
my.lengths <- unlist(lapply(Overlap, function (x) { length(unlist(x))}))
my.values <- as.character(unlist(c(do.call(rbind, lapply(Overlap, as.data.frame)))))
my.matrix <- matrix(NA, nrow = max(my.lengths), ncol = length(my.lengths))
my.cumsum <- cumsum(my.lengths)
mm <- 1
for(i in 1:length(my.lengths)) {
  
  my.matrix[1:my.lengths[i],i] <- my.values[mm:my.cumsum[i]]
  
  mm <- my.cumsum[i]+1
  
}
my.gene.list <- as.data.frame(my.matrix)
write.xlsx(my.gene.list,file = "Overlapping_up_DEGs.xlsx",rowNames=F,colNames=T)

#Common up GO terms
gene.list19<-my.gene.list$V1
eg19 = bitr(gene.list19, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg19)
iEC.CCM_vs_HD01.common.up.ego <- enrichGO(gene= eg19$ENTREZID,
                                         OrgDb         = org.Hs.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE)

iEC.CCM_vs_HD01.common.up.ego.res<-iEC.CCM_vs_HD01.common.up.ego@result
iEC.CCM_vs_HD01.common.up.ego.res.f<-subset(iEC.CCM_vs_HD01.common.up.ego.res,iEC.CCM_vs_HD01.common.up.ego.res$p.adjust<0.05)
iEC.CCM_vs_HD01.common.up.ego.res.f$log2foldchange<-apply(BS.HD01.vs.iEC.HD01.down.ego.res.f,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(log2(GeneRatio/BgRatio),2)
  enrichment_fold
})

tiff("genes_up-regulated_in_all.tiff",width = 2000,height = 3000,compression = "none",res = 400)
pheatmap(counts.norm[Overlap[[1]],c(7:18)],cellwidth = 10,cellheight = 10,cluster_rows = T,cluster_cols = F,scale = "row")
dev.off()
tiff("genes_up-regulated_in_CCM01&02.tiff",width = 2400,height = 4200,compression = "none",res = 400)
pheatmap(counts.norm[Overlap[[2]],c(7:18)],cellwidth = 10,cellheight = 10,cluster_rows = T,cluster_cols = F,scale = "row")
dev.off()
tiff("genes_up-regulated_in_CCM01&03.tiff",width = 2000,height = 3600,compression = "none",res = 400)
pheatmap(counts.norm[Overlap[[3]],c(7:18)],cellwidth = 10,cellheight = 10,cluster_rows = T,cluster_cols = F,scale = "row")
dev.off()
tiff("genes_up-regulated_in_CCM02&03.tiff",width = 2000,height = 3600,compression = "none",res = 400)
pheatmap(counts.norm[Overlap[[4]],c(7:18)],cellwidth = 10,cellheight = 10,cluster_rows = T,cluster_cols = F,scale = "row")
dev.off()

####################################################################Visualization###########################################################################
library(viridis)
library(ggplot2)
mycol1<-c("#6BA5CE","#F5AA5F")
cmap<-c("viridis","magma","inferno","plasma","cividis","rocket","mako","turbo")
##############enrichGO plot
#common BS.vs.iECs up GO
BS.HD01.vs.iEC.HD01.up.ego.plot<-BS.HD01.vs.iEC.HD01.up.ego.res.f[c(2,5,7,10),]
BS.HD01.vs.iEC.HD01.up.ego.plot<-BS.HD01.vs.iEC.HD01.up.ego.plot[order(BS.HD01.vs.iEC.HD01.up.ego.plot$log2foldchange,decreasing = F),]
BS.HD01.vs.iEC.HD01.up.ego.plot$Description<-factor(BS.HD01.vs.iEC.HD01.up.ego.plot$Description,levels = BS.HD01.vs.iEC.HD01.up.ego.plot$Description)

BS.CCM01.vs.iEC.CCM01.up.ego.plot<-BS.CCM01.vs.iEC.CCM01.up.ego.res.f[BS.CCM01.vs.iEC.CCM01.up.ego.res.f$Description%in%BS.HD01.vs.iEC.HD01.up.ego.plot$Description,]

BS.iECs.common<-rbind(BS.HD01.vs.iEC.HD01.up.ego.plot,BS.CCM01.vs.iEC.CCM01.up.ego.plot)
BS.iECs.common$Sample<-factor(c(rep("HD01",4),rep("CCM01",4)),levels = c("HD01","CCM01"))

#p1<-ggplot(data = BS.HD01.vs.iEC.HD01.up.ego.plot,aes(x=log2foldchange,y=Description,fill=-log10(qvalue)))+
#  geom_bar(width = 0.5,stat = "identity",position=position_dodge())+
#  theme_classic()+
#  scale_x_continuous(expand = c(0,0))+
#  scale_fill_viridis(option = cmap[1],begin = 0.5,end = 1,direction = -1)
#p1<-p1+theme(axis.text.y = element_blank(),
#             axis.ticks.y = element_blank())+
#  geom_text(data =  BS.HD01.vs.iEC.HD01.up.ego.plot,aes(x=0.1,y=Description,label=Description),
#             size=4,hjust=0.01,vjust=-0.1)
#p1<-p1+geom_text(data =  BS.HD01.vs.iEC.HD01.up.ego.plot,aes(x=0.1,y=Description,label=geneID),
#             size=2.5,fontface="italic",hjust=0.01,vjust=1.5)



#p2<-ggplot(data = BS.iECs.common,aes(x=log2foldchange,y=Description,fill=Sample))+
#  geom_bar(width = 0.5,stat = "identity",position=position_dodge(width = 0.5))+
#  theme_classic()+
#  scale_x_continuous(expand = c(0,0))+
#  scale_fill_manual(values = alpha(mycol1,0.66))

#p2<-p2+theme(axis.text.y = element_blank(),
#             axis.ticks.y = element_blank())+
#  geom_text(data =   BS.iECs.common,aes(x=0.1,y=Description,label=Description),
#            size=4,hjust=0.01,vjust=-2.5)
# 
#p2<-p2+geom_text(data =  BS.iECs.common[1:4,],aes(x=0.1,y=Description,label=geneID,color=-log10(qvalue)),
#                 size=2.5,fontface="italic",hjust=0.01,vjust=-1.4)+
#  scale_color_viridis(option = cmap[1],begin = 0.5,end = 1,direction = -1)

#p2<-p2+geom_text(data =  BS.iECs.common[5:8,],aes(x=0.1,y=Description,label=geneID,color=-log10(qvalue)),
#                 size=2.5,fontface="italic",hjust=0.01,vjust=2)+
#  scale_color_viridis(option = cmap[1],begin = 0.5,end = 1,direction = -1)

#p2

ggplot(BS.HD01.vs.iEC.HD01.up.ego.plot,aes(y=log2foldchange,x=Description))+coord_flip()+
  geom_point(aes(size=Count,color=qvalue),shape=16)+
  scale_size_continuous(range = c(6,11))+
  scale_color_viridis(option = cmap[1],begin = 0.5,end = 1,direction = -1)+
  scale_y_continuous(limits = c(min(BS.HD01.vs.iEC.HD01.up.ego.plot$log2foldchange)-0.05,max(BS.HD01.vs.iEC.HD01.up.ego.plot$log2foldchange)+0.02))+
  labs(x=NULL,y="Log2Foldchange")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,linewidth = 1),
        axis.text = element_text(color = "black",size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))



library(pheatmap)

norm.count<-log2(counts(dds, normalized=TRUE) + 1)
iEC.CCM01.vs.iEC.HD01.DEG.f<-iEC.CCM01.vs.iEC.HD01.DEG.f[order(iEC.CCM01.vs.iEC.HD01.DEG.f$log2FoldChange,decreasing = T),]
gene<-rownames(iEC.CCM01.vs.iEC.HD01.DEG.f)[1:20]
pheatmap(norm.count[gene,],cluster_rows = F,cluster_cols = F,scale = "row")

iEC.CCM02.vs.iEC.HD01.DEG.f<-iEC.CCM02.vs.iEC.HD01.DEG.f[order(iEC.CCM02.vs.iEC.HD01.DEG.f$log2FoldChange,decreasing = T),]
gene<-rownames(iEC.CCM02.vs.iEC.HD01.DEG.f)[1:20]
pheatmap(norm.count[gene,],cluster_rows = F,cluster_cols = F,scale = "row")

iEC.CCM03.vs.iEC.HD01.DEG.f<-iEC.CCM03.vs.iEC.HD01.DEG.f[order(iEC.CCM03.vs.iEC.HD01.DEG.f$log2FoldChange,decreasing = T),]
gene<-rownames(iEC.CCM03.vs.iEC.HD01.DEG.f)[1:20]
pheatmap(norm.count[gene,],cluster_rows = F,cluster_cols = F,scale = "row")

iEC.CCM02.vs.iEC.CCM03.DEG.f<-iEC.CCM02.vs.iEC.CCM03.DEG.f[order(iEC.CCM02.vs.iEC.CCM03.DEG.f$log2FoldChange,decreasing = T),]
gene<-rownames(iEC.CCM02.vs.iEC.CCM03.DEG.f)[1:20]
pheatmap(norm.count[gene,],cluster_rows = F,cluster_cols = F,scale = "row")

BS.CCM01.vs.iEC.CCM01.DEG.f<-BS.CCM01.vs.iEC.CCM01.DEG.f[order(BS.CCM01.vs.iEC.CCM01.DEG.f$log2FoldChange,decreasing = T),]
gene<-rownames(BS.CCM01.vs.iEC.CCM01.DEG.f)[1:20]
pheatmap(norm.count[gene,],cluster_rows = F,cluster_cols = F,scale = "row")



#Selected BBB related gene expression
g.select=c("ABCB1","ABCA2","ABCG2","ABCC1","ABCC2","ABCC3","ABCC4","ABCC5","SLC2A1","SLC5A1","SLC2A13","SLC5A3","SLC7A1","SLC7A3",
           "SLC7A5","SLC7A6","SLC38A1","SLC38A2","SLC38A3","SLC38A5","SLC1A1","SLC1A2","SLC1A3","SLC1A4","SLC1A15","SLC6A6","SLC6A9",
           "SLC16A1","SLC16A2","SLC16A7","SLC27A1","SLC27A4","SLC28A2","SLC29A1","SLC29A2","SLC22A8","SLCO1A4","SLCO2B1","SLCO1C1",
           "SLC22A1","SLC22A2","SLC22A3","SLC22A5","SLC29A4","SLC44A1","SLC5A6","AVPR1A","TFRC","TFR2","LEPR","INSR","LRP1","LRP2",
           "AGER","MFSD2A","OCLN","CLDN1","CLDN3","CLDN5","CLDN12","TJP1","TJP2","TJP3","CDH5","PECAM1","F11R","JAM2","JAM3","ESAM",
           "GJB6","GJA1","DMD","PECAM1","PDGFRB")

counts.norm<-counts(dds, normalized=T)
g.select.f<-intersect(g.select,rownames(counts.norm))
tiff("BBB_related_genes_heatmap_without_CCM02_03.tiff",width = 2400,height = 4600,compression = "none",res = 400)
pheatmap(counts.norm[g.select.f,-c(10:15)],cellwidth = 10,cellheight = 10,cluster_rows = T,cluster_cols = T,scale = "row")
dev.off()

barplot(counts.norm["PDGFRB",])


#Selected CCM gene expression

g.select.1<-c("MAP3K3", "KLF2", "KLF4", "CLDN5", "CDH5", "RHOA", "RHOC", "CAV1", "VWF", "TM", "BMI1", "ADAMTS4", "KRIT1", "PDCD10", "CCM2")
counts.norm<-counts(dds, normalized=T)
g.select1.f<-intersect(g.select.1,rownames(counts.norm))

tiff("Selected_genes_expression_heatmap_in_iECs.tiff",width = 2000,height = 3000,compression = "none",res = 400)
pheatmap(counts.norm[g.select1.f,c(7:18)],cellwidth = 10,cellheight = 10,cluster_rows = T,cluster_cols = F,scale = "row")
dev.off()

tiff("Selected_genes_expression_heatmap_in_BS.tiff",width = 2000,height = 3000,compression = "none",res = 400)
pheatmap(counts.norm[g.select1.f,c(1:6)],cellwidth = 10,cellheight = 10,cluster_rows = T,cluster_cols = F,scale = "row")
dev.off()

#Zonation genes heatmap
zonation.genes<-read.xlsx("Hs_brain_EC_zonation_genes.xlsx",rowNames = F,colNames = T,sheet = 2)
artery.genes<-intersect(rownames(counts.norm),zonation.genes$Artery)
arterial.genes<-intersect(rownames(counts.norm),zonation.genes$Arteriole)
capillary.genes<-intersect(rownames(counts.norm),zonation.genes$Capillary)
venous.genes<-intersect(rownames(counts.norm),zonation.genes$Venule)
vein.genes<-intersect(rownames(counts.norm),zonation.genes$Vein)
                        

zonation.genes.f<-data.frame("gene"=c("SOX17","DLL4","HEY2",unique(c(artery.genes,arterial.genes)),capillary.genes,unique(c(venous.genes,vein.genes))),
                             "cluster"=c(rep("Arterial",29),rep("Capillary",19),rep("Venous",26)))
zonation.genes.f<-zonation.genes.f[-c(6,42,63,68),]

annotation_row<-data.frame("Cluster"=zonation.genes.f$cluster,
                           row.names = zonation.genes.f$gene)

color =colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50)
legend_breaks <- seq(-1, 1, length.out = 50) 

pheatmap(counts.norm[zonation.genes.f$gene,],cluster_rows = F,cluster_cols = T,scale = "column",show_rownames = T,
         annotation_row = annotation_row,color = color,breaks = legend_breaks)


pheatmap(counts.norm[zonation.genes.f$gene,],cluster_rows = F,cluster_cols = T,scale = "row",show_rownames = T,
         annotation_row = annotation_row,color = color,breaks = legend_breaks)

tiff("EC_zonation_gene_heatmap.tiff",width = 2000,height = 6000,compression = "none",res = 400)
pheatmap(counts.logcpm[zonation.genes.f$gene,],cluster_rows = F,cluster_cols = T,scale = "column",show_rownames = T,
         annotation_row = annotation_row,color = color)
dev.off()


tiff("BS_EC_zonation_gene_heatmap.tiff",width = 1400,height = 6000,compression = "none",res = 400)
pheatmap(counts.logcpm[zonation.genes.f$gene,1:6],cluster_rows = F,cluster_cols = T,scale = "column",show_rownames = T,
         annotation_row = annotation_row,color = color)
dev.off()

tiff("iECs_EC_zonation_gene_heatmap.tiff",width = 1900,height = 6000,compression = "none",res = 400)
pheatmap(counts.logcpm[zonation.genes.f$gene,7:18],cluster_rows = F,cluster_cols = T,scale = "column",show_rownames = T,
         annotation_row = annotation_row,color = color)
dev.off()

arterial.genes.counts<-counts.norm[c("SOX17","DLL4","HEY2",unique(c(artery.genes,arterial.genes))),]
capillary.genes.counts<-counts.norm[capillary.genes,]
Vein.gene.counts<-counts.norm[unique(c(venous.genes,vein.genes)),]

arterial.avg.exp <- colMeans(arterial.genes.counts)
capillary.avg.exp <- colMeans(capillary.genes.counts)
vein.avg.exp <- colMeans(Vein.gene.counts)

avg.exp <- rbind(arterial.avg.exp,capillary.avg.exp)
avg.exp <- rbind(avg.exp, vein.avg.exp)

pheatmap(avg.exp,scale = "column",cluster_cols = F,cluster_rows = F)


#ssGSEA
zonation.list<-list("Artery"=c("SOX17","DLL4","HEY2",unique(c(artery.genes,arterial.genes))),
                    "Capillary"=capillary.genes,
                    "Vein"=unique(c(venous.genes,vein.genes)))

counts<-counts(dds, normalized=F)
counts.logcpm<-log2(cpm(counts)+1)

gsva_mat <- gsva(expr=counts.logcpm, 
                 gset.idx.list=zonation.list,
                 method="ssgsea",
                 ssgsea.norm = TRUE,
                 kcdf="Gaussian", #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())
write.xlsx(gsva_mat,"EC_zonation_ssGSEA_score.xlsx",rowNames=T,colNames=T)

tiff("BS_EC_zonation_ssGSEA_heatmap.tiff",width = 1200,height = 1200,compression = "none",res = 400)
pheatmap(gsva_mat[,1:6],cluster_cols = F,scale = "none",cluster_rows = F,breaks = seq(0.5, 1.3, length.out = 100),color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
dev.off()

tiff("iECs_EC_zonation_ssGSEA_heatmap.tiff",width = 1800,height = 1200,compression = "none",res = 400)
pheatmap(gsva_mat[,7:18],cluster_cols = F,scale = "none",cluster_rows = F)
dev.off()


> sessionInfo()
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 15.3.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base  

