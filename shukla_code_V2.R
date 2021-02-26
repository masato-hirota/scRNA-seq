library(Seurat)
library(cowplot)
library(patchwork)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

data <- Read10X(data.dir = "/Volumes/bucket/IshikawaU/Masato/shukla_RNA/Integrated/outs/filtered_feature_bc_matrix")
data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)

Tag9 <- read.table("Tag9.csv", sep =',', header = TRUE) 
Tag10 <- read.table ("Tag10.csv", sep = ',', header= TRUE) 

Tag9$genotype <- rep("WT", 392)
Tag10$genotype <- rep("KO", 335)

Tag9 <- Tag9 [-2]
Tag10 <- Tag10 [-2]

Tags <- rbind(Tag9,Tag10)
data <- subset(data, cells = Tags$Barcode)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

data <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 7.5)

VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

data@meta.data$Barcode <- rownames(data@meta.data)

data@meta.data <- left_join(data@meta.data, Tags, by="Barcode")
rownames(data@meta.data) <- data@meta.data$Barcode
data@meta.data <- dplyr::select(data@meta.data,-Barcode)

#Dimentional reduction
# Run the standard workflow for visualization and clustering
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(data)
JunB.combined <- ScaleData(data, features = all.genes, verbose = FALSE)
JunB.combined <- RunPCA(JunB.combined, npcs = 30, verbose = FALSE)
ElbowPlot(JunB.combined,ndims = 30)
# t-SNE and Clustering
JunB.combined <- RunTSNE(JunB.combined, reduction = "pca", dims = 1:15, check_duplicates = FALSE)
JunB.combined <- FindNeighbors(JunB.combined, reduction = "pca", dims = 1:15)
JunB.combined <- FindClusters(JunB.combined,resolution = 0.7)

p <- DimPlot(JunB.combined, reduction = "tsne", split.by = "genotype", pt.size =0.3)
p$data$genotype <- factor(x = p$data$genotype, levels = c("WT", "KO")) # change the order of the factor levels
plot(p)
ggsave("split.by.genotype.png", plot=p, dpi = 500, width = 6, height = 5)

p1 <- DimPlot(JunB.combined, reduction = "tsne", group.by = "genotype", pt.size =0.3) + labs(title = "")
plot(p1)
ggsave("group.by.genotype.png", plot=p1, dpi = 500, width = 6, height = 5)

p2 <- DimPlot(JunB.combined, reduction = "tsne", pt.size =0.3,)
plot(p2)
ggsave("combined.png", plot=p2, dpi = 500, width = 6, height = 5)

saveRDS(JunB.combined, file = "JunB.combined.rds")
all.markers <- FindAllMarkers(object = JunB.combined)

FeaturePlot(JunB.combined, features = c("Sell"),  
            reduction = "tsne", cols = c("lightgrey", "darkred"), ncol = 3,min.cutoff = 0, 
            split.by ="genotype") & theme(plot.title = element_text(size = 10))

Cluster0 <- subset(JunB.combined, subset = seurat_clusters == 0)
Cluster1 <- subset(JunB.combined, subset = seurat_clusters == 1)
Cluster2 <- subset(JunB.combined, subset = seurat_clusters == 2)
Cluster3 <- subset(JunB.combined, subset = seurat_clusters == 3)
Cluster4 <- subset(JunB.combined, subset = seurat_clusters == 4)
Cluster5 <- subset(JunB.combined, subset = seurat_clusters == 5)


C0.WT <- nrow(Cluster0@meta.data[Cluster0@meta.data$genotype=='WT',])/nrow(Cluster0@meta.data)*100
C0.KO <- 100-C0.WT

C1.WT <- nrow(Cluster1@meta.data[Cluster1@meta.data$genotype=='WT',])/nrow(Cluster1@meta.data)*100
C1.KO <- 100-C1.WT

C2.WT <- nrow(Cluster2@meta.data[Cluster2@meta.data$genotype=='WT',])/nrow(Cluster2@meta.data)*100
C2.KO <- 100-C2.WT

C3.WT <- nrow(Cluster3@meta.data[Cluster3@meta.data$genotype=='WT',])/nrow(Cluster3@meta.data)*100
C3.KO <- 100-C3.WT

C4.WT <- nrow(Cluster4@meta.data[Cluster4@meta.data$genotype=='WT',])/nrow(Cluster4@meta.data)*100
C4.KO <- 100-C4.WT

C5.WT <- nrow(Cluster5@meta.data[Cluster5@meta.data$genotype=='WT',])/nrow(Cluster5@meta.data)*100
C5.KO <- 100-C5.WT

df1 <- cbind(C0.WT,C1.WT,C2.WT,C3.WT,C4.WT,C5.WT)
df2 <- cbind(C0.KO,C1.KO,C2.KO,C3.KO,C4.KO,C5.KO)

df <- rbind(df1,df2)
rownames(df) <- c("WT","KO")
colnames (df) <- c("0","1","2","3","4","5")
library(reshape2)
df <- melt(df)

g <- ggplot(df, aes(x = Var2  , y = value, fill = Var1)) + 
  geom_bar(stat = 'identity',width=0.7) + 
  ggtitle('') + 
  xlab('Clusters') + 
  ylab('Cell percentage in each cluster') + 
  theme_minimal()
g <- g + theme_bw()
g <- g + scale_y_continuous(expand = c(0,0), limits = c(0,100))
g <- g + theme(axis.text=element_text(colour="black"))
g <- g + theme(legend.title = element_blank())
g <- g + scale_x_continuous(breaks=seq(0,5,1))
g <- g + theme(axis.text=element_text(size=20, color = 'black', family = 'Helvetica'),
               axis.title=element_text(size=20,color = 'black', family = 'Helvetica'))
g <- g + theme(legend.text =element_text(size=15,color = 'black', family = 'Helvetica') )
g <- g + scale_fill_manual(values = c("dodgerblue4","firebrick"))
g

ggsave(file = "Percentage.png", plot = g, dpi = 500, width = 5, height = 5)

features <- c("Ccr7", "Sell", "Bcl2", "Stab1", "Tcf7", "Il7r","Cd27","Ctla2a","Cd28","Cxcr3", "Cxcr6", "Ly6a","Zeb2","Klrg1","Gzma","Gzmb","Klrk1","Stmn1","Cks1b","Birc5","Pcna","Ptma","Mcm5")
features <- as.character(features)

averages <- AverageExpression(JunB.combined, return.seurat=TRUE, add.ident = c("genotype","seurat_clusters"))

p3 <- DoHeatmap(averages, features = features, size = 5.5,hjust = 0.5, angle = 0) +
  scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"),na.value = "white")
ggsave(file = "heatmap.png", plot = p3, dpi = 300, width = 10, height = 5)

p4 <- DoHeatmap(JunB.combined, features = features, size = 5.5,hjust = 0.5, angle = 0) +
  scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"),na.value = "white") 
ggsave(file = "heatmap_each_cell.png", plot = p4, dpi = 300, width = 10, height = 5)

top20 <- FindAllMarkers(object = JunB.combined, slot = "data") %>% group_by(cluster) %>% top_n(20, avg_logFC)

p5 <- DoHeatmap(JunB.combined, features = top20$gene, size = 5 ,hjust = 0.5, angle = 0) +
  scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"),na.value = "white") + theme(text = element_text(size = 8), legend.key.size = unit(1, 'cm'), 
                                                                                             legend.text = element_text(size=10),
                                                                                             legend.title = element_text(size=14),
                                                                                             axis.text.y = element_text(face = 'bold'))
                                                                                             
ggsave(file = "heatmap_each_cell.top20.png", plot = p5, dpi = 300, width = 8, height = 10)

p6 <- DoHeatmap(averages, features = top20$gene, size = 5.5,hjust = 0.5, angle = 0) +
  scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"),na.value = "white") + theme(text = element_text(size = 8), legend.key.size = unit(1, 'cm'), 
                                                                                             legend.text = element_text(size=10),
                                                                                             legend.title = element_text(size=14),
                                                                                             axis.text.y = element_text(face = 'bold'))
ggsave(file = "heatmap.top20.png", plot = p6, dpi = 500, width = 8, height = 10)


###get_gene_count_data
Normalized.count <- GetAssayData(JunB.combined, slot = "data")
write.csv(Normalized.count, "LogNormalized.csv")
Zscore <- GetAssayData(JunB.combined, slot ="scale.data")
write.csv(Zscore, "scaled.csv")
Raw.count <- GetAssayData(JunB.combined@assays$RNA, slot ="counts")
write.csv(Raw.count, "Raw.csv")

###Vlnplot
VlnPlot(JunB.combined,features = "Junb", slot ="data", group.by = 'genotype') 
VlnPlot(JunB.combined,features = "Junb", slot ="counts", group.by = 'genotype') 

###Pseudo_bulk_analysis
rawdata <- read.csv("Raw.csv", header = TRUE, row.names = 1)
col <- colnames(rawdata)

library(stringr)
Tag9.Barcode <- str_replace(Tag9$Barcode, pattern="-", replacement=".")
Tag10.Barcode <- str_replace(Tag10$Barcode, pattern="-", replacement=".")

WT.count <- dplyr::select(rawdata, matches(Tag9.Barcode))
KO.count <- dplyr::select(rawdata, matches(Tag10.Barcode))

WT.bulk <- apply(WT.count, 1, sum)
KO.bulk <- apply(KO.count, 1, sum)

PB <- cbind(WT.bulk, KO.bulk)
colnames(PB) <- c("WT","KO")
#remove lowly expressed genes (you can modify this parameter)
write.csv(PB, "Pseudobulk.csv")

library(TCC)
group <- c("WT", "KO")
tcc <- new("TCC", PB, group)
tcc <- filterLowCountGenes(tcc, low.count = 10)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors

tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)

write.csv(result, "result.csv")

DEGs <- result[result$estimatedDEG == 1,]
write.csv(DEGs, "DEGs.csv", row.names = TRUE)

library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)

gene <- DEGs$gene_id

my.symbols <- as.character(gene)
all.symbols <- as.character(result$gene_id)

EnzID.my.symbols <- AnnotationDbi::select(org.Mm.eg.db, 
                      keys = my.symbols,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")

EnzID.all.symbols <- AnnotationDbi::select(org.Mm.eg.db, 
                                          keys = all.symbols,
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "SYMBOL")


gene <- EnzID.my.symbols$ENTREZID
#all.gene <- EnzID.all.symbols$ENTREZID

Gsea <- enrichGO(gene, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01, qvalueCutoff = 0.05, readable = TRUE)

head(as.data.frame(Gsea))

#Gsea.simple<-simplify(Gsea)
#head(as.data.frame(Gsea.simple))
barplot(Gsea, drop=TRUE, showCategory=15, title="Biological Process")
clusterProfiler::dotplot(Gsea)
clusterProfiler::emapplot(Gsea)

clusterProfiler::cnetplot(Gsea, categorySize="pvalue", foldChange=gene)

result$entrez <- EnzID.all.symbols$ENTREZID
result.fc <- dplyr::select(result, c(entrez, m.value))
result.fc <- result.fc[order(result.fc$m.value,decreasing = TRUE),]
result.fc <- na.omit(result.fc)
geneList = result.fc[,2]
names(geneList) = as.character(result.fc[,1])

gse_result<- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "BP",nPerm = 1000, minGSSize = 120, verbose = FALSE, pvalueCutoff = 0.05)
head(as.data.frame(gse_result))

gseaplot(gse_result, geneSetID = "GO:0051301", title = "cell division")


gse_result@geneSets


##########################
apop <- gse_result@geneSets$'GO:0051301' # you can serch GO term of interest in google and put the ID here *GO:0006915 is apoptotic process

apop.fc <- NULL
for (i in 1:length(apop)){
  if (apop[i] %in% names(geneList))
  apop.fc <- rbind(apop.fc,c(geneList[names(geneList) == apop[i]], apop[i]))
                   }

colnames(apop.fc) <- c("log2FC", "ENTREZID")
merged <- merge(apop.fc, EnzID.all.symbols, by = "ENTREZID")
colnames(merged)[3] <- "gene_id"
merged.pval <- merge(merged, result, by = "gene_id")
View(merged.pval)


##############################################
#Idents(JunB.combined) <- "seurat_clusters"
Idents(JunB.combined) <- "genotype"
fc <- FindMarkers(JunB.combined, ident.1 = "KO", ident.2 = "WT", verbose = FALSE, logfc.threshold = 0)
fc <- na.omit(fc)
scDEGs <- subset(fc, p_val_adj < 0.05 & fc$avg_logFC > 0.25 | fc$avg_logFC < -0.25)
#write.csv(scDEGs, "scDEGs.csv")

my.symbols <- as.character(rownames(scDEGs))
all.symbols <- as.character(rownames(fc))


EnzID.my.symbols <- AnnotationDbi::select(org.Mm.eg.db, 
                                          keys = my.symbols,
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "SYMBOL")

EnzID.all.symbols <- AnnotationDbi::select(org.Mm.eg.db, 
                                           keys = all.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")



gene <- EnzID.my.symbols$ENTREZID
all.gene <- EnzID.all.symbols$ENTREZID

Gsea <- enrichGO(gene, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01, qvalueCutoff = 0.05, readable = TRUE)

head(as.data.frame(Gsea))

#Gsea.simple<-simplify(Gsea)
#head(as.data.frame(Gsea.simple))

barplot(Gsea, drop=TRUE, showCategory=15, title="Biological Process")
clusterProfiler::dotplot(Gsea)
clusterProfiler::emapplot(Gsea)
clusterProfiler::cnetplot(Gsea, categorySize="pvalue", foldChange=gene)

fc$entrez <- EnzID.all.symbols$ENTREZID
result.fc <- dplyr::select(fc, c(entrez, avg_logFC))
result.fc <- result.fc[order(result.fc$avg_logFC, decreasing = TRUE),]
result.fc <- na.omit(result.fc)
geneList = result.fc[,2]
names(geneList) = as.character(result.fc[,1])

gse_result<- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "BP",nPerm = 1000, minGSSize = 120, verbose = FALSE, pvalueCutoff = 0.05)
head(as.data.frame(gse_result))

gseaplot(gse_result, geneSetID = "GO:0043066", title = "negative regulation of apoptotic signaling pathway")

apop <- gse_result@geneSets$'GO:0043066' # you can serch GO term of interest in google and put the ID here *GO:0006915 is apoptotic process

apop.fc <- NULL
for (i in 1:length(apop)){
  if (apop[i] %in% names(geneList))
    apop.fc <- rbind(apop.fc,c(geneList[names(geneList) == apop[i]], apop[i]))
}

colnames(apop.fc) <- c("log2FC", "ENTREZID")
merged <- merge(apop.fc, EnzID.all.symbols, by = "ENTREZID")
colnames(merged)[3] <- "gene_id"
result.fc$gene_id <- rownames(result.fc)
merged.pval <- merge(merged, result.fc, by = "gene_id")
View(merged.pval)

