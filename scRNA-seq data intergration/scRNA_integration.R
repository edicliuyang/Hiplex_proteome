#Note : This is not a full tutorial but an annotated command list
#For a more detailed tutorial and intro to Seurat, please see our tutorial of the 3k dataset

library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)

dir <- "E:/Desktop/sample0509/skin/"  
setwd(dir)
#change filename1 to name of txt file you want to load
data1 <- read.table("Filtered_matrix_correct.tsv", header = TRUE, sep = "\t", row.names = 1)
data2 <- t(data1)
sample1.name <- "DBiT_seq"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
DBiT_seq <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
DBiT_seq$tech <-"DBiT_seq"
DBiT_seq$celltype <-"DBiT_seq"
# 
# temp <- Read10X(data.dir = "E:/Desktop/sample0509/skin_ref/filtered_feature_bc_matrix")
# temp <- CreateSeuratObject(counts = temp, min.cells = 3, min.features = 200)
# temp[["percent.mt"]] <- PercentageFeatureSet(temp, pattern = "^MT-")
# #temp <- subset(temp, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# scRNA_seq <- temp
# scRNA_seq$tech <- "scRNA_seq"

scRNA_seq <- readRDS("E:/Desktop/sample0509/skin_ref/ref.RDS")
scRNA_seq$tech <- "scRNA_seq"

#VlnPlot(HNLN2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pancreas.list <- list(DBiT_seq,scRNA_seq)
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = FALSE)
}

pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 2000)
options(future.globals.maxSize= 3791289600)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = FALSE)

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)

immune.combined <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                 verbose = FALSE)

saveRDS(immune.combined, file = "E:/Desktop/sample0509/skin_ref/Skin_int.RDS")


immune.combined = readRDS("E:/Desktop/sample0509/skin_ref/Skin_int.RDS")

immune.combined <- RunPCA(immune.combined, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:10)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:10)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p3 <- DimPlot(immune.combined, reduction = "umap",group.by = "tech", order =c("scRNA-seq","DBiT_seq") , cols = c("lightblue","red"))
p3

p5 <- DimPlot(immune.combined, reduction = "umap",group.by = "celltype")
p5

p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,,pt.size = 0.01)
p3+p4+p5

pbmc.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)

immune.combined <-ScaleData(object=immune.combined, features = rownames(immune.combined))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(immune.combined, features = top10$gene) + scale_fill_gradientn(colors = c("red", "black","green")) 


features <- c("CD3", "CD4", "CD8", "CD19", "CD20", "LAG3", "CXCR6", "CXCR5", "ICOS","CD278","FOXP3", "IL2RA", "CD25", "CD69", "TNFRSF4","CD134","CD185","CD186" )

FeaturePlot(immune.combined, features = features)

id <- Idents(immune.combined)

write.table(Idents(immune.combined), "E:/Desktop/sample0509/skin_ref/idents.txt",sep="\t",quote = F)


df <- data.frame(V1 = names(id),V2=id )

test <- df %>% separate(V1, c("A", "B"),  sep = "x")
test$V2 = as.factor(test$V2)
dir <- "E:/Desktop/sample0509/skin_ref/"  
setwd(dir)

imported_raster=OpenImageR::readImage("E:/Desktop/sample0509/skin/skin.jpg")
g1 <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)


g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_main_color_integrate_30.pdf",sep =""), width=12.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=V2)) +
  #scale_color_manual(values=c("#F68282", "#31C53F", "#1FA195",'#B95FBB','#D4D915','#28CECA','#ff9a36','#2FF18B')) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape = 15, size = 3)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()



pancreas.list <- list(pbmc_rna, DBiT_seq )


for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = FALSE)
}

# reference.list <- list(pancreas.list[[1]], pancreas.list[[3]] )
# features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 2000)
# reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = features)
# 
# immune.anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", 
#                                          anchor.features = features)
# 
# pancreas.integrated <- IntegrateData(anchorset = immune.anchors,normalization.method = "SCT",  dims = 1:30)
# 
# DefaultAssay(pancreas.integrated) <- "integrated"
# # Run the standard workflow for visualization and clustering
# #pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
# pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
# pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
# p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
# p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
#   NoLegend()
# p1 + p2

pancreas.query <- pancreas.list[[1]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.list[[2]], query = pancreas.query, normalization.method ="SCT",
                                        dims = 1:10)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.list[[2]]$celltype, 
                            dims = 1:10)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

pancreas.query$predicted.id

count(pancreas.query$predicted.id)

saveRDS(pancreas.query, file = "E:/Desktop/lymphoma_ref/patient/FL1/predicted.RDS")

pancreas.query = readRDS("E:/Desktop/lymphoma_ref/lymph/FL4/predicted.RDS")

pancreas.query <- SCTransform(pancreas.query, vars.to.regress = "percent.mt", verbose = FALSE)
immune.combined <- RunPCA(pancreas.query, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:10)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:10)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)
DimPlot(immune.combined, group.by = "predicted.id", label = TRUE) + NoLegend()


dir <- "E:/Desktop/lymphoma_ref/patient/FL1/"  
setwd(dir)

pancreas.query$celltype = pancreas.query$predicted.id

names(pancreas.query$celltype)

count(pancreas.query$celltype)$x

df <- data.frame(V1 = names(pancreas.query$celltype),V2=pancreas.query$celltype )

test <- df %>% separate(V1, c("A", "B"),  sep = "x")
test$V2 = as.factor(test$V2)


imported_raster=OpenImageR::readImage("E:/Desktop/sample0509/skin/skin.jpg")
g1 <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)


g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_predict_10.pdf",sep =""), width=12.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=V2)) +
  scale_color_manual(values=c("#F68282", "#31C53F", "#1FA195",'#B95FBB','#D4D915','#28CECA','#ff9a36','#2FF18B','#aeadb3','#faf4cf','#CCB1F1','#25aff5','#A4DFF2','#4B4BF7','#AC8F14','#E6C122', "grey","blue")) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape = 15, size = 3)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()





