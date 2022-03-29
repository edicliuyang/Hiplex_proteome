library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(SingleR)
library("monocle")
library(cowplot)
require(scales)

dir <- "C:/Users/EDIC-THINKPAD/Desktop/FFPEs/FFPE-1/"  
setwd(dir)

#change filename1 to name of txt file you want to load
data1 <- read.table("Filtered_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)
data2 <- t(data1)
sample1.name <- "FFPE-1"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
FFPE1 <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
FFPE1$tech <-"FFPE-1"

dir <- "C:/Users/EDIC-THINKPAD/Desktop/FFPEs/FFPE-2/"  
setwd(dir)

#change filename1 to name of txt file you want to load
data1 <- read.table("Filtered_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)
data2 <- t(data1)
sample1.name <- "FFPE-2"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
FFPE2 <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
FFPE2$tech <-"FFPE-2"




#read ref data from RDS file

dir <- "C:/Users/EDIC-THINKPAD/Desktop/nature2019"  
setwd("C:/Users/EDIC-THINKPAD/Desktop/nature2019")
fData(data)
pData(data)
data = readRDS("cds_cleaned_sampled_100k.RDS")
df = data.frame(gene = fData(data)$gene_short_name)
idname = names(table(pData(data)$id))
sc = data.frame(pData(data)$id, pData(data)$day)
sc <- unique(sc)
sc

s1 = row.names(subset(pData(data), id == 4))
s1 = data[,s1]
x <- exprs(s1)
row.names(x) <- fData(s1)$gene_short_name

s2 = row.names(subset(pData(data), id == 5))
s2 = data[,s2]
y <- exprs(s2)
row.names(y) <- fData(s2)$gene_short_name
z = cbind2(x,y)
row.names(z) <- fData(s2)$gene_short_name


ref <- SummarizedExperiment(
  list(counts=z),
  colData=DataFrame(label.main=c(s1@phenoData@data$Main_trajectory,s2@phenoData@data$Main_trajectory),
                    label.sub=c(s1@phenoData@data$Sub_trajectory_name,s2@phenoData@data$Sub_trajectory_name) )
)
ref <- scater::logNormCounts(ref)

idname = c("6",
           "44",
           "24",
           "5",
           "60",
           "59",
           "4",
           "43",
           "26",
           "25"
)
temp = c()
name1 = c()
name2 = c()
for (i in idname){
  s1 = row.names(subset(pData(data), id == i))
  s1 = data[,s1]
  x <- exprs(s1)
  temp = cbind2(temp,x)
  name1 = c(name1, s1@phenoData@data$Main_trajectory)
  name2 = c(name2, s1@phenoData@data$Sub_trajectory_name)
}

row.names(temp) <- fData(s2)$gene_short_name

#change filename1 to name of txt file you want to load
sample1.name <- "E10.5"
matrix1.data <- Matrix(temp, sparse = TRUE)
E10 <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
E10$tech <-"E10.5"

idname = c("10",
           "50",
           "49",
           "8",
           "7",
           "29"
          
           
)
temp = c()
name1 = c()
name3 = c()
for (i in idname){
  s1 = row.names(subset(pData(data), id == i))
  s1 = data[,s1]
  x <- exprs(s1)
  temp = cbind2(temp,x)
  name1 = c(name1, s1@phenoData@data$Main_trajectory)
  name3 = c(name3, s1@phenoData@data$Sub_trajectory_name)
}

row.names(temp) <- fData(s2)$gene_short_name
names = c(name3,name2)

na1 = rep("E11 brain", 1840)
na2 = rep("E11 tail", 2109)
names = c(na1,na2,names)

#change filename1 to name of txt file you want to load
sample1.name <- "E9.5"
matrix1.data <- Matrix(temp, sparse = TRUE)
E9 <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
E9$tech <-"E11.5"




pancreas.list <- list(FFPE1,FFPE2, E10)
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = FALSE)
}

pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
options(future.globals.maxSize= 3791289600)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = FALSE)

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)


immune.combined <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                 verbose = FALSE)

immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)
immune.combined$annotation = names
p3 <- DimPlot(immune.combined, reduction = "umap",group.by = "tech", order = c("E11 brain 10um","E11 tail 10um","E10.5","E11.5") ,pt.size = 0.01,cols = c("lightblue", "pink", "red","blue"))
p3
p4
p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,,pt.size = 0.01)
plot_grid(p3,p4)

p3 <- DimPlot(immune.combined, reduction = "umap",group.by = "annotation", order = c("E10 Eye"))
p3

sort(table(names))
names_other = replace(names, names %in% c("Len epithelial trajectory"，
                                     "Placenta endodermal trajectory"，
                                     "Satellite glia trajectory"，
                                     "Primordial germ cell trajectory"，
                                     "Lymphatic endothelial trajectory"，
                                     "Liver endothelial trajectory"，
                                     "Retina epithelial trajectory"，
                                     "Stomach epithelial trajectory"，
                                     "Lung epithelial trajectory"，
                                     "Enteric neuron trajectory 2"，
                                     "Olfactory ensheathing cell trajectory"，
                                     "Megakaryocyte trajectory"，
                                     "Melanocyte trajectory"，
                                     "Olfactory sensory neuron trajectory"，
                                     "Enteric neuron trajectory 1"，
                                     "Retinal fibroblast trajectory"，
                                     "Schwann cell and Enteric glia trajectory 2"，
                                     "Urothelium trajectory"，
                                     "Apical ectodermal ridge trajectory"，
                                     "Cardiac muscle trajectory"，
                                     "Endocardium trajectory"，
                                     "Renal epithelial trajectory"，
                                     "Shisa6 positive neuron trajectory"，
                                     "Midgut/Hindgut epithelial trajectory"，
                                     "Pericardium trajectory"，
                                     "Arterial endothelial trajectory"，
                                     "Retina trajectory"，
                                     "Schwann cell and Enteric glia trajectory 1"，
                                     "Brain endothelial trajectory"，
                                     "Auditory epithelial trajectory"，
                                     "Hepatocyte trajectory"，
                                     "White blood cell trajectory"，
                                     "Ependymal cell trajectory"，
                                     "Olfactory epithelial trajectory"，
                                     "Definitive erythroid trajectory"
), "other")
sort(table(names_other))


immune.combined$annotation = names_other

my_color_palette <- c(hue_pal()(20),"black","blue")

p3 <- DimPlot(immune.combined, reduction = "umap",group.by = "annotation", order = c("E11 brain","E11 tail"), label = TRUE) +
  scale_color_manual(values = my_color_palette)
p3

immune.combined$annotation = names


my_color_palette <- c(hue_pal()(53),"black")

p3 <- DimPlot(immune.combined, reduction = "umap",group.by = "annotation", order = c("E10 Eye"), label = TRUE) + NoLegend() +
  scale_color_manual(values = my_color_palette)
p3



write.table(sort(table(names)), "annotation.txt",sep="\t",quote = F)

p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p4
plot_grid(p3,p4)



Idents(immune.combined)

write.table(Idents(immune.combined), "ident.txt",sep="\t",quote = F)



pbmc.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + scale_fill_gradientn(colors = c("red", "black","green")) 


go <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.table(go, "go_E11.txt",sep="\t")






#build ref data

dir <- "C:/Users/EDIC-THINKPAD/Desktop/nature2019"  
setwd("C:/Users/EDIC-THINKPAD/Desktop/nature2019")

data = readRDS("cds_cleaned_sampled_100k.RDS")
df = data.frame(gene = fData(data)$gene_short_name)
idname = names(table(pData(data)$id))
sc = data.frame(pData(data)$id, pData(data)$day)
sc <- unique(sc)
sc

s1 = row.names(subset(pData(data), id == 4))
s1 = data[,s1]
x <- exprs(s1)
row.names(x) <- fData(s1)$gene_short_name

s2 = row.names(subset(pData(data), id == 5))
s2 = data[,s2]
y <- exprs(s2)
row.names(y) <- fData(s2)$gene_short_name
z = cbind2(x,y)
row.names(z) <- fData(s2)$gene_short_name


ref <- SummarizedExperiment(
  list(counts=z),
  colData=DataFrame(label.main=c(s1@phenoData@data$Main_trajectory,s2@phenoData@data$Main_trajectory),
                    label.sub=c(s1@phenoData@data$Sub_trajectory_name,s2@phenoData@data$Sub_trajectory_name) )
)
ref <- scater::logNormCounts(ref)

idname = c("6",
           "44",
           "24",
           "5",
           "60",
           "59",
           "4",
           "43",
           "26",
           "25"
)
temp = c()
name1 = c()
name2 = c()
for (i in idname){
  s1 = row.names(subset(pData(data), id == i))
  s1 = data[,s1]
  x <- exprs(s1)
  temp = cbind2(temp,x)
  name1 = c(name1, s1@phenoData@data$Main_trajectory)
  name2 = c(name2, s1@phenoData@data$Sub_trajectory_name)
}

row.names(temp) <- fData(s2)$gene_short_name


ref <- SummarizedExperiment(
  list(counts=temp),
  colData=DataFrame(label.main=name1,
                    label.sub=name2 )
)

#run from here!!!!

hESCs = GetAssayData(pbmc, slot = "counts")
common <- intersect(rownames(hESCs), rownames(ref))
ref <- ref[common,]
hESCs <- hESCs[common,]


ref <- scater::logNormCounts(ref)


pred.hpca <- SingleR(test = hESCs, ref = ref, 
                     labels = ref$label.main)

pred.hpca1 <- SingleR(test = hESCs, ref = ref, 
                      labels = ref$label.sub)


model1 <- pred.hpca 
model2 <- pred.hpca1

pred.hpca <- model1
pred.hpca <- model2

dir <- "C:/Users/EDIC-THINKPAD/Desktop/pseudo_bulk/0713cL"  
setwd("C:/Users/EDIC-THINKPAD/Desktop/pseudo_bulk/0713cL")
df = data.frame(name = pred.hpca@rownames, type = pred.hpca@listData$first.labels)

pred.hpca$labels
pbmc[["SingleR.cluster.labels"]] <- 
  pred.hpca$labels
DimPlot(pbmc, label = TRUE,group.by = "SingleR.cluster.labels") + NoLegend()

sort(table(pbmc$SingleR.cluster.labels))
dimnames(sort(table(pbmc$SingleR.cluster.labels)))

df <- data.frame(pbmc$SingleR.cluster.labels)
df1 <-data.frame(X =row.names(df), count= df$pbmc.SingleR.cluster.labels)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

#imported_raster=OpenImageR::readImage("FFPE-2.jpg")
#g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_main.pdf",sep =""), width=12.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape = 15, size = 3)+
  expand_limits(x = 0, y = 0)  +
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



pred.hpca <- model2


df = data.frame(name = pred.hpca@rownames, type = pred.hpca@listData$first.labels)

pred.hpca$labels
pbmc[["SingleR.cluster.labels"]] <- 
  pred.hpca$labels
DimPlot(pbmc, label = TRUE,group.by = "SingleR.cluster.labels")

pbmc$SingleR.cluster.labels

df <- data.frame(pbmc$SingleR.cluster.labels)
df1 <-data.frame(X =row.names(df), count= df$pbmc.SingleR.cluster.labels)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

imported_raster=OpenImageR::readImage("FFPE-2.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_main.pdf",sep =""), width=12.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
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








setwd("C:/Users/EDIC-THINKPAD/Desktop/pseudo_bulk/0713cL")
#plot sub catagory of nature paper
pred.hpca <- model2

xx <- pred.hpca$labels
lesscells <- c("Liver endothelial trajectory",
               "Retinal fibroblast trajectory",
               "Branchial arch epithelial trajectory",
               "Neural epithelial trajectory",
               "Renal epithelial trajectory",
               "Apical ectodermal ridge trajectory",
               "Epidermis trajectory",
               "Cholinergic neuron trajectory",
               "Stomach epithelial trajectory",
               "Placenta endodermal trajectory",
               "Osteoblast trajectory",
               "Chondrocyte trajectory",
               "Schwann cell and Enteric glia trajectory 1",
               "Megakaryocyte trajectory",
               "Pericardium trajectory",
               "Schwann cell and Enteric glia trajectory 2",
               "PNS glia precursor cell trajectory",
               "Shisa6 positive neuron trajectory",
               "Neuron progenitor trajectory",
               "Olfactory sensory neuron trajectory",
               "Olfactory epithelial trajectory",
               "Urothelium trajectory",
               "Olfactory ensheathing cell trajectory",
               "Lung epithelial trajectory",
               "Midgut/Hindgut epithelial trajectory",
               "Brain endothelial trajectory",
               "Ependymal cell trajectory",
               "Arterial endothelial trajectory",
               "Retina trajectory",
               "Auditory epithelial trajectory",
               "Inhibitory neuron trajectory",
               "Enteric neuron trajectory 1",
               "Endocardium trajectory"
               
               
               
               
               
)
xxx <- replace(xx, xx %in% lesscells, "other")
xxx

pbmc[["SingleR.cluster.others"]] <- 
  xxx

DimPlot(pbmc, label = TRUE,group.by = "SingleR.cluster.others") 

df = data.frame(name = pred.hpca@rownames, type = pred.hpca@listData$first.labels)


df <- data.frame(pbmc$SingleR.cluster.others)
df1 <-data.frame(X =row.names(df), count= df$pbmc.SingleR.cluster.others)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

imported_raster=OpenImageR::readImage("FFPE-2.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_sub-1.pdf",sep =""), width=12.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
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




plotScoreHeatmap(pred.hpca)

some <- table(pred.hpca$labels)
x <- dimnames(sort(some,decreasing=T))
x















head(Idents(pbmc), 5)
Idents(pbmc)
ident <- Idents(pbmc)
df <- data.frame(ident[])
df1 <-data.frame(X =row.names(df), count= df$ident..)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

imported_raster=OpenImageR::readImage("FFPE-2.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_10_08_nobackground.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
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


imported_raster=OpenImageR::readImage("FFPE-2.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_10_08_biggerpoints.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape = 15, size = 4)+
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





ref = MouseRNAseqData()
hESCs = GetAssayData(pbmc, slot = "counts")
common <- intersect(rownames(hESCs), rownames(ref))
ref <- ref[common,]
hESCs <- hESCs[common,]

pred.hpca <- SingleR(test = hESCs, ref = ref, 
                     labels = ref$label.main)
table(pred.hpca$labels)


df = data.frame(name = pred.hpca@rownames, type = pred.hpca@listData$first.labels)

pred.hpca$labels
pbmc[["SingleR.cluster.labels"]] <- 
  pred.hpca$labels
DimPlot(pbmc, label = TRUE,group.by = "SingleR.cluster.labels")

pbmc$SingleR.cluster.labels

df <- data.frame(pbmc$SingleR.cluster.labels)
df1 <-data.frame(X =row.names(df), count= df$pbmc.SingleR.cluster.labels)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

imported_raster=OpenImageR::readImage("ventricle.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_cell type_main_with tissue.pdf",sep =""), width=12.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
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


plotScoreHeatmap(pred.hpca)



