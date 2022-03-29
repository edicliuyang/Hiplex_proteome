library(Seurat)
library(SeuratData)
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

#change filename1 to name of txt file you want to load
data1 <- read.table("Filtered_matrix_correct.tsv", header = TRUE, sep = "\t", row.names = 1)

data3 <- data.frame(X = rownames(data1), data1)
temp1 <- data3 %>% separate(X, c("A", "B"),  sep = "x")

#repair the strips of the data, here is the x= 34 and 46, y = 20. 

temp1$A = NULL
temp1$B = NULL
temp1$unmapped = NULL
data2 <- t(temp1)
sample1.name <- "atrium"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
sample          <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)

#clustering based on SCT
sample <- PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")
sample <- SCTransform(sample, vars.to.regress = "percent.mt", verbose = FALSE)
sample <- RunPCA(sample, verbose = FALSE)
sample <- RunUMAP(sample, dims = 2:20, verbose = FALSE)
sample <- FindNeighbors(sample, dims = 2:20, verbose = FALSE)
sample <- FindClusters(sample, resolution=0.8, verbose = FALSE)
DimPlot(sample, label = TRUE) + NoLegend()
ident <- Idents(sample)
df <- data.frame(ident[])
df1 <-data.frame(X =row.names(df), count= df$ident..)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

#dot size = 3
imported_raster=OpenImageR::readImage("*.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_SCT_new.pdf",sep =""), width=8.6, height=8.6)
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
  theme(plot.title = element_text(hjust = 0.8, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()


#clustering based on CLR
sample <- CreateSeuratObject(matrix1.data, min.cells = 20, project = sample1.name)
sample <- NormalizeData(sample, normalization.method = 'CLR', margin = 2) 
sample <- ScaleData(sample)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 1000)
sample <- RunPCA(sample, verbose = FALSE)
sample <- RunUMAP(sample, dims = 2:20, verbose = FALSE)

sample <- FindNeighbors(sample, dims = 2:20, verbose = FALSE)
sample <- FindClusters(sample, resolution=0.8, verbose = FALSE)

ident <- Idents(sample)
df <- data.frame(ident[])
df1 <-data.frame(X =row.names(df), count= df$ident..)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

pdf(file = paste("clustering-CLR_new.pdf",sep =""), width=8.6, height=8.6)
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






