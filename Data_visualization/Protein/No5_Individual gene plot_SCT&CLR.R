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

dir = getwd()

#read expression matrix
data1 <- read.table("Filtered_matrix_correct.tsv", header = TRUE, sep = "\t", row.names = 1)
data3 <- data.frame(X = rownames(data1), data1)
temp1 <- data3 %>% separate(X, c("A", "B"),  sep = "x")

temp1$A = NULL
temp1$B = NULL
temp1$unmapped = NULL
data2 <- t(temp1)
sample1.name <- "atrium"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
sample          <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)

#use sctransform to normalize
sample <- PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")
sample <- SCTransform(sample, vars.to.regress = "percent.mt", verbose = FALSE)

imported_raster=OpenImageR::readImage("tonsil.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
proteindata <- sample[["RNA"]]@data
proteindata <- t(proteindata)
protein <-as.data.frame(as.matrix(proteindata))
protein$X = row.names(protein)
protein$X
test <- protein %>% separate(X, c("A", "B"),  sep = "x")
test$A

subdir= paste(dir, "/single_protein_SCT",sep ="")
dir.create(subdir, showWarnings = FALSE)
setwd(subdir)
protein = colnames(test)


for (i in 1:(length(colnames(test))-2)) {
 # test1 = data.frame(A = test$A, B = test$B, C =test[,protein[i]])
 test2 = test[test[protein[i]]>0, ]
  
 p <- ggplot(test2, aes(x = as.numeric(test2$A), y = as.numeric(test2$B), colour=test2[, protein[i]])) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  scale_color_gradientn(colours = c("blue","green", "red"),
                        oob = scales::squish) +
  ggtitle(protein[i]) +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 3)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(colour="black",size=30),
        axis.title=element_text(colour="black",size=30,face="bold"),
        legend.text=element_text(colour="black",size=30),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size =1, fill = NA))
 pdf(file = paste(protein[i], ".pdf",sep =""), width=8.6, height=8.6)
 print(p)
dev.off()
}

#normalize with CLR
sample          <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
sample <- NormalizeData(sample, normalization.method = "CLR", margin=2)
proteindata <- sample[["RNA"]]@data
proteindata <- t(proteindata)
protein <-as.data.frame(as.matrix(proteindata))
protein$X = row.names(protein)
protein$X
test <- protein %>% separate(X, c("A", "B"),  sep = "x")
test$A

dir.create(paste(dir, "/single_protein_CLR_asinteger",sep =""), showWarnings = FALSE)

setwd(paste(dir, "/single_protein_CLR_asinteger",sep =""))
protein = colnames(test)

for (i in 1:(length(colnames(test))-2)) {
  # test1 = data.frame(A = test$A, B = test$B, C =test[,protein[i]])
  test2 = test[test[protein[i]]>0, ]
  
  p <- ggplot(test2, aes(x = as.numeric(test2$A), y = as.numeric(test2$B), colour=test2[, protein[i]])) +
    #scale_color_gradientn(colours = c("black", "green")) + 
    scale_color_gradientn(colours = c("blue","green", "red"),
                          oob = scales::squish) +
    ggtitle(protein[i]) +
    #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
    geom_point(shape = 15, size = 3)+
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text=element_text(colour="black",size=30),
          axis.title=element_text(colour="black",size=30,face="bold"),
          legend.text=element_text(colour="black",size=30),
          legend.title = element_blank(),
          #legend.title = element_text(colour="black", size=15, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(size =1, fill = NA))
  pdf(file = paste(protein[i], ".pdf",sep =""), width=8.6, height=8.6)
  print(p)
  dev.off()
}




