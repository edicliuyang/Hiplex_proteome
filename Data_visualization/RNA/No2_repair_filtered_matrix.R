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

#replace YOUR_EXPRESSION_MATRIX.tsv with your own .tsv file
data1 <- read.table(file = 'YOUR_EXPRESSION_MATRIX.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)

#extract the coordinates of each pixel
temp1 <- data1 %>% separate(X, c("A", "B"),  sep = "x")

#repair the strips of the data, if not you can skip this part..
#the repair was done by averaging the two neighboring columns.
col1 = 34
col2 = 28
row = 26
row2 = 34
row3 = 48
row4 = 44


#repair Row
for (i in 1:50) {
  temp1[(temp1$A==i&temp1$B==row),] = as.integer((as.integer(temp1[(temp1$A==i&temp1$B==row+1),]) + as.integer(temp1[(temp1$A==i&temp1$B==row-1),])) * 0.5)
  #temp1[(temp1$A==i&temp1$B==row2),] = as.integer((as.integer(temp1[(temp1$A==i&temp1$B==row2+1),]) + as.integer(temp1[(temp1$A==i&temp1$B==row2-1),])) *0.5)
  #temp1[(temp1$A==i&temp1$B==row3),] = as.integer((as.integer(temp1[(temp1$A==i&temp1$B==row3+1),]) + as.integer(temp1[(temp1$A==i&temp1$B==row3-1),])) *0.5)
  #temp1[(temp1$A==i&temp1$B==row4),] = as.integer((as.integer(temp1[(temp1$A==i&temp1$B==row4+1),]) + as.integer(temp1[(temp1$A==i&temp1$B==row4-1),])) *0.5)
}

#repair column
for (i in 1:50) {
  temp1[(temp1$A==col1&temp1$B==i),] = as.integer((as.integer(temp1[(temp1$A==(col1-1)&temp1$B==i),]) + as.integer(temp1[(temp1$A==(col1+1)&temp1$B==i),])) *0.5)
  #temp1[(temp1$A==col2&temp1$B==i),] = as.integer((as.integer(temp1[(temp1$A==(col2-1)&temp1$B==i),]) + as.integer(temp1[(temp1$A==(col2+1)&temp1$B==i),])) *0.5)
}




temp1 = data.frame(X=paste0(temp1$A, "x", temp1$B), temp1)

temp1$A = NULL
temp1$B = NULL

location <- read.table("position.txt", sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])
x = x[-1]

data_filtered <- temp1[temp1$X %in% x,]
write.table(data_filtered, file = 'Filtered_matrix_correct.tsv', sep = '\t',col.names=TRUE, row.names = FALSE,quote = FALSE)