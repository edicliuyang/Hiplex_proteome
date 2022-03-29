library(ggplot2)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(dplyr)


##read expression matrix
my_data <- read.table(file = 'Filtered_matrix_correct.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)
data_filtered <- my_data

##Count UMIs per pixel
UMI_count <- rowSums(data_filtered[,2:(ncol(data_filtered)-1)])

#Count genes per pixel
data_filtered_binary <- data_filtered[,2:ncol(data_filtered)] %>% mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)

range_UMI <- 2000  #change the x axis maxium of UMI
range_gene = 200 #change the x axis maxium of gene

##UMI Count
test <- data_filtered %>% separate(X, c("A", "B"),  sep = "x")
df <- data.frame(number=1, c=UMI_count)
pdf(file = paste("UMI_prerun.pdf",sep =""), width=8.6, height=8.6)
ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=..density..),binwidth=range_UMI/20, color="black", fill="white",size=1)+ 
  geom_density(alpha=.2, fill="#FF6666",size=1,color ="red") +
  scale_x_continuous(name="UMI",limits = c(0,range_UMI)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
  #xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(colour="black",size=20),
        axis.title=element_text(colour="black",size=25,face="bold"),
        legend.text=element_text(colour="black",size=20),
        legend.title = element_text(colour="black", size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
dev.off()

##gene Count
df <- data.frame(number=1, c=gene_count)
pdf(file = paste("gene_prerun.pdf",sep =""), width=8.6, height=8.6)
ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=..density..),binwidth=range_gene/20, color="black", fill="white",size=1)+ 
  geom_density(alpha=.2, fill="#FF6666",size=1,color ="red") +
  scale_x_continuous(name="gene",limits = c(0,range_gene)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
  #xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(colour="black",size=20),
        axis.title=element_text(colour="black",size=25,face="bold"),
        legend.text=element_text(colour="black",size=20),
        legend.title = element_text(colour="black", size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
dev.off()

imported_raster=OpenImageR::readImage("*.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)

#UMI heatmap
pdf(file = paste("UMI_heatmap_prerun.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=UMI_count)) +
  scale_color_gradientn(colours = c("blue","green", "red"),limits=c(0,range_UMI),
                        oob = scales::squish) +
  ggtitle("UMI") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
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

#gene heatmap
pdf(file = paste("gene_heatmap_prerun.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=gene_count)) +
  scale_color_gradientn(colours = c("blue","green", "red"),limits=c(0,range_gene),
                        oob = scales::squish) +
  ggtitle("gene") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
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


