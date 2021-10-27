## ##Analysis of ANALYSED flow data - T cells####
setwd("~/Google Drive/Shared drives/Okkengroup/Experiments/Julius/Experiments/Phenotyping_Cyp11a1_E1020K_Summary/All_Populations_Excel/T_cell_MASTER")

####Load all required packages####
library(tidyverse)
library(ggfortify)
library(readxl)
library(geaom)
library(ggbeeswarm)
library(RColorBrewer)
library(viridis)
library(ggthemes)
library(ggforce)
library(car)
library(umap)


####Load dataset####
T_cell_master <- read_excel('Populations_accross_Experiments_T_cells_MASTER_R.xlsx') %>% as.data.frame()

####Modify dataset####
T_cell_master_numeric <- T_cell_master[1:96,12:28]
T_cell_master_labels <- T_cell_master[1:96, 1:12]


head(T_cell_master)
str(T_cell_master)

####Generate PCA Plot####
pca_res <- prcomp(T_cell_master_numeric, scale. = TRUE)
PCA_experiment_T <- autoplot(pca_res, data = T_cell_master, colour = "Experiment", size = 2.5) +
  theme_bw() +
  ggtitle("Principle Component Analysis - Date of Experiment", subtitle ="No batch effect based on experiments") +
  theme(plot.title = element_text(size = 20,face = "bold"),
       plot.subtitle = element_text(size = 15),
       axis.title.x = element_text(size = 15),
       axis.title.y = element_text(size = 15),
       plot.margin = unit (c(1, 1, 1, 1), "cm"))


ggsave("PCA_experiment_T.png", width = 15, height = 10)

####Generate UMAP Plot####
umap.T_cell_master_numeric = umap(T_cell_master_numeric)
umap.T_cell_master_numeric <- data.frame(x = umap.T_cell_master_numeric$layout[ ,1],
                                         y = umap.T_cell_master_numeric$layout[ ,2],
                                         Age = T_cell_master[ ,4])

ggplot(umap.T_cell_master_numeric, aes (x, y, colour = Age)) + geom_point() + theme_bw()

####Generate tSNE Plot####

