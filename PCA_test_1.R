#Analysis of ANALYSED flow data - B cells
setwd("~/Desktop/Flow Cytometry Analysis in R/PCA")

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

####Load required dataset####
B_MASTER_RAW <- read_excel('B_cells_MASTER_for_R_copy.xlsx') %>% as.data.frame()
B_MASTER <-B_MASTER_RAW[1:90, ]
B_MASTER_N <- B_MASTER_RAW[1:90,12:28]
B_MASTER_LABELS <- B_MASTER_RAW[1:90, 1:11]
B_MASTER_ALL <- B_MASTER_RAW[1:96,12:28]

####Modify dataset####
X_axis_order <- c("WT", "Mb1 E1020K", "Mb1 Cyp11a1 KO", "Mb1 E1020K Cyp11a1 KO")
X_axis_order_2 <- c("Young WT", "Old WT", "Young Mb1 E1020K", "Old Mb1 E1020K", "Young Mb1 Cyp11a1 KO", "Old Mb1 Cyp11a1 KO", "Young Mb1 E1020K Cyp11a1 KO", "Old Mb1 E1020K Cyp11a1 KO")
X_axis_order_3 <- c("Young WT", "Old WT", "Young Mb1 E1020K", "Old Mb1 E1020K", "Young Mb1 Cyp11a1 KO", "Old Mb1 Cyp11a1 KO", "Young Mb1 E1020K Cyp11a1 KO", "Old Mb1 E1020K Cyp11a1 KO", "Young E1020K", "Old E1020K")
X_axis_order_4 <- c("B cells (Total Lymphocyts)", "T cells (Total Lymphocytes)", "B220- B cells (Total B cells)", "GC B cells (CD38) (Total B cells)", "GC B cells (GL-7) (Total B cells)", "B1 B cells (Total B cells)", "FO B cells (Total B cells)", "MZ B cells (Total B cells)", "IgM- IgD- B cells (Total B cells)", "T1 B cells (Tota B cells)", "T2 B cells (Total B cells)", "T3 B cells (Total B cells)", "Dark Zone B cells (Total B cells)", "Light Zone B cells (Total B cells)")
X_axis_order_4_od <- sort(X_axis_order_4)


####Generate UMAP Plot####
umap.B_MASTER_N = umap(B_MASTER_N)
umap.B_MASTER_N
head(umap.B_MASTER_N$layout, 3)
plot(umap.B_MASTER_N, B_MASTER_LABELS)

####Generate PCA Plot####
pca_res <- prcomp(B_MASTER_ALL, scale. = TRUE)
autoplot(pca_res)

autoplot(pca_res, data = B_MASTER_RAW, colour = "Genotype", shape = "Age", size = 5, loadings = TRUE, loadings.colour = "grey", loadings.label = TRUE, loadings.label.size = 3, laodings.label.vjust = 20) +
  ggtitle("Immature and class-swicthed B cells\nare altered in Old Mb1 Cyp11a1 KO and Old Mb1 E1020K mice", subtitle = "Cell populations associated with differences: Germinal Center (GC), B1-like, IgM- IgD-, Dark Zone and Light Zone B cells") +
  theme_bw() +
  scale_shape_manual(values = c(4, 5)) +
  theme(plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
  
ggsave("plot_PCA.png", width = 15, height = 10)

####Plot multiple variables on same plot####
#Extract column names from datafram
B_MASTER_Col <- colnames(B_MASTER)
str(B_MASTER_Col)

#Place data into long format and plot
plot_overview_7_8 <- B_MASTER %>% select(Genotype, Genotype_2, Gender, Age, Age_days, Age_weeks, B_cells, T_cells, B220_neg, GC_B_cells_A, GC_B_cells_B, B1_B_cells, FO_B_cells, MZ_B_cells, IgMneg_IgDneg_B_cells, T1_B_cells, T2_B_cells, T3_B_cells, Dark_Zone_B_cells, Light_Zone_B_cells) %>%
  pivot_longer(., c(B_cells, T_cells, B220_neg, GC_B_cells_A, GC_B_cells_B, B1_B_cells, FO_B_cells, MZ_B_cells, IgMneg_IgDneg_B_cells, T1_B_cells, T2_B_cells, T3_B_cells, Dark_Zone_B_cells, Light_Zone_B_cells), names_to = "Var", values_to = "Val") %>%
  ggplot(aes(x = Var, y = Val, color = Gender)) +
  ggtitle("Phenotyping study of Cyp11a1 knockout and activated PI3Kdelta (E1020K) mice", subtitle = "Overview accross all genotypes, separated by age") +
  labs(y = "", x = "") +
  scale_x_discrete(label = X_axis_order_4_od) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L), breaks = scales::pretty_breaks(n = 10)) +
  scale_colour_colorblind() +
 geom_boxplot() + geom_jitter(width = 0.15) +
  facet_wrap_paginate(~fct_relevel(Genotype_2, "Young WT", "Old WT", "Young Mb1 E1020K", "Old Mb1 E1020K", "Young Mb1 Cyp11a1 KO", "Old Mb1 Cyp11a1 KO", "Young Mb1 E1020K Cyp11a1 KO", "Old Mb1 E1020K Cyp11a1 KO"), ncol = 1, nrow = 2, page = 4) +
  theme_bw() + theme(strip.text.x = element_text(size = 10, color = "black", face = "bold"),
                     strip.background = element_rect(color = "black", fill = "white", size = 1.5, linetype = "solid"),
                     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                     axis.text.y = element_text(size = 10),
                     plot.title = element_text(size = 25, face = "bold.italic"),
                     plot.subtitle = element_text(size = 15),
                     plot.margin = unit (c(2, 2, 2, 2), "cm"))
plot_overview_1
plot_overview_2
plot_overview_3
plot_overview_4
plot_overview_5
plot_overview_6
plot_overview_7
plot_overview_8

plot_overview_7_8
ggsave("plot_overview_7_8.png", width = 15, height = 10)

scale_x_discrete(label = X_axis_order_4_od) +
####Generate Dot Plot####
#Simple plot
ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype, level = X_axis_order), y = B1_B_cells, fill = Gender, shape = Gender)) +
  labs(x = "Genotypes", y = "% B1 B cells/Lymphocytes") +
  scale_y_continuous(labels = scales::percent) +
  geom_dotplot(binaxis = "y", stackdir = "center", method = "histodot", binwidth = 0.02, alpha = 0.8) +
  theme_bw()

#All information plot
plot_1 <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype, level = X_axis_order), y = B1_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "Genotypes", y = "% B1 B cells/Lymphocytes", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(0,1)) +
  scale_colour_gradient(low = "#56B1F7", high = "#132B43") +
  scale_size_continuous(range = c(3, 12)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "red", size = 20, shape = 95, alpha = 0.8) +
  geom_beeswarm(dodge.width = 0.02) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15,  angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, face = "bold", vjust = -0.15),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 2),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))

plot_1
ggsave("plot_1.png", width = 15, height = 10)

plot_2 <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = B1_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "Genotypes", y = "% B1 B cells/Lymphocytes", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#56B1F7", high = "#132B43") +
  scale_size_continuous(range = c(5, 12)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "red", size = 20, shape = 95, alpha = 0.8) +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, face = "bold", vjust = -0.15),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 2),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
plot_2
ggsave("plot_2.png", width = 15, height = 10)


####Boxplot####
ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype, level = X_axis_order), y = B1_B_cells, fill = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "Genotypes", y = "% B1 B cells/Lymphocytes") +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_gradient(low = "#56B1F7", high = "#132B43") +
  geom_boxplot() +
  theme_bw()


####Individual plots####
#Spleen weight
plot_spleen_weight <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = Spleen_Weight_mg, color = Age_weeks, shape = Gender, size = 2)) +
  labs(x = "", y = "Spleen wieght (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("Spleen weight (in mg) across genotypes", subtitle = "Spleenweight is elevated in some old Mb1 Cyp11a1 KO mice") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))

plot_spleen_weight
ggsave("plot_spleen_weight.png", width = 15, height = 10)

#B220- B cells
plot_B220neg <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = B220_neg, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% B220- B cells/Lymphocytes", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("B220- B cells", subtitle = "B220- B cells are elevated upon PI3Kdelta activation, but decrease with age") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))

plot_B220neg
ggsave("plot_B220neg.png", width = 15, height = 10)


#B cells
plot_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% B cells/Lymphocytes", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("B cells", subtitle = "% B cells of total lymphocytes are affected by inreasing age") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))

plot_B_cells
ggsave("plot_B_cells.png", width = 15, height = 10)

#B220- B cells across all genotypes
plot_all_B220neg <- ggplot(data = B_MASTER_RAW, mapping = aes(x = factor(Genotype_2, level = X_axis_order_3), y = B220_neg, color = Age_weeks, shape = Gender)) +
  labs(x = "", y = "% B220- B cells/Lymphocytes", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("B220- B cells across all genotypes", subtitle = "Old germline E1020K appear to deacrease B220- B cell propotions while old Mb1 E1020K do not") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_all_B220neg
ggsave("plot_all_B220negs.png", width = 15, height = 10)

#MZ B cells
plot_MZ_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = MZ_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% MZ B cells/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("Marginal Zone (MZ) B cells", subtitle = "% MZ B cells of total B cells is only affected by activated PI3Kdelta and not Cyp11a1 loss") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))

plot_MZ_B_cells
ggsave("plot_MZ_B_cells.png", width = 15, height = 10)

#GC B cells (1)
plot_GC_1_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = GC_B_cells_A, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% GC B cells (CD38)/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("Germinal Center (GC) (CD38) B cells", subtitle = "% GC B cells of total B cells is significantly elevated in mice with splenomegaly") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))


plot_GC_1_B_cells
ggsave("plot_GC_1_B_cells.png", width = 15, height = 10)

#GC B cells (2)
plot_GC_2_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = GC_B_cells_B, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% GC B cells (GL-7)/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("Germinal Center (GC) (GL-7) B cells", subtitle = "% GC B cells of total B cells is significantly elevated in mice with splenomegaly") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))

plot_GC_2_B_cells
ggsave("plot_GC_2_B_cells.png", width = 15, height = 10)

#IgM- IgD- B cells
plot_IgMneg_IgDneg_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = IgMneg_IgDneg_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% IgM- IgD-/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("IgM- IgD- B cells", subtitle = "% IgM- IgD- B cells of total B cells is significantly elevated in old Mb1 E1020K and old Cyp11a1 KO mice") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_IgMneg_IgDneg_B_cells
ggsave("plot_IgMneg_IgDneg_B_cells.png", width = 15, height = 10)

#B1 B cells
plot_B1_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = B1_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% B1 B cells/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("B1-like B cells ", subtitle = "% B1-like B cells of total B cells is significantly elevated in old Mb1 E1020K and old Cyp11a1 KO mice") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_B1_B_cells
ggsave("plot_B1_B_cells.png", width = 15, height = 10)

#Dark Zone B cells
plot_DZ_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = Dark_Zone_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% Dark Zone B cells/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("Dark Zone (DZ) B cells ", subtitle = "% DZ B cells of total B cells is  elevated in old Mb1 E1020K and old Cyp11a1 KO mice") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_DZ_B_cells
ggsave("plot_DZ_B_cells.png", width = 15, height = 10)


#Light Zone B cells
plot_LZ_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = Light_Zone_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "% Dark Zone B cells/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("Light Zone (LZ) B cells ", subtitle = "% LZ B cells of total B cells is  elevated in old Mb1 E1020K and old Cyp11a1 KO mice") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_LZ_B_cells
ggsave("plot_LZ_B_cells.png", width = 15, height = 10)

#T cells
plot_T_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = T_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "T cells/Lymphocytes", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("T cells ", subtitle = "% T cells of total lymphocytes is  similar across genotypes") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_T_cells
ggsave("plot_T_cells.png", width = 15, height = 10)


#T1 B cells
plot_T1_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = T1_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "T1 B cells/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("T1 B cells ", subtitle = "") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_T1_B_cells
ggsave("plot_T1_B_cells.png", width = 15, height = 10)

#T2 B cells
plot_T2_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = T2_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "T2 B cells/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("T2 B cells ", subtitle = "") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_T2_B_cells
ggsave("plot_T2_B_cells.png", width = 15, height = 10)

#T3 B cells
plot_T3_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = T3_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "T3 B cells/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("T3 B cells ", subtitle = "") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_T3_B_cells
ggsave("plot_T3_B_cells.png", width = 15, height = 10)

#FO B cells
plot_FO_B_cells <- ggplot(data = B_MASTER, mapping = aes(x = factor(Genotype_2, level = X_axis_order_2), y = FO_B_cells, color = Age_weeks, shape = Gender, size = Spleen_Weight_mg)) +
  labs(x = "", y = "FO B cells/B cells", size = "Spleen weight (in mg)", color = "Age (in weeks)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 10)) +
  scale_x_discrete(expand = c(-0.1,1)) +
  scale_colour_gradient(low = "#E69F00", high = "#332288") +
  scale_size_continuous(range = c(0.5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 20, shape = 95, alpha = 1.2) +
  ggtitle("FO B cells ", subtitle = "") +
  geom_beeswarm() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, face = "bold.italic"),
        plot.subtitle = element_text(size = 15),
        plot.margin = unit (c(1, 1, 1, 1), "cm"))
plot_FO_B_cells
ggsave("plot_FO_B_cells.png", width = 15, height = 10)

####Statistical analysis####
##B_cells
lm.B_cells <- lm(B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
anova_B_cells <- aov(lm.B_cells)
summary(anova_B_cells)
step(lm.B_cells)
TukeyHSD(anova_B_cells)

#Diagnostics B_cells
leveneTest(B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)#does not work
diagnostics_B_cell <- plot(lm.B_cells)
par(mfrow=c(2,2))
diagnostics_B_cell

#Simplified
lm.B_cells_2 <- lm(B_cells~Genotype_2, data = B_MASTER)
plot(lm.B_cells_2)
anova_B_cells_2 <- aov(lm.B_cells_2)
TukeyHSD(anova_B_cells_2)

##B1_B_cells
lm.B1_B_cells <- lm(B1_B_cells~Genotype_2, data = B_MASTER)
plot(lm.B1_B_cells)
anova_B1_B_cells <- aov(lm.B1_B_cells)
tukey_B1_B_cells <- TukeyHSD(anova_B1_B_cells)
plot(tukey_B1_B_cells)

lm.B1_B_cells <- lm(B1_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.B1_B_cells)

##B220- B cells
lm.B220neg_B_cells <- lm(B220_neg~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.B220neg_B_cells)
lm.B220_neg_cells <- lm(B220_neg~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.B220_neg_cells)

##MZ B cells
lm.MZ_B_cells <- lm(MZ_B_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.MZ_B_cells)
lm.MZ_B_cells <- lm(MZ_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.MZ_B_cells)

##IgM- IgD- B cells
lm.IgMneg_IgDneg_B_cells <- lm(IgMneg_IgDneg_B_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.IgMneg_IgDneg_B_cells)
lm.IgMneg_IgDneg_B_cells <- lm(IgMneg_IgDneg_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.MZ_B_cells)

##GC B cells (1)
lm.GC_B_cells_A <- lm(GC_B_cells_A~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.GC_B_cells_A)
lm.GC_B_cells_A <- lm(GC_B_cells_A~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.GC_B_cells_A)

##GC B cells (2)
lm.GC_B_cells_B <- lm(GC_B_cells_B~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.GC_B_cells_B)
lm.GC_B_cells_B <- lm(GC_B_cells_A~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.GC_B_cells_B)

##DZ B cells
lm.DZ_B_cells <- lm(Dark_Zone_B_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.DZ_B_cells)
lm.DZ_B_cells <- lm(Dark_Zone_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.DZ_B_cells)

##LZ B cells
lm.LZ_B_cells <- lm(Light_Zone_B_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.LZ_B_cells)
lm.LZ_B_cells <- lm(Light_Zone_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.LZ_B_cells)

##T cells
lm.T_cells <- lm(T_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.T_cells)
lm.T_cells <- lm(T_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.T_cells)

##T1 B cells
lm.T1_B_cells <- lm(T1_B_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.T1_B_cells)
lm.T1_B_cells <- lm(T1_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.T1_B_cells)

##T2 B cells
lm.T2_B_cells <- lm(T2_B_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.T2_B_cells)
lm.T2_B_cells <- lm(T2_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.T2_B_cells)

##T3 B cells
lm.T3_B_cells <- lm(T3_B_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.T3_B_cells)
lm.T3_B_cells <- lm(T3_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.T3_B_cells)

##FO B cells
lm.FO_B_cells <- lm(FO_B_cells~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.FO_B_cells)
lm.FO_B_cells <- lm(FO_B_cells~Genotype*Gender*Spleen_Weight_mg*Age_weeks, data = B_MASTER)
step(lm.FO_B_cells)

##Spleen weight
lm.Spleen_weight <- lm(Spleen_Weight_mg~Genotype_2, data = B_MASTER)
par(mfrow=c(2,2))
plot(lm.Spleen_weight)
lm.Spleen_weight <- lm(Spleen_Weight_mg~Genotype*Gender*Age_weeks, data = B_MASTER)
step(lm.Spleen_weight)
