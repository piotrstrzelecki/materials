# 1. Required packages and functions --------------------------------------
# Install libraries; remove '#' to execute the line
# install.packages("dplyr")
# install.packages("psych")
# install.packages("qgraph")
# install.packages("corrplot")
# install.packages("ggplot2")
# install.packages("factoextra")
# install.packages("ggrepel")
# install.packages("tidyverse")
# install.packages("viridis")
# install.packages("plotly")
# install.packages("reshape")
# install.packages("Morpho")
# install.packages("cluster")
# install.packages("gridExtra")
# install.packages("ggpmisc")
# install.packages("openair")

# Load libraries
library(psych)
library(corrplot)
library(ggplot2)
library(qgraph)
library(factoextra)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(Morpho)
library(viridis)
library(cluster)
library(ggpmisc)
library(openair)

# 2. Load data ------------------------------------------------------------
# Set your working directory
setwd("C:/...") # COMPLETE!
data2 <- read.csv("https://raw.githubusercontent.com/piotrstrzelecki/materials/main/dataset_mm.csv", header = T, sep = ",", dec = ".", check.names = F) # load the data
res <- 0.012 # set the resolution
unit <- "mm" # set the unit

## 2.1 Column names preparation
# Prepare the names of the variables: 3 columns including original names, numerical id, names with units
{
  coln <- cbind(original = colnames(data), id = c(1:ncol(data)), ggplot = c(1:ncol(data)))
  coln[1:6, 3] <- paste0(coln[1:6, 1], " (pixel value)")
  coln[9:12, 3] <- paste0(coln[9:12, 1], " (", unit, "^3)")
  coln[13, 3] <- paste0(coln[13, 1], " (", unit, "^2)")
  coln[c(14:24, 36:38), 3] <- paste0(coln[c(14:24, 36:38), 1], " (", unit, ")")
  coln[c(7:8, 25:35), 3] <- paste0(coln[c(7:8, 25:35), 1], " (over(,phantom(0)))")
  coln[39:44, 3] <- paste0(coln[39:44, 1], " (degree)")
  coln[, 3] <- gsub(" ", "~", coln[, 3])
}

# 3.0 The data set overview - basic statistics -----------------------------
# Basic statistics
describeBy(data)

# "min intensity" shows a zero variance. It was excluded from further analyses
data <- data[, -4]
coln <- coln[-4, ]

# 4.0 Correlation between variables ---------------------------------------

colnames(data) <- coln[, 2] # switch to numerical id of variables names
# Figure 5a
corrplot(cor(data), method = "color", tl.col = "black") # corrplot::
# Figure 5b
qgraph(cor(data), posCol = "blue3", negCol = "brown3", minimum = 0, cut = 1, vsize = 5,
       legend.cex = .5, label.scale.equal = T, legend = T, borders = T, palette = "ggplot2",
       GLratio = 5, groups = list(`chemcial composition` = 1:7, `size` = 8:23, `shape` = 24:34,
                                  `spatial arrangement` = 35:43)) # qgraph::

# 4.1 Principal component analysis ----------------------------------------
# variables from the spatial arrangement category were excluded
data.pca <- prcomp(data[, 1:34], center = T, scale. = T) 
summary(data.pca)

# 4.2 Number of principal components determination ------------------------
# Figure 6a
fviz_eig(data.pca, addlabels = T, ggtheme = theme_classic()) # factoextra::

# 4.3 Factor analysis -----------------------------------------------------

data.fa <- principal(data[, 1:34], rotate = "none", nfactors = 4, scores = TRUE) #psych::

# data manipulation
fa.loadings <- data.fa$loadings[1:34, ]
fa.loadings <- data.frame(fa.loadings)
fa.loadings$var <- factor(rownames(fa.loadings), levels = rev(rownames(fa.loadings)))
loadings.m <- melt(fa.loadings, id.vars = "var")

# Figure 6b
ggplot(loadings.m, aes(x = value, y = var, fill = value)) +
  facet_wrap(~variable, nrow = 1, scales = "fixed") + # place the factors in separate facets
  geom_bar(stat = "identity") +
  scale_fill_gradient2(name = "Loading", high = "blue", mid = "white", low = "red",
                       midpoint = 0, guide = F) +
  theme_bw(base_size = 10) +
  xlab(label = "Standardized loadings") +
  ylab(label = "")

# statistics for the first and second PCs
# data manipulation
data.fa <- principal(data[, 1:34], rotate = "none", nfactors = 2, scores = TRUE) 
fa.loadings <- data.frame(data.fa$uniquenesses)
colnames(fa.loadings) <- "uniqueness"
fa.loadings$communality <- data.fa$communality
fa.loadings$complexity <- data.fa$complexity
fa.loadings$var <- factor(rownames(fa.loadings), levels = rev(rownames(fa.loadings)))

# Figure 6c
loadings.m <- melt(fa.loadings, id.vars = "var")
ggplot(loadings.m, aes(x = value, y = var, fill = variable)) +
  geom_bar(stat = "identity", fill = "grey") +
  facet_wrap(~variable, nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(axis.title.y = element_blank())

# 4.4 Data in PC coordinates ----------------------------------------------

pca.xy <- data.frame(data.fa$scores[, 1:2])

# Sample of the data set for visual assessment 
set.seed(1)
sample <- pca.xy[mcNNindex(target = as.matrix(pca.xy), query = as.matrix(kmeans(pca.xy, 
                                                                                centers = 10,nstart = 25, iter.max = 100)$centers), k = 1), ]

# Figure 7a
ggplot(pca.xy, aes(x = PC1, y = PC2)) +
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") +
  geom_point(sample, mapping = aes(x = PC1, y = PC2), size = 1, color = "red") +
  geom_label_repel(sample[order(sample[, 1]), 1:2], mapping = aes(x = PC1, y = PC2, 
                                                                  label = 1:nrow(sample)), max.overlaps = 15, box.padding = 0.4, segment.color = "red") +
  theme(legend.position = c(0.8, 0.3))

# 5.0 Cluster analysis ----------------------------------------------------
# 5.1 Optimal number of clusters determination

# Figure 8a
fviz_nbclust(scale(data[, 1:34]),  FUNcluster = clara, 
             diss = dist(scale(data[, 1:34]), method = "manhattan"), method = "gap_stat")

# 5.2 Data clustering 

clust <- clara(scale(data[, 1:34]), k = 3, metric = "manhattan", samples = 50, sampsize = 500)

# Figure 9a
ggplot() +
  geom_point(pca.xy, mapping = aes(x = PC1, y = PC2, color = as.factor(clust$clustering))) +
  labs(col = "cluster") + xlim(-2, 10) + ylim(-2.5, 5) + theme(legend.position = c(0.9, 0.2))

write.csv(clust$clustering, "clusters.csv", row.names = FALSE)

## 5.3 Quality of the clustering – silhouette statists

# Figure 8b
clust <- pam(scale(data[, 1:34]), k = length(clust$i.med),
             metric = "manhattan", medoids = clust$i.med, do.swap = F)
fviz_silhouette(clust)

# 6.0 Microstructural assessment  ---------------------------------------
# 6.1 Chemical composition, size and shape of components in clusters

data.ggplot <- data[, 1:34]
colnames(data.ggplot) <- coln[1:34, 3]
data.ggplot$cluster <- as.factor(clust$clustering)
data.ggplot <- melt(data.ggplot, id = "cluster")

# Figure 10
ggplot(data.ggplot, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin() +
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0,fill = "white") +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~variable, ncol = 5, scale = "free", labeller = label_parsed) +
  theme(legend.position = c(0.9, 0.05), axis.title.x = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())

# 6.2 Spatial arrangement of components

# Radial distance calculation
colnames(data) <- coln[, 2]
rad.dist <- data.frame(sqrt((data$`36` - ceiling(max(data$`36`)) / 2)^2 +
                              (data$`37` - ceiling(max(data$`37`)) / 2)^2 +
                              (data$`38` - ceiling(max(data$`38`)) / 2)^2))
rad.dist$cluster <- as.factor(clust$clustering)

# Radial distribution function g(r)
dr <- res #
gr <- matrix(ncol = length(unique(clust$clustering)), nrow = ceiling(max(data$`38`) / 2 / res))

for (i in 1:length(unique(clust$clustering))) {
  rd <- rad.dist[clust$clustering == i, 1]
  for (j in 1:nrow(gr)) {
    gr[j, i] <- length(rd[rd > (j * res) & rd < (j * res + dr)]) / (4 * pi * (j * res)^2 * dr)
  }
  rm(list = "rd")
}

# Smoothed radial distribution function g(r) calculation
gr.sm <- gr
for (i in 1:ncol(gr)) {
  gr.sm[, i] <- smooth.spline(gr[, i], spar = 0.2)$y
}

# Periodicity of peaks in g(r) determination
gr.his <- list()
for (i in 1:(ncol(gr.sm))) {
  gr.his[[i]] <- ((1:nrow(gr)) * res)[ggpmisc:::find_peaks(gr.sm[, i])]
  for (j in seq_along(gr.his[[i]])) {
    gr.his[[i]][j] <- gr.his[[i]][j + 1] - gr.his[[i]][j]
  }
}
for (i in 1:length(gr.his)) {
  names(gr.his)[i] <- paste0("cluster ", i)
}

# Visualization and data manipulation
gr.vis <- data.frame(gr)
for (i in 1:ncol(gr.vis)) {
  colnames(gr.vis)[i] <- paste0("cluster ", i)
}
gr.vis$x <- (1:nrow(gr)) * res
gr.vis <- melt(gr.vis, id = "x")
gr.vis$smooth <- melt(gr.sm)$value
gr.his <- melt(gr.his)

# Figure 11
r1 <- ggplot(gr.vis[gr.vis$x > 0.7, ]) +
  geom_line(mapping = aes(y = value, x = x, color = variable)) +
  theme(legend.position = "none") +
  labs(y = "g(r)") +
  labs(x = paste0("r ", "(", unit, ")")) +
  facet_wrap(~variable, ncol = 1, scale = "free_y")
r2 <- ggplot(gr.vis[gr.vis$x > 0.7, ]) +
  geom_line(mapping = aes(y = smooth, x = x, color = variable)) +
  stat_peaks(mapping = aes(y = smooth, x = x), span = 11) +
  theme(legend.position = "none") +
  labs(x = paste0("r ", "(", unit, ")"), y = "smoothed g(r) & peaks") +
  facet_wrap(~variable, ncol = 1, scale = "free_y")
r3 <- ggplot(gr.his, aes(y = L1, x = value, colour = L1)) +
  geom_boxplot() +
  theme(legend.position = "none",axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = paste0("distance between successive peaks ", "(", unit, ")"), y = "peaks distribution") +
  facet_wrap(~L1, ncol = 1, scales = "free_y")
grid.arrange(r1, r2, r3, ncol = 3)

## 6.3 Spatial orientation of components ----

# data manipulation
data.dir <- data[, 38:43]
colnames(data.dir) <- coln[38:43, 1]
data.dir$cluster <- as.factor(paste0("cluster ", clust$clustering))

# Figure 12
windRose(data.dir,
         type = "cluster", ws = "semi-axis a plunge", wd = "semi-axis a trend", 
         layout = c(nlevels(data.dir$cluster), 1), angle = 10, 
         annotate = F, breaks = seq(0, 90, 15), paddle = F,
         grid.line = 1, offset = 0, cols = rev("inferno"(4)[2:4]), 
         border = 8, key.header = "semi-axis a plunge (degrees)",
         key.footer = "", key.position = "right")

windRose(data.dir,
         type = "cluster", ws = "semi-axis b plunge", wd = "semi-axis b trend",
         layout = c(nlevels(data.dir$cluster), 1), angle = 10,
         annotate = F, breaks = seq(0, 90, 15), paddle = F,
         grid.line = 1, offset = 0, cols = rev("inferno"(4)[2:4]),
         border = 8, key.header = "semi-axis a plunge (degrees)",
         key.footer = "", key.position = "right")

windRose(data.dir,
         type = "cluster", ws = "semi-axis c plunge", wd = "semi-axis c trend",
         layout = c(nlevels(data.dir$cluster), 1), angle = 10, 
         annotate = F, breaks = seq(0, 90, 15), paddle = F,
         grid.line = 1, offset = 0, cols = rev("inferno"(4)[2:4]), 
         border = 8, key.header = "semi-axis c plunge (degrees)",
         key.footer = "", key.position = "right")
