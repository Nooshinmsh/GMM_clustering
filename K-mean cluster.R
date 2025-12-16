library(tidyverse)
library(dplyr)
library(ggplot2)
library(tibble)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(viridis)
library(ggthemes)
library(stringr)
library(ragg)
library(future.apply)
library(terra)
library(cluster)
library(bigstatsr)
kmax = 30
ks = 2: kmax
#read raster data
img1 = rast("satimg/landsat_norm_bandready/img11_norm.tif")

#calculate the best score for band 1-25
bands = img1[[1:25]]

#elbow method

# Extract pixel values (matrix: n_pixels × 20)
vals <- values(bands)
vals <- vals[complete.cases(vals), ]   # remove rows with any NA

#  Sample pixels (important)
set.seed(123)
n_sample <- 200000
if (nrow(vals) > n_sample) {
  vals <- vals[sample(nrow(vals), n_sample), ]
}

# Elbow curve: total within-cluster sum of squares (WSS)

wss <- sapply(ks, function(k) {
  kmeans(vals, centers = k, nstart = 30)$tot.withinss
})

elbow_df <- data.frame(k = ks, wss = wss)

# Plot elbow
ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Elbow Method Bands 1 to 25",
    x = "Number of clusters",
    y = "Total within cluster"
  )
