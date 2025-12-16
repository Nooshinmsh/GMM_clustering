# Sequential GMM Clustering of Landsat Multi-Band Images

This repository contains R code for **unsupervised clustering** of **multi-band Landsat imagery** using **Gaussian Mixture Models (GMM)**.  
The implementation is designed to handle **large raster datasets** by processing the image sequentially in chunks, making it memory efficient for high-dimensional data (e.g., >90 bands).

The sequential_gmm function processes a raster image in chunks and applies a GMM for clustering. It supports automatic selection of the optimal number of clusters (G) using the **Bayesian Information Criterion (BIC)**.

## Features
- Unsupervised clustering using **Gaussian Mixture Models (GMM)**
- Automatic selection of the optimal number of clusters using **BIC**
- Chunk-wise processing for large raster images
- Supports high-dimensional Landsat predictor stacks
- Built using `mclust` and `terra`
