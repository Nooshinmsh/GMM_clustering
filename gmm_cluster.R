library(mclust)
library(terra)

sequential_gmm <- 
function(raster_path, chunk_size = 1000, 
                          n_clusters = NULL, max_G = 15,
                          overlap_pct = 0.1) {
                            img1 <- rast(raster_path)
                            # Get raster dimensions
                            nrows <- nrow(img1)
                            ncols <- ncol(img1)
                            total_cells <- ncell(img1)
                            n_chunks <- ceiling(total_cells / chunk_size)
                            # Initialize variables
                            all_data <- NULL
                            gmm_model <- NULL
                            n_clusters_final <- n_clusters
                            
                            # Track model evolution
                            models_history <- list()
                            cat("Processing", total_cells, "cells in", n_chunks, "chunks\n")
                            for(i in 1:n_chunks) {
                               cat("\nProcessing chunk", i, "/", n_chunks, "...")
                               # Calculate chunk boundaries with overlap
                               start_cell <- max(1, (i-1) * chunk_size +1)
                               end_cell <- min(total_cells, i * chunk_size)
                                start_row <- ceiling(start_cell / ncols)
                                 end_row <- ceiling(end_cell / ncols)
                                 nrows_chunk <- end_row - start_row + 1
                                 # Read chunk
                                 chunk_vals <- tryCatch({
                                   values(img1, row = start_row, nrows = nrows_chunk)
                                   }, error = function(e) {
                                     cat("Error reading chunk:", e$message, "\n")
                                     return(NULL)
                                     })
                                     # Remove NA rows
                                     if(!is.null(chunk_vals)) {
                                       na_rows <- rowSums(is.na(chunk_vals)) > 0
                                       chunk_vals <- chunk_vals[!na_rows, , drop = FALSE]
                                       }
                                       # Skip if no data
                                        if(is.null(chunk_vals) || nrow(chunk_vals) < 100) {
                                          cat("Skipping chunk (insufficient data)")
                                           next
                                            }
                                             # If first chunk, create initial model
                                              if(is.null(gmm_model)) {
                                                 # Determine optimal G if not specified
                                                 if(is.null(n_clusters_final)) {
                                                   # Use BIC to find optimal G on first substantial chunk
                                                    cat("\nFinding optimal number of clusters...")
                                                    test_gmm <- Mclust(chunk_vals, G = 1:max_G, verbose = FALSE)
                                                    n_clusters_final <- test_gmm$G
                                                    cat(" Optimal G =", n_clusters_final)
                                                    }
                                                    # Train initial model
                                                    gmm_model <- Mclust(chunk_vals, G = n_clusters_final, 
                                                     modelNames = "VVV", verbose = FALSE)
                                                     } else 
                                                     # Update model with new data
                                                     gmm_model <- Mclust(chunk_vals, G = n_clusters_final,
                                                      modelNames = "VVV",
                                                       verbose = FALSE)
                                                        }
                                                        # Predict clusters for current chunk
                                                        predictions <- predict(gmm_model, newdata = chunk_vals)
                                                        # Calculate cell numbers for this chunk
                                                        cells_in_chunk <- seq(start_cell, end_cell)
                                                        cells_in_chunk <- cells_in_chunk[!na_rows]  # Remove NA cells
                                                        # Store results
                                                        chunk_result <- data.frame(
                                                          cell = cells_in_chunk,
                                                          cluster = predictions$classification)
                                                          all_data <- rbind(all_data, chunk_result)
                                                          models_history[[i]] <- gmm_model
                                                          cat(" Done. Clusters:", n_clusters_final, "BIC:", round(gmm_model$bic, 0))
                                                          return(list(
                                                            model = gmm_model,
                                                            clusters = all_data,
                                                             history = models_history,
                                                             optimal_G = n_clusters_final))
                          }
                                                              
 # nolint
result <- sequential_gmm(
  raster_path = "satimg/landsat_norm_bandready/img11_norm.tif",
  chunk_size = 100000,          # Smaller chunks for 91 bands
  max_G = 15,                 # Test up to 15 clusters
  overlap_pct = 0.05          # 5% overlap between chunks
)
