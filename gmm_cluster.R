library(mclust)
library(terra)

sequential_gmm <- function(raster_path, chunk_size = 1000, n_clusters = NULL, 
                           max_G = 15, overlap_pct = 0.1) {
  
  # Load raster
  img <- rast(raster_path)
  
  # Get raster dimensions
  nrows <- nrow(img)
  ncols <- ncol(img)
  total_cells <- ncell(img)
  n_chunks <- ceiling(total_cells / chunk_size)
  
  # Initialize
  models_history <- list()
  chunk_results <- vector("list", n_chunks)
  gmm_model <- NULL
  n_clusters_final <- n_clusters
  all_chunk_data <- NULL
  
  cat("Processing", total_cells, "cells in", n_chunks, "chunks\n")
  
  for (i in 1:n_chunks) {
    cat("\nProcessing chunk", i, "/", n_chunks, "...")
    
    start_cell <- max(1, (i - 1) * chunk_size + 1)
    end_cell <- min(total_cells, i * chunk_size)
    all_cells_in_chunk <- seq(start_cell, end_cell)
    
    start_row <- ceiling(start_cell / ncols)
    end_row <- ceiling(end_cell / ncols)
    nrows_chunk <- end_row - start_row + 1
    
    # Read chunk values
    chunk_vals <- tryCatch({
      values(img, row = start_row, nrows = nrows_chunk)
    }, error = function(e) {
      cat(" Error reading chunk:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(chunk_vals) || !is.matrix(chunk_vals)) {
      cat(" Invalid chunk data\n")
      chunk_results[[i]] <- data.frame(cell = all_cells_in_chunk, cluster = NA)
      next
    }
    
    valid_rows <- rowSums(is.na(chunk_vals)) == 0
    chunk_vals_clean <- chunk_vals[valid_rows, , drop = FALSE]
    valid_cells <- all_cells_in_chunk[valid_rows]
    
    if (nrow(chunk_vals_clean) < 50) {
      cat(" Skipping chunk (insufficient valid data:", nrow(chunk_vals_clean), "samples)\n")
      chunk_results[[i]] <- data.frame(cell = all_cells_in_chunk, cluster = NA)
      next
    }
    
    cat(" (", nrow(chunk_vals_clean), "valid samples of", length(all_cells_in_chunk), "total)")
    
    if (is.null(gmm_model)) {
      # First chunk - determine optimal G
      if (is.null(n_clusters_final)) {
        cat("\n  Finding optimal number of clusters...")
        test_gmm <- tryCatch({
          Mclust(chunk_vals_clean, G = 1:max_G, verbose = FALSE)
        }, error = function(e) {
          cat(" Model fitting failed:", e$message)
          return(NULL)
        })
        
        if (!is.null(test_gmm)) {
          n_clusters_final <- test_gmm$G
          cat(" Optimal G =", n_clusters_final)
        } else {
          cat("\n  Using default G = 5")
          n_clusters_final <- 5
        }
      }
      
      # Train initial model
      gmm_model <- tryCatch({
        Mclust(chunk_vals_clean, G = n_clusters_final, modelNames = "VVV", verbose = FALSE)
      }, error = function(e) {
        cat(" Initial model failed:", e$message)
        cat("\n  Trying simpler model (EII)...")
        tryCatch({
          Mclust(chunk_vals_clean, G = n_clusters_final, modelNames = "EII", verbose = FALSE)
        }, error = function(e2) {
          cat(" Fallback model also failed")
          return(NULL)
        })
      })
      
      if (!is.null(gmm_model)) {
        all_chunk_data <- chunk_vals_clean
        predictions <- predict(gmm_model, newdata = chunk_vals_clean)
        predictions_class <- predictions$classification
      } else {
        cat("\n  Assigning random clusters")
        predictions_class <- sample(1:n_clusters_final, nrow(chunk_vals_clean), replace = TRUE)
      }
      
    } else {
      # Update model
      combined_data <- rbind(all_chunk_data, chunk_vals_clean)
      max_training_samples <- 50000
      if (nrow(combined_data) > max_training_samples) {
        set.seed(123)
        training_data <- combined_data[sample(nrow(combined_data), max_training_samples), , drop = FALSE]
      } else {
        training_data <- combined_data
      }
      
      cat("\n  Updating model with", nrow(training_data), "samples...")
      gmm_model <- tryCatch({
        Mclust(training_data, G = n_clusters_final, modelNames = "VVV", verbose = FALSE)
      }, error = function(e) {
        cat(" Model update failed:", e$message)
        return(gmm_model)
      })
      
      all_chunk_data <- combined_data
      if (nrow(all_chunk_data) > 2 * max_training_samples) {
        all_chunk_data <- tail(all_chunk_data, max_training_samples)
      }
      
      if (!is.null(gmm_model)) {
        predictions <- predict(gmm_model, newdata = chunk_vals_clean)
        predictions_class <- predictions$classification
      } else {
        cat("\n  Assigning fallback random clusters")
        predictions_class <- sample(1:n_clusters_final, nrow(chunk_vals_clean), replace = TRUE)
      }
    }
    
    # Fill all cells (with NA for invalid ones)
    all_clusters <- rep(NA, length(all_cells_in_chunk))
    all_clusters[valid_rows] <- predictions_class
    
    chunk_results[[i]] <- data.frame(
      cell = all_cells_in_chunk,
      cluster = all_clusters
    )
    
    models_history[[i]] <- gmm_model
    
    if (!is.null(gmm_model)) {
      cat("\n  Done. Clusters:", n_clusters_final, 
          " BIC:", ifelse(!is.null(gmm_model$bic), round(gmm_model$bic, 0), "NA"))
    } else {
      cat("\n  Done.")
    }
  }
  
  # Combine results
  all_data <- do.call(rbind, chunk_results)
  all_data <- all_data[order(all_data$cell), ]
  all_data <- all_data[!is.na(all_data$cell), ]
  
  return(list(
    model = gmm_model,
    clusters = all_data,
    history = models_history,
    optimal_G = n_clusters_final
  ))
}

# Test 
result <- sequential_gmm(
  raster_path = "satimg/landsat_norm_bandready/img33_norm.tif",
  chunk_size = 10000,
  n_clusters = 5,
  max_G = 10
)
