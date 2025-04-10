# sc_pipeline/r_scripts/run_soupx.R

# This script is called via the RBridge, which passes args and data via Pickle files

# Load required libraries
suppressPackageStartupMessages({
  library(SoupX)
  library(Matrix)
  library(DropletUtils)
  library(Seurat)
  library(dplyr)
})

# Get input and output arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Load input data
input_data <- readRDS(input_file)
args <- input_data$args
data <- input_data$data

# Load the data
raw_data <- Read10X_h5(args$raw_counts_path)
filtered_data <- Read10X_h5(args$filtered_counts_path)

# Check if data was multimodal
if(typeof(raw_data) == "list") {
  raw_data <- raw_data$`Gene Expression`
}
if(typeof(filtered_data) == "list") {
  filtered_data <- filtered_data$`Gene Expression`
}

# Ensure gene dimensions match by finding common genes
common_genes <- intersect(rownames(raw_data), rownames(filtered_data))
message("Found ", length(common_genes), " common genes between raw and filtered data")

if (length(common_genes) == 0) {
  # Prepare error output
  results <- list(
    success = FALSE,
    error = "No common genes found between raw and filtered data",
    messages = c("Found 0 common genes between raw and filtered data")
  )
  saveRDS(results, output_file)
  stop("No common genes found between raw and filtered data")
}

# Subset both matrices to the common genes
raw_gene_idx <- match(common_genes, rownames(raw_data))
filtered_gene_idx <- match(common_genes, rownames(filtered_data))

raw_counts_subset <- raw_data[raw_gene_idx, ]
filtered_counts_subset <- filtered_data[filtered_gene_idx, ]

# Prepare SoupX input
message("Creating SoupX channel...")
sc <- SoupX::SoupChannel(raw_counts_subset, filtered_counts_subset)

# Collect messages for output
all_messages <- c()

# Generate clusters using Seurat if not provided
if (is.null(data$clusters)) {
  msg <- "No clusters provided. Generating clusters with Seurat..."
  message(msg)
  all_messages <- c(all_messages, msg)
  
  # Create Seurat object
  msg <- "Creating Seurat object..."
  message(msg)
  all_messages <- c(all_messages, msg)
  seurat_obj <- CreateSeuratObject(counts = filtered_counts_subset)
  
  # Standard Seurat workflow
  msg <- "Regularizing data via SCTransform..."
  message(msg)
  all_messages <- c(all_messages, msg)
  seurat_obj <- SCTransform(seurat_obj, verbose = F)
  
  
  msg <- "Running PCA..."
  message(msg)
  all_messages <- c(all_messages, msg)
  seurat_obj <- RunPCA(seurat_obj)
  
  # Determine number of PCs to use
  ndims <- min(args$ndims %||% 30, ncol(seurat_obj@reductions$pca@cell.embeddings))
  
msg <- "Running UMAP..."
  message(msg)
  all_messages <- c(all_messages, msg)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:ndims, verbose = F)


  msg <- paste0("Finding neighbors using...")
  message(msg)
  all_messages <- c(all_messages, msg)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:ndims, verbose = F)
  
  msg <- "Finding clusters..."
  message(msg)
  all_messages <- c(all_messages, msg)
  resolution <- args$resolution %||% 0.8
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  # Extract clusters
  meta <- seurat_obj@meta.data
  
  # Set clusters in SoupX
  msg <- "Setting clusters in SoupX..."
  message(msg)
  all_messages <- c(all_messages, msg)
  sc <- SoupX::setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
} else {
  msg <- "Using provided clusters..."
  message(msg)
  all_messages <- c(all_messages, msg)
  sc <- SoupX::setClusters(sc, unlist(data$clusters))
}

# Run SoupX
msg <- "Running SoupX contamination estimation..."
message(msg)
all_messages <- c(all_messages, msg)

# Determine if we should use tfidf based on dataset size
use_tfidf <- TRUE
if (ncol(filtered_counts_subset) < 500) {
  msg <- "Small dataset detected, disabling tfidf"
  message(msg)
  all_messages <- c(all_messages, msg)
  use_tfidf <- FALSE
}

sc <- SoupX::autoEstCont(sc, tfidf = use_tfidf)
contamination <- round(sc$fit$prop.cont, 4)

msg <- paste0("Estimated contamination fraction: ", contamination)
message(msg)
all_messages <- c(all_messages, msg)

# Adjust counts
msg <- "Correcting counts..."
message(msg)
all_messages <- c(all_messages, msg)
out <- SoupX::adjustCounts(sc)

# Prepare results
results <- list(
  success = TRUE,
  corrected_counts = out,
  contamination_fraction = contamination,
  messages = all_messages
)

# Save results
saveRDS(results, output_file)
msg <- "SoupX processing complete."
message(msg)

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x