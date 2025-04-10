# sc_pipeline/r_scripts/run_soupx.R

# Load required libraries
suppressPackageStartupMessages({
  library(SoupX)
  library(Matrix)
  library(DropletUtils)
  library(Seurat)
  library(jsonlite)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# First argument is the workspace directory
workspace_dir <- args[1]
args <- args[-1]  # Remove workspace_dir from args

# Parse remaining arguments
arg_pairs <- seq(1, length(args), by=2)
arg_names <- gsub("^--", "", args[arg_pairs])
arg_values <- args[arg_pairs + 1]
named_args <- setNames(arg_values, arg_names)

# Define paths for output files
output_matrix_path <- file.path(workspace_dir, named_args$output_matrix)
output_metadata_path <- file.path(workspace_dir, named_args$output_metadata)

# Load the data
raw_data <- Read10X_h5(named_args$raw_counts_path)
filtered_data <- Read10X_h5(named_args$filtered_counts_path)

# Collect messages for output
all_messages <- c()

# Check if data was multimodal
if(typeof(raw_data) == "list") {
  raw_data <- raw_data$`Gene Expression`
}
if(typeof(filtered_data) == "list") {
  filtered_data <- filtered_data$`Gene Expression`
}

# Ensure gene dimensions match by finding common genes
common_genes <- intersect(rownames(raw_data), rownames(filtered_data))
msg <- paste0("Found ", length(common_genes), " common genes between raw and filtered data")
all_messages <- c(all_messages, msg)

if (length(common_genes) == 0) {
  # Prepare error output
  metadata <- list(
    success = FALSE,
    error = "No common genes found between raw and filtered data",
    messages = all_messages
  )
  write_json(metadata, output_metadata_path)
  stop("No common genes found between raw and filtered data")
}

# Subset both matrices to the common genes
raw_gene_idx <- match(common_genes, rownames(raw_data))
filtered_gene_idx <- match(common_genes, rownames(filtered_data))

raw_counts_subset <- raw_data[raw_gene_idx, ]
filtered_counts_subset <- filtered_data[filtered_gene_idx, ]

# Prepare SoupX input
msg <- "Creating SoupX channel..."
all_messages <- c(all_messages, msg)
sc <- SoupX::SoupChannel(raw_counts_subset, filtered_counts_subset)

# START CLUSTERING WITH SEURAT TO PROVIDE CLUSTERING TO SOUPX
msg <- "Creating Seurat object..."
all_messages <- c(all_messages, msg)
seurat_obj <- CreateSeuratObject(counts = filtered_counts_subset)

# Standard Seurat workflow with SCTransform
msg <- "Regularizing data via SCTransform..."
all_messages <- c(all_messages, msg)
seurat_obj <- SCTransform(seurat_obj, verbose = F)

msg <- "Running PCA..."
all_messages <- c(all_messages, msg)
seurat_obj <- RunPCA(seurat_obj)

# Determine number of PCs to use
ndims <- min(args$ndims %||% 30, ncol(seurat_obj@reductions$pca@cell.embeddings))

msg <- "Running UMAP..."
all_messages <- c(all_messages, msg)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:ndims, verbose = F)

msg <- paste0("Finding neighbors using...")
all_messages <- c(all_messages, msg)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:ndims, verbose = F)

msg <- "Finding clusters..."
all_messages <- c(all_messages, msg)
resolution <- args$resolution %||% 0.8
seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

# Extract clusters
meta <- seurat_obj@meta.data

# Set clusters in SoupX
msg <- "Setting clusters in SoupX..."
all_messages <- c(all_messages, msg)
sc <- SoupX::setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))

# Run SoupX
msg <- "Running SoupX contamination estimation..."
all_messages <- c(all_messages, msg)

# Determine if we should use tfidf based on dataset size. See docs for info about smaller datasets causing issues
use_tfidf <- TRUE
if (ncol(filtered_counts_subset) < 500) {
  msg <- "Small dataset detected, disabling tfidf"
  all_messages <- c(all_messages, msg)
  use_tfidf <- FALSE
}

sc <- SoupX::autoEstCont(sc, tfidf = use_tfidf)
contamination <- round(sc$fit$prop.cont, 4)

msg <- paste0("Estimated contamination fraction: ", contamination)
all_messages <- c(all_messages, msg)

# Adjust counts
msg <- "Correcting counts..."
all_messages <- c(all_messages, msg)
out <- SoupX::adjustCounts(sc)

# Save results to matrix output file
writeMM(out, output_matrix_path)
# Save row and column names for checking before adding back to AnnData
writeLines(rownames(out), file.path(workspace_dir, "features.txt"))
writeLines(colnames(out), file.path(workspace_dir, "barcodes.txt"))
msg <- "SoupX processing complete."
all_messages <- c(all_messages, msg)

# Prepare results
metadata <- list(
  success = TRUE,
  contamination_fraction = contamination,
  messages = all_messages
)

write_json(metadata, output_metadata_path, auto_unbox = TRUE)

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x