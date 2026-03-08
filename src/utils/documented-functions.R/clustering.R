#' @title Run harmony pipeline on a Seurat Object
#' @description This function performs normalization, feature selection, scaling,
#' PCA, batch correction using Harmony, and optional UMAP and clustering on a Seurat object.
#' It is useful for integrating data across batches or technical sources.
#'
#' @param SeuratObject A Seurat object with raw or preprocessed RNA data.
#' @param batch Character string specifying the metadata column that defines batches (default: `"orig.file"`).
#' @param nPCs Integer. Number of principal components to use for Harmony and downstream steps (default: `20`).
#' @param cluster_resolution Numeric. Resolution parameter for `FindClusters()` (default: `0.5`).
#' @param k.param Integer. Number of nearest neighbors for `FindNeighbors()` (default: `20`).
#' @param show.plots Logical. Whether to plot PCA and Harmony reduction colored by batch (default: `TRUE`).
#' @param run.umap Logical. Whether to run UMAP embedding (default: `TRUE`). If `FALSE`, UMAP is skipped.
#' @param save.clusters Logical. Whether to save current clusters as a new metadata column `old_seurat_clusters` before re-clustering (default: `FALSE`).
#'
#' @return A Seurat object with batch-corrected dimensionality reductions (`harmony`), optional UMAP embedding,
#' and new clustering assignments. If `save.clusters = TRUE`, a new column `old_seurat_clusters` is added.
#'
#' @details This function assumes the input data is stored in the `"RNA"` assay. It performs:
#' \itemize{
#'   \item Normalization via `NormalizeData()`
#'   \item Variable feature selection using VST
#'   \item Scaling and PCA
#'   \item Batch correction using Harmony
#'   \item Visualization of batch effects (if `show.plots = TRUE`)
#'   \item Optional UMAP and clustering
#' }
#'
#' @examples
#' \dontrun{
#' harmonized_obj <- Harmonize(seurat_obj, batch = "sample_id", nPCs = 30)
#' }
#'
Harmonize = function(SeuratObject, batch = "orig.file", nPCs = 20, cluster_resolution = 0.5, k.param = 20, show.plots = TRUE, run.umap = TRUE, save.clusters = FALSE){
  
  # Save old clusters
  if(save.clusters) SeuratObject$old_seurat_clusters = SeuratObject$seurat_clusters
  
  DefaultAssay(SeuratObject) = "RNA"
  
  # Pre-processing
  SeuratObject <- SeuratObject %>% 
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(nPCs = nPCs, verbose = FALSE)
  
  # Run harmony
  harmony <- SeuratObject %>% RunHarmony(batch, plot_convergence = {if(show.plots) TRUE else FALSE})
  
  # Check PCA for mitigation of batche effects
  if(show.plots){
    print(plot_grid(
      DimPlot(SeuratObject, reduction = "pca", pt.size = .1, group.by = batch) + NoLegend(),
      DimPlot(harmony, reduction = "harmony", pt.size = .1, group.by = batch) + NoLegend()
    ))
  }
  
  # Downstream analysis
  if(!run.umap){
    harmony <- harmony %>%
      FindNeighbors(reduction = "harmony", dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution) %>%
      identity()
  } else {
    harmony <- harmony %>%
      RunUMAP(reduction = "harmony", dims = 1:nPCs) %>%
      FindNeighbors(reduction = "harmony", dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution) %>%
      identity()
  }
  
  return(harmony)
}


#' @title Preprocess, Integrate, and Cluster a Seurat Object
#' @description This function performs normalization, feature selection, optional batch integration,
#' dimensionality reduction, and clustering on a Seurat object. Optionally includes UMAP and t-SNE embedding.
#'
#' @param object A Seurat object to be processed.
#' @param nfeatures Integer. Number of variable features to select (default: `2000`).
#' @param numPCs Integer. Number of principal components to use for dimensionality reduction and clustering (default: `20`).
#' @param normalization.method Character. Method for normalization (default: `"LogNormalize"`).
#' @param scale.factor Numeric. Scale factor for normalization (default: `10000`).
#' @param selection.method Character. Method for selecting variable features (default: `"vst"`).
#' @param scale.all Logical. Whether to scale all genes, not just variable features (default: `FALSE`).
#' @param cluster_resolution Numeric. Resolution parameter for clustering (default: `0.5`).
#' @param elbow Logical. Whether to show an elbow plot and prompt for manual selection of number of PCs (default: `FALSE`).
#' @param integrate.by Character. Metadata column to use for batch integration. If `NULL`, no integration is performed (default: `NULL`).
#' @param k.weight Integer. Weight parameter for `IntegrateData()` (default: `100`).
#' @param do.tsne Logical. Whether to run t-SNE for visualization (default: `FALSE`).
#' @param do.umap Logical. Whether to run UMAP for visualization (default: `TRUE`).
#' @param remove.batches Integer. Minimum number of cells required to retain a batch during integration (default: `60`).
#' @param reduction Character. Reduction method for `FindIntegrationAnchors()` (default: `"cca"`). Use `"rpca"` for reciprocal PCA.
#' @param k.anchor Integer. Number of anchors to use when finding integration anchors (default: `5`).
#'
#' @return A Seurat object that has undergone normalization, scaling, PCA, and clustering, with optional batch correction and UMAP/t-SNE embedding.
#'
#' @details 
#' If `integrate.by` is provided, the function will:
#' \itemize{
#'   \item Remove underpowered batches (fewer than `remove.batches` cells),
#'   \item Split the object by batch,
#'   \item Normalize and identify variable features within each batch,
#'   \item Integrate using Seurat’s integration pipeline (CCA or RPCA),
#'   \item Restore previous `misc` slot contents,
#'   \item Ensure original batch factor levels are preserved.
#' }
#' If `elbow = TRUE`, an elbow plot is generated and the user is prompted to input the number of PCs.
#'
#' @examples
#' \dontrun{
#' clustered_obj <- ClusterSeurat(seurat_obj, integrate.by = "sample", do.tsne = TRUE)
#' }
ClusterSeurat = function(object, nfeatures = 2000, numPCs = 20, normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst",  
                         scale.all = FALSE,  cluster_resolution = .5, elbow = FALSE, integrate.by = NULL, k.weight = 100, do.tsne = FALSE, do.umap = TRUE, 
                         remove.batches = 60, reduction = "cca", k.anchor = 5){
  if(!is.null(integrate.by)){
    
    if(is.factor(object@meta.data[[integrate.by]])){
      factor.levels = levels(object@meta.data[[integrate.by]])
    } else {
      factor.levels = NULL
    }
    
    # Remove batches smaller than remove.batches
    object@meta.data[,integrate.by] = as.character(object@meta.data[,integrate.by])
    batches.remove = names(which(table(object@meta.data[,integrate.by]) < remove.batches))
    message("Removing batches: ", paste0(batches.remove, collapse = ", "))
    object = object[,!object@meta.data[,integrate.by] %in% batches.remove] #subset(object, eval(parse(text = integrate.by)) %in% batches.remove, invert = TRUE)
    misc = object@misc
    smallest_ident <- min(table(object@meta.data[,integrate.by]))
    obj.list <- SplitObject(object, split.by = integrate.by)
    rm(object)
    for (i in 1:length(obj.list)) {
      obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
      obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = selection.method, nfeatures = nfeatures, verbose = FALSE)
    }
    
    # If using RPCA
    if(reduction == "rpca"){
      features <- SelectIntegrationFeatures(object.list = obj.list)
      obj.list <- lapply(X = obj.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
    }
    
    obj.anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = reduction, k.anchor = k.anchor)
    
    if(smallest_ident < k.weight){
      print("Smallest ident less than k.weight, setting k.weight to size of smallest ident.")
      k.weight = smallest_ident
    }
    object <- IntegrateData(anchorset = obj.anchors, k.weight = k.weight)
    object@misc = misc
    DefaultAssay(object) <- "integrated"
    
    # Make integrate.by a factor again if it was one previously
    if(!is.null(factor.levels)) object@meta.data[[integrate.by]] = factor(object@meta.data[[integrate.by]], levels = factor.levels)
    
  }
  else{
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object, selection.method = selection.method, nfeatures = nfeatures, verbose = TRUE)
  }
  
  if(scale.all){
    all.genes <- rownames(object)
    object <- ScaleData(object, features = all.genes)
  }
  else{
    object <- ScaleData(object)
  }
  
  object <- RunPCA(object)
  if(elbow){
    ElbowPlot(object, ndims = 50)
    numPCs <- readline(prompt = "Enter the number of PCs desired for clustering analysis: ")
  }
  object <- FindNeighbors(object, dims = 1:numPCs)
  object <- FindClusters(object, resolution = cluster_resolution)
  if(do.tsne){
    object <- RunTSNE(object, dims = 1:numPCs)
  }
  if(do.umap){
    object <- RunUMAP(object, dims = 1:numPCs)
  }
  
  DefaultAssay(object) = "RNA"
  
  return(object)
}
