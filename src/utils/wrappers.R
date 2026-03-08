# Wrappers for functions


TrainTestModels_wrapper = function(object, 
                                   assay = 'RNA', 
                                   plot = FALSE, 
                                   expressed.genes = NULL, 
                                   proportion = 0.6, 
                                   LENGTH.GENE.SET = 30, 
                                   size = 100, 
                                   ...){
  
  # conf.obj = readRDS('../../patchseq_data/seurat_objects/rgc-pch-sc-int-v1.rds')
  # conf.obj = readRDS('../../patchseq_data/seurat_objects/rgc-pch-sc-int-5-100.rds')
  
  # Train on the scRNA-seq data
  scrna = subset(object, orig.ident == 'scRNA-seq')
  pch.seq = subset(object, orig.ident == 'patch-seq')
  
  # Use top 30 DEGs per cluster
  de = TopNDEGs(scrna, n = LENGTH.GENE.SET, group.by = 'type', assay = assay)
  if(is.null(expressed.genes)) expressed.genes = rownames(scrna) #names(which(PercentageExpressed(pchseq.resampled, features = rownames(pchseq.obj)) > 0))
  features.use = sort(unique(intersect(de$gene, expressed.genes)))
  
  # Integrated features
  # features.use = conf.obj@assays$integrated@var.features
  
  # Scale features
  scrna = ScaleData(scrna, features = features.use, assay = assay)
  pch.seq = ScaleData(pch.seq, features = features.use, assay = assay)
  
  # Train and test set
  object.train = GetTrainingSet(scrna, 
                                group.by = 'type', 
                                max.size = size, 
                                proportion = proportion, 
                                )
  
  object.test = subset(scrna, cells = Cells(object.train), invert = TRUE)
  
  # Train on integrated data
  system.time({xgb.model = XGBoost_train2(object.train, 
                                          object.test, 
                                          predict.by = 'type', 
                                          assay = assay, 
                                          ...)})
  # saveRDS(xgb.model, "../../Species_Reference/TranRGC_int_v2_xgb_int_30_cells.rds")
  
  system.time({cv.fit = TrainLogistic(scrna, 
                                      proportion = proportion,
                                      group.by = 'type', 
                                      size = size, 
                                      alpha = 1, 
                                      parallel = TRUE, 
                                      assay = assay)})
  # saveRDS(cv.fit, '../../Species_Reference/TranRGC_int_v2_glmnet_int_30_cells.rds')
    
  # conf.obj = readRDS('../../patchseq_data/seurat_objects/rgc-pch-sc-int-v1.rds')
  # xgb.model = readRDS("../../Species_Reference/TranRGC_int_v2_xgb_int_30_cells.rds")
  # cv.fit = readRDS('../../Species_Reference/TranRGC_int_v2_glmnet_int_30_cells.rds')
  
  # Predict logistic
  features.use = rownames(cv.fit$glmnet.fit$beta[[1]])
  # features.use = rownames(coef(cv.fit, s = 0.01))[-1]
  X_test = t(pch.seq@assays[[assay]]@scale.data)[,features.use]
  pch.seq$glm.class = predict(cv.fit, newx = X_test, s = 'lambda.min', type = "class")
  prob.matrix = as.data.frame(predict(cv.fit, newx = X_test, s = 'lambda.min', type = "response")[,,1])
  prob.matrix$type = ExtractString(rownames(prob.matrix), after = '\\.')
  avg.matrix.glm = prob.matrix %>% 
    group_by(type) %>% 
    summarise(across(.cols = where(is.numeric), .fns = mean)) %>% 
    column_to_rownames('type') %>%
    as.matrix()
  avg.matrix.glm = avg.matrix.glm[(unique(pch.seq$type)),]
  
  # Predict xgboost
  pch.seq = PredictLabels2(pch.seq, xgb.model, assay = assay, verbose = FALSE)
  prob.matrix = as.data.frame(t(PredictLabels2(pch.seq, xgb.model, assay = assay, verbose = FALSE, return.prob.matrix = TRUE)))
  prob.matrix$type = ExtractString(rownames(prob.matrix), after = '\\.')
  avg.matrix.xgb = prob.matrix %>% 
    group_by(type) %>% 
    summarise(across(.cols = where(is.numeric), .fns = mean)) %>% 
    column_to_rownames('type') %>%
    as.matrix()
  avg.matrix.xgb = avg.matrix.xgb[unique(pch.seq$type),]
  
  if(plot){
    print(
      ggHeatmap(avg.matrix.glm, title = 'Logistic (RGC type)', legend.name = 'Prob', xlab = 'Class\n', ylab = NULL, label = FALSE, border = T) + NoLegend() |
      ggHeatmap(avg.matrix.xgb, title = 'XGBoost (RGC type)', legend.name = 'Prob', xlab = 'Class\n', ylab = NULL, label = FALSE, border = T) + theme(axis.text.y = element_blank()) + NoLegend() |
      ggHeatmap((avg.matrix.glm*avg.matrix.xgb), title = 'Product', legend.name = 'Prob', xlab = 'Class\n', ylab = NULL, label = FALSE, border = T, max.value = 1) + theme(axis.text.y = element_blank())
    )
  }
  
  pchseq.obj = DownsampleSeurat(pch.seq, group.by = 'type', size = 1)
  
  pchseq.obj$xgb.anno.rgc.1 = apply(avg.matrix.xgb, 1, function(row) names(sort(row, decreasing = TRUE))[[1]])
  pchseq.obj$xgb.anno.rgc.2 = apply(avg.matrix.xgb, 1, function(row) names(sort(row, decreasing = TRUE))[[2]])
  pchseq.obj$xgb.anno.rgc.3 = apply(avg.matrix.xgb, 1, function(row) names(sort(row, decreasing = TRUE))[[3]])
  pchseq.obj$xgb.prob.rgc.1 = round(apply(avg.matrix.xgb, 1, function(row) (sort(row, decreasing = TRUE))[[1]]), 2)
  pchseq.obj$xgb.prob.rgc.2 = round(apply(avg.matrix.xgb, 1, function(row) (sort(row, decreasing = TRUE))[[2]]), 2)
  pchseq.obj$xgb.prob.rgc.3 = round(apply(avg.matrix.xgb, 1, function(row) (sort(row, decreasing = TRUE))[[3]]), 2)
  
  pchseq.obj$glm.anno.rgc.1 = apply(avg.matrix.glm, 1, function(row) names(sort(row, decreasing = TRUE))[[1]])
  pchseq.obj$glm.anno.rgc.2 = apply(avg.matrix.glm, 1, function(row) names(sort(row, decreasing = TRUE))[[2]])
  pchseq.obj$glm.anno.rgc.3 = apply(avg.matrix.glm, 1, function(row) names(sort(row, decreasing = TRUE))[[3]])
  pchseq.obj$glm.prob.rgc.1 = round(apply(avg.matrix.glm, 1, function(row) (sort(row, decreasing = TRUE))[[1]]), 2)
  pchseq.obj$glm.prob.rgc.2 = round(apply(avg.matrix.glm, 1, function(row) (sort(row, decreasing = TRUE))[[2]]), 2)
  pchseq.obj$glm.prob.rgc.3 = round(apply(avg.matrix.glm, 1, function(row) (sort(row, decreasing = TRUE))[[3]]), 2)
  
  pchseq.obj$goetz.assignment.clean = ExtractString(pchseq.obj$goetz.assignment, after = '_')
  pchseq.obj$xgb.anno.rgc.1.clean = gsub('C', '', ExtractString(pchseq.obj$xgb.anno.rgc.1, after = '_'))
  pchseq.obj$glm.anno.rgc.1.clean = gsub('C', '', ExtractString(pchseq.obj$glm.anno.rgc.1, after = '_'))
  pchseq.obj$xgb.anno.rgc.2.clean = gsub('C', '', ExtractString(pchseq.obj$xgb.anno.rgc.2, after = '_'))
  pchseq.obj$glm.anno.rgc.2.clean = gsub('C', '', ExtractString(pchseq.obj$glm.anno.rgc.2, after = '_'))
  pchseq.obj$xgb.anno.rgc.3.clean = gsub('C', '', ExtractString(pchseq.obj$xgb.anno.rgc.3, after = '_'))
  pchseq.obj$glm.anno.rgc.4.clean = gsub('C', '', ExtractString(pchseq.obj$glm.anno.rgc.3, after = '_'))
  
  knitr::kable(subset(pchseq.obj, high.conf == 1)@meta.data[,c('goetz.assignment.clean', 
                                                               'xgb.anno.rgc.1.clean', 
                                                               'xgb.prob.rgc.1', 
                                                               'glm.anno.rgc.1.clean', 
                                                               'glm.prob.rgc.1')
  ], align = "c")
  
  message('XGBoost classifications with high concordance: ', 
          mean(subset(pchseq.obj, high.conf == 1)@meta.data$high.conf.assignment == 
                 subset(pchseq.obj, high.conf == 1)@meta.data$xgb.anno.rgc.1.clean))
  message('Glmnet classifications with high concordance: ', 
          mean(subset(pchseq.obj, high.conf == 1)@meta.data$high.conf.assignment == 
                 subset(pchseq.obj, high.conf == 1)@meta.data$glm.anno.rgc.1.clean))
  
  # message('XGBoost classifications recovered from Goetz: ', 
  #         mean(subset(pchseq.obj@meta.data, !is.na(goetz.assignment))$goetz.assignment.clean == 
  #                subset(pchseq.obj@meta.data, !is.na(goetz.assignment))$xgb.anno.rgc.1.clean))
  # message('Glmnet classifications recovered from Goetz: ', 
  #         mean(subset(pchseq.obj@meta.data, !is.na(goetz.assignment))$goetz.assignment.clean == 
  #                subset(pchseq.obj@meta.data, !is.na(goetz.assignment))$glm.anno.rgc.1.clean))
  
  pchseq.obj
}

build_omim_gene_disease <- function(
    morbidmap,
    mim2gene = NULL,
    keep_only_genes = TRUE
) {
  stopifnot(file.exists(morbidmap))
  if (!is.null(mim2gene)) stopifnot(file.exists(mim2gene))
  
  # --- read morbidmap (gene–disease relationships) ---
  morbid <- read.delim(
    morbidmap,
    header = FALSE,
    comment.char = "#",
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  colnames(morbid) <- c(
    "disease",
    "gene_symbols",
    "disease_mim",
    "chromosome"
  )
  
  # expand multiple genes per disease
  gene_disease <- tidyr::separate_rows(
    morbid,
    gene_symbols,
    sep = ",\\s*"
  )
  colnames(gene_disease)[colnames(gene_disease) == "gene_symbols"] <- "gene"
  
  # --- optionally add gene identifiers from mim2gene ---
  if (!is.null(mim2gene)) {
    m2g <- read.delim(
      mim2gene,
      comment.char = "#",
      sep = '\t',
      stringsAsFactors = FALSE
    )
    
    if (keep_only_genes && "type" %in% colnames(m2g)) {
      m2g <- m2g[m2g$type == "gene", ]
    }
    
    gene_disease <- dplyr::left_join(
      gene_disease,
      m2g,
      by = c("gene" = "gene_symbol")
    )
  }
  
  return(gene_disease)
}


#' Interactively refine and annotate a Seurat cell class
#'
#' This function subsets a Seurat object to a specified cell class,
#' harmonizes batches, and provides an interactive workflow for
#' removing and merging clusters before renumbering and annotating them.
#'
#' During execution, the user is prompted to:
#' \itemize{
#'   \item Select clusters to remove
#'   \item Iteratively specify clusters to merge
#' }
#'
#' Diagnostic plots (DimPlot + DotPlot) are shown before and after
#' refinement to guide decisions.
#'
#' @param seurat_obj A \code{Seurat} object containing clustering results.
#' @param cell_class Character scalar specifying the cell class to subset
#'   (matched against \code{seurat_obj$cell_class}).
#' @param keep_extra_clusters Optional numeric vector of cluster IDs to
#'   retain even if they are not labeled as the specified cell class.
#'   Useful for rescuing borderline clusters.
#' @param batch Character scalar giving the metadata column used for
#'   batch harmonization (passed to \code{Harmonize}).
#'   Default is \code{"orig.file"}.
#' @param annotation_prefix Character prefix used when constructing the
#'   annotated cluster labels (e.g. \code{"BC"} → \code{"BC-0"}).
#'   Defaults to \code{cell_class}.
#' @param annotation_genes A gene annotation object or gene set passed to
#'   \code{Genes()} for DotPlot visualization.
#'
#' @return A refined \code{Seurat} object with:
#' \itemize{
#'   \item Updated clustering after merges
#'   \item Renumbered \code{seurat_clusters}
#'   \item A new \code{annotated} metadata column
#' }
#'
#' @details
#' Cluster removal is performed using \code{subset(..., invert = TRUE)}.
#' Cluster merging is handled by repeated calls to \code{MergeClusters}.
#' Final cluster IDs are renumbered contiguously using
#' \code{RenumberClusters}.
#'
#' This function is intended for exploratory, analyst-in-the-loop
#' cell type refinement rather than fully automated pipelines.
#'
#' @examples
#' \dontrun{
#' BC <- ClusterCellClass(
#'   initial,
#'   cell_class = "BC",
#'   keep_extra_clusters = 14,
#'   annotation_prefix = "BC",
#'   annotation_genes = major_annotation
#' )
#' }
#'
#' @seealso
#' \code{\link{subset}}, \code{\link{MergeClusters}},
#' \code{\link{RenumberClusters}}, \code{\link{DimPlot}}
#'
#' @export
ClusterCellClass <- function(
    seurat_obj,
    class,
    keep_extra_clusters = NULL,
    batch = "orig.file",
    use = c('harmony', 'seurat'),
    annotation_prefix = NULL,
    annotation_genes = major_annotation, 
    ...
) {
  
  message("▶ Subsetting cell class: ", class)
  
  obj <- subset(
    seurat_obj,
    cell_class %in% class |
      seurat_clusters %in% keep_extra_clusters
  )
  
  use = match.arg(use)
  if(use == 'harmony'){
    obj <- Harmonize(obj, batch = batch, ...)
  } else {
    obj <- ClusterSeurat(obj, integrate.by = batch, ...)
  }
  
  print(
    DimPlotLabeled(obj, group.by = "seurat_clusters") |
      DotPlot3(
        obj,
        features = Genes(annotation_genes),
        group.by = "seurat_clusters"
      )
  )
  
  ## -------------------------
  ## Interactive: remove clusters
  ## -------------------------
  remove <- readline(
    "Clusters to REMOVE (comma-separated, empty = none): "
  )
  
  if (nzchar(remove)) {
    remove <- as.numeric(strsplit(remove, ",")[[1]])
    obj <- subset(obj, seurat_clusters %in% remove, invert = TRUE)
  }
  
  ## -------------------------
  ## Interactive: merge clusters
  ## -------------------------
  repeat {
    merge <- readline(
      "Clusters to MERGE (e.g. 0,1 or empty to stop): "
    )
    
    if (!nzchar(merge)) break
    
    merge <- as.numeric(strsplit(merge, ",")[[1]])
    obj <- MergeClusters(obj, merge)
  }
  
  ## -------------------------
  ## Renumber + annotate
  ## -------------------------
  obj <- RenumberClusters(obj)
  
  if (is.null(annotation_prefix)) {
    annotation_prefix <- class
  }
  
  obj$annotated <- paste0(
    annotation_prefix,
    "-",
    as.numeric(obj$seurat_clusters)
  )
  
  print(
    DimPlotLabeled(obj, group.by = "annotated") |
      DotPlot3(
        obj,
        features = Genes(annotation_genes),
        group.by = "annotated"
      )
  )
  
  return(obj)
}


MakeBestMatchDF = function(BIGSEURAT, species_cluster = 'species_cluster', orthotype = 'type'){
  
  speciesList = as.character(unique(BIGSEURAT$species))
  
  message('Subsetting objects...')
  speciesACList = lapply(speciesList, function(currentSpecies) subset(BIGSEURAT, species == currentSpecies))
  names(speciesACList) = speciesList
  
  # Generate jaccard matrices
  message('Computing jaccard...')
  jaccardList = lapply(seq_along(speciesACList), function(index) JSMatrix(t(table(speciesACList[[index]]@meta.data[[species_cluster]], speciesACList[[index]]@meta.data[[orthotype]]))))
  names(jaccardList) = speciesList
  # saveRDS(jaccardList, '../../Ortho_Objects/jaccardList.rds')
  
  jaccardList = lapply(jaccardList, as.data.frame.matrix)
  
  # Get the best match in each species 
  message('Finding best matches...')
  bestMatch = lapply(jaccardList, function(matrix) apply(matrix, 1, function(x) colnames(matrix)[which.max(x)]))
  # reducedList = Reduce(function(x, y) merge(x, y, by = "row.names", all = T), jaccardList[1:3])
  bestMatchDf = as.data.frame(t(do.call(rbind, bestMatch)))
  colnames(bestMatchDf) = speciesList
  # bestMatchDf$OrthoType = rownames(bestMatchDf)
  # saveRDS(bestMatchDf, '../../Ortho_Objects/best-match.rds')
  
  message('Finding best jaccard value...')
  bestJaccard = lapply(jaccardList, function(matrix) apply(matrix, 1, max))
  # reducedList = Reduce(function(x, y) merge(x, y, by = "row.names", all = T), jaccardList[1:3])
  bestJaccardDf = as.data.frame(t(do.call(rbind, bestJaccard)))
  colnames(bestJaccardDf) = speciesList
  # bestJaccardDf$OrthoType = rownames(bestJaccardDf)
  # saveRDS(bestJaccardDf, '../../Ortho_Objects/best-jaccard.rds')
  
  list(bestMatchDf, bestJaccardDf)
}


RunLDA = function(x_train, x_test, class){
  
  df <- as.data.frame(x_train)
  df$class = class
  
  # -----------------------------------------------------------
  # 2. Fit Linear Discriminant Analysis
  # -----------------------------------------------------------
  lda_model <- lda(class ~ ., data = df)
  
  # Variance explained by each LDA axis
  prop_var <- lda_model$svd^2 / sum(lda_model$svd^2)
  round(prop_var, 3)
  
  # Project data into LDA space
  lda_result = predict(lda_model)
  EvaluateModel(table(lda_result$class, df$class))
  lda_proj <- lda_result$x
  lda_df <- data.frame(lda_proj, class = df$class)
  
  # -----------------------------------------------------------
  # 3. Plot LDA projection
  # -----------------------------------------------------------
  
  if('LD2' %in% colnames(lda_df)){
    p_lda <- ggplot(lda_df, aes(LD1, LD2, color = class)) +
      geom_point(size = 3, alpha = 0.8) +
      theme_minimal(base_size = 14) +
      labs(
        title = "LDA projection",
        x = paste0("LD1 (", round(prop_var[1] * 100, 1), "% var)"),
        y = paste0("LD2 (", round(prop_var[2] * 100, 1), "% var)")
      ) +
      scale_color_brewer(palette = "Dark2")
  } else {
    p_lda <- ggplot(lda_df, aes(LD1, color = class)) +
      geom_density(alpha = 0.5) +
      # geom_point(size = 3, alpha = 0.8) +
      theme_minimal(base_size = 14) +
      labs(
        title = "LDA projection",
        x = paste0("LD1 (", round(prop_var[1] * 100, 1), "% var)")
        # y = paste0("LD2 (", round(prop_var[2] * 100, 1), "% var)")
      ) +
      scale_color_brewer(palette = "Dark2")
  }
  
  
  # -----------------------------------------------------------
  # 4. Compare with PCA projection
  # -----------------------------------------------------------
  # pca_res <- prcomp(x_train, scale. = TRUE)
  pca_df <- data.frame(x_train, class = df$class)
  
  p_pca <- ggplot(pca_df, aes(PC_1, PC_2, color = class)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(
      title = "PCA projection",
      x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "% var)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "% var)")
    ) +
    scale_color_brewer(palette = "Dark2")
  
  # -----------------------------------------------------------
  # 5. Display side-by-side
  # -----------------------------------------------------------
  print(gridExtra::grid.arrange(p_pca, p_lda, ncol = 2))
  
  # Return prediction for test set
  lda_test_result = predict(lda_model, newdata = as.data.frame(x_test))
  
  list(lda_result, lda_test_result)
}

plot_metrics = function(data, add = c('mean_sd'), pr_bc = TRUE){
  # (PrettyBarplot(data, x = 'species', y = 'gaba.rgc', add = add, title = 'RGC-AC distance', ylab = 'Relative distance') |
  #  PrettyBarplot(data, x = 'species', y = 'gaba.rgc.norm.rgc.bc', add = add, title = 'RGC-AC distance / RGC-BC distance', ylab = 'Relative distance', jitter.params = list(alpha = 0.4, width = 0.2)) |
  #  PrettyBarplot(data, x = 'species', y = 'gaba.rgc.norm.gaba.bc', add = add, title = 'RGC-AC distance / AC-BC distance', ylab = 'Relative distance', jitter.params = list(alpha = 0.4, width = 0.2)))
  # PrettyBarplot(data, x = 'species', y = 'pr.bc', add = c('mean_sd', 'jitter'), title = 'PR-BC distance', ylab = 'Relative distance'))
  
  
  if(pr_bc){
    (PrettyBarplot(data, x = 'species', y = 'gaba.rgc.norm.gaba.bc', add = add, title = 'RGC-gabaAC / gabaAC-BC', ylab = 'Relative distance') |
       PrettyBarplot(data, x = 'species', y = 'gly.rgc.norm.gaba.bc', add = add, title = 'RGC-glyAC / gabaAC-BC', ylab = 'Relative distance') |
       PrettyBarplot(data, x = 'species', y = 'pr.bc.norm.gaba.bc', add = add, title = 'PR-BC / gabaAC-BC', ylab = 'Relative distance'))
  } else {
    (PrettyBarplot(data, x = 'species', y = 'gaba.rgc.norm.gaba.bc', add = add, title = 'RGC-gabaAC / gabaAC-BC', ylab = 'Relative distance') |
       PrettyBarplot(data, x = 'species', y = 'gly.rgc.norm.gaba.bc', add = add, title = 'RGC-glyAC / gabaAC-BC', ylab = 'Relative distance'))
    # PrettyBarplot(data, x = 'species', y = 'pr.bc.norm.gaba.bc', add = add, title = 'PR-BC / gabaAC-BC', ylab = 'Relative distance'))
  }
  
}

collect_metrics2 = function(avg.pc.list, FUN = `/`){
  
  do.call(rbind, lapply(names(avg.pc.list), function(species){
    # lapply(avg.pc.list[[species]], function(this.mean){
    this.mean = avg.pc.list[[species]]
    gaba.rgc.norm.rgc.bc = bootstrap_operator(subset(this.mean, sample1 == 'RGC' & sample2 == 'gabaAC')[['similarity']], 
                                              subset(this.mean, sample1 == 'RGC' & sample2 == 'BC')[['similarity']], 
                                              FUN)
    gaba.rgc.norm.gaba.bc = bootstrap_operator(subset(this.mean, sample1 == 'RGC' & sample2 == 'gabaAC')[['similarity']], 
                                               subset(this.mean, sample1 == 'gabaAC' & sample2 == 'BC')[['similarity']], 
                                               FUN)
    gly.rgc.norm.gaba.bc = bootstrap_operator(subset(this.mean, sample1 == 'RGC' & sample2 == 'glyAC')[['similarity']], 
                                              subset(this.mean, sample1 == 'gabaAC' & sample2 == 'BC')[['similarity']], 
                                              FUN)
    gly.rgc.norm.rgc.bc = bootstrap_operator(subset(this.mean, sample1 == 'RGC' & sample2 == 'glyAC')[['similarity']], 
                                             subset(this.mean, sample1 == 'RGC' & sample2 == 'BC')[['similarity']], 
                                             FUN)
    gaba.rgc.frac = bootstrap_operator(subset(this.mean, sample1 == 'RGC' & sample2 == 'gabaAC')[['similarity']], 
                                       c(0), 
                                       `+`)
    pr.bc = bootstrap_operator(subset(this.mean, sample1 == 'PR' & sample2 == 'BC')[['similarity']], 
                               c(0), 
                               `+`)
    pr.bc.norm.gaba.bc = bootstrap_operator(subset(this.mean, sample1 == 'PR' & sample2 == 'BC')[['similarity']], 
                                            subset(this.mean, sample1 == 'gabaAC' & sample2 == 'BC')[['similarity']], 
                                            FUN)
    
    data.frame(species = species, 
               gaba.rgc.norm.rgc.bc = gaba.rgc.norm.rgc.bc, 
               gaba.rgc.norm.gaba.bc = gaba.rgc.norm.gaba.bc, 
               gly.rgc.norm.gaba.bc = gly.rgc.norm.gaba.bc,
               gly.rgc.norm.rgc.bc = gly.rgc.norm.rgc.bc, 
               gaba.rgc = gaba.rgc.frac, 
               pr.bc = pr.bc, 
               pr.bc.norm.gaba.bc = pr.bc.norm.gaba.bc)
    # })
  }) #%>% unlist(recursive = FALSE)
  )
}

ComputeDistancesFromDistList = function(dist.list, ...){
  
  # Compute distances in PC space
  pc.mats.all = dapply(names(dist.list), function(species) {
    message(species)
    metadata = metadata.list[[species]]
    order = levels(objectList[[species]]$annotated)
    dist.mat = dist.list[[species]]
    dist.mat = dist.mat[order,order]
    SmartMatrix(dist.mat, 
                metadata[match(colnames(dist.mat), metadata$annotated),], 
                metadata[match(colnames(dist.mat), metadata$annotated),])
  })
  
  # Average across cell classes
  avg.pc.list = dapply(names(pc.mats.all), function(species){
    CollectSmartMatrix(pc.mats.all[[species]], 'cell_class2')
  })
  
  # Retrieve distances of interest
  metrics.list = collect_metrics2(avg.pc.list, ...) #lapply(avg.pc.list, collect_metrics2)
  
  # Return objects
  list(pc.mat = pc.mats.all, metrics = metrics.list)
}

ComputePCDistance = function(nPCs){
  
  # Retrieve binarized matrices
  bin.mats = lapply(speciesList.cp %>% setNames(speciesList.cp), function(species){
    message(species)
    basename = paste0('celltypephylo/', species, '_v33/')
    bin.mat.files = list.files(path = basename, pattern = "\\-bin.mat.rds$", include.dirs = TRUE)
    lapply(paste0(basename, bin.mat.files), function(x) {
      bin.mat = readRDS(x)
      order = gsub(' ', '-', levels(factor(objectList[[species]]$annotated)))
      bin.mat = bin.mat[,order]
    })
  })
  
  # Compute distances in PC space
  pc.mats.all = lapply(names(bin.mats), function(species) {
    message(species)
    lapply(bin.mats[[species]], function(bin.mat){
      metadata = metadata.list[[species]]
      dist.mat = PseudoBulkPCA(bin.mat, label = TRUE, verbose = FALSE, nPCs = nPCs, return.plot = FALSE)$dist %>% as.matrix()
      # dist.mat = dist.mat[order,order]
      SmartMatrix(dist.mat, 
                  metadata[match(colnames(dist.mat), metadata$annotated),], 
                  metadata[match(colnames(dist.mat), metadata$annotated),])
    })
  }) %>% setNames(speciesList.cp)
  
  # Average across cell classes
  avg.pc.list = lapply(names(pc.mats.all), function(species){
    lapply(pc.mats.all[[species]], function(pc.mat){
      AverageSmartMatrix(pc.mat, 'cell_class2')
    })
  }) %>% setNames(speciesList.cp)
  
  # Retrieve distances of interest
  do.call(rbind, lapply(names(avg.pc.list), function(species){
    lapply(avg.pc.list[[species]], function(this.mean){
      gaba.rgc.frac = subset(this.mean, sample1 == 'RGC' & sample2 == 'gabaAC')[['mean_similarity']]
      pr.bc = subset(this.mean, sample1 == 'PR' & sample2 == 'BC')[['mean_similarity']]
      gaba.rgc.norm.rgc.bc = subset(this.mean, sample1 == 'RGC' & sample2 == 'gabaAC')[['mean_similarity']]/
        subset(this.mean, sample1 == 'RGC' & sample2 == 'BC')[['mean_similarity']]
      gaba.rgc.norm.gaba.bc = subset(this.mean, sample1 == 'RGC' & sample2 == 'gabaAC')[['mean_similarity']]/
        subset(this.mean, sample1 == 'gabaAC' & sample2 == 'BC')[['mean_similarity']]
      gly.rgc.norm.gaba.bc = subset(this.mean, sample1 == 'RGC' & sample2 == 'glyAC')[['mean_similarity']]/
        subset(this.mean, sample1 == 'gabaAC' & sample2 == 'BC')[['mean_similarity']]
      pr.bc.norm.gaba.bc = subset(this.mean, sample1 == 'PR' & sample2 == 'BC')[['mean_similarity']]/
        subset(this.mean, sample1 == 'gabaAC' & sample2 == 'BC')[['mean_similarity']]
      
      data.frame(species = species, 
                 gaba.rgc.norm.rgc.bc = gaba.rgc.norm.rgc.bc, 
                 gaba.rgc.norm.gaba.bc = gaba.rgc.norm.gaba.bc, 
                 gly.rgc.norm.gaba.bc = gly.rgc.norm.gaba.bc,
                 gaba.rgc = gaba.rgc.frac, 
                 pr.bc = pr.bc, 
                 pr.bc.norm.gaba.bc = pr.bc.norm.gaba.bc)
    })
  }) %>% unlist(recursive = FALSE))
  
  # browser()
  # collect_metrics2(avg.pc.list)
  
}

ComputePCDistance2 = function(objectList, metadata.list, nPCs = 2:5, show_example = NULL, downsample = TRUE, seed = 2, size = 100){
  
  # Compute distances in PC space
  pc.mats.all = dapply(speciesList.cp, function(species) {
    message(species)
    dapply(seq_len(n_iter), function(i){
      metadata = metadata.list[[species]]
      if(downsample) { 
        object = DownsampleSeurat(objectList[[species]], 
                                  group.by = 'annotated', 
                                  size = size, 
                                  seed = seed) 
      } else {
        object = objectList[[species]]
      }
      pca.res = PseudoBulkPCA(object, 
                              group.by = 'annotated', 
                              nPCs = nPCs, 
                              return.dist = TRUE, 
                              return.plot = FALSE)
      dist.mats = pca.res$dist %>% setNames(nPCs)
      lapply(dist.mats, function(dist.mat){
        order = levels(objectList[[species]]$annotated)
        dist.mat = dist.mat[order,order]
        SmartMatrix(dist.mat, 
                    metadata[match(colnames(dist.mat), metadata$annotated),], 
                    metadata[match(colnames(dist.mat), metadata$annotated),], 
                    misc = list(var_exp = pca.res$var_exp))
      })
    })
  })
  
  # Reorder such that seeds are the bottom rung
  pc.mats.all = lapply(pc.mats.all, function(species.mats) {
    dapply(nPCs, function(nPC){
      
      seed.list = lapply(species.mats, function(this.iter) {
        this.iter[[as.character(nPC)]]
      })
      
      # browser()
      
      # Average across seeds
      elementwise_mean(seed.list)
    })
  })
  
  
  if(!is.null(show_example)){
    print(SmartHeatmap(pc.mats.all[[show_example]]$`2`, 
                       annotate.by = 'cell_class2', 
                       colors = list('cell_class2' = cell_class2_colors),
                       title = show_example))
    print(SmartHeatmap(pc.mats.all[[show_example]]$`3`, 
                       annotate.by = 'cell_class2', 
                       colors = list('cell_class2' = cell_class2_colors),
                       title = show_example))
    print(SmartHeatmap(pc.mats.all[[show_example]]$`4`, 
                       annotate.by = 'cell_class2', 
                       colors = list('cell_class2' = cell_class2_colors),
                       title = show_example))
  }
  
  # Average across cell classes
  avg.pc.list = dapply(names(pc.mats.all), function(species){
    print(species)
    lapply(pc.mats.all[[species]], function(pc.mat){
      # lapply(n_iter, function(i){
      # lapply(nPCs, function(nPC){
      # pc.mat[[nPC]]
      # if(species == 'ZebrafishLyu') browser()
      CollectSmartMatrix(pc.mat, 'cell_class2')
      # AverageSmartMatrix(pc.mat, 'cell_class2')
      # })
      # })
    })
  })
  
  # Reorder by nPCs at top layer
  npc.list = dapply(nPCs, function(nPC) {
    lapply(avg.pc.list, function(avg.pc){
      avg.pc[[as.character(nPC)]]
    })
  })
  
  # Retrieve distances of interest
  metrics.list = lapply(npc.list, collect_metrics2)
  
  # Return objects
  list(pc.mat = pc.mats.all, metrics = metrics.list)
}

ReadSAMapObjects = function(new.lamprey = FALSE, 
                            nonmammals = TRUE, 
                            mammals = FALSE, 
                            mouse = FALSE, 
                            lyu = FALSE, 
                            new = FALSE, 
                            metadata = FALSE, 
                            lighten = FALSE){
  
  objectList = list()
  
  # Read in AC object for matching cells with AC object
  message('Reading AC object...')
  vertebrateAC = readRDS('../../Ortho_Objects/vertebrateAC_unintegrated.rds')
  
  # Annotate on and OFF BCs
  message('Reading full objects...')
  
  if(metadata) {
    FUN = function(object){
      readRDS(object)@meta.data
    }
    OnOffBcClassification = function(object, ...){
      object
    }
    MatchACs = function(object, ...){
      object
    }
  } else {
    FUN = readRDS
  }
  
  if(nonmammals){
    message('Reading non-mammal objects...')
    if(new.lamprey) objectList$Lamprey = FUN('../../Full_Objects/Lamprey_Wang_full_v2.rds') # Already processed through on off BC functions
    if(!new.lamprey) objectList$Lamprey = OnOffBcClassification(FUN('../../Full_Objects/Lamprey_Wang_full_v1.rds'), 'LOC116937345', plot = FALSE, cutoff = 20)
    objectList$Shark = OnOffBcClassification(FUN('../../Full_Objects/Shark_full_v5.rds'), 'LOC119962784', cutoff = 20)
    objectList$Killifish = OnOffBcClassification(FUN('../../Species_Objects/Killifish_ncbi_initial_v5.rds'), 'isl1a', cutoff = 5)
    objectList$Goldfish = OnOffBcClassification(FUN('../../Full_Objects/Goldfish_full_v2.rds'), 'isl1', cutoff = 5)
    if(!lyu) objectList$Zebrafish = OnOffBcClassification(FUN('../../Full_Objects/Zebrafish_full_v3.rds'), 'isl1', cutoff = 5)
    if(lyu) objectList$ZebrafishLyu = OnOffBcClassification(FUN('../../Full_Objects/Zebrafish_lyu_full_v2.rds'), 'isl1', cutoff = 5)
    objectList$Newt = OnOffBcClassification(FUN('../../Full_Objects/Newt_full_v3.rds'), 'ISL1', plot = FALSE, cutoff = 5)
    objectList$Axolotl = OnOffBcClassification(FUN('../../Full_Objects/Axolotl_full_v3.rds'), 'ISL1', plot = FALSE, cutoff = 5)
    objectList$Chicken = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Chicken_full_v3.rds')))
    objectList$Lizard = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Lizard_full_v7.rds')), cutoff = 10)
  }
  
  if(lighten){
    message('Lightening objects by removing scaled data')
    objectList = lapply(objectList, function(object) {
      if(is(object@assays$RNA@scale.data, "matrix")){
        object@assays$RNA@scale.data = matrix()
      }
      
      return(object)
    }) 
  }
    
  if(mammals){
    message('Reading mammal objects...')
    objectList$Opossum = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Opossum_full_v2.rds')), cutoff = 10) 
    objectList$Sheep = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Sheep_full_v1.rds')), cutoff = 10) 
    objectList$Squirrel = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Squirrel_full_v2.rds')), cutoff = 10)
    objectList$Rat = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Rat_full_v2.rds')), cutoff = 20) 
    if(mouse) objectList$Mouse = OnOffBcClassification((FUN('../../Full_Objects/Mouse_full_v3.rds')), cutoff = 10) # Don't match ACs as this will fail for mouse do to barcode mismatch
    objectList$TreeShrew = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/TreeShrew_full_v2.rds')), cutoff = 10) 
    objectList$MouseLemur = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/MouseLemur_full_v2.rds')), cutoff = 10) 
    objectList$Marmoset = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Marmoset_full_v2.rds')), cutoff = 10) 
    objectList$Macaque = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Macaque_full_v2.rds')), cutoff = 10) 
    objectList$Human = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Human_full_v2.rds')), cutoff = 10) 
  }
  
  if(new){
    objectList$Axolotl = OnOffBcClassification(FUN('../../Full_Objects/Axolotl_full_v3.rds'), 'ISL1', plot = FALSE, cutoff = 5)
    objectList$Newt = OnOffBcClassification(FUN('../../Full_Objects/Newt_full_v3.rds'), 'ISL1', plot = FALSE, cutoff = 5)
    objectList$Rat = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/Rat_full_v2.rds')), cutoff = 20) 
    objectList$MouseLemur = OnOffBcClassification(MatchACs(FUN('../../Full_Objects/MouseLemur_full_v2.rds')), cutoff = 10) 
  }
  
  if(lighten){
    message('Lightening objects by removing scaled data')
    objectList = lapply(objectList, function(object) {
      if(is(object@assays$RNA@scale.data, "matrix")){
        object@assays$RNA@scale.data = matrix()
      }
      
      return(object)
    })
  }
  
  message('Finished reading objects...')
  
  if(!metadata){
    message('Adding factor levels to annotated...')
    objectList = lapply(objectList, OrderAnnotated)
  }
  
  objectList
}

OrderAnnotated = function(object){
  object$annotated = factor(object$annotated,
                            levels = (Metadata(object, 'annotated', 'cell_class2') %>%
                                        arrange(factor(cell_class2, levels = unique(Annotation(major_annotation_pr_gaba_gly_on_off)))))$annotated)
  # object$annotated = gsub(' ', '-', object$annotated)
  Idents(object) = 'annotated'
  object
}

ReadFullObjects = function(){
  
    nonmammals = data.frame(species = c("Lamprey", 'Zebrafish', "Chicken", "Lizard"), 
                            path = c("Lamprey_kPetMar1_initial_v2", 'Zebrafish_full_v2', "Chicken_full_v1", "Lizard_full_v1"))
    
    # Read in objects
    objectList = lapply(nonmammals$path, function(path) readRDS(paste0('../../Full_Objects/', path, '.rds')))
    names(objectList) = nonmammals$species
    
    # Add species ident
    objectList = lapply(names(objectList), function(species) {
      objectList[[species]]$species = species
      return(objectList[[species]])
    })
    names(objectList) = nonmammals$species
    
    # Lamprey cell class
    objectList$Lamprey$cell_class = objectList$Lamprey$peng_annotation
    
    # Make lamprey AC clusters glycinergic
    objectList$Lamprey$annotated = objectList$Lamprey$peng_clusters
    objectList$Lamprey$annotated[objectList$Lamprey$cell_class == 'AC'] = paste0('gly', objectList$Lamprey$peng_clusters[objectList$Lamprey$cell_class == 'AC'])
    
    # Make lamprey SACs as AC-0
    objectList$Lamprey$cell_class[objectList$Lamprey$peng_clusters == 'SAC'] = 'AC'
    objectList$Lamprey$annotated[objectList$Lamprey$peng_clusters == 'SAC'] = 'gabaAC-0'
    objectList$Lamprey$annotated = factor(objectList$Lamprey$annotated)
    
    # Coerce dual/ngng classifications into GABA/gly
    objectList$Chicken$classification[objectList$Chicken$annotated %in% c('AC-23', 'AC-62')] = 'GABA'
    objectList$Lizard$classification[objectList$Lizard$annotated %in% c('AC-16', 'AC-11')] = 'GABA'
    objectList$Zebrafish$classification[objectList$Zebrafish$annotated %in% c('AC-15_SEG', 'AC-9')] = 'Gly'
    
    # Assign ACs their subclass
    objectList[c('Zebrafish', 'Chicken', 'Lizard')] = lapply(objectList[c('Zebrafish', 'Chicken', 'Lizard')], function(object) {
      object$annotated = as.character(object$annotated)
      object$annotated[object$cell_class == 'AC'] = paste0(gsub('both', 'dual', tolower(object$classification[object$cell_class == 'AC'])), object$annotated[object$cell_class == 'AC'])
      return(object)
    })
    
    pm2ze = readRDS('../cones/orthology_graphs/pm2ze.rds')
    pm2ch = readRDS('../cones/orthology_graphs/pm2ch.rds')
    pm2li = readRDS('../cones/orthology_graphs/pm2li.rds')
    ze2li = readRDS('../cones/orthology_graphs/ze2li.rds')
    ze2ch= readRDS('../cones/orthology_graphs/ze2ch.rds')
    li2ch = readRDS('../cones/orthology_graphs/li2ch.rds')
    
    # Collapse chicken HCs and MGs
    objectList$Chicken$annotated[objectList$Chicken$cell_class == 'HC'] = 'HC'
    objectList$Chicken$annotated[objectList$Chicken$cell_class == 'MG'] = 'MG'
    
    # Make annotated the ident
    objectList = lapply(objectList, function(object) {
      
      # Add cell class that distinguishes gaba and gly AC
      object$cell_class2 = ExtractString(object$annotated, after = '-')
      
      # Sort annotated
      object$annotated = factor(object$annotated, levels = Metadata(object, 'cell_class2', 'annotated')$annotated)
      
      Idents(object) = 'annotated'
      
      return(object)
    })
    
    lapply(objectList, function(object) table(object$annotated))
    
    # Normalize zebrafish data (which is a chimeric object)
    objectList$Zebrafish = NormalizeData(objectList$Zebrafish)
    
    # Zebrafish all from Lyu
    objectList$ZebrafishLyu = readRDS('../../Full_Objects/Zebrafish_lyu_full_v1.rds')
    objectList$ZebrafishLyu$cell_class2 = objectList$ZebrafishLyu$cell_class
    objectList$ZebrafishLyu$cell_class2[objectList$ZebrafishLyu$cell_class == 'AC'] = paste0(tolower(objectList$ZebrafishLyu$classification[objectList$ZebrafishLyu$cell_class == 'AC']), 'AC')
    objectList$ZebrafishLyu = AddGabaGly(objectList$ZebrafishLyu)
    
    # Control for PR subclustering
    objectList$ZebrafishLyuPR = objectList$ZebrafishLyu
    objectList$ZebrafishLyuPR$annotated[objectList$ZebrafishLyuPR$cell_class2 == 'PR'] = 'PR-1'
    
    # Add cat shark
    objectList$Shark = readRDS('../../Full_Objects/Shark_full_v3.rds')
    objectList$Shark = subset(objectList$Shark, annotated %in% c('PR-Precursor', 'AC-15'), invert = TRUE)
    objectList$Shark = subset(objectList$Shark, annotated %in% c('RGC-1', 'PR-Rod2'), invert = TRUE)
    objectList$Shark$annotated[objectList$Shark$cell_class2 == 'gabaAC'] = paste0('gaba', objectList$Shark$annotated[objectList$Shark$cell_class2 == 'gabaAC'])
    objectList$Shark$annotated[objectList$Shark$cell_class2 == 'glyAC'] = paste0('gly', objectList$Shark$annotated[objectList$Shark$cell_class2 == 'glyAC'])
    objectList$Shark$annotated[objectList$Shark$annotated == 'AC-10'] = paste0('gly', objectList$Shark$annotated[objectList$Shark$annotated == 'AC-10'])
    objectList$Shark$cell_class2 = ExtractString(objectList$Shark$annotated, after = '-')
    # saveRDS(objectList$Shark, '../../Full_Objects/Shark_full_v4.rds')
    
    # Killifish v3
    objectList$Killifish = readRDS('../../Full_Objects/Killifish_full_v3.rds')
    objectList$Killifish = subset(objectList$Killifish, annotated %in% c('AC-21', 'AC-9'), invert = TRUE) # AC-9 is the BHLHE22+, we will keep it moving forward
    objectList$Killifish$annotated[objectList$Killifish$cell_class2 == 'gabaAC'] = paste0('gaba', objectList$Killifish$annotated[objectList$Killifish$cell_class2 == 'gabaAC'])
    objectList$Killifish$annotated[objectList$Killifish$cell_class2 == 'glyAC'] = paste0('gly', objectList$Killifish$annotated[objectList$Killifish$cell_class2 == 'glyAC'])

    # Killifish ncbi v5
    # objectList$KillifishNCBI = readRDS('../../Full_Objects/Killifish_ncbi_initial_v5.rds')
    
    # Goldfish
    objectList$Goldfish = readRDS('../../Full_Objects/Goldfish_full_v2.rds')
    objectList$Goldfish = subset(objectList$Goldfish, annotated %in% c('AC-1', 'AC-3', 'RGC-10', 'RGC-3'), invert = TRUE)
    
    # Axolotl
    objectList$Axolotl = readRDS('../../Full_Objects/Axolotl_full_v2.rds')
    
    # Opossum
    op.full = readRDS('../../Full_Objects/Opossum_full_v1.rds')
    op.full$cell_class2 = as.character(op.full$cell_class)
    op.full$cell_class2[op.full$classification == 'Gly'] = 'glyAC'
    op.full$cell_class2[op.full$classification != 'Gly'] = 'gabaAC'
    objectList$Opossum = op.full
    rm(op.full)
    
    # Rat
    rat.full = readRDS('../../Full_Objects/Rat_full.rds')
    rat.full$cell_class = rat.full$cell_class2
    rat.full$cell_class2 = as.character(rat.full$cell_class)
    rat.full$cell_class2[rat.full$classification == 'GABA'] = 'gabaAC'
    rat.full$cell_class2[rat.full$classification != 'GABA'] = 'glyAC'
    rat.full$cell_class2[rat.full$annotated == 'AC-13'] = 'gabaAC'
    rat.full$cell_class2[rat.full$cell_class2 %in% c('Rod', 'Cone')] = 'PR'
    rat.full$cell_class = replace_factor(rat.full$cell_class, c('Rod', 'Cone'), 'PR')
    # saveRDS(rat.full, '../../Full_Objects/Rat_full_v1.rds')
    objectList$Rat = rat.full
    rm(rat.full)
    
    # Squirrel
    objectList$Squirrel = readRDS('../../Full_Objects/Squirrel_full_v1.rds')
    objectList$Squirrel$cell_class2 = as.character(objectList$Squirrel$cell_class)
    objectList$Squirrel$cell_class2[objectList$Squirrel$classification == 'Gly'] = 'glyAC'
    objectList$Squirrel$cell_class2[objectList$Squirrel$classification != 'Gly'] = 'gabaAC'
    
    # Tree Shrew
    objectList$TreeShrew = readRDS('../../Full_Objects/TreeShrew_full_v1.rds')
    objectList$TreeShrew = subset(objectList$TreeShrew, annotated %in% c('BC-13', 'BC-14', 'BC-10'), invert = TRUE)
    
    # Remove LINC and/or antisense transcripts
    # objectList$TreeShrew = SubsetSeuratGenes(objectList$TreeShrew, 
    #                                          features = rownames(objectList$TreeShrew)[!(startsWith(rownames(objectList$TreeShrew), 'LINC') | 
    #                                                                                             grepl('-AS-', rownames(objectList$TreeShrew)))])
    
    # Mouse lemur
    objectList$MouseLemur = readRDS('../../Full_Objects/MouseLemur_full_v1.rds')
    objectList$MouseLemur$cell_class2 = as.character(objectList$MouseLemur$cell_class)
    objectList$MouseLemur$cell_class2[objectList$MouseLemur$classification == 'Gly'] = 'glyAC'
    objectList$MouseLemur$cell_class2[objectList$MouseLemur$classification != 'Gly'] = 'gabaAC'
    objectList$MouseLemur$cell_class2[objectList$MouseLemur$cell_class2 %in% c('Rod', 'Cone')] = 'PR'
    objectList$MouseLemur$cell_class[objectList$MouseLemur$cell_class %in% c('Rod', 'Cone')] = 'PR'
    
    # Marmoset
    objectList$Marmoset = readRDS('../../Full_Objects/Marmoset_full_v1.rds')
    objectList$Marmoset$cell_class2 = as.character(objectList$Marmoset$cell_class)
    objectList$Marmoset$cell_class2[objectList$Marmoset$classification == 'Gly'] = 'glyAC'
    objectList$Marmoset$cell_class2[objectList$Marmoset$classification != 'Gly'] = 'gabaAC'
    
    # Macaque
    objectList$Macaque = readRDS('../../Full_Objects/Macaque_full_v2.rds')
    
    # Human
    human.full = readRDS('../../Full_Objects/Human_full_v1.rds')
    human.full$cell_class2 = human.full$cell_class
    human.full$cell_class2[human.full$classification == 'GABA'] = 'gabaAC'
    human.full$cell_class2[human.full$classification == 'Gly'] = 'glyAC'
    objectList$Human = human.full
    rm(human.full)
    
    # Clean up lamprey
    lamprey = readRDS('../../Full_Objects/Lamprey_kPetMar1_initial_v3_converted.rds')
    # lamprey$annotated2 = lamprey$annotated
    # lamprey$annotated = lamprey$code
    lamprey = subset(lamprey, annotated %in% c('BC-4', paste0('RGC-', c(28, 32, 33, 34, 35, 37))), invert = TRUE)
    lamprey$cell_class2[lamprey$cell_class == 'RGC'] = 'gabaAC'
    lamprey$cell_class2[lamprey$annotated == 'RGC-15'] = 'RGC'
    lamprey$cell_class2[lamprey$annotated == 'RGC-17'] = 'glyAC'
    objectList$Lamprey = lamprey
    rm(lamprey)
    # lamprey = readRDS('../../Full_Objects/Lamprey_Wang_full_v1.rds')
    
    objectList = lapply(objectList, function(object) {
      object$annotated = factor(object$annotated,
                                levels = (Metadata(object, 'annotated', 'cell_class2') %>%
                                  arrange(factor(cell_class2, levels = unique(Annotation(major_annotation_pr_gaba_gly)))))$annotated)
      object$annotated = gsub(' ', '-', object$annotated)
      Idents(object) = 'annotated'
      object
    })
    
    objectList = objectList[speciesList.cp]
    objectList
}

SearchContrainedTrees = function(object, orig.dir, options = ''){
  
  message('Searching all possible 5-tip trees...')
  all_possible_trees = allTrees(5, rooted = FALSE)
  
  metadata = Metadata(object, 'annotated', 'cell_class2', 'lit_type')
  PR = paste0('<', paste0(metadata$annotated[metadata$cell_class2 == 'PR'], collapse = '|'), '>')
  BC = paste0('<', paste0(metadata$annotated[metadata$cell_class2 == 'BC'], collapse = '|'), '>')
  glyAC = paste0('<', paste0(metadata$annotated[metadata$cell_class2 == 'glyAC'], collapse = '|'), '>')
  gabaAC = paste0('<', paste0(metadata$annotated[metadata$cell_class2 == 'gabaAC'], collapse = '|'), '>')
  # gabaAC = paste0(metadata$annotated[which(metadata$cell_class2 == 'gabaAC' & (metadata$lit_type != 'SAC' | is.na(metadata$lit_type)))], collapse = '|')
  # gabaAC = gsub('\\|AC-17', '', gabaAC)
  RGC = paste0('<', paste0(metadata$annotated[metadata$cell_class2 == 'RGC'], collapse = '|'), '>')
  HC = paste0('<', paste0(metadata$annotated[metadata$cell_class2 == 'HC'], collapse = '|'), '>')
  MG = paste0('<', paste0(metadata$annotated[metadata$cell_class2 == 'MG'], collapse = '|'), '>')
  
  # Rename the tip labels in each tree
  skeletons <- lapply(all_possible_trees, function(tr) {
    tr$tip.label <- c('PR', 'BC', 'glyAC', 'gabaAC', 'RGC')
    tr
  })
  
  trees_renamed <- lapply(all_possible_trees, function(tr) {
    tr$tip.label <- c(PR, BC, glyAC, gabaAC, RGC) #, HC, MG)
    tr
  })
  
  # Change class
  class(trees_renamed) = 'multiPhylo'
  
  # Check first tree
  # lapply(trees_renamed, plot)
  
  # Save each topology to separate file
  lapply(seq_along(trees_renamed), function(i){
    
    filename = paste0('constraint_tests/', orig.dir, '/tip5-topologies-', i, '.tre')
    write.tree(trees_renamed[[i]], filename)
    
    # Replace pipe with comma (separator between types)
    system(paste0("sed -i '' 's/|/,/g' ", filename))
    
    # Replace < with ( (separator between classes)
    system(paste0("sed -i '' 's/<//g' ", filename))
    
    # Replace < with ( (separator between classes)
    system(paste0("sed -i '' 's/>//g' ", filename))
    
  })
  
  ### Test constrained topologies
  
  # Loop over all the permutations for this species
  
  # Loop over constrained trees and run iqtree for each constraint
  original.file = paste0('celltypephylo/', orig.dir, '/', orig.dir, '-30-1e-10--none-k-2000-1')
  bin.mat = readRDS(paste0(original.file, '-bin-mat.rds'))
  ml.treefile = paste0(original.file, '.treefile')
  message('Analyzing file: ', original.file)
  
  # Get iqtree model
  iqtree.model = ExtractString(system(paste0("grep 'Best-fit model according to BIC: ' ", original.file, ".iqtree"),
                                      intern = TRUE),
                               before = 'BIC: ')
  # iqtree.model = 'MFP'
  message('Using model: ', iqtree.model)
  
  contrained.res = lapply(seq_along(trees_renamed), function(i){
    
    message('Trying topology ', i, '!')
    
    # Get data
    filename = paste0('constraint_tests/', orig.dir, '/tip5-topologies-', i, '.tre')
    constraint.tree = read.tree(filename)
    phydat = BinaryPhyDat(t(bin.mat))
    output.file = filename
    
    # Check names match
    all(constraint.tree$tip.label %in% colnames(bin.mat))
    setdiff(constraint.tree$tip.label, colnames(bin.mat))
    setdiff(colnames(bin.mat), constraint.tree$tip.label)
    
    # Run iqtree
    system.time({tree = RunIQTree2(phydat, 
                                   prefix = output.file, 
                                   iqtree.model = iqtree.model,
                                   constraint = filename,
                                   options = options
    )})
    
    # Get log-likelihood
    log.lik = as.numeric(ExtractString(system(paste0("grep 'Log-likelihood of the tree' ", output.file, ".iqtree"), 
                                              intern = TRUE), 
                                       before = 'tree: ', 
                                       after = ' \\('))
    
    list(tree = tree, log.lik = log.lik, output = output.file)
  })
  
  # Run tree topology tests
  all.trees.file = paste0('constraint_tests/', orig.dir, '/tip5-topologies.trees')
  tree.files = paste0(c(paste0(sapply(contrained.res, function(x) x$output), '.treefile'), 
                        ml.treefile), # the original tree too
                      collapse = ' ')
  system(paste0('cat ', tree.files, ' > ', all.trees.file))
  
  # Compare
  test.call = paste0('iqtree2 -s ', original.file, '.phy -z ', all.trees.file, ' -m ', iqtree.model, ' -au -zb 10000 -redo', 
                     ' -te ', ml.treefile,
                     ' -pre ', all.trees.file)
  message('Running call: ', test.call)
  system(test.call)
  
  contrained.res
}

PlotMeanHistograms = function(version, by = 5){
  
  # Load libraries
  library(magick)
  library(cowplot)
  
  # Set path and pattern
  v12_dirs <- list.dirs(path = "celltypephylo", recursive = TRUE, full.names = TRUE)
  v12_dirs <- v12_dirs[grepl(paste0(version, "$"), v12_dirs)]
  
  files_in_v12_dirs <- unlist(lapply(v12_dirs, function(dir) {
    list.files(path = dir, pattern = "\\-mean-cv.pdf$", full.names = TRUE)
  }))
  
  # Extract page 2 from each PDF
  images <- lapply(files_in_v12_dirs, function(f) {
    img <- image_read_pdf(f, density = 150)
    if (length(img) >= 2) img[2] else NULL  # skip if no second page
  })
  
  # Remove NULLs (PDFs without a second page)
  images <- Filter(Negate(is.null), images)
  images = FlipOrder(images, by)
  titles = sapply(seq_along(files_in_v12_dirs), function(i) ExtractString(basename(files_in_v12_dirs[[i]]), after = '-'))
  titles = FlipOrder(titles, by)
  
  # Wrap images into ggdraw plots
  plots <- lapply(seq_along(images), function(i) {
    img = images[[i]]
    ggdraw() + 
      draw_image(img)+
      draw_label(label = titles[[i]], y = 1, vjust = 1.2, fontface = "bold", size = 10)
  })
  
  # Combine and display
  plot_grid(plotlist = plots, nrow = by)  # change nrow/ncol if needed
  
  # Optional: save the plot
  # ggsave("combined_page2_plot.pdf", combined, width = 4 * length(plots), height = 4)
}

asrAnalaysis = function(object, path){
  bin.mat = readRDS(paste0(path, '-bin-mat.rds'))
  # genes = readRDS(paste0(path, '-genes.rds'))$genes.use
  asr = read.table(paste0(path, '.state'), header=TRUE)
  tree = read.tree(paste0(path, '.treefile'))
  tree$node.label = gsub('Node', '', tree$node.label)
  
  asr_matrix <- asr %>%
    mutate(State = ifelse(State == "-", 0.5, State)) %>% 
    mutate(State = as.numeric(State)) %>% 
    dplyr::select(Node, Site, State) %>%
    pivot_wider(names_from = Site, values_from = State) %>%
    column_to_rownames('Node')
  # arrange(Node)  # Optional: sort by node name
  
  rownames(asr_matrix) = gsub('Node', '', rownames(asr_matrix))
  colnames(asr_matrix) = rownames(bin.mat)
  
  # merge with tip data
  full_matrix = rbind(t(bin.mat), asr_matrix)
  full_matrix[1:10,1:10]
  node.list = c(tree$tip.label, tree$node.label)
  full_matrix = full_matrix[node.list,]
  
  # # Find the edge leading to PR
  # tip.of.interest = which(node.list == 'PR')
  # node.of.interest = tree$edge[tree$edge[,2] == tip.of.interest]
  # node.of.interest
  # 
  # # Difference matrix
  # row = as.data.frame(t(full_matrix[node.list[node.of.interest],]))
  # row$diff = row[,1] - row[,2]
  # row
  
  # diff.gene = c(rownames(row)[row$diff > 0], rownames(row)[row$diff < 0])
  # DotPlot3(objectList$Chicken, features = colnames(row)[row > 0], group.by = 'annotated')
  
  plot(tree)
  nodelabels()
  edgelabels(round(tree$edge.length, 2))
  
  # plot_tree(tree, Metadata(objectList$Chicken, 'annotated', 'cell_class', 'cell_class2'), max.bl = 1.5)
  tree.dat = plot_tree(tree, Metadata(object, 'annotated', 'cell_class', 'cell_class2'), return.df = TRUE)
  tree.dat
  
  tree.dat$from = node.list[tree.dat$parent]
  tree.dat$to = node.list[tree.dat$node]
  
  gene.switches = apply(tree.dat, 1, function(row){
    diff =  full_matrix[row[['to']],] - full_matrix[row[['from']],]
    gained = sort(colnames(diff)[which(diff == 1)])
    lost = sort(colnames(diff)[which(diff == -1)])
    list(gained = gained, lost = lost)
  })
  
  parent_node <- getMRCA(tree, tip = subset(tree.dat, cell_class2 == 'RGC')$label)
  switch.rgc = gene.switches[[which(tree.dat$node == parent_node)]]
  message('Found ', length(do.call(c, switch.rgc)), ' switches at the RGC MRCA!')
  
  # DotPlot3(objectList$Chicken, features = c('RBPMS'), group.by = 'annotated')
  DotPlot3(object, 
           features = do.call(c, switch.rgc), 
           group.by = 'annotated', 
           binarization = bin.mat[do.call(c, switch.rgc),])
}

SNFBinarize = function(gene, avg.expr, pct.expr, k = 20, sigma = 0.5){
  
  data1 = t(avg.expr[gene, ])
  data2 = t(pct.expr[gene, ])
  
  # Use pairwise squared Euclidean distance
  dist1 <- dist2(as.matrix(data1), as.matrix(data1))
  dist2 <- dist2(as.matrix(data2), as.matrix(data2))
  
  # Convert to affinity matrices using a Gaussian kernel
  W1 <- affinityMatrix(dist1, K = k, sigma = sigma)
  W2 <- affinityMatrix(dist2, K = k, sigma = sigma)
  
  # 3. Fuse the networks
  W_fused <- SNF(list(W1, W2), K = k)
  
  # 4. Perform spectral clustering on the fused network
  clusters <- spectralClustering(W_fused, K = 2)
  
  # 5. Plot the fused similarity matrix and clusters
  # print("Cluster assignments:")
  names(clusters) = colnames(W_fused)
  clusters.ordered = sort(clusters)
  order = names(clusters.ordered)
  
  # heatmap(W1[order,order], Rowv = NA, Colv = NA, symm = TRUE, main = "W1")
  # heatmap(W2[order,order], Rowv = NA, Colv = NA, symm = TRUE, main = "W2")
  # heatmap(W_fused[order,order], Rowv = NA, Colv = NA, symm = TRUE, main = "Fused Similarity Matrix")
  
  means = tapply(data1, clusters, mean)
  
  if(means[[1]] > means[[2]]) {
    -1*clusters+2
  } else {
    clusters-1
  }
}

PlotNTrees = function(version, n = 1, max.bl = NULL){
  lamprey.distance = ComputeDistancesFromTrees2(objectList$Lamprey, paste0('celltypephylo/Lamprey_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  shark.distance = ComputeDistancesFromTrees2(objectList$Shark, paste0('celltypephylo/Shark_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  kf.distance = ComputeDistancesFromTrees2(objectList$Killifish, paste0('celltypephylo/Killifish_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  ze.distance = ComputeDistancesFromTrees2(objectList$ZebrafishLyu, paste0('celltypephylo/ZebrafishLyu_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  ch.distance = ComputeDistancesFromTrees2(objectList$Chicken, paste0('celltypephylo/Chicken_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  li.distance = ComputeDistancesFromTrees2(objectList$Lizard, paste0('celltypephylo/Lizard_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  op.distance = ComputeDistancesFromTrees2(objectList$Opossum, paste0('celltypephylo/Opossum_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  rat.distance = ComputeDistancesFromTrees2(objectList$Rat, paste0('celltypephylo/Rat_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  ml.distance = ComputeDistancesFromTrees2(objectList$MouseLemur, paste0('celltypephylo/MouseLemur_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  mar.distance = ComputeDistancesFromTrees2(objectList$Marmoset, paste0('celltypephylo/Marmoset_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  macaque.distance = ComputeDistancesFromTrees2(objectList$Macaque, paste0('celltypephylo/Macaque_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
  human.distance = ComputeDistancesFromTrees2(objectList$Human, paste0('celltypephylo/Human_', version, '/'), plot = TRUE, n = n, max.bl = max.bl)
}

PlotDivergenceFinal = function(version = 'v16', ...){
  lamprey.distance = ComputeDistancesFromTrees2(objectList$Lamprey, paste0('celltypephylo/Lamprey_', version, '/'), plot = FALSE, ...)
  # shark.distance = ComputeDistancesFromTrees2(objectList$Shark, paste0('celltypephylo/Shark_', version, '/'), plot = FALSE, ...)
  # gf.distance = ComputeDistancesFromTrees2(objectList$Goldfish, paste0('celltypephylo/Goldfish_', version, '/'), plot = FALSE, ...)
  kf.distance = ComputeDistancesFromTrees2(objectList$Killifish, paste0('celltypephylo/Killifish_', version, '/'), plot = FALSE, ...)
  ze.distance = ComputeDistancesFromTrees2(objectList$ZebrafishLyu, paste0('celltypephylo/ZebrafishLyu_', version, '/'), plot = FALSE, ...)
  ch.distance = ComputeDistancesFromTrees2(objectList$Chicken, paste0('celltypephylo/Chicken_', version, '/'), plot = FALSE, ...)
  li.distance = ComputeDistancesFromTrees2(objectList$Lizard, paste0('celltypephylo/Lizard_', version, '/'), plot = FALSE, ...)
  op.distance = ComputeDistancesFromTrees2(objectList$Opossum, paste0('celltypephylo/Opossum_', version, '/'), plot = FALSE, ...)
  sq.distance = ComputeDistancesFromTrees2(objectList$Squirrel, paste0('celltypephylo/Squirrel_', version, '/'), plot = FALSE, ...)
  rat.distance = ComputeDistancesFromTrees2(objectList$Rat, paste0('celltypephylo/Rat_', version, '/'), plot = FALSE, ...)
  ml.distance = ComputeDistancesFromTrees2(objectList$MouseLemur, paste0('celltypephylo/MouseLemur_', version, '/'), plot = FALSE, ...)
  mar.distance = ComputeDistancesFromTrees2(objectList$Marmoset, paste0('celltypephylo/Marmoset_', version, '/'), plot = FALSE, ...)
  macaque.distance = ComputeDistancesFromTrees2(objectList$Macaque, paste0('celltypephylo/Macaque_', version, '/'), plot = FALSE, ...)
  human.distance = ComputeDistancesFromTrees2(objectList$Human, paste0('celltypephylo/Human_', version, '/'), plot = FALSE, ...)
  
  all.distance = rbind(lamprey.distance %>% mutate(species = 'Lamprey'),
                       # shark.distance %>% mutate(species = 'Shark'),
                       # gf.distance %>% mutate(species = 'Goldfish'),
                       kf.distance %>% mutate(species = 'Killifish'), 
                       ze.distance %>% mutate(species = 'Zebrafish'), 
                       ch.distance %>% mutate(species = 'Chicken'), 
                       li.distance %>% mutate(species = 'Lizard'),
                       op.distance %>% mutate(species = 'Opossum'),
                       sq.distance %>% mutate(species = 'Squirrel'),
                       rat.distance %>% mutate(species = 'Rat'),
                       ml.distance %>% mutate(species = 'MouseLemur'),
                       mar.distance %>% mutate(species = 'Marmoset'),
                       macaque.distance %>% mutate(species = 'Macaque'), 
                       human.distance %>% mutate(species = 'Human'))
  
  all.distance$gaba.rgc.norm.mg = all.distance$gaba.rgc / all.distance$rgc.mg
  all.distance$gly.rgc.norm.mg = all.distance$gly.rgc / all.distance$rgc.mg
  all.distance$gaba.rgc.norm.tree = all.distance$gaba.rgc / all.distance$tree.length
  all.distance$gly.rgc.norm.tree = all.distance$gly.rgc / all.distance$tree.length
  all.distance$gaba.rgc.norm.bb = all.distance$gaba.rgc / all.distance$backbone.length
  all.distance$gly.rgc.norm.bb = all.distance$gly.rgc / all.distance$backbone.length
  all.distance$gaba.rgc.norm.bc = all.distance$gaba.rgc / all.distance$rgc.bc
  all.distance$gly.rgc.norm.bc = all.distance$gly.rgc / all.distance$rgc.bc
  all.distance$pr.bc.norm.mg = all.distance$bc.pr / all.distance$bc.mg
  all.distance$pr.bc.norm.pr.mg = all.distance$bc.pr / all.distance$pr.mg
  all.distance$gaba.gly.norm.gaba.bc = all.distance$gaba.gly / all.distance$gaba.bc
  all.distance$gaba.bc.norm.gaba.mg = all.distance$gaba.bc / all.distance$gaba.mg
  
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.bb')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.bb')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt4 = (PrettyBarplot(summary.df, 
                        x = 'Species', 
                        y = 'Divergence', 
                        fill = 'Subclass',
                        color = 'black', 
                        add = c('mean_sd', 'jitter'), 
                        position = position_dodge(0.7), 
                        ylab = 'Normalized RGC distance',
                        add.params = list(fill = "Subclass", shape = 21, color = 'black'))
  ) + NoLegend()+
    #coord_cartesian(ylim = c(0,1.1)) + 
    RotatedAxis()
  
  plt4
}


ComputeAndPlotAllDistances2 <- function(species.list, version, ...) {
  
  distance.list = lapply(species.list, function(this.species){
    this.distance = ComputeDistancesFromTrees2(objectList[[this.species]], paste0('celltypephylo/', this.species, '_', version, '/'), plot = FALSE, ...)
    this.distance %>% mutate(species = this.species)
  })
  
  all.distance = do.call(rbind, distance.list)
  
  all.distance$gaba.rgc.norm.mg = all.distance$gaba.rgc / all.distance$rgc.mg
  all.distance$gly.rgc.norm.mg = all.distance$gly.rgc / all.distance$rgc.mg
  all.distance$gaba.rgc.norm.tree = all.distance$gaba.rgc / all.distance$tree.length
  all.distance$gly.rgc.norm.tree = all.distance$gly.rgc / all.distance$tree.length
  all.distance$gaba.rgc.norm.bb = all.distance$gaba.rgc / all.distance$backbone.length
  all.distance$gly.rgc.norm.bb = all.distance$gly.rgc / all.distance$backbone.length
  all.distance$gaba.rgc.norm.bc = all.distance$gaba.rgc / all.distance$rgc.bc
  all.distance$gly.rgc.norm.bc = all.distance$gly.rgc / all.distance$rgc.bc
  all.distance$pr.bc.norm.bb = all.distance$bc.pr / all.distance$backbone.length
  all.distance$pr.bc.norm.pr.mg = all.distance$bc.pr / all.distance$pr.mg
  all.distance$gaba.gly.norm.gaba.bc = all.distance$gaba.gly / all.distance$gaba.bc
  all.distance$gaba.bc.norm.gaba.mg = all.distance$gaba.bc / all.distance$gaba.mg
  
  # Raw score (numerator)
  print(PrettyBarplot(all.distance, x = 'species', y = 'gaba.rgc', add = c('mean_sd', 'jitter'), 
                      ylab = 'GABA AC-RGC distance') | 
          
          # Raw score (numerator)
          PrettyBarplot(all.distance, x = 'species', y = 'gly.rgc', add = c('mean_sd', 'jitter'), 
                        ylab = 'Gly AC-RGC distance'))
  
  # Normalizer 1 (denominator)
  print(PrettyBarplot(all.distance, x = 'species', y = 'rgc.mg', add = c('mean_sd', 'jitter'), 
                      ylab = 'RGC-MG distance') | 
          
          # Normalizer 2 (denominator)
          PrettyBarplot(all.distance, x = 'species', y = 'rgc.bc', add = c('mean_sd', 'jitter'), 
                        ylab = 'RGC-BC distance') |
          
          # Normalizer 3: tree length
          PrettyBarplot(all.distance, x = 'species', y = 'tree.length', add = c('mean_sd', 'jitter'),
                        ylab = 'Tree length') |
          
          # Back bone
          PrettyBarplot(all.distance, x = 'species', y = 'backbone.length', add = c('mean_sd', 'jitter'),
                        ylab = 'Backbone length'))
  
  # Foil 1: gaba-gly / gaba-bc
  print(PrettyBarplot(all.distance, x = 'species', y = 'gaba.gly.norm.gaba.bc', add = c('mean_sd', 'jitter'), 
                      ylab = 'GABA-Gly / GABA-BC') | 
          
          # Foil 2: gaba-gly / gaba-bc
          PrettyBarplot(all.distance, x = 'species', y = 'gaba.bc.norm.gaba.mg', add = c('mean_sd', 'jitter'), 
                        ylab = 'GABA-BC / GABA-MG'))
  
  # BC-PR / BC-MG
  print(PrettyBarplot(all.distance, x = 'species', y = 'pr.bc.norm.pr.mg', add = c('mean_sd', 'jitter'),
                      ylab = 'PR-BC distance normalized to MG') |
          # coord_cartesian(ylim = c(0,1))
          
          # BC-PR / PR-MG
          PrettyBarplot(all.distance, x = 'species', y = 'pr.bc.norm.bb', add = c('mean_sd', 'jitter'),
                        ylab = 'PR-BC distance normalized to backbone'))
  # coord_cartesian(ylim = c(0,1.5))
  
  
  # GABA and gly together
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.mg')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.mg')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt1 = (PrettyBarplot(summary.df, 
                        x = 'Species', 
                        y = 'Divergence', 
                        fill = 'Subclass',
                        color = 'black', 
                        add = c('mean_sd', 'jitter'), 
                        position = position_dodge(0.7), 
                        ylab = 'RGC divergence normalized to MG', 
                        add.params = list(fill = "Subclass", shape = 21, color = 'black')
  ) + 
    # coord_cartesian(ylim = c(0,1)) + 
    RotatedAxis()+
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey', alpha = 0.7))
  # ggsave('../../figures/my_figs/rgc-divergence-norm-mg.pdf', height = 4, width = 6)
  
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.bc')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.bc')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt2 = (PrettyBarplot(summary.df, 
                        x = 'Species', 
                        y = 'Divergence', 
                        fill = 'Subclass',
                        color = 'black', 
                        add = c('mean_sd', 'jitter'), 
                        position = position_dodge(0.7), 
                        ylab = 'RGC divergence normalized to BC',
                        add.params = list(fill = "Subclass", shape = 21, color = 'black')
  ) + 
    #coord_cartesian(ylim = c(0,1.1)) + 
    RotatedAxis() + 
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey', alpha = 0.7))
  
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.tree')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.tree')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt3 = (PrettyBarplot(summary.df, 
                        x = 'Species', 
                        y = 'Divergence', 
                        fill = 'Subclass',
                        color = 'black', 
                        add = c('mean_sd', 'jitter'), 
                        position = position_dodge(0.7), 
                        ylab = 'RGC divergence normalized to tree length',
                        add.params = list(fill = "Subclass", shape = 21, color = 'black'))
  ) + 
    #coord_cartesian(ylim = c(0,1.1)) + 
    RotatedAxis() 
  # geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey', alpha = 0.7)
  
  # ggsave('../../figures/my_figs/rgc-divergence-norm-bc.pdf', height = 4, width = 6)  
  
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.bb')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.bb')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt4 = (PrettyBarplot(summary.df, 
                        x = 'Species', 
                        y = 'Divergence', 
                        fill = 'Subclass',
                        color = 'black', 
                        add = c('mean_sd', 'jitter'), 
                        position = position_dodge(0.7), 
                        ylab = 'RGC divergence normalized to backbone',
                        add.params = list(fill = "Subclass", shape = 21, color = 'black'))
  ) + 
    #coord_cartesian(ylim = c(0,1.1)) + 
    RotatedAxis() 
  
  print(plt1 | plt2 | plt3 | plt4)
  
  return(summary.df)
}

ComputeAndPlotAllDistances <- function(version, ...) {

  lamprey.distance = ComputeDistancesFromTrees2(objectList$Lamprey, paste0('celltypephylo/Lamprey_', version, '/'), plot = FALSE, ...)
  shark.distance = ComputeDistancesFromTrees2(objectList$Shark, paste0('celltypephylo/Shark_', version, '/'), plot = FALSE, ...)
  gf.distance = ComputeDistancesFromTrees2(objectList$Goldfish, paste0('celltypephylo/Goldfish_', version, '/'), plot = FALSE, ...)
  kf.distance = ComputeDistancesFromTrees2(objectList$Killifish, paste0('celltypephylo/Killifish_', version, '/'), plot = FALSE, ...)
  ze.distance = ComputeDistancesFromTrees2(objectList$ZebrafishLyu, paste0('celltypephylo/ZebrafishLyu_', version, '/'), plot = FALSE, ...)
  zepr.distance = ComputeDistancesFromTrees2(objectList$ZebrafishLyu, paste0('celltypephylo/ZebrafishLyuPR_', version, '/'), plot = FALSE, ...)
  ch.distance = ComputeDistancesFromTrees2(objectList$Chicken, paste0('celltypephylo/Chicken_', version, '/'), plot = FALSE, ...)
  li.distance = ComputeDistancesFromTrees2(objectList$Lizard, paste0('celltypephylo/Lizard_', version, '/'), plot = FALSE, ...)
  op.distance = ComputeDistancesFromTrees2(objectList$Opossum, paste0('celltypephylo/Opossum_', version, '/'), plot = FALSE, ...)
  sq.distance = ComputeDistancesFromTrees2(objectList$Squirrel, paste0('celltypephylo/Squirrel_', version, '/'), plot = FALSE, ...)
  rat.distance = ComputeDistancesFromTrees2(objectList$Rat, paste0('celltypephylo/Rat_', version, '/'), plot = FALSE, ...)
  ts.distance = ComputeDistancesFromTrees2(objectList$TreeShrew, paste0('celltypephylo/TreeShrew_', version, '/'), plot = FALSE, ...)
  ml.distance = ComputeDistancesFromTrees2(objectList$MouseLemur, paste0('celltypephylo/MouseLemur_', version, '/'), plot = FALSE, ...)
  mar.distance = ComputeDistancesFromTrees2(objectList$Marmoset, paste0('celltypephylo/Marmoset_', version, '/'), plot = FALSE, ...)
  macaque.distance = ComputeDistancesFromTrees2(objectList$Macaque, paste0('celltypephylo/Macaque_', version, '/'), plot = FALSE, ...)
  human.distance = ComputeDistancesFromTrees2(objectList$Human, paste0('celltypephylo/Human_', version, '/'), plot = FALSE, ...)
  
  all.distance = rbind(lamprey.distance %>% mutate(species = 'Lamprey'),
                       shark.distance %>% mutate(species = 'Shark'),
                       gf.distance %>% mutate(species = 'Goldfish'),
                       kf.distance %>% mutate(species = 'Killifish'), 
                       ze.distance %>% mutate(species = 'Zebrafish'),
                       zepr.distance %>% mutate(species = 'ZebrafishPR'),
                       ch.distance %>% mutate(species = 'Chicken'), 
                       li.distance %>% mutate(species = 'Lizard'),
                       op.distance %>% mutate(species = 'Opossum'),
                       sq.distance %>% mutate(species = 'Squirrel'),
                       rat.distance %>% mutate(species = 'Rat'),
                       ts.distance %>% mutate(species = 'TreeShrew'),
                       ml.distance %>% mutate(species = 'MouseLemur'),
                       mar.distance %>% mutate(species = 'Marmoset'),
                       macaque.distance %>% mutate(species = 'Macaque'), 
                       human.distance %>% mutate(species = 'Human'))
  
  all.distance$gaba.rgc.norm.mg = all.distance$gaba.rgc / all.distance$rgc.mg
  all.distance$gly.rgc.norm.mg = all.distance$gly.rgc / all.distance$rgc.mg
  all.distance$gaba.rgc.norm.tree = all.distance$gaba.rgc / all.distance$tree.length
  all.distance$gly.rgc.norm.tree = all.distance$gly.rgc / all.distance$tree.length
  all.distance$gaba.rgc.norm.bb = all.distance$gaba.rgc / all.distance$backbone.length
  all.distance$gly.rgc.norm.bb = all.distance$gly.rgc / all.distance$backbone.length
  all.distance$gaba.rgc.norm.bc = all.distance$gaba.rgc / all.distance$rgc.bc
  all.distance$gly.rgc.norm.bc = all.distance$gly.rgc / all.distance$rgc.bc
  all.distance$pr.bc.norm.bb = all.distance$bc.pr / all.distance$backbone.length
  all.distance$pr.bc.norm.pr.bb = all.distance$bc.pr / all.distance$backbone.length
  all.distance$gaba.gly.norm.gaba.bc = all.distance$gaba.gly / all.distance$gaba.bc
  all.distance$gaba.bc.norm.gaba.mg = all.distance$gaba.bc / all.distance$gaba.mg
  
  # Raw score (numerator)
  print(PrettyBarplot(all.distance, x = 'species', y = 'gaba.rgc', add = c('mean_sd', 'jitter'), 
                  ylab = 'GABA AC-RGC distance') | 
    
  # Raw score (numerator)
  PrettyBarplot(all.distance, x = 'species', y = 'gly.rgc', add = c('mean_sd', 'jitter'), 
                  ylab = 'Gly AC-RGC distance'))
  
  # Normalizer 1 (denominator)
  print(PrettyBarplot(all.distance, x = 'species', y = 'rgc.mg', add = c('mean_sd', 'jitter'), 
                  ylab = 'RGC-MG distance') | 
    
  # Normalizer 2 (denominator)
  PrettyBarplot(all.distance, x = 'species', y = 'rgc.bc', add = c('mean_sd', 'jitter'), 
                  ylab = 'RGC-BC distance') |
  
  # Normalizer 3: tree length
  PrettyBarplot(all.distance, x = 'species', y = 'tree.length', add = c('mean_sd', 'jitter'),
                  ylab = 'Tree length') |
    
  # Back bone
  PrettyBarplot(all.distance, x = 'species', y = 'backbone.length', add = c('mean_sd', 'jitter'),
                  ylab = 'Backbone length'))

  # Foil 1: gaba-gly / gaba-bc
  print(PrettyBarplot(all.distance, x = 'species', y = 'gaba.gly.norm.gaba.bc', add = c('mean_sd', 'jitter'), 
                  ylab = 'GABA-Gly / GABA-BC') | 
    
  # Foil 2: gaba-gly / gaba-bc
  PrettyBarplot(all.distance, x = 'species', y = 'gaba.bc.norm.gaba.mg', add = c('mean_sd', 'jitter'), 
                  ylab = 'GABA-BC / GABA-MG'))
  
  # Another possible normalizer
  # PrettyBarplot(all.distance, x = 'species', y = 'gaba.bc', add = c('mean_sd', 'jitter'), 
                # ylab = 'GABA-BC / GABA-MG'))
  
  # # Raw / MG
  # PrettyBarplot(all.distance, x = 'species', y = 'gaba.rgc.norm.mg', add = c('mean_sd', 'jitter')) + 
  #   coord_cartesian(ylim = c(0,1))
  # 
  # PrettyBarplot(all.distance, x = 'species', y = 'gly.rgc.norm.mg', add = c('mean_sd', 'jitter'))+
  #   coord_cartesian(ylim = c(0,1))
  # 
  # # Raw / total tree length
  # PrettyBarplot(all.distance, x = 'species', y = 'gaba.rgc.norm.tree', add = c('mean_sd', 'jitter'))+
  #   coord_cartesian(ylim = c(0,1))
  # 
  # # Raw / BC
  # PrettyBarplot(all.distance, x = 'species', y = 'gaba.rgc.norm.bc', add = c('mean_sd', 'jitter'))+
  #   coord_cartesian(ylim = c(0,1))
  
  # BC-PR / BC-MG
  print(PrettyBarplot(all.distance, x = 'species', y = 'pr.bc.norm.bb', add = c('mean_sd', 'jitter'),
                ylab = 'PR-BC distance normalized to backbone') |
  # coord_cartesian(ylim = c(0,1))
  
  # BC-PR / PR-MG
  PrettyBarplot(all.distance, x = 'species', y = 'pr.bc.norm.pr.bb', add = c('mean_sd', 'jitter'),
                ylab = 'PR-BC distance normalized to backbone'))
  # coord_cartesian(ylim = c(0,1.5))
  
  
  # GABA and gly together
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.mg')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.mg')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt1 = (PrettyBarplot(summary.df, 
                x = 'Species', 
                y = 'Divergence', 
                fill = 'Subclass',
                color = 'black', 
                add = c('mean_sd', 'jitter'), 
                position = position_dodge(0.7), 
                ylab = 'RGC divergence normalized to MG', 
                add.params = list(fill = "Subclass", shape = 21, color = 'black')
  ) + 
    # coord_cartesian(ylim = c(0,1)) + 
    RotatedAxis()+
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey', alpha = 0.7))
  # ggsave('../../figures/my_figs/rgc-divergence-norm-mg.pdf', height = 4, width = 6)
  
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.bc')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.bc')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt2 = (PrettyBarplot(summary.df, 
                x = 'Species', 
                y = 'Divergence', 
                fill = 'Subclass',
                color = 'black', 
                add = c('mean_sd', 'jitter'), 
                position = position_dodge(0.7), 
                ylab = 'RGC divergence normalized to BC',
                add.params = list(fill = "Subclass", shape = 21, color = 'black')
  ) + 
    #coord_cartesian(ylim = c(0,1.1)) + 
    RotatedAxis() + 
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey', alpha = 0.7))
  
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.tree')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.tree')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt3 = (PrettyBarplot(summary.df, 
                x = 'Species', 
                y = 'Divergence', 
                fill = 'Subclass',
                color = 'black', 
                add = c('mean_sd', 'jitter'), 
                position = position_dodge(0.7), 
                ylab = 'RGC divergence normalized to tree length',
                add.params = list(fill = "Subclass", shape = 21, color = 'black'))
  ) + 
    #coord_cartesian(ylim = c(0,1.1)) + 
    RotatedAxis() 
  # geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey', alpha = 0.7)
  
  # ggsave('../../figures/my_figs/rgc-divergence-norm-bc.pdf', height = 4, width = 6)  
  
  summary.df = rbind(all.distance[, c('species', 'gaba.rgc.norm.bb')] %>% mutate(Subclass = 'gabaAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')),
                     all.distance[, c('species', 'gly.rgc.norm.bb')] %>% mutate(Subclass = 'glyAC') %>%
                       setNames(c('Species', 'Divergence', 'Subclass')))
  
  plt4 = (PrettyBarplot(summary.df, 
                        x = 'Species', 
                        y = 'Divergence', 
                        fill = 'Subclass',
                        color = 'black', 
                        add = c('mean_sd', 'jitter'), 
                        position = position_dodge(0.7), 
                        ylab = 'RGC divergence normalized to backbone',
                        add.params = list(fill = "Subclass", shape = 21, color = 'black'))
  ) + 
    #coord_cartesian(ylim = c(0,1.1)) + 
    RotatedAxis() 
  
  print(plt1 | plt2 | plt3 | plt4)
  
  return(summary.df)
}

ReadTrees = function(basename){
  files = list.files(path = basename, pattern = "\\.treefile$", include.dirs = TRUE)
  lapply(paste0(basename, files), read.tree)
}

ComputeDistancesFromTrees2 = function(object, basename, plot = FALSE, n = NULL, select = NULL, normalize.by.backbone = FALSE, winsorize = FALSE, ...){
  
  files = list.files(path = basename, pattern = "\\.treefile$", include.dirs = TRUE)
  if(!is.null(n)) files = files[1:n]
  if(!is.null(select)) files = files[select]
  print(files)
  do.call(rbind, lapply(paste0(basename, files), function(file){
    metadata = Metadata(object, 'annotated', 'cell_class', 'cell_class2')
    
    # tree = read.tree(paste0(basename, '-', seed, '.treefile'))
    tree = read.tree(file)
    
    # Plot trees
    if(plot) {
      return(plot_tree(tree, metadata, title = basename(file), ...))
    }
    
    # Compute length of tree backbone (non-tip length)
    backbone.length = FindBackboneLength(tree, winsorize = winsorize)
    if(normalize.by.backbone) {
      tree$edge.length = tree$edge.length/backbone.length
    }
    
    # Compute distances 
    dist.mat = ComputeCopheneticDistance(tree, metadata)
    
    df = data.frame(gaba.rgc = dist.mat['gabaAC', 'RGC'], 
                   gaba.gly = dist.mat['gabaAC', 'glyAC'], 
                   gly.rgc = dist.mat['glyAC', 'RGC'],
                   rgc.mg = dist.mat['MG', 'RGC'], 
                   rgc.bc = dist.mat['BC', 'RGC'], 
                   bc.pr = dist.mat['BC', 'PR'], 
                   bc.mg = dist.mat['BC', 'MG'], 
                   pr.mg = dist.mat['PR', 'MG'], 
                   gaba.bc = dist.mat['gabaAC', 'BC'], 
                   gaba.mg = dist.mat['gabaAC', 'MG'], 
                   tree.length = sum(tree$edge.length), 
                   backbone.length = backbone.length)
    
  }))
}

ComputeDistancesFromTrees = function(object, basename, plot = TRUE, n = 10, ...){
  
  do.call(rbind, lapply(seq_len(n), function(seed){
    metadata = Metadata(object, 'annotated', 'cell_class', 'cell_class2')
    
    tree = read.tree(paste0(basename, '-', seed, '-2000.treefile'))
    
    # Plot trees
    if(plot) {
      print(plot_tree(tree, metadata, ...))
    }
    
    # Compute distances 
    dist.mat = ComputeCopheneticDistance(tree, 
                                         metadata)
    
    data.frame(gaba.rgc = dist.mat['gabaAC', 'RGC'], 
               gaba.bc = dist.mat['gabaAC', 'BC'],
               gaba.gly = dist.mat['gabaAC', 'glyAC'],
               gaba.mg = dist.mat['gabaAC', 'MG'],
               gly.rgc = dist.mat['glyAC', 'RGC'],
               rgc.mg = dist.mat['MG', 'RGC'], 
               rgc.bc = dist.mat['BC', 'RGC'], 
               bc.pr = dist.mat['BC', 'PR'], 
               bc.mg = dist.mat['BC', 'MG'], 
               pr.mg = dist.mat['PR', 'MG'], 
               tree.length = sum(tree$edge.length))
  }))
}

GenerateBinaryMatrix3 = function(object, 
                                 species, 
                                 pct.exp.thres = 30, 
                                 pval.thres = 1e-10, 
                                 noise.sd = NULL, 
                                 seed = 1, 
                                 top.n.genes = NULL, 
                                 verbose = FALSE, 
                                 nthreads = 1, 
                                 sort.genes.by = c('none', 'r.squared', 'LRT'),
                                 selection.method = 'BIC',
                                 r.squared.threshold = 0.6, 
                                 iqtree.model = 'MFP+ASC', 
                                 normalize.pct.expr = FALSE, 
                                 use = c('mv', 'mv3', 'expr', 'pct', 'best', 'k'),
                                 use.norm.exp = FALSE, 
                                 use.expr.score = FALSE,
                                 use.mv = TRUE,
                                 n.bootstraps = 500, 
                                 bootstrap = TRUE,
                                 compute.var.genes = TRUE, 
                                 nfeatures = 3000,
                                 noise.alpha = 0.25, 
                                 asr = FALSE, 
                                 expansion = 10, 
                                 fast = FALSE){
  
  # withr::with_seed(seed, {
  
  # Step 0: file paths
  dir.create(paste0('celltypephylo/', species), showWarnings = FALSE)
  superprefix = paste0('celltypephylo/', species, '/', species, '-', pct.exp.thres, '-', pval.thres, '-', noise.sd, '-', sort.genes.by, '-', use, '-', nfeatures)
  fafile = paste0(superprefix, '-', seed, '.fa')
  phyfile = gsub(".fa", ".phy", fafile)
  prefix = gsub(".fa", "", fafile)
  treefile = gsub("fa", "treefile", fafile)
  binfile = paste0(prefix, '-binarization.pdf')
  genefile = paste0(prefix, '-genes.rds')
  plotfile = paste0(prefix, '-plot.pdf')
  cvfile = paste0(prefix, '-mean-cv.pdf')
  ratefile = paste0(prefix, '-rate-sorted.txt')
  datsumfile = paste0(prefix, '-datsum.rds')
  datfile = paste0(prefix, '-raw-input.rds')
  binmatfile = paste0(prefix, '-bin-mat.rds')
  message("Output files are named: ", prefix)
  
  # Step 1: identify genes that fluctuate a lot across clusters (high variance in terms of percentage expressed)
  object.ds = DownsampleSeurat(object, group.by = 'annotated', size = 100, seed = seed)
  message('Step 1: calculating percentage expressed for clusters...')
  
  # Establish type order
  object.ds$annotated = factor(object.ds$annotated) # Remove unused factors
  if(is.factor(object.ds$annotated)) {
    type.order = gsub(' ', '-', levels(object.ds$annotated))
  } else {
    type.order = gsub(' ', '-', Metadata(object.ds, 'cell_class2', 'annotated')$annotated)
  }
  object.ds$annotated = factor(gsub(' ', '-', as.character(object.ds$annotated)), levels = type.order)
  
  avg.expr = LogAvgExpr(object.ds, 
                        group.by = 'annotated', 
                        features = rownames(object), 
                        assay = 'RNA')[,type.order]
  
  # Scale to a maximum value of 100
  # avg.expr = avg.expr*100/max(avg.expr) # no need
  pct.expr = PercentageExpressed2(object.ds, 
                                  group.by = 'annotated', 
                                  features = rownames(object))[,type.order]
  
  if(compute.var.genes){
    message('Computing ', nfeatures, ' highly variable features')
    # var.genes = FindVariableFeatures(object.ds, nfeatures = nfeatures)@assays$RNA@var.features
    var.genes = SelectFeatures(object.ds, nfeatures = nfeatures, group.by = 'annotated', min.pct.expressed = 20)
    # var.genes = tail(FindVariableFeatures(object.ds, nfeatures = nfeatures)@assays$RNA@var.features, 100)
    
    # Ensure that they vary by at least X percent
    # vars = apply(pct.expr[var.genes,], 1, function(x) max(x) - min(x))
    # var.genes = names(which(vars > 10))
    
  } else {
    vars = apply(pct.expr, 1, function(x) max(x) - min(x))
    var.genes = names(which(vars > pct.exp.thres))
  }
  
  # Step 2: identify which of these genes are bimodal in their distribution
  message('Step 2: testing for bimodality at a LRT p-value < ', pval.thres, ' on ', length(var.genes), ' genes!')
  gmm_res = pbapply::pblapply(var.genes, function(gene) {
    BinarizeExpressionMV2(gene, 
                          avg.expr,
                          pct.expr,
                          seed = seed, 
                          verbose = verbose, 
                          noise.alpha = noise.alpha, 
                          expansion = expansion, 
                          fast = fast)
  }) %>% setNames(var.genes)
  
  gmm_res.clean = gmm_res[sapply(gmm_res, function(x) length(x$mv.vector) > 1)]
  message('GMM succeeded on ', length(gmm_res.clean), ' genes')
  pvals.clean = p.adjust(as.numeric(sapply(gmm_res.clean, function(x) x$p_value)) %>% setNames(names(gmm_res.clean)), method = 'fdr')
  message('p-values range: ', paste0(range(pvals.clean), collapse = ' - '))
  bic.delta = as.numeric(sapply(gmm_res.clean, function(x) x$delta.bic)) %>% setNames(names(gmm_res.clean))
  message('delta BIC range: ', paste0(range(bic.delta), collapse = ' - '))
  r.squared.values.mv = as.numeric(sapply(gmm_res.clean, function(x) x$r.squared.mv)) %>% setNames(names(gmm_res.clean))
  r.squared.values.mv3 = as.numeric(sapply(gmm_res.clean, function(x) x$r.squared.mv3)) %>% setNames(names(gmm_res.clean))
  r.squared.values.expr = as.numeric(sapply(gmm_res.clean, function(x) x$r.squared.expr)) %>% setNames(names(gmm_res.clean))
  r.squared.values.pct = as.numeric(sapply(gmm_res.clean, function(x) x$r.squared.pct)) %>% setNames(names(gmm_res.clean))
  r.squared.values.k = as.numeric(sapply(gmm_res.clean, function(x) x$r.squared.k)) %>% setNames(names(gmm_res.clean))
  r.squared.values.pct.expr = as.numeric(sapply(gmm_res.clean, function(x) x$r.squared.pct.expr)) %>% setNames(names(gmm_res.clean))
  
  r.squared.df = data.frame(r.squared.values.mv, 
                            r.squared.values.mv3,
                            r.squared.values.expr,
                            r.squared.values.pct, 
                            r.squared.values.k)
  r.squared.best = (t(apply(r.squared.df, 1, function(x) x == max(x))))
  print(colSums(r.squared.best))
  # r.squared.df$best = apply(r.squared.df, 1, function(x) paste0(c('mv', 'expr', 'pct')[which(x == max(x))], collapse = '-'))
  # table(r.squared.values > r.squared.values.expr)
  # table(r.squared.values > r.squared.values.pct)
  # table(r.squared.values.pct > r.squared.values.expr)
  
  message('median r^2 value mv: ', median(r.squared.values.mv))
  message('median r^2 value mv3: ', median(r.squared.values.mv3))
  message('median r^2 value expr: ', median(r.squared.values.expr))
  message('median r^2 value pct: ', median(r.squared.values.pct))
  message('median r^2 value k-means: ', median(r.squared.values.k))
  # message('R^2 range: ', paste0(range(r.squared.values), collapse = '-'))
  
  high.pvalue = names(which(pvals.clean > pval.thres))
  message('Removing ', length(high.pvalue), ' genes with LRT p-value > ', pval.thres)
  high.bic = names(which(bic.delta > 0))
  message('Removing ', length(high.bic), ' genes with BIC > 0')
  low.consistency = names(which(r.squared.values.pct.expr < 0.5))
  message('Removing ', length(low.consistency), ' genes with low expr/pct consistency')
  bimodal.genes = setdiff(names(pvals.clean), c(high.pvalue, high.bic, low.consistency))
  message('Found ', length(bimodal.genes), ' bimodal genes!')
  
  # Step 3: binarize their expression
  message('Step 3: binarizing expression...')
  dat.use = list(avg.expr[bimodal.genes,], pct.expr[bimodal.genes,])
  saveRDS(dat.use, datfile)
  use = match.arg(use)
  message('Using binarization method ', use, '!')
  if(use == 'mv'){
    bin.mat = do.call(rbind, lapply(gmm_res.clean, function(x) t(x$mv.vector)))
    r.squared.values = r.squared.values.mv
  } else if(use == 'mv3'){
    bin.mat = do.call(rbind, lapply(gmm_res.clean, function(x) t(x$mv3.vector)))
    r.squared.values = r.squared.values.mv3
  } else if(use == 'pct'){
    bin.mat = do.call(rbind, lapply(gmm_res.clean, function(x) t(x$pct.vector)))
    r.squared.values = r.squared.values.pct
  } else if(use == 'expr'){
    bin.mat = do.call(rbind, lapply(gmm_res.clean, function(x) t(x$expr.vector)))
    r.squared.values = r.squared.values.expr
  } else if(use == 'k'){
    bin.mat = do.call(rbind, lapply(gmm_res.clean, function(x) t(x$k.vector)))
    r.squared.values = r.squared.values.k
  } else if(use == 'best'){
    pick = apply(r.squared.df, 1, function(x) c('mv.vector', 'mv3.vector', 'expr.vector', 'pct.vector', 'k.vector')[which.max(x)])
    bin.mat = do.call(rbind, lapply(seq_along(pick), function(i) (gmm_res.clean[[i]][[ pick[[i]] ]])))
    r.squared.values = apply(r.squared.df, 1, function(x) x[which.max(x)])
  }
  
  colnames(bin.mat) = colnames(dat.use[[1]]) #gsub(' ', '-', colnames(dat.use[[1]]))
  rownames(bin.mat) = names(gmm_res.clean)
  bin.mat = bin.mat[bimodal.genes,]
  
  # Step 4.1: Compute mean variance explained and filter out poorly represented genes
  message('Step 4: checking quality of binarization...')
  gene.pool = intersect(bimodal.genes, names(which(r.squared.values > r.squared.threshold)))
  saveRDS(bin.mat[gene.pool,], binmatfile)
  message('Found ', length(gene.pool), ' genes with r^2 greater than ', r.squared.threshold)
  
  sort.genes.by <- match.arg(sort.genes.by)
  switch(sort.genes.by,
         none = message("Not sorting genes; will sample randomly..."),
         r.squared = message("Sorting genes by R^2 value..."),
         LRT = message("Sorting genes by LRT p-value..."))
  
  # Step 4.2: Randomly sample genes from the pool
  if(sort.genes.by == 'none'){
    if(length(gene.pool) > top.n.genes){
      genes.use = sample(gene.pool, top.n.genes)
      message('Randomly sampling ', top.n.genes, ' genes from this pool...')
    } else {
      genes.use = gene.pool
    }
  } else if(sort.genes.by == 'r.squared'){
    r.squared.values.sorted = sort(r.squared.values, decreasing = TRUE)
    genes.use = head(names(r.squared.values.sorted), top.n.genes)
    message('Picking top ', top.n.genes, ' genes by r^2...')
  } else if(sort.genes.by == 'LRT'){
    pvals.clean.filt.sorted = sort(pvals.clean[names(pvals.clean) %in% gene.pool])
    genes.use = head(names(pvals.clean.filt.sorted), top.n.genes)
    message('Picking top ', top.n.genes, ' genes by LRT p-value...')
  }
  
  message('Mean r^2 value: ', mean(r.squared.values[genes.use]))
  message('95% CI: ', quantile(r.squared.values[genes.use], probs = c(0.025, 0.975))[[1]], 
          '-', quantile(r.squared.values[genes.use], probs = c(0.025, 0.975))[[2]])
  bin.mat = bin.mat[genes.use,]
  message('Using ', length(genes.use), ' genes across ', ncol(bin.mat), ' clusters!')
  
  # SNF binarization
  # message('Running SNF binarization...')
  # bin.mat = do.call(rbind, pbapply::pblapply(genes.use, SNFBinarize, avg.expr, pct.expr))
  # rownames(bin.mat) = genes.use
  
  # Step 5: Save binarization
  col_fun = colorRamp2(c(0, 1), c("white", "red"))
  hc_rows <- hclust(dist(bin.mat))
  row_order <- hc_rows$order
  ordered_row_names <- rownames(bin.mat)[row_order]
  ordered_row_names_rsquared = paste0(ordered_row_names, ' (', round(r.squared.values[ordered_row_names], 2), ')')
  # pdf(binfile, width = 0.2*length(ordered_row_names) + 1, height = 12.5) #family = 'ArialMT')
  DotPlot3(object, 
           features = ordered_row_names, 
           coord.flip = FALSE, 
           binarization = bin.mat[ordered_row_names,type.order], 
           group.by = 'annotated')
  ggsave(binfile, width = 0.25*length(ordered_row_names) + 1, height = 0.25*length(type.order)+1, limitsize = FALSE)

  # dev.off()
  # print(DotPlot3(object, features = ordered_row_names, coord.flip = FALSE))
  # print(ggHeatmap(t(bin.mat[ordered_row_names,unique(object$annotated)]), 
  #                 high_color = '#584B9FFF', 
  #                 legend.name = 'Scaled\nexpression') + 
  #         theme_cowplot() + 
  #         RotatedAxis())
  # dev.off()
  
  # print(Heatmap2((apply(avg.expr[ordered_row_names,], 1, min_max_norm)), 
  #                col = col_fun, 
  #                row.font.size = 5, 
  #                column_labels = ordered_row_names_rsquared))
  # print(Heatmap2((apply(pct.expr[ordered_row_names,], 1, min_max_norm)), 
  #                col = col_fun, 
  #                row.font.size = 5, 
  #                column_labels = ordered_row_names_rsquared))
  # print(Heatmap2(t(bin.mat[ordered_row_names,]), 
  #                col = col_fun, 
  #                row.font.size = 5, 
  #                column_labels = ordered_row_names_rsquared))
  # dev.off()
  
  saveRDS(list(gene.pool = gene.pool, 
               genes.use = ordered_row_names, 
               bimodal = bimodal.genes, 
               var.genes = var.genes), 
          genefile)
  
  # Step 5.1: check mean/cv relationship
  means = rowMeans(bin.mat)
  sds = apply(bin.mat, 1, function(x) sqrt(sum((x - mean(x))^2)/(length(x))))
  datsum = data.frame(mean = means, 
                      sd = sds, 
                      cv = sds/means)
  saveRDS(datsum, datsumfile)
  
  pdf(cvfile, height = 5, width = 5)
  print(ScatterPlot(datsum, x = 'mean', y = 'cv', logX = TRUE, logY = TRUE, lm = FALSE, r_pval = FALSE)+
          stat_function(fun = function(x) 1/sqrt(x), color = "red", size = 1.2, linetype = "dashed"))
  print(PrettyHistogram(datsum$mean, xlab = 'Mean'))
  dev.off()
  
  # Step 6: Construct phylogeny for binary trait data using IQTree
  library(phangorn)
  phydat = phyDat(t(bin.mat), type="USER", levels=c("0", "1"), compress = FALSE)
  write.phyDat(phydat, phyfile, format = "phylip")
  if(asr) {
    args = c('-s', phyfile, '-redo', '-pre', prefix, '-nt', nthreads, '-wsr', '-m', iqtree.model, '--ancestral')
  } else {
    args = c('-s', phyfile, '-redo', '-pre', prefix, '-nt', nthreads, '-wsr', '-m', iqtree.model)
  }
  system2('iqtree2', 
          args, 
          stdout = ifelse(verbose, "", FALSE))
  
  tree = read.tree(treefile)
  # tree = root(tree, outgroup = 'MG-1')
  pdf(plotfile, height = 50, width = 30)
  plot(tree)
  edgelabels(round(tree$edge.length, 2))
  dev.off()
  
  # Step 7: rank genes by gene-specific substitution rate
  tab = read.table(paste0(prefix, '.rate'), header=TRUE)
  if(nrow(tab) > 0){ # if best model does not include rate heterogeneity, this file will be blank
    tab$gene = rownames(bin.mat)
    tab.sorted = tab %>% arrange(Rate)
    write.table(tab.sorted, ratefile)
    message('Top conserved genes: ', paste(head(tab.sorted$gene), collapse = ', '))
  }
  
  superprefix
  
  # })
}

GenerateBinaryMatrix2 = function(object, 
                                 species, 
                                 pct.exp.thres = 30, 
                                 pval.thres = 1e-6, 
                                 noise.sd = NULL, 
                                 seed = 1, 
                                top.n.genes = NULL, 
                                verbose = FALSE, 
                                nthreads = 1, 
                                sort.genes.by = c('none', 'r.squared', 'LRT'),
                                r.squared.threshold = 0.6, 
                                iqtree.model = 'MFP+ASC', 
                                normalize.pct.expr = FALSE, 
                                use.norm.exp = TRUE, 
                                use.expr.score = FALSE,
                                n.bootstraps = 500, 
                                bootstrap = TRUE,
                                compute.var.genes = FALSE, 
                                noise.alpha = 0.05, 
                                asr = FALSE){
  # withr::with_seed(seed, {
  
  # Step 0: file paths
  dir.create(paste0('celltypephylo/', species), showWarnings = FALSE)
  superprefix = paste0('celltypephylo/', species, '/', species, '-', pct.exp.thres, '-', pval.thres, '-', noise.sd, '-', sort.genes.by, '-', top.n.genes)
  fafile = paste0(superprefix, '-', seed, '.fa')
  phyfile = gsub(".fa", ".phy", fafile)
  prefix = gsub(".fa", "", fafile)
  treefile = gsub("fa", "treefile", fafile)
  binfile = paste0(prefix, '-binarization.pdf')
  genefile = paste0(prefix, '-genes.rds')
  plotfile = paste0(prefix, '-plot.pdf')
  cvfile = paste0(prefix, '-mean-cv.pdf')
  ratefile = paste0(prefix, '-rate-sorted.txt')
  datsumfile = paste0(prefix, '-datsum.rds')
  datfile = paste0(prefix, '-raw-input.rds')
  binmatfile = paste0(prefix, '-bin-mat.rds')
  message("Output files are named: ", prefix)
  
  # Step 1: identify genes that fluctuate a lot across clusters (high variance in terms of percentage expressed)
  object.ds = DownsampleSeurat(object, group.by = 'annotated', size = 100, seed = seed)
  message('Step 1: calculating percentage expressed for clusters...')
  
  avg.expr = LogAvgExpr(object.ds, group.by = 'annotated', 
                        features = rownames(object), assay = 'RNA')
  
  # Scale to a maximum value of 100
  # avg.expr = avg.expr*100/max(avg.expr) # no need
  pct.expr = PercentageExpressed2(object.ds, 
                                  group.by = 'annotated', 
                                  features = rownames(object))

    
  if(normalize.pct.expr) {
    message('Normalizing percentage expressed by size factors...')
    norm.factors = median(colMeans(pct.expr))/colMeans(pct.expr)
    pct.expr = as.data.frame(t(t(pct.expr) * norm.factors))
  }
  
  if(use.norm.exp) {
    message('Using log-normalized expression instead of percentage expression')
    dat.use = avg.expr
  } else if(use.expr.score){
    # expression score = pct.expr/100 * avg.expr
    message('Using expression score (pct.expr/100 * avg.expr)')
    dat.use = pct.expr/100 * avg.expr
  } else {
    dat.use = pct.expr
  }
  
  if(compute.var.genes){
    var.genes = FindVariableFeatures(object.ds, nfeatures = 8000)@assays$RNA@var.features
  } else {
    vars = apply(pct.expr, 1, function(x) max(x) - min(x))
    var.genes = names(which(vars > pct.exp.thres))
  }
  
  # Step 2: identify which of these genes are bimodal in their distribution
  message('Step 2: testing for bimodality at a LRT p-value < ', pval.thres, ' on ', length(var.genes), ' genes!')
  gmm_res = pbapply::pbapply(dat.use[var.genes,], 1, function(expr) {
    BimodalityTest2(expr,
                    seed = seed, 
                    noise.sd = noise.sd, 
                    verbose = verbose, 
                    n.bootstraps = n.bootstraps, 
                    noise.alpha = noise.alpha, 
                    bootstrap = bootstrap)
    })
  
  gmm_res.clean = gmm_res[sapply(gmm_res, function(x) length(x[[1]]) > 1)]
  message('GMM succeeded on ', length(gmm_res.clean), ' genes')
  pvals.clean = as.numeric(sapply(gmm_res.clean, function(x) x[[2]])) %>% setNames(names(gmm_res.clean))
  bimodal.genes = names(which(pvals.clean < pval.thres))
  message('Found ', length(bimodal.genes), ' bimodal genes!')
  
  # Step 3: binarize their expression
  message('Step 3: binarizing expression...')
  dat.use = dat.use[bimodal.genes,]
  saveRDS(dat.use, datfile)
  bin.mat = do.call(rbind, lapply(gmm_res.clean, function(x) t(x[[1]])))
  colnames(bin.mat) = gsub(' ', '-', colnames(dat.use))
  rownames(bin.mat) = names(gmm_res.clean)
  bin.mat = bin.mat[bimodal.genes,]
  
  # Step 4.1: Compute mean variance explained and filter out poorly represented genes
  message('Step 4: checking quality of binarization...')
  r.squared.values = sapply(seq_len(nrow(bin.mat)), function(i){
    # print(i)
    summary(lm(as.numeric(dat.use[i,]) ~ as.numeric(bin.mat[i,])))$r.squared
  }) %>% setNames(bimodal.genes)
  gene.pool = names(which(r.squared.values > r.squared.threshold))
  saveRDS(bin.mat[gene.pool,], binmatfile)
  message('Found ', length(gene.pool), ' genes with r^2 greater than ', r.squared.threshold)
  
  sort.genes.by <- match.arg(sort.genes.by)
  switch(sort.genes.by,
         none = print("Not sorting genes; will sample randomly..."),
         r.squared = print("Sorting genes by R^2 value..."),
         LRT = print("Sorting genes by LRT p-value..."))
  
  # Step 4.2: Randomly sample genes from the pool
  if(sort.genes.by == 'none'){
    if(length(gene.pool) > top.n.genes){
      genes.use = sample(gene.pool, top.n.genes)
      message('Randomly sampling ', top.n.genes, ' genes from this pool...')
    } else {
      genes.use = gene.pool
    }
  } else if(sort.genes.by == 'r.squared'){
    r.squared.values.sorted = sort(r.squared.values, decreasing = TRUE)
    genes.use = head(names(r.squared.values.sorted), top.n.genes)
    message('Picking top ', top.n.genes, ' genes by r^2...')
  } else if(sort.genes.by == 'LRT'){
    pvals.clean.filt.sorted = sort(pvals.clean[names(pvals.clean) %in% gene.pool])
    genes.use = head(names(pvals.clean.filt.sorted), top.n.genes)
    message('Picking top ', top.n.genes, ' genes by LRT p-value...')
  }
  
  message('Mean r^2 value: ', mean(r.squared.values[genes.use]))
  message('95% CI: ', quantile(r.squared.values[genes.use], probs = c(0.025, 0.975))[[1]], 
          '-', quantile(r.squared.values[genes.use], probs = c(0.025, 0.975))[[2]])
  bin.mat = bin.mat[genes.use,]
  message('Using ', length(genes.use), ' genes across ', ncol(bin.mat), ' clusters!')
  
  # SNF binarization
  # message('Running SNF binarization...')
  # bin.mat = do.call(rbind, pbapply::pblapply(genes.use, SNFBinarize, avg.expr, pct.expr))
  # rownames(bin.mat) = genes.use
  
  # Step 5: Save binarization
  col_fun = colorRamp2(c(0, 1), c("white", "red"))
  hc_rows <- hclust(dist(bin.mat))
  row_order <- hc_rows$order
  ordered_row_names <- rownames(bin.mat)[row_order]
  ordered_row_names_rsquared = paste0(ordered_row_names, ' (', round(r.squared.values[ordered_row_names], 2),')')
  pdf(binfile, height = 110, width = 12.5)
  print(Heatmap2(t(apply(pct.expr[ordered_row_names,], 1, min_max_norm)), 
                 col = col_fun, 
                 row.font.size = 5, 
                 row_labels = ordered_row_names_rsquared))
  print(Heatmap2(t(apply(dat.use[ordered_row_names,], 1, min_max_norm)), 
                 col = col_fun, 
                 row.font.size = 5, 
                 row_labels = ordered_row_names_rsquared))
  print(Heatmap2(bin.mat[ordered_row_names,], 
                 col = col_fun, 
                 row.font.size = 5, 
                 row_labels = ordered_row_names_rsquared))
  dev.off()
  
  saveRDS(list(gene.pool = gene.pool, 
               genes.use = ordered_row_names, 
               bimodal = bimodal.genes, 
               var.genes = var.genes), 
          genefile)
  
  # Step 5.1: check mean/cv relationship
  means = rowMeans(bin.mat)
  sds = apply(bin.mat, 1, function(x) sqrt(sum((x - mean(x))^2)/(length(x))))
  datsum = data.frame(mean = means, 
                      sd = sds, 
                      cv = sds/means)
  saveRDS(datsum, datsumfile)
  
  pdf(cvfile, height = 5, width = 5)
  print(ScatterPlot(datsum, x = 'mean', y = 'cv', logX = TRUE, logY = TRUE, lm = FALSE, r_pval = FALSE)+
    stat_function(fun = function(x) 1/sqrt(x), color = "red", size = 1.2, linetype = "dashed"))
  print(PrettyHistogram(datsum$mean, xlab = 'Mean'))
  dev.off()

  # Step 6: Construct phylogeny for binary trait data using IQTree
  library(phangorn)
  phydat = phyDat(t(bin.mat), type="USER", levels=c("0", "1"), compress = FALSE)
  write.phyDat(phydat, phyfile, format = "phylip")
  if(asr) {
    args = c('-s', phyfile, '-redo', '-pre', prefix, '-nt', nthreads, '-wsr', '-m', iqtree.model, '--ancestral')
  } else {
    args = c('-s', phyfile, '-redo', '-pre', prefix, '-nt', nthreads, '-wsr', '-m', iqtree.model)
  }
  system2('iqtree2', 
          args, 
          stdout = ifelse(verbose, "", FALSE))
  
  tree = read.tree(treefile)
  # tree = root(tree, outgroup = 'MG-1')
  pdf(plotfile, height = 50, width = 30)
  plot(tree)
  edgelabels(round(tree$edge.length, 2))
  dev.off()
  
  # Step 7: rank genes by gene-specific substitution rate
  tab = read.table(paste0(prefix, '.rate'), header=TRUE)
  if(nrow(tab) > 0){ # if best model does not include rate heterogeneity, this file will be blank
    tab$gene = rownames(bin.mat)
    tab.sorted = tab %>% arrange(Rate)
    write.table(tab.sorted, ratefile)
    message('Top conserved genes: ', paste(head(tab.sorted$gene), collapse = ', '))
  }
  
  superprefix
  
  # })
}





### OLD
GenerateBinaryMatrix = function(object, species, pct.exp.thres, pval.thres, noise.sd = 2, seed = 1, 
                                top.n.genes = NULL, verbose = FALSE, nthreads = 1){
  
  # Step 0: file paths
  dir.create(paste0('celltypephylo/', species), showWarnings = FALSE)
  fafile = paste0('celltypephylo/', species, '/', species, '-', pct.exp.thres, '-', pval.thres, '-', noise.sd, '-', seed, '-', top.n.genes, '.fa')
  phyfile = gsub(".fa", ".phy", fafile)
  prefix = gsub(".fa", "", fafile)
  treefile = gsub("fa", "treefile", fafile)
  binfile = paste0(prefix, '-binarization.pdf')
  genefile = paste0(prefix, '-genes.rds')
  plotfile = paste0(prefix, '-plot.pdf')
  ratefile = paste0(prefix, '-rate-sorted.txt')
  
  # Step 1: identify genes that fluctuate a lot across clusters (high variance in terms of percentage expressed)
  message('Step 1: calculating percentage expressed for clusters...')
  pct.expr = PercentageExpressed2(DownsampleSeurat(object, group.by = 'annotated', size = 100, seed = seed), 
                                  group.by = 'annotated', features = rownames(object))
  vars = apply(pct.expr, 1, function(x) max(x) - min(x))
  # print(PrettyHistogram(vars, vline = pct.exp.thres))
  var.genes = names(which(apply(pct.expr, 1, function(x) max(x) - min(x) > pct.exp.thres)))
  
  # Step 2: identify which of these genes are bimodal in their distribution
  message('Step 2: identifying bimodal genes...')
  gmm_res = apply(pct.expr[var.genes,], 1, BimodalityTest2, seed = seed, noise.sd = noise.sd, verbose = verbose)
  gmm_res.clean = gmm_res[sapply(gmm_res, function(x) length(x[[1]]) > 0)]
  
  pvals.clean = sapply(gmm_res.clean, function(x) x[[2]])
  # pvals.clean = setNames(as.numeric(pvals[which(pvals != 'error')]), names(pvals[which(pvals != 'error')]))
  # print(PrettyHistogram(-log10(pvals.clean), vline = 2))
  
  if(!is.null(top.n.genes)) bimodal.genes = head(names(sort(pvals.clean[pvals.clean < 1e-2])), top.n.genes) else bimodal.genes = names(which(pvals.clean < 1e-2))
  
  # Step 3: binarize their expression
  message('Step 3: binarizing expression...')
  dat.use = pct.expr[bimodal.genes,]
  # bin.mat = as.data.frame(t(apply(dat.use, 1, BimodalityTest, seed = seed, alpha = 1, return.pval = FALSE, noise.sd = noise.sd, verbose = verbose))) %>% setNames(colnames(pct.expr))
  bin.mat = do.call(rbind, lapply(gmm_res.clean, function(x) t(x[[1]])))
  colnames(bin.mat) = gsub(' ', '-', colnames(pct.expr))
  rownames(bin.mat) = names(gmm_res.clean)
  bin.mat = bin.mat[bimodal.genes,]
  
  # Step 4: Compute mean variance explained and filter out poorly represented genes
  message('Step 4: checking quality of binarization...')
  r.squared.values = sapply(seq_len(nrow(bin.mat)), function(i){
    summary(lm(as.numeric(dat.use[i,]) ~ as.numeric(bin.mat[i,])))$r.squared
  }) %>% setNames(bimodal.genes)
  
  message('Mean r^2 value: ', mean(r.squared.values))
  message('95% CI: ', quantile(r.squared.values, probs = c(0.025, 0.975))[[1]], '-', quantile(r.squared.values, probs = c(0.025, 0.975))[[2]])
  genes.remove = names(which(r.squared.values < 0.3))
  message('Removing ', length(genes.remove), ' genes with r.squared < 0.3...')
  bin.mat = bin.mat[!rownames(bin.mat) %in% genes.remove,]
  genes.use = rownames(bin.mat)
  
  message('Using ', nrow(bin.mat), ' genes and ', ncol(bin.mat), ' clusters!')
  
  # Step 5: Save binarization
  col_fun = colorRamp2(c(0, 1), c("white", "red"))
  hc_rows <- hclust(dist(bin.mat))
  row_order <- hc_rows$order
  ordered_row_names <- rownames(mat)[row_order]
  pdf(binfile, height = 110, width = 12.5)
  print(Heatmap2(t(apply(pct.expr[ordered_row_names,], 1, min_max_norm)), col = col_fun, row.font.size = 5))
  print(Heatmap2(bin.mat[ordered_row_names,], col = col_fun, row.font.size = 5))
  dev.off()
  saveRDS(ordered_row_names, genefile)
  
  # Step 6: Construct phylogeny for binary trait data using IQTree
  library(phangorn)
  phydat = phyDat(t(bin.mat), type="USER", levels=c("0", "1"), compress = FALSE)
  write.phyDat(phydat, phyfile, format = "phylip")
  system2('iqtree2', 
          c('-s', phyfile, '-redo', '-pre', prefix, '-nt', nthreads, '-wsr'), 
          stdout = ifelse(verbose, "", NULL))
  
  tree = read.tree(treefile)
  # tree = root(tree, outgroup = 'MG-1')
  pdf(plotfile, height = 50, width = 30)
  plot(tree)
  edgelabels(round(tree$edge.length, 2))
  dev.off()
  
  # Step 7: rank genes by gene-specific substitution rate
  tab=read.table(paste0(prefix, '.rate'),header=TRUE)
  tab$gene = rownames(bin.mat)
  tab.sorted = tab %>% arrange(Rate)
  write.table(tab.sorted, ratefile)
  message('Top conserved genes: ', paste(head(tab.sorted$gene), collapse = ', '))
  
  tree
}

TrainLogistic = function(object, proportion = 0.6, group.by = 'seurat_clusters', 
                         size = 100, seed = 12345, cv = TRUE, assay = 'RNA',...){
  
  library(glmnet)
  
  # Classes
  class_levels = levels(object@meta.data[[group.by]])
  clusters = object@meta.data[[group.by]]
  
  # Get 60% of cluster cells
  Idents(object) = group.by
  # training_cells = unlist(sapply(seq_along(levels(clusters)), function(index) 
  #   WhichCells(object, 
  #              idents = class_levels[[index]], 
  #              downsample = floor(table(clusters)[[index]])*proportion)))
  # object_train = object[,training_cells]
  
  # Ensure that every cluster is no more than n size cells
  # object_train = DownsampleSeurat(object_train, group.by = group.by, size = size, seed = seed)
  
  # The rest are training cells
  # object_test = subset(object, cells = Cells(object_train), invert = TRUE)
  
  # Edit: 1/17/26
  object_train = GetTrainingSet(object, 
                                group.by = group.by, 
                                proportion = proportion, 
                                max.size = size)
  object_test = subset(object, cells = Cells(object_train), invert = TRUE)
  
  # Print message
  message('Training model on ', length(Cells(object_train)), ' cells')
  message('Testing model on ', length(Cells(object_test)), ' cells')
  # message("Using ", length(Cells(object_train)), " cells to train classifier")
  # message("Using ", length(Cells(object_test)), " cells to test classifier")
  
  # Prepare training and testing data
  X_train = t(object_train@assays[[assay]]@scale.data)
  stopifnot(nrow(X_train) > 0)
  # Y_train = model.matrix(~ seurat_clusters, data = object_train@meta.data)
  # colnames(Y_train) = gsub('seurat_clusters', '', colnames(Y_train))
  Y_train = object_train@meta.data[[group.by]]
  X_test = t(object_test@assays[[assay]]@scale.data)
  Y_test = object_test@meta.data[[group.by]]
  
  # Cross-validation or not? 
  if(cv){
    message("Starting cross-validation...")
    cv.fit = cv.glmnet(X_train, Y_train, family = 'multinomial', ...)
    message("Cross-validation complete!")
    lambda_min <- cv.fit$lambda.min
  } else {
    message("Training glmnet (skipping cross-validation)...")
    cv.fit = glmnet(X_train, Y_train, family = 'multinomial', ...)
    lambda_min = exp(-4)
  }
  
  message('Using lambda: ', lambda_min)
  train_predictions_class = predict(cv.fit, newx = X_train, s = lambda_min, type = "class")
  message('Training accuracy: ', mean(train_predictions_class == Y_train))
  predictions_class = predict(cv.fit, newx = X_test, s = lambda_min, type = "class")
  message('Test accuracy: ', mean(predictions_class == Y_test))
  
  cv.fit
}


# train.preds <- class_levels[apply(predict(cv.fit, newx = X_train, s = "lambda.min", type = "response"), 1, which.max)]
# train.accuracy <- mean(train.preds == clusters_train)
# message('Training accuracy:', train.accuracy)
# 
# # plot(cv.fit) 
# 
# test.preds <- class_levels[apply(predict(cv.fit, newx = X_test, s = "lambda.min", type = "response"), 1, which.max)]
# test.accuracy <- mean(test.preds == clusters_test)
# message("Test accuracy: ", test.accuracy)

RunEve = function(object, ...){
  
  # Tree
  phylo = readRDS('../../Ortho_Objects/phylo.hs.liz.rds')
  # phylo
  
  phylo$tip.label = c("Lizard", "Chicken", "Opossum", "Ferret", "Pig", "Cow", "Sheep", "Squirrel", "Peromyscus" , "Rat", "Mouse", "Rhabdomys", "TreeShrew" ,  "MouseLemur",  "Marmoset",  "Human", "Macaque")
  # ape::plot.phylo(phylo)
  
  library('evemodel')
  exprMat = GetAssayData(object)
  res <- betaSharedTest(tree = phylo, 
                        gene.data = exprMat, 
                        colSpecies = object$species, 
                        ...)
  
  # Beta change
  parameter.summary = as.data.frame(res$indivBetaRes$par)
  parameter.summary$gene = rownames(object)
  parameter.summary$beta.shared = res$sharedBeta
  parameter.summary$pval = pchisq(res$LRT, df = 1, lower.tail = F)
  parameter.summary$padj = p.adjust(parameter.summary$pval, method = 'fdr')
  parameter.summary$log10.padj = -log10(parameter.summary$padj)
  parameter.summary$beta.diff = res$indivBetaRes$par[,'beta']-res$sharedBeta
  head(parameter.summary)
  
  # Low beta: high species divergence
  parameter.summary.sorted = arrange(parameter.summary, padj)
  
  proportion.divergent.high = nrow(subset(parameter.summary.sorted, padj < 0.05 & beta.diff < 0))/nrow(parameter.summary.sorted)
  proportion.divergent.low = nrow(subset(parameter.summary.sorted, padj < 0.001 & beta.diff < 0))/nrow(parameter.summary.sorted)
  
  list(res = parameter.summary.sorted, proportion.divergent = c(proportion.divergent.high, proportion.divergent.low))
}

PrepareExpressionData = function(object, this.type){
  object = subset(object, type %in% c(this.type))
  # object = ac.ortho.present.species
  object$type = factor(object$type)
  object$species_animal = factor(paste0(object$species, '-', object$animal), levels = unique(paste0(object$species, '-', object$animal)))
  table(object$type, object$species_animal)
  quant.data = MakeSmartMatrix(table(object$type, object$species_animal))
  quant.data
  
  samples.keep = names(which(apply(quant.data@matrix, 2, function(x) all(x > 20))))
  object.filt = subset(object, species_animal %in% samples.keep)
  object.filt$species_animal = factor(object.filt$species_animal)
  object.filt
  table(object.filt$type, object.filt$species_animal)
  
  # Downsample all bio reps to 30 cells
  object.filt$species_animal_type = CombineAndRefactor(object.filt$species_animal, object.filt$type) #paste0(object.filt$species_animal, '-', object.filt$type)
  object.filt$species_animal_type_by_species = CombineAndRefactor(object.filt$species_animal, object.filt$type, reverse = TRUE)
  object.filt.ds = DownsampleSeurat(object.filt, group.by = 'species_animal_type', size = 30)
  table(object.filt.ds$type, object.filt.ds$species_animal)
  
  # Get genes
  n.species.expressed = readRDS('../../Ortho_Objects/n.species.expressed.rds')
  expressed.genes = names(which(n.species.expressed == 17))
  
  # Avg expression
  avg.expr.species.type.by.type = LogAvgExpr(object.filt.ds, group.by = 'species_animal_type', assay = 'RNA', features = expressed.genes)
  avg.expr.species.type.by.species = LogAvgExpr(object.filt.ds, group.by = 'species_animal_type_by_species', assay = 'RNA', features = expressed.genes)
  
  # Build matrix
  species.anno = str_split_fixed(colnames(avg.expr.species.type.by.type), '-', 3)[,1]
  animal.anno = paste0(str_split_fixed(colnames(avg.expr.species.type.by.type), '-', 3)[,1], '-',
                       str_split_fixed(colnames(avg.expr.species.type.by.type), '-', 3)[,2])
  type.anno = str_split_fixed(colnames(avg.expr.species.type.by.type), '-', 3)[,3]
  avg.expr.mat = SmartMatrix(as.matrix(avg.expr.species.type.by.type), 
                             row.data = data.frame(gene = rownames(avg.expr.species.type.by.type)),
                             col.data = data.frame(
                               species = species.anno,
                               technology = species_seq_method$Method[match(species.anno, species_seq_method$CommonName)],
                               type = type.anno, 
                               animal = animal.anno))
  
  avg.expr.mat
}

ComputeDivergenceScore = function(smatrix, this.type, use.median = FALSE){
  resid.var.mat = SmartMatrix(1-cor(smatrix@matrix)^2, 
                              row.data = smatrix@col.data, 
                              col.data = smatrix@col.data)
  
  # Diagonal to NA
  diag(resid.var.mat@matrix) = NA
  
  heatmap = grid.grabExpr(draw(SmartHeatmap(resid.var.mat, 
                                            annotation_names = c('species', 'animal', 'type', 'technology'), 
                                            cluster_columns = FALSE, 
                                            cluster_rows = FALSE, 
                                            rect_gp = gpar(col = NA))
                                            # legend.name = '1-R^2')
  ))
  # print(plot_grid(heatmap))
  
  # Extract residual variances
  within.rvar = na.omit(resid.var.mat@matrix[outer(resid.var.mat@row.data$species, resid.var.mat@row.data$species, function(x,y) x == y)])
  within.rvar
  across.rvar = resid.var.mat@matrix[outer(resid.var.mat@row.data$species, resid.var.mat@row.data$species, function(x,y) x != y)]
  across.rvar
  
  if(use.median) {
    stat.across = median(across.rvar)
    stat.within = median(within.rvar, na.rm = TRUE) 
  } else {
    stat.across = mean(across.rvar)
    stat.within = mean(within.rvar, na.rm = TRUE) 
  } 
  
  histogram = PrettyHistogram2(within.rvar, across.rvar, xlab = 'Residual variance', name1 = 'within-species', name2 = 'across-species', 
                         vline = c(stat.across, stat.within), title = this.type)

  # print(histogram)
  print(plot_grid(heatmap, histogram, nrow = 2, rel_heights = c(5, 1)))
  
  stat.diff = stat.across - stat.within
  
  data.frame(type = this.type, 
             stat.across = stat.across, 
             stat.across.sd = sd(across.rvar),
             stat.within = stat.within, 
             stat.within.sd = sd(within.rvar), 
             stat.diff = stat.diff,
             stat.diff.se = sqrt(sd(across.rvar)^2/length(across.rvar) + sd(within.rvar)^2/length(within.rvar)))
}

PrimateRodentAnalysis = function(object, downsample.size, features.use, type = NULL, n_genes = 8){
  
  # Filter to specific type if necessary 
  # if(!is.null(type)) object = subset(object, expr)
  # 10% cutoff calculation
  # cutoff = quantile(table(BIGSEURAT$species), probs = 0.1)
  # seurat_filt
  
  object = DownsampleSeurat(object, group.by = 'species', size = downsample.size)
  object = SubsetSeuratGenes(object, features = features.use)
  object = ScaleData(object, features = features.use)
  
  object$clade = ifelse(object$species %in% rodents, 'rodent', ifelse(object$species %in% primates, 'primate', 'other'))
  de_primate_rodent = FindMarkersFast(object, ident.1 = 'primate', ident.2 = 'rodent', group.by = 'clade', p_val_adj_cutoff = 1, avg_log2FC_cutoff = 0)
  
  # Volcano plot
  print(volcanoPlot(de_primate_rodent, labels = TRUE, max_fdr = 1e-300))
  
  # Significant
  signif.primate.genes = de_primate_rodent %>% filter(avg_log2FC > 0.5 & p_val_adj < 0.01) %>% pull(gene) 
  signif.rodent.genes = de_primate_rodent %>% filter(avg_log2FC < -0.5 & p_val_adj < 0.01) %>% pull(gene) 
  print(SeuratHeatmap(DownsampleSeurat(object, group.by = 'species', size = 100), 
                      features = c(signif.primate.genes, signif.rodent.genes), 
                      group.by = 'species', 
                      annotate.by = 'species', 
                      show_column_names = FALSE, 
                      row_names_gp = gpar(fontsize = 2)))
  
  # Top N
  top.primate.genes = de_primate_rodent %>% filter(avg_log2FC > 0) %>% head(n_genes) %>% pull(gene) 
  top.rodent.genes = de_primate_rodent %>% filter(avg_log2FC < 0) %>% head(n_genes) %>% pull(gene) 
  print(VlnPlot(object, features = c(top.primate.genes, top.rodent.genes), group.by = 'species', stack = TRUE))
  
  return(list(summary = de_primate_rodent, primate = signif.primate.genes, rodent = signif.rodent.genes))
}

InstallMyPackages = function(libraries = LIBRARIES){
  
  lapply(libraries, install_package)
  
}

RunOU = function(df, expression){
  options(warn = 2)
  
  tree <- with(df,ouchtree(node,ancestor,Time/max(Time),spcode))
  
  features = colnames(expression)
  results = do.call(rbind, lapply(features, function(gene) {
    
    # Add gene to df
    df$gene = as.numeric(expression[,gene])
    message('working on gene ', gene)
    
    # h1 <- brown(tree_edges_full_final['calb1'],tree)
    h2 <- tryCatch(hansen(df['gene'],tree,(df['OU.1']),sqrt.alpha=1,sigma=1), error = function(e) NULL)
    # h3 <- hansen(tree_edges_full_final['calb1'],tree,tree_edges_full_final['OU.2'],sqrt.alpha=1,sigma=1)
    h4 <- tryCatch(hansen(df['gene'],tree,(df['OU.3']),sqrt.alpha=1,sigma=1), error = function(e) NULL)
    h5 <- tryCatch(hansen(df['gene'],tree,(df['OU.4']),sqrt.alpha=1,sigma=1), error = function(e) NULL)
    
    # Check for non-convergence
    if(is.null(h2) | is.null(h4) | is.null(h5)){
      return(data.frame(gene = gene, 
                        logLik_null = NA, 
                        logLik_primate = NA, 
                        logLik_rodent = NA,
                        lrt_primate = NA, 
                        lrt_rodent = NA,
                        p_value_primate = NA, 
                        p_value_rodent = NA))
    }
    
    # Get the log-likelihoods of both models
    logLik_null <- logLik(h2)
    logLik_primate <- logLik(h4)
    logLik_rodent <- logLik(h5)
    
    # Perform the Likelihood Ratio Test
    lrt_primate <- 2 * (logLik_primate - logLik_null)
    lrt_rodent <- 2 * (logLik_rodent - logLik_null)
    
    # Calculate the p-value using chi-square distribution
    p_value_primate <- pchisq(lrt_primate, df = summary(h4)$dof - summary(h2)$dof, lower.tail = FALSE)
    p_value_rodent <- pchisq(lrt_rodent, df = summary(h5)$dof - summary(h2)$dof, lower.tail = FALSE)
    
    return(data.frame(gene = gene, 
                      logLik_null = logLik_null, 
                      logLik_primate = logLik_primate, 
                      logLik_rodent = logLik_rodent,
                      lrt_primate = lrt_primate, 
                      lrt_rodent = lrt_rodent,
                      p_value_primate = p_value_primate, 
                      p_value_rodent = p_value_rodent))
  }))
  
  return(results)
}

GlaucomaDETest = function(this.object, group.by, ident.1, ident.2, filename, overwrite = FALSE, top.n = 2, 
                          cols = conditions_palette6, volcano.plot = TRUE, method = 'seurat', 
                          avg_log2FC_cutoff = 0.25, p_val_adj_cutoff = 0.05, downsample = NULL){
  
  filepath = dirname(filename)
  prefix = gsub('.xlsx', '', basename(filename))
  types.use = names(which(!apply(as.data.frame.matrix(table(this.object$type, this.object@meta.data[[group.by]])), 1, function(x) any(x == 0))))
  
  # delete volcano plot directory if already exists
  dir_path = paste0(filepath, '/volcano_plots/', prefix)
  if (dir.exists(dir_path) & volcano.plot) {
    unlink(dir_path, recursive = TRUE)  # Delete the directory and its contents
  }
  
  if(overwrite | !file.exists(filename)){
    if(method == 'seurat'){
      
      de_table = do.call(rbind, lapply(types.use, function(this.type){
        
        message("Working on type: ", this.type)
          
          # tryCatch to avoid cases where there are no DEGs returned
          tryCatch({
            object = subset(this.object, type == this.type)
            if(!is.null(downsample)) object = DownsampleSeurat(object, group.by = 'sample_core', size = downsample)
            de = FindMarkersFast(object, group.by = group.by, ident.1 = ident.1, ident.2 = ident.2, p_val_adj_cutoff = 1, avg_log2FC_cutoff = 0)
            de$type = factor(this.type, levels = ALLTYPES)
            
            if(volcano.plot){
              suppressWarnings(dir.create(dir_path))
              pdf(paste0(filepath, '/volcano_plots/', prefix, '/', this.type, '.pdf'))
              print(volcanoPlot(de, max_fdr = 1e-100, labels = TRUE, max.overlaps = 20, size = 3))
              dev.off()
            }
            
            return(de %>% filter(p_val_adj <= p_val_adj_cutoff & abs(avg_log2FC) > avg_log2FC_cutoff))
          }, error = function(e) {
            print(e$message)
            return(NULL)
          })
      }))
  } else if(method == 'DESeq2'){
    library(DESeq2)
    
    pseudobulk <- AggregateExpression(this.object, assays = "RNA", slot = 'counts', return.seurat = TRUE, group.by = c(group.by, "type", 'sample_core'))
    pseudobulk[['condition']] = str_split_fixed(Cells(pseudobulk), '_', 4)[,1]
    pseudobulk[['number']] = str_split_fixed(Cells(pseudobulk), '_', 4)[,2]
    pseudobulk[['lit_type']] = str_split_fixed(Cells(pseudobulk), '_', 4)[,3]
    pseudobulk[['sample_core']] = str_split_fixed(Cells(pseudobulk), '_', 4)[,4]
    pseudobulk[['type']] = paste0(pseudobulk$number, '_', pseudobulk$lit_type)
    
    de_table = do.call(rbind, lapply(types.use, function(this.type){
      
      # tryCatch to avoid cases where there are no DEGs returned
      tryCatch({
        de = RunDESeq2(subset(pseudobulk, type == this.type), ident.1, ident.2)
        de$type = factor(this.type, levels = ALLTYPES)
        
        if(volcano.plot){
          suppressWarnings(dir.create(paste0(filepath, '/volcano_plots/', prefix)))
          pdf(paste0(filepath, '/volcano_plots/', prefix, '/', this.type, '.pdf'))
          print(volcanoPlot(de, fc_cutoff = 1, max_fdr = 1e-50, labels = TRUE, max.overlaps = 20, size = 3, max_fc = 6))
          dev.off()
        }
        
        return(de)
      }, error = function(e) {
        return(NULL)
      })
    }))
    
  } else {
    stop('method should be either seurat or DESeq2!')
  }
    de_table = subset(de_table, abs(avg_log2FC) > avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff)
    write.xlsx(de_table, file = filename) # saveRDS(de_table, filename)
  } 
  
  de_table = read.xlsx(filename) # readRDS(filename)
  de_table$type = factor(de_table$type, levels = ALLTYPES)
  
  de_up = as.data.frame(table(subset(de_table, avg_log2FC > avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff)$type)) %>% mutate(Direction = ident.1)
  de_down = as.data.frame(table(subset(de_table, avg_log2FC < avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff)$type)) %>% mutate(Direction = ident.2)
  
  p1 = ggbarplot(rbind(de_up, de_down) %>% setNames(c('Type', 'Frequency', 'Direction')), 
                 x = 'Type', y = 'Frequency', fill = 'Direction', position = position_stack()) +
    scale_fill_manual(values = cols)+
    # NoLegend()+
    RotatedAxis()
  
  # Get top genes
  de_table %>% 
    group_by(type) %>%
    arrange(-abs(avg_log2FC)) %>% 
    slice_head(n = top.n) %>%
    ungroup() -> top.genes
  
  # p2 = VlnPlot(this.object,
  #         features = top.genes$gene,
  #         group.by = 'type',
  #         split.by = group.by,
  #         pt.size = 0,
  #         stack = TRUE,
  #         flip = TRUE,
  #         cols = conditions_palette6[match(c(ident.1, ident.2), names(conditions_palette6))])
  
  # top DEG
  top_gene = (de_table %>% arrange(p_val_adj, avg_log2FC))[1,]
  p2 = TitlePlot(
       VlnPlot2(this.object,
                features = top_gene$gene,
                split.by = 'sample_core',
                group.by = group.by,
                idents = top_gene$type,
                cols = conditions_palette6,
                stack = FALSE,
                combine= TRUE,
                pt.size = 0.1), 
       title = paste0(top_gene$gene, ' in ', top_gene$type))

  return(list(p1, p2))
  
}

SplitJaccardMatrix = function(object, group.by = 'species', ident.1 = 'seurat_clusters', ident.2 = 'type', nrow = 1, ncol = NULL, 
                              return.list = FALSE, widths = 1, remove.y.axis = FALSE, FUN = NULL, remove.x.axis = FALSE, ari = FALSE, 
                              width.adjust = 10, ...){
  
  objectList = SplitObject(object, split.by = group.by)
  
  # Order as original 
  if(is.factor(object@meta.data[[group.by]])) objectList = objectList[levels(object@meta.data[[group.by]])]
  
  heatmapList = lapply(seq_along(objectList), function(index) {
    
    if(ari){
      JSHeatmap2(objectList[[index]]@meta.data[[ident.1]], objectList[[index]]@meta.data[[ident.2]], 
                 title = paste0(names(objectList)[[index]], '\nARI = ', 
                                round(adj.rand.index(objectList[[index]]@meta.data[[ident.1]], objectList[[index]]@meta.data[[ident.2]]), 2), ''),
                 ...) + 
        theme(plot.background = element_rect(fill = "transparent", colour = NA)) # make transparent so doesn't block titles
    } else {
      JSHeatmap2(objectList[[index]]@meta.data[[ident.1]], objectList[[index]]@meta.data[[ident.2]], 
                 title = names(objectList)[[index]], 
                 ...) + 
        theme(plot.background = element_rect(fill = "transparent", colour = NA))
    }
    
  })
  
  # Remove y axis text? 
  if(remove.y.axis){
    heatmapList[2:length(heatmapList)] = lapply(heatmapList[2:length(heatmapList)], function(plt) plt + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()))
  }
  
  # Remove x axis text? 
  if(remove.x.axis){
    heatmapList = lapply(heatmapList, function(plt) plt + 
                           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  }
  
  # Assign names
  names(heatmapList) = names(objectList)
  
  if(!is.null(FUN)){
    heatmapList = lapply(heatmapList, function(plt) plt + FUN)
  }
  
  # Number of columns
  if(is.null(ncol)) ncol = length(objectList)
  if(return.list) return(heatmapList)
  
  #print(length(heatmapList))
  
  if(length(widths) > 1){
  # Add small pseudocount to widths to account for y-axis names
    widths[[1]] = widths[[1]] + width.adjust
  }

  # Plot
  ggarrange(plotlist = heatmapList,
            ncol = ncol,
            nrow = nrow,
            common.legend = TRUE,
            legend = "right",
            align = 'h',
            widths = widths)
  
  # Patchwork solution
  # plot_grid(plotlist = heatmapList, 
  #           rel_widths = widths)
  # Reduce(`+`, heatmapList) +
  #   plot_layout(ncol = ncol, 
  #               nrow = nrow,
  #               widths = widths,
  #               guides = "collect") &
  #   theme(legend.position = "right")
  
  # wrap_plots(heatmapList, 
  #            ncol = ncol, 
  #            nrow = nrow,
  #            widths = widths) 
    # plot_layout(guides = "collect")
  
  # wrap_plots(heatmapList, 
  #            ncol = ncol, 
  #            nrow = nrow,
  #            widths = widths) +
  #   plot_layout(guides = "collect") &  # Collect legends into one
  #   theme(legend.position = "right", plot.title = element_text())   # Place legend on right

}

BinaryTree = function(genes = NULL, masks = NULL, return.masks = FALSE, ...){
  
  if(!is.null(genes)){
    full_masks = lapply(matrix.list[genes], BinarizeExpression, return.vector = FALSE)
    names(full_masks) = genes
    masks = lapply(full_masks, function(full_mask) {
      new = apply(full_mask, 2, mean)
      new[new > 0.5] = 1
      new[new <= 0.5] = 0
      return(new)
    })
    names(masks) = genes
    
    # Remove failures and TFs with no variance
    variance = sapply(masks, var)
    filtered.masks = masks[!is.na(masks) & variance != 0]
    filtered.full_masks = full_masks[!is.na(masks) & variance != 0]
    
    if(return.masks) return(list(filtered.full_masks, filtered.masks))
    
  } else if(!is.null(masks)){
    filtered.masks = masks
  } else{
    stop("One of genes or masks must be defined!")
  }
  
  nontf.expr.mask = do.call(rbind, filtered.masks)
  colnames(nontf.expr.mask) = colnames(matrix.list[[1]])
  nontf.order = rownames(nontf.expr.mask)
  nontf.expr.mask = nontf.expr.mask[nontf.order, ]
  colnames(nontf.expr.mask) = Metadata(bc_ortho, feature.1 = "seurat_clusters", feature.2 = "NOG")$NOG
  bc.orthotypes = colnames(nontf.expr.mask)
  
  # Heatmap(t(nontf.expr.mask), 
  #         rect_gp = gpar(col = "grey", lwd = 1), col = c("white", "black"), 
  #         show_row_dend = TRUE, show_column_dend = TRUE, 
  #         clustering_distance_rows = "manhattan", 
  #         clustering_method_rows = "average", 
  #         clustering_distance_columns = "manhattan", 
  #         clustering_method_columns = "average")
  
  nontf.trees = ConstructTrees(t(nontf.expr.mask), ...)
  
  return(nontf.trees)
  
}
  
ConstructTree = function(matrix, bootstrap = TRUE, seed = 42, root = NULL, nIter = 1000, mc.cores = 16){
  library(phangorn)
  library(TreeTools)
  
  # Set seed
  set.seed(seed)
  
  # Make phyDat object
  phydat = phyDat(matrix, type="USER", levels=c("0", "1"), compress = FALSE)
  
  # Distance methods
  treeUPGMA <- upgma(dist(PhyDatToMatrix(phydat), method = "manhattan"))
  treeNJ <- NJ(dist(PhyDatToMatrix(phydat), method = "manhattan"))
  
  # Maximum parsimony
  treeMP <- pratchet(phydat, trace = 0, minit=100)
  
  # Assign branch lengths
  treeMP  <- acctran(treeMP, phydat)
  treeMP  <- di2multi(treeMP, tol = 0) # Keep as binary tree
  if(inherits(treeMP, "multiPhylo")){
    treeMP <- unique(treeMP)
  }
  
  # Labeling internal nodes
  treeMP$node.label = (length(treeMP$tip.label)+1):((length(treeMP$tip.label)+treeMP$Nnode))
  
  if(!is.null(root)) {
    treeMP = root(phy = treeMP, node = root)
    # treeMP = midpoint(treeMP)
  } 
  
  return(list(phydat = phydat, UPGMA = treeUPGMA, NJ = treeNJ, MP = treeMP))
}

BootstrapTrees = function(matrix, method = 'binary', seed = 42, nIter = 100, mc.cores = 1, 
                          run = c('UPGMA', 'NJ', 'MP')){
  
  # Set seed
  set.seed(seed)
  
  bs_upgma = NULL
  bs_nj = NULL
  bs_mp = NULL
  bs_ml = NULL
  
  message('UPGMA bootstrap...')
  bs_upgma = as.multiPhylo(lapply(seq_len(nIter), function(iter) 
    upgma(dist(matrix[,sample(1:ncol(matrix), replace = TRUE)], 
               method = method))))
  
  message('NJ bootstrap...')
  bs_nj = as.multiPhylo(lapply(seq_len(nIter), function(iter) 
    NJ(dist(matrix[,sample(1:ncol(matrix), replace = TRUE)], 
            method = method))))
  
  message('MP bootstrap...')
  bs_mp = as.multiPhylo(mclapply(seq_len(nIter), function(iter) {
    phydat = phyDat(matrix[,sample(1:ncol(matrix), replace = TRUE)], 
                    type="USER", levels=c("0", "1"), compress = FALSE)
    pratchet(phydat, trace = 0, minit = 100)
  }, mc.cores = mc.cores))
  
  if('ML' %in% run){
    message('ML bootstrap...')
    bs_ml = as.multiPhylo(mclapply(seq_len(nIter), function(iter) {
      phydat = phyDat(matrix[,sample(1:ncol(matrix), replace = TRUE)], 
                      type="USER", levels=c("0", "1"), compress = FALSE)
      RunIQTree2(phydat)
    }, mc.cores = mc.cores))
  } 
  
  return(list(UPGMA = bs_upgma, NJ = bs_nj, MP = bs_mp, ML = bs_ml))
}

BinaryPhyDat = function(matrix){
  phyDat(matrix, type="USER", levels=c("0", "1"), compress = FALSE)
}

ConstructTrees = function(matrix, method = 'binary', bootstrap = FALSE, seed = 42, root = NULL, nIter = 1000, mc.cores = 1, consensus = TRUE, minit = 1000){
  
  library(phangorn)
  library(TreeTools)
  
  # Set seed
  set.seed(seed)
  
  # Make phyDat object
  phydat = phyDat(matrix, type="USER", levels=c("0", "1"), compress = FALSE)
  
  # Distance methods
  treeUPGMA <- upgma(dist(PhyDatToMatrix(phydat), method = method))
  treeNJ <- NJ(dist(PhyDatToMatrix(phydat), method = method))
  
  # Maximum parsimony
  treeMPs <- pratchet(phydat, trace = 0, minit = minit, all = TRUE)
  
  # Check for multiple trees
  if(inherits(treeMPs, 'multiPhylo')){
    
    message('Found ', length(treeMPs), ' equally good parsimony trees!')
    
    # Make consensus tree out of all the equal parsimony trees?
    if(consensus){
      
      message('Making consensus tree for MP...')
      treeMP <- consensus(treeMPs, p = 0.5, check.labels = TRUE)
      
    } else {
      
      # Get the first tree
      treeMP = treeMPs[[1]]
      
    }
    
  } else {
    
    # Get the only tree
    treeMP = treeMPs
    
  }
  
  # Assign branch lengths
  if(!consensus){
    treeMP  <- acctran(treeMP, phydat)
    treeMP  <- di2multi(treeMP, tol = 0) # Keep as binary tree
  }
  
  if(inherits(treeMP, "multiPhylo")){
    treeMP <- unique(treeMP)
  }
  
  # Labeling internal nodes
  treeMP$node.label = (length(treeMP$tip.label)+1):((length(treeMP$tip.label)+treeMP$Nnode))
  
  if(!is.null(root)) {
    treeMP = root(phy = treeMP, node = root)
    # treeMP = midpoint(treeMP)
  } 
  
  # ML tree 
  treeML = RunIQTree2(phydat)
  
  # Bootstrap
  if(bootstrap){
    
    message('UPGMA bootstrap...')
    bs_upgma = as.multiPhylo(lapply(seq_len(nIter), function(iter) 
      upgma(dist(PhyDatToMatrix(phydat)[,sample(1:ncol(matrix), replace = TRUE)], 
                 method = method))))
    treeUPGMA = plotBS(treeUPGMA, bs_upgma, main="UPGMA", type = "none")
    treeUPGMA$bs.support = treeUPGMA$node.label
    treeUPGMA$bs.support[is.na(treeUPGMA$bs.support)] = 0
    
    message('NJ bootstrap...')
    bs_nj = as.multiPhylo(lapply(seq_len(nIter), function(iter) 
      NJ(dist(PhyDatToMatrix(phydat)[,sample(1:ncol(matrix), replace = TRUE)], 
              method = method))))
    treeNJ = plotBS(treeNJ, bs_nj, main="NJ", type = "none")
    treeNJ$bs.support = treeNJ$node.label
    treeNJ$bs.support[is.na(treeNJ$bs.support)] = 0
    
    message('MP bootstrap...')
    bs_mp = as.multiPhylo(mclapply(seq_len(nIter), function(iter) {
      phydat = phyDat(PhyDatToMatrix(phydat)[,sample(1:ncol(matrix), replace = TRUE)], 
                      type="USER", levels=c("0", "1"), compress = FALSE)
      pratchet(phydat, trace = 0, minit = 100)
    }, mc.cores = mc.cores))
    treeBS = plotBS(treeMP, bs_mp, main="MP", type = "none")
    treeMP$bs.support = treeBS$node.label
    treeMP$bs.support[is.na(treeMP$bs.support)] = 0
    
    # message('ML bootstrap...')
    # bs_ml = as.multiPhylo(mclapply(seq_len(nIter), function(iter) {
    #   phydat = phyDat(PhyDatToMatrix(phydat)[,sample(1:ncol(matrix), replace = TRUE)], 
    #                   type="USER", levels=c("0", "1"), compress = FALSE)
    #   RunIQTree2(phydat)
    # }, mc.cores = mc.cores))
    # treeBS = plotBS(treeML, bs_ml, main="ML", type = "none")
    # treeML$bs.support = treeBS$node.label
    # treeML$bs.support[is.na(treeML$bs.support)] = 0
    
  }
  
  # Label internal nodes
  treeUPGMA$node.labels = seq(length(treeUPGMA$tip.label) + 1, length(treeUPGMA$tip.label) + treeUPGMA$Nnode)
  treeNJ$node.labels = seq(length(treeNJ$tip.label) + 1, length(treeNJ$tip.label) + treeNJ$Nnode)
  treeMP$node.labels = seq(length(treeMP$tip.label) + 1, length(treeMP$tip.label) + treeMP$Nnode)
  treeML$node.labels = seq(length(treeML$tip.label) + 1, length(treeML$tip.label) + treeML$Nnode)
  
  return(list(phydat = phydat, 
              UPGMA = treeUPGMA, 
              NJ = treeNJ, 
              MP = treeMP, 
              ML = treeML))
}

RunSoupx = function(object, feature.low = 800){
  
  # Find preliminary clusters for SoupX cleaning
  cells=ClusterSeurat(subset(object, nFeature_RNA >= feature.low))
  droplet_matrix = subset(object, nFeature_RNA < feature.low)@assays$RNA@counts
  
  # Make SoupChannel object
  sc = SoupChannel(droplet_matrix, # droplets
                   cells@assays$RNA@counts) # cells
  
  # Set clusters
  sc = setClusters(sc, cells$RNA_snn_res.0.5)
  
  # Automated method
  par(mfrow=c(1,3))
  sc = autoEstCont(sc)
  
  # Evidence of ambient RNA contamination in 1363, setting a higher contamination rate
  # sc.1363 = setContaminationFraction(sc.1363, 0.2)
  
  # save adjustment
  adjust = adjustCounts(sc)
  
  # Quantify corrected counts per cell
  cells$soupx.correction = colSums(cells@assays$RNA@counts-adjust)
  
  # Sanity check: verify that the appropriate fractions of counts were adjusted
  message(sum(cells$soupx.correction)/sum(cells@assays$RNA@counts))
  
  # Adjust counts
  cells@assays$RNA@counts=adjust
  
  # Combine
  cells = ClusterSeurat(cells)
  
  # Check SoupX correction
  umap.line = DimPlot(cells, reduction = "umap", group.by = "orig.ident")
  umap.soupx = FeaturePlot(cells, reduction = "umap", feature="soupx.correction")
  # df = data.frame(corrected.counts = cells$soupx.correction, orig.ident = cells$orig.ident)
  # violin = ggplot(df, aes(orig.ident, corrected.counts)) + geom_violin() + theme_bw()
  violin = VlnPlot(cells, group.by = 'orig.ident', features = 'soupx.correction')
  
  print(umap.line | umap.soupx | violin)
  
  return(cells)
}

OrthotypeHeatmaps2 = function(object){
  
  speciesACList = lapply(unique(object$species), function(currentSpecies) subset(object, species == currentSpecies))
  names(speciesACList) = unique(object$species)
  
  # Generate jaccard matrices
  jaccardList = lapply(seq_along(speciesACList), function(index) JSMatrix(t(table(speciesACList[[index]]$annotated, speciesACList[[index]]$seurat_clusters))))
  jaccardList = lapply(jaccardList, as.data.frame.matrix)
  
  # Get the best match in each species 
  bestMatch = lapply(jaccardList, function(matrix) apply(matrix, 1, max))
  bestMatchDf = as.data.frame(t(do.call(rbind, bestMatch)))
  colnames(bestMatchDf) = unique(object$species)
  bestMatchDf$OrthoType = rownames(bestMatchDf)
  
  # Melt and plot
  melted = reshape2::melt(bestMatchDf)
  
  ## Confusion matrices, sorted by OrthoType conservation
  OT_order = levels(reorder(melted$OrthoType, melted$value, FUN = median, decreasing = TRUE))
  heatmapList = lapply(seq_along(speciesACList), function(index) {
    JSHeatmap(JSMatrix(t(table(speciesACList[[index]]$annotated, speciesACList[[index]]$seurat_clusters))), 
              heatmap = TRUE, 
              title = unique(object$species)[[index]], 
              row.order = OT_order, 
              stagger.threshold = 0.25)
  })
  
  # Plot
  ggarrange(plotlist = heatmapList, 
            ncol = length(heatmapList), 
            nrow = 1, 
            common.legend = TRUE, 
            legend = "right")
}


OrthotypeAnalysis2 = function(reference, 
                              objectList, 
                              homologyList, 
                              group.by = "annotated",
                              types.use = NULL, 
                              downsample = 100, 
                              cluster_resolution = 0.4, 
                              sample.basis = NULL, # for sampling the features
                              use = c('seurat', 'harmony'),
                              ...){
  
  # Assign species names
  # objectList = dapply(names(objectList), function(name) {
  #   objectList[[name]]$species = name
  #   objectList[[name]]
  # })
  
  # First dataset must be reference
  # reference = objectList[[1]]
  
  # Subset homologyList to genes present in each object
  homologyList = lapply(seq_along(homologyList), function(index){
    homology = homologyList[[index]][,intersect(colnames(homologyList[[index]]), rownames(objectList[[index]]))]
    return(homology[rowSums(homology) > 0,])
  })
  
  # Generate a basis (1:1 orthologs across all species)
  basis = Reduce(intersect, lapply(homologyList, function(x) rownames(x)))
  if(!is.null(sample.basis)) basis = sample(basis, size = sample.basis)
  message("Using a basis of ", length(basis), " genes!")
  
  # Subset to genes and types in ref
  if(length(types.use) > 0){
    seurat1 = DownsampleSeurat(subset(reference, annotated %in% types.use), group.by = group.by, size = downsample)
  } else {
    seurat1 = DownsampleSeurat(reference, group.by = group.by, size = downsample)
  }
  seurat1 = SubsetSeuratGenes(seurat1, basis)
  
  # Subset to genes and types in each species
  objectList = lapply(seq_along(objectList), function(index){
    
    # Downsample object
    if(length(types.use) > 0){
      seurat2 = DownsampleSeurat(subset(objectList[[index]], annotated %in% types.use), group.by = group.by, size = downsample)
    } else {
      seurat2 = DownsampleSeurat(objectList[[index]], group.by = group.by, size = downsample)
    }
    
    # Transform expression
    homology = homologyList[[index]][basis,]
    seurat2.converted = CreateSeuratObject(homology %*% seurat2@assays$RNA@counts[colnames(homologyList[[index]]),])
    
    # Transfer metadata
    seurat2.converted = TransferMetadata(seurat2, seurat2.converted)
    
    return(seurat2.converted)
  })
  
  # Merge
  ortho = merge(seurat1, objectList)
  
  message('Removing objects to save space...')
  rm(objectList)
  rm(seurat1)
  rm(homologyList)
  gc()
  
  # Without integration
  # ortho = ClusterSeurat(ortho)
  # print(DimPlot(ortho, group.by = c("species", "annotated")))
  
  # With integration
  use = match.arg(use)
  if(use == 'seurat'){
    ortho = ClusterSeurat(ortho, integrate.by = "species", cluster_resolution = cluster_resolution, ...)
  } else if(use == 'harmony'){
    ortho = Harmonize(ortho, batch = 'species', cluster_resolution = cluster_resolution, ...)
  } 
  
  # print(DimPlot(ortho, group.by = "species") | DimPlot(ortho, label = TRUE) + NoLegend() | DimPlot(ortho, group.by = "annotated", label = TRUE) + NoLegend())
  # print(OrthotypeHeatmaps(ortho))
  
  return(ortho)
}

OrthotypeHeatmaps = function(object, stagger.threshold = 0.2){
  
  OrthoObjectList = SplitObject(object, split.by = "species")
  heatmapList = lapply(seq_along(OrthoObjectList), function(index) JSHeatmap(JSMatrix((table(OrthoObjectList[[index]]$seurat_clusters, 
                                                                                             OrthoObjectList[[index]]$annotated))), 
                                                                             heatmap = TRUE, 
                                                                             title = names(OrthoObjectList)[[index]], 
                                                                             stagger.threshold = stagger.threshold))
  
  # Plot
  ggarrange(plotlist = heatmapList, 
            ncol = length(OrthoObjectList), 
            nrow = 1, 
            common.legend = TRUE, 
            legend = "right", 
            align = "h")
}

OrthotypeAnalysis = function(objectList, orthology_key, downsample = 200, types.remove = NULL){
  
  # Generate the ortholog seurat object for each species
  speciesOrthoList = lapply(objectList, function(object) OrthologSeurat(object, 
                                                                        orthology_key = orthology_key, 
                                                                        common.genes = TRUE, 
                                                                        mart_filepath = "../../Orthology/martRefChicken.csv", 
                                                                        reference_species = "Chicken"))
  
  # Subset each matrix to common genes
  common.genes = Reduce(intersect, sapply(speciesOrthoList, function(x) rownames(x)))
  message("Found ", length(common.genes), " common genes")
  speciesOrthoList = lapply(speciesOrthoList, function(x) SubsetSeuratGenes(x, features = common.genes))
  
  # Downsample each species type to equal number of cells
  speciesSubsetList = lapply(speciesOrthoList, function(object) {
    Idents(object) = "annotated"
    object = subset(object, cells = WhichCells(object, downsample = downsample, seed = 12345))
    return(object)
  })
  
  # Merge
  OrthoObject <- merge(speciesSubsetList[[1]], y = speciesSubsetList[2:length(speciesSubsetList)])
  
  # Keep features in at least 3 cells
  # OrthoObject <- CreateSeuratObject(GetAssayData(OrthoObject), min.cells = 3)
  
  # Add species name
  # OrthoObject$species = rep(c("Chicken", "Lizard", "Opossum"), unlist(lapply(speciesSubsetList, ncol)))
  
  # Add celltype information
  OrthoObject$species_class = as.character(unlist(lapply(speciesSubsetList, function(speciesAC) speciesAC$cell_class)))
  OrthoObject$species_cluster = as.character(unlist(lapply(speciesSubsetList, function(speciesAC) speciesAC$annotated)))
  
  # Remove certain types
  if(!is.null(types.remove)) {
    OrthoObject = subset(OrthoObject, species_cluster %in% types.remove, invert = TRUE)
    message("Removing ", paste0(types.remove, collapse = ", "))
  }
  
  message("Using the following types: ")
  print(table(OrthoObject$species_class))
  print(table(OrthoObject$species_class, OrthoObject$species))
  
  # Save object
  # saveRDS(OrthoObject, output_file_v1)
  
  # Tabulate species
  message("# of cells from each species: ")
  print(t(t(table(OrthoObject[["species"]]))))
  
  # Integrate by species
  OrthoObject = ClusterSeurat(OrthoObject, integrate.by = "species", cluster_resolution = 0.5)
  
  return(OrthoObject)
}

ClusterEnrichmentComparison = function(object, cell.threshold = 1000){
  marker = ifelse(any(grep("NEUN", object$enrichment, ignore.case = T)), "RBFOX3", ifelse(grep("CD90", object$enrichment, ignore.case = T), "THY1", stop("Couldn't find marker")))
  
  # Save original clustering 
  object$orig.clusters = object$seurat_clusters
  objectList = SplitObject(object, split.by = "enrichment")
  
  # Include only enrichment batches with more than N cells
  objectList = objectList[sapply(objectList, ncol) > cell.threshold]
  
  # Downsample to same number of cells
  min_cells = min(sapply(objectList, function(object) length(Cells(object))))
  message(paste0("Subsetting down to ", min_cells, " cells!"))
  downsampledList = lapply(objectList, function(object) object[,sample(colnames(object), min_cells, replace = FALSE)])
  
  # Cluster each one separately
  downsampledList = lapply(downsampledList, function(obj) {
    if(length(unique(obj$animal)) > 1) {
      Harmonize(obj, batch = "animal", cluster_resolution = 1.5, run.umap = FALSE, show.plots = FALSE)
    } else {
      ClusterSeurat(obj, cluster_resolution = 1.5, do.umap = FALSE)
    }
  })
  
  # Order based on marker expression
  marker.order = rev(names(sort(as.matrix(AvgExpr(object, features = marker, assay = "RNA", group.by = "seurat_clusters"))[1,])))
  
  # Comparison to original clustering
  p.list = lapply(seq_along(downsampledList), function(index){
    object = downsampledList[[index]]
    name = names(downsampledList)[index]
    
    # How many original clusters were retrieved above given threshold
    stats = OverlapStatistics(table(object$orig.clusters, object$seurat_clusters))
    nRetrieved = length(unique(subset(stats, log.p >= 10 & overlap >= 30)$ident1))
    message("Retrieved ", nRetrieved, " clusters from ", name)
    
    JSHeatmap(JSMatrix((table(object$orig.clusters, object$seurat_clusters))), 
              heatmap = TRUE, 
              title = paste0(name, " subset"),
              row.order = marker.order,
              stagger.threshold = 0.25) + NoLegend()
  })
  
  
  object$ordered = factor(object$seurat_clusters, levels = rev(marker.order))
  
  return(plot_grid(plotlist = c(p.list,
                         list(VlnPlot(object, group.by = "ordered", features = marker, pt.size = 0) + coord_flip() + NoLegend() + NoAxes())), 
            ncol = length(p.list)+1, align = "h", axis = "bt",  rel_widths = c(rep(1, length(p.list)),0.4)))
  
}

ModifyRetinalAtlases = function(species){
  
  # Initial filepath
  initial_filepath = paste0("../../Species_Objects/", species, "_initial.rds")
  
  # Read in full object
  species_initial = readRDS(initial_filepath)
  
  # Print object to show number of cells and features
  message(paste0(capture.output(species_initial), collapse = "\n"))
  
  if(species %in% c("Macaque", "Marmoset", "Peromyscus", "Squirrel",
                    "Opossum", "Cow", "Sheep", "Lizard", "Mouse", 
                    "Pig", "Chicken", "Ferret", "Zebrafish", 
                    "MouseLemur", "Rat", "Mouse", "Lamprey")){
    species_initial = ConvertGeneSymbols(species_initial, "../../Orthology/martMergeRefHuman.txt", species) 
    
    # Check that original feature names were saved correctly and that no genes were lost
    stopifnot(length(rownames(species_initial@assays$RNA@counts)) == length(species_initial@misc$orig.features))
    
  } else if(species %in% c("Goldfish")){
    # goldfish_key = fread("../../Orthology/ZF_LA_SB_symbol_for_Seurat.txt")
    # new_genes = goldfish_key$V3[match(rownames(object), goldfish_key$V7)]
    # object@misc$orig.features = rownames(object@assays$RNA@counts)
    # rownames(object@assays$RNA@counts) = new_genes
    # rownames(object@assays$RNA@data) = new_genes
    
    # Goldfish has whole genome duplication so all symbols are duplicated
    species_initial = ConvertGeneSymbols(species_initial, "../../Orthology/martMergeRefHuman.txt", species, make.unique = FALSE) 
    metadata = species_initial@meta.data
    
    # Aggregate duplicated genes
    mydat = as.data.frame(as.matrix(species_initial@assays$RNA@counts))
    mydat$gene = rownames(species_initial)
    mydat.sum <- aggregate(. ~ gene, data = mydat, sum)
    
    # Check that it worked
    mydat[mydat$gene == "rbfox1",30:40]
    mydat.sum[mydat.sum$gene == "rbfox1",30:40]
    
    # Transfer metadata
    rownames(mydat.sum) = mydat.sum$gene
    species_initial2 <- CreateSeuratObject(mydat.sum[,colnames(mydat.sum) != "gene"])
    species_initial2 = TransferMetadata(from = species_initial, to = species_initial2)
    species_initial = species_initial2
    rm(species_initial2)
    
    # species_initial$orig.file = species_initial$orig.ident
    
    # Process for downstream analysis
    species_initial = ClusterSeurat(species_initial, cluster_resolution = 0.5)
    
    # Change meta features or will cause issues downstream
    # species_initial[["RNA"]]@meta.features <- data.frame(row.names = rownames(species_initial[["RNA"]]))
  }
  
  if(species %in% c("Macaque")) species_initial$seurat_clusters = species_initial$annotated
  Idents(species_initial) = "seurat_clusters"
  species_initial = UpperCase_genes(species_initial)
  # DotPlot(species_initial, features = rev(c("PAX6", "TFAP2A", "TFAP2B", "CHAT", "SLC5A7"))) + coord_flip() + RotatedAxis()
  
  # Recover mislabeled amacrines
  if(species == "Human") species_initial$cell_class[species_initial$seurat_clusters %in% c(5, 14, 24)] = "AC"
  if(species == "Marmoset") species_initial$cell_class[species_initial$seurat_clusters %in% c(27)] = "Rod"
  if(species == "Peromyscus") species_initial$cell_class[species_initial$seurat_clusters %in% c(14,48)] = "AC"
  if(species == "Cow") species_initial$cell_class[species_initial$seurat_clusters %in% c(8,9,20)] = "AC"
  if(species == "Opossum") species_initial$cell_class[species_initial$seurat_clusters %in% c(1,29,30)] = "AC"
  if(species == "Lizard") species_initial$cell_class[species_initial$seurat_clusters %in% c(18)] = "AC"
  if(species == "Zebrafish") species_initial$cell_class[species_initial$seurat_clusters %in% c(9,18,24)] = "AC"
  
  # Fails for chicken and zebrafish because seurat_clusters are not always in the same cell class
  if(species == "Chicken"){
    
    # Convert to major cell class
    species_initial$cell_class = as.character(species_initial$cell_class)
    species_initial$cell_class[species_initial$cell_class == "GabaAC" | species_initial$cell_class == "GlyAC"] = "AC"
    species_initial$cell_class[species_initial$cell_class == "BP"] = "BC"
    # species_initial$cell_class[!species_initial$cell_class %in% Annotation(major_annotation)] = "Other"
    species_initial$cell_class[species_initial$cell_class == "MicroG"] = "Other"
    
    # species_initial$cell_class = ifelse(startsWith(as.character(species_initial$annotated), "AC-"), "AC", "Other")
    DotPlot(species_initial, features = rev(Genes(major_annotation))) + coord_flip()
  } else if(species %in% c("Zebrafish", "Goldfish")){
    DotPlot(species_initial, features = rev(Genes(major_annotation))) + coord_flip()
  } else {
    
    # Convert to major cell class
    species_initial$cell_class = as.character(species_initial$cell_class)
    species_initial$cell_class[species_initial$cell_class == "GabaAC" | species_initial$cell_class == "GlyAC"] = "AC"
    species_initial$cell_class[species_initial$cell_class == "BP"] = "BC"
    # species_initial$cell_class[!species_initial$cell_class %in% Annotation(major_annotation)] = "Other"
    species_initial$cell_class[species_initial$cell_class == "MicroG"] = "Other"
    
    # Check if cones and rods are in same cluster
    if(max(apply(table(species_initial$cell_class, species_initial$seurat_clusters), 2, function(x) length(which(x != 0)))) > 1){
      species_initial$cell_class[species_initial$cell_class == "Rod" | species_initial$cell_class == "Cone"] = "PR"
    }
    
    species_initial$cell_class = factor(species_initial$cell_class, levels = unique(major_annotation@annotation))
    
    # Default palette for DimPlot
    # palette = c("cyan", "chartreuse2", "red", "gold", "magenta", "grey", "blue", "darkred")
    # if(!"Cone" %in% unique(species_initial$cell_class)) palette = c("cyan", "chartreuse2", "gold", "magenta", "grey", "blue", "darkred")
    
    # Plot
    print(plot_grid(
    {if(species != "Macaque") DimPlot(species_initial, group.by = "cell_class", cols = myPalette(7))}, # c("darkred", "red", "gold", "chartreuse2", "cyan", "blue", "magenta", "grey")
    AnnotatedDotPlot(species_initial, major_annotation_custom, group.by = "seurat_clusters", color.clusters.by = "cell_class",
                     features = rev(Genes(major_annotation_custom)), color.genes = FALSE) + coord_flip() + RotatedAxis(),
    rel_widths = c(1,2.5)))
  }
  
  return(species_initial)
}

RunSubsampling = function(object, nPermutations = 50, fraction.use = 0.8, shuffle.types = FALSE, nCores = 1, method, group.by){
  
  # Sample cells
  subsample.size = floor(length(Cells(object))*fraction.use)
  rand.cells.use = lapply(seq_len(nPermutations), function(x) sample(colnames(object), size = subsample.size, replace=F))
  
  # Sample parameters
  rand.nPCs = sample(c(15:30), nPermutations, replace = TRUE)
  rand.k = sample(c(16:24), nPermutations, replace = TRUE)
  rand.res = sample(seq(1, 2, by = 0.1), nPermutations, replace = TRUE)
  
  clusters = lapply(seq_len(nPermutations), function(iteration){
    
    message("Running iteration #", iteration, " with the following params: \n nPCs: ", rand.nPCs[[iteration]], "\n k.param: ", rand.k[[iteration]], "\n res: ", rand.res[[iteration]])
    
    sample = suppressMessages(SubsampleCluster(object, 
                                               cells.use = rand.cells.use[[iteration]],
                                               nPCs = rand.nPCs[[iteration]], 
                                               k.param = rand.k[[iteration]], 
                                               resolution = rand.res[[iteration]], 
                                               shuffle.types = shuffle.types, 
                                               method = method, 
                                               group.by = group.by))
    
    return(list(original = sample$old_seurat_clusters, permuted = sample$seurat_clusters))
  })#, mc.cores = nCores)
  
  # Name by permutation
  names(clusters) = seq_len(nPermutations)
  
  return(list(clusters = clusters, 
              rand.cells.use = rand.cells.use, 
              rand.nPCs = rand.nPCs, 
              rand.k = rand.k, 
              rand.res = rand.res))
}

SubsampleCluster = function(seurat, cells.use, nPCs, k.param, resolution, recompute.var.genes = FALSE, shuffle.types = FALSE, method = "seurat", group.by = "animal") {
  
  # Subsample
  message("Subsampling to ", length(cells.use), " cells!")
  subsample = seurat[, cells.use]
  
  # Re-process
  if(method == "harmony"){
    subsample = Harmonize(subsample, 
                          batch = group.by, 
                          nPCs = nPCs, 
                          k.param = k.param, 
                          cluster_resolution = resolution, 
                          show.plots = FALSE, 
                          run.umap = FALSE)
  } else if(method == "seurat"){
    subsample = ReprocessIntegrated(subsample, 
                                    nPCs = nPCs,
                                    k.param = k.param,
                                    cluster_resolution = resolution,
                                    recompute.var.genes = recompute.var.genes,
                                    run.umap = FALSE,
                                    verbose = FALSE,
                                    method = method)
  }
  
  return(subsample)
}

# Species-level amacrine cell analysis
render_report = function(species, 
                         initial = TRUE, 
                         batch_int = FALSE, 
                         integrate_by = "animal",
                         harmony = FALSE, 
                         contamination_threshold = 6, 
                         contamination = NULL, 
                         nFeature_threshold = -2, 
                         doublet_finder = FALSE, 
                         manual_annotation = NULL, 
                         sac_annotation = FALSE, 
                         de_expression = TRUE, 
                         save = TRUE, 
                         output_file = NULL){
  
  if(!is.null(output_file)) {
    output_file = output_file
  } else {
    output_file = paste0("html_reports/Report-", species, ".html")
  }
  
  rmarkdown::render("Species_AC_analysis_v11.Rmd", 
                    params = list(
                      species = species,
                      initial = initial,
                      batch_int = batch_int,
                      integrate_by = integrate_by,
                      harmony = harmony,
                      contamination_threshold = contamination_threshold,
                      contamination = contamination,
                      nFeature_threshold = nFeature_threshold,
                      doublet_finder = doublet_finder,
                      manual_annotation = manual_annotation,
                      sac_annotation = sac_annotation,
                      de_expression = de_expression,
                      save = save
                    ),
                    output_file = output_file
  )
}

ScTypeAnnotation = function(object, gs_list){
  library(HGNChelper)
  
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  # DB file
  # db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  # tissue = "Eye" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = object[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  
  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
  # In case Seurat is used, it is either seurat[["RNA"]]@scale.data (default), seurat[["SCT"]]@scale.data, in case sctransform is used for normalization,
  # or seurat[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
  
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(object@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(object@meta.data[object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(object@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Other"
  
  object@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    object@meta.data$customclassif[object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  return(object)
}

AnnotatedUmap = function(object, annotation, group.by = "seurat_clusters", 
                         umap.group.by = NULL, color.clusters.by = "cell_class", 
                         rel_widths = c(1,2.5), plot.umap = TRUE, 
                         plot.dotplot = TRUE, color.genes = FALSE, 
                         title = NULL, bar_cols = NULL, 
                         plot.proportions = FALSE, umap.legend = FALSE, 
                         pretty.umap = FALSE, return.list = FALSE, ...){
  
  Idents(object) = group.by
  orig.object = object # AnnotatedDotplot function is sensitive to the factor levels of color.clusters.by
  
  object@meta.data[,color.clusters.by] = factor(factor(object@meta.data[,color.clusters.by], levels = unique(Annotation(annotation))))
  colors = GetAnnotationColors(levels(object@meta.data[,color.clusters.by]), annotation)
  message("Using the following colors: ", paste0(colors, collapse = ", "))
  
  if(is.null(umap.group.by)) umap.group.by = group.by
  
  if(plot.umap){
    if(pretty.umap){
      umap = PrettyUmap2(object, group.by = color.clusters.by, cols = colors, show.legend = ifelse(umap.legend, TRUE, FALSE), ...) 
    } else {
      umap = TitlePlot(theme_umap(ClusterBatchPlot(object, batch = color.clusters.by, cols = colors, shuffle = TRUE, group.by = umap.group.by, ...)), title = title) + 
        {if(!umap.legend) NoLegend()}
    }
  }
  
  if(plot.dotplot & plot.proportions){
    if(is.null(bar_cols)) bar_cols = AnnotationColors(annotation)[ Metadata(object, group.by, color.clusters.by)[[color.clusters.by]] ]
    print(bar_cols)
    
    plt = plot_grid(
      # UMAP
      umap,
      
      # Dotplot
      AnnotatedDotPlot(orig.object, annotation, group.by = group.by, color.clusters.by = color.clusters.by,
                       features = Genes(annotation), color.genes = color.genes) + 
        RotatedAxis() + 
        NoLegend() +
        coord_flip()+
        # {if(coord_flip) coord_flip()} +
        theme(axis.title = element_blank()),
      
      # Barplot
      CelltypeProportionBarplot(object, x = group.by, y = 'animal', show.all = FALSE) + 
        scale_fill_manual(values = bar_cols) + 
        theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
        coord_flip(),
      rel_widths = rel_widths,
      ncol = 3,
      align = "h", 
      axis = "bt"
    )
  } else if(plot.dotplot & !plot.umap){
    plt = AnnotatedDotPlot(orig.object, annotation, group.by = group.by, color.clusters.by = color.clusters.by,
                       features = (Genes(annotation)), color.genes = color.genes) + 
        RotatedAxis() + 
        # {if(coord_flip) coord_flip()} +
        theme(axis.title = element_blank())
    
  } else if(plot.dotplot){
    dp = AnnotatedDotPlot(orig.object, annotation, group.by = group.by, color.clusters.by = color.clusters.by,
                          features = (Genes(annotation)), color.genes = color.genes) + 
      RotatedAxis() + 
      # {if(coord_flip) coord_flip()} +
      theme(axis.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
    
    # Return as a list? 
    if(return.list) return(list(umap, dp))
    
    plt = plot_grid(
      umap,
      dp,
      rel_widths = rel_widths, 
      align = "h", 
      axis = "bt"
    )
  } else {
    plt = umap
  }
  
  return(plt)
}

PlotAmacrineTypes = function(object, Gaba.gene.1 = "GAD1", Gaba.gene.2 = "GAD2", group.by = "seurat_clusters", annotation = NULL){
  
  # annotation = eval(parse(text = annotation))
  
  order = arrange(Metadata(object, feature.1 = group.by, feature.2 = "classification", feature.3 = "lit_type"), 
                  classification, 
                  factor(lit_type, levels = unique(Annotation(annotation))))[[group.by]]
  
  plt = AnnotatePlot(
          OrderedDotPlot(object, annotation, group.by = group.by, color.clusters.by = group.by, order = order, 
                         features = rev(c(Gaba.gene.1, Gaba.gene.2, "SLC6A9", setdiff(Genes(annotation), "SLC6A9")))), 
          annotation, object_y = object, y_annotation = "lit_type"
          ) + coord_flip() + theme(axis.title = element_blank())
  
  return(plt)
}

GabaGlyClassification2 = function(object, Gaba.gene.1 = "GAD1", Gaba.gene.2 = "GAD2", Gly.gene = "SLC6A9", 
                                  mus1 = c(0.01,0.5), mus2 = c(0.01,0.5), sigmas1 = c(0.01,0.1), sigmas2 = c(0.01,0.1), 
                                  log.transform = FALSE, zero.base.clusters = TRUE){
  
  Idents(object) = "seurat_clusters"
  GabaGly.df = if(log.transform){
    t(log1p(AverageExpression(object, features = c(Gaba.gene.1, Gaba.gene.2, Gly.gene), assays = c("RNA"), slot = "data")$RNA))
  } else {
    t((AverageExpression(object, features = c(Gaba.gene.1, Gaba.gene.2, Gly.gene), assays = c("RNA"), slot = "data")$RNA))
  }
  
  GabaGly.df = as.data.frame(apply(GabaGly.df, 2, rescale))
  GabaGly.df$cluster = rownames(GabaGly.df)
  
  # Do this in the GADPlot function
  # if(zero.base.clusters){
  #   GabaGly.df$cluster = rownames(GabaGly.df)
  # } else {
  #   GabaGly.df$cluster = as.numeric(as.character(rownames(GabaGly.df)))+1
  # }
  
  # Average of GAD1 and GAD2 scaled, normalized expression
  if(Gaba.gene.1 %in% colnames(GabaGly.df) & Gaba.gene.2 %in% colnames(GabaGly.df)) {
    Gaba.gene = paste0(Gaba.gene.1, "n", Gaba.gene.2)
    GabaGly.df[[Gaba.gene]] = rescale((GabaGly.df[[Gaba.gene.1]] + GabaGly.df[[Gaba.gene.2]])/2)
    # GabaGly.df[[Gaba.gene]] = rescale(GabaGly.df[[Gaba.gene.1]] + GabaGly.df[[Gaba.gene.2]]) # rescale(log2(rescale((GabaGly.df$GAD1 + GabaGly.df$GAD2), to = c(1,10))))
  } else if(Gaba.gene.2 %in% colnames(GabaGly.df)) {
    Gaba.gene = Gaba.gene.2
    GabaGly.df[[Gaba.gene]] = GabaGly.df[[Gaba.gene.2]]
  } else if(Gaba.gene.1 %in% colnames(GabaGly.df)) {
    Gaba.gene = Gaba.gene.1
    GabaGly.df[[Gaba.gene]] = GabaGly.df[[Gaba.gene.1]]
  }
  
  GAD.model <- mixtools::normalmixEM(GabaGly.df[[Gaba.gene]], lambda=c(0.5,0.5), mu=mus1, sigma=sigmas1) # The starting values of sigma are important for the GAD model
  Glyt1.model <- mixtools::normalmixEM(GabaGly.df[[Gly.gene]], lambda=c(0.5,0.5), mu=mus2, sigma=sigmas2)
  
  GABA.df = cbind(GabaGly.df[,c(Gaba.gene, Gly.gene, "cluster")], ComputePosteriors(GabaGly.df[[Gaba.gene]], GAD.model))
  Glyt1.df = cbind(GabaGly.df[,c(Gaba.gene, Gly.gene, "cluster")], ComputePosteriors(GabaGly.df[[Gly.gene]], Glyt1.model))         
  
  # Find cutpoints
  cutoff.GABA = FindCutpoint(GAD.model)
  cutoff.Gly = FindCutpoint(Glyt1.model)
  
  # Plot
  plt = plot_grid(
    MixtureHistogram(GABA.df, Gaba.gene, GAD.model, intersects = FindCutpoint(GAD.model), legend.name = "p_2_xi", xlab = paste0(Gaba.gene, " expression"), plot.points = TRUE) + 
      NoLegend() + rremove("xlab") + rremove("x.text"), NULL, NULL,
    NULL, NULL, NULL,
    GADPlot(GabaGly.df, Gaba.gene = Gaba.gene, Gly.gene = Gly.gene, cutoff.Gly = cutoff.Gly, cutoff.GABA = cutoff.GABA, 
            zero.base.clusters = zero.base.clusters) + xlab(Gaba.gene), NULL, 
    MixtureHistogram(Glyt1.df, Gly.gene, Glyt1.model, intersects = FindCutpoint(Glyt1.model), legend.name = "p_2_xi", xlab = paste0(Gly.gene, " expression"), plot.points = TRUE) + 
      coord_flip() + NoLegend() + rremove("ylab") + rremove("y.text"), 
    ncol = 3, nrow = 3, rel_widths = c(2, 0, 1), rel_heights = c(1, 0, 2), align = "hv", axis = "bt")
  
  nGnGs = subset(GabaGly.df, eval(parse(text = Gaba.gene)) < cutoff.GABA & 
                   eval(parse(text = Gly.gene)) < cutoff.Gly)$cluster
  Dual.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) >= cutoff.GABA & 
                        eval(parse(text = Gly.gene)) >= cutoff.Gly)$cluster
  Gly.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) < cutoff.GABA & 
                       eval(parse(text = Gly.gene)) >= cutoff.Gly)$cluster
  Gaba.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) >= cutoff.GABA & 
                        eval(parse(text = Gly.gene)) < cutoff.Gly)$cluster
  
  message(length(nGnGs), " nGnGs types found: ", paste0(nGnGs, collapse = ", "))
  message(length(Dual.types), " dual types found: ", paste0(Dual.types, collapse = ", "))
  message(length(Gly.types), " glycinergic types found: ", paste0(Gly.types, collapse = ", "))
  message(length(Gaba.types), " GABAergic types found: ", paste0(Gaba.types, collapse = ", "))
  
  object$classification = NA
  object$classification[object$seurat_clusters %in% nGnGs] = "nGnG"
  object$classification[object$seurat_clusters %in% Dual.types] = "Both"
  object$classification[object$seurat_clusters %in% Gly.types] = "Gly"
  object$classification[object$seurat_clusters %in% Gaba.types] = "GABA"
  
  # Add posterior probabilities to seurat object
  object$GABA_posterior = GABA.df$p_2_xi[match(object$seurat_clusters,GABA.df$cluster)]
  object$Gly_posterior = Glyt1.df$p_2_xi[match(object$seurat_clusters,Glyt1.df$cluster)]
  object$GABA_score = GabaGly.df[[Gaba.gene]][match(object$seurat_clusters,GabaGly.df$cluster)]
  object$Gly_score = GabaGly.df[[Gly.gene]][match(object$seurat_clusters,GabaGly.df$cluster)]
  object
  
  # Order based on group
  object$ac.order = factor(object$seurat_clusters, 
                           levels = arrange(Metadata(object, feature.1 = "seurat_clusters", feature.2 = "classification", feature.3 = "lit_type"), classification, lit_type)[["seurat_clusters"]])
  
  return(list(object = object, plt = plt))
}

GabaGlyClassification = function(object, group.by = 'seurat_clusters', cutoff.Gly = 0.15, cutoff.GABA = 0.15, Gaba.gene.1 = "GAD1", Gaba.gene.2 = "GAD2", Gly.gene = "SLC6A9", ac.order = FALSE){
  
  Idents(object) = group.by
  GabaGly.df = t(AverageExpression(object, features = c(Gaba.gene.1, Gaba.gene.2, Gly.gene), assays = c("RNA"), slot = "data")$RNA)
  GabaGly.df = as.data.frame(apply(GabaGly.df, 2, rescale))
  GabaGly.df$cluster = rownames(GabaGly.df)
  
  # Sum of GAD1 and GAD2 average, scaled, normalized expression
  if(Gaba.gene.1 %in% colnames(GabaGly.df) & Gaba.gene.2 %in% colnames(GabaGly.df)) {
    Gaba.gene = paste0(Gaba.gene.1, "n", Gaba.gene.2)
    GabaGly.df[[Gaba.gene]] = rescale(GabaGly.df[[Gaba.gene.1]] + GabaGly.df[[Gaba.gene.2]]) # rescale(log2(rescale((GabaGly.df$GAD1 + GabaGly.df$GAD2), to = c(1,10))))
  } else if(Gaba.gene.2 %in% colnames(GabaGly.df)) {
    Gaba.gene = Gaba.gene.2
    GabaGly.df[[Gaba.gene]] = GabaGly.df[[Gaba.gene.2]]
  } else if(Gaba.gene.1 %in% colnames(GabaGly.df)) {
    Gaba.gene = Gaba.gene.1
    GabaGly.df[[Gaba.gene]] = GabaGly.df[[Gaba.gene.1]]
  }
  
  nGnGs = subset(GabaGly.df, eval(parse(text = Gaba.gene)) < cutoff.GABA & 
                   eval(parse(text = Gly.gene)) < cutoff.Gly)$cluster
  Dual.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) >= cutoff.GABA & 
                        eval(parse(text = Gly.gene)) >= cutoff.Gly)$cluster
  Gly.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) < cutoff.GABA & 
                       eval(parse(text = Gly.gene)) >= cutoff.Gly)$cluster
  Gaba.types = subset(GabaGly.df, eval(parse(text = Gaba.gene)) >= cutoff.GABA & 
                        eval(parse(text = Gly.gene)) < cutoff.Gly)$cluster
  
  message(length(nGnGs), " nGnGs types found: ", paste0(nGnGs, collapse = ", "))
  message(length(Dual.types), " dual types found: ", paste0(Dual.types, collapse = ", "))
  message(length(Gly.types), " glycinergic types found: ", paste0(Gly.types, collapse = ", "))
  message(length(Gaba.types), " GABAergic types found: ", paste0(Gaba.types, collapse = ", "))
  
  print(plot_grid({if(Gaba.gene.1 %in% colnames(GabaGly.df)) GADPlot(GabaGly.df, Gaba.gene = Gaba.gene.1, Gly.gene = Gly.gene)}, 
                  {if(Gaba.gene.2 %in% colnames(GabaGly.df)) GADPlot(GabaGly.df, Gaba.gene = Gaba.gene.2, Gly.gene = Gly.gene)}, 
                  GADPlot(GabaGly.df, Gaba.gene = Gaba.gene, Gly.gene = Gly.gene, cutoff.Gly = cutoff.Gly, cutoff.GABA = cutoff.GABA), 
                  ncol = 3))
  
  object$classification = NA
  object$classification[object@meta.data[[group.by]] %in% nGnGs] = "nGnG"
  object$classification[object@meta.data[[group.by]] %in% Dual.types] = "Both"
  object$classification[object@meta.data[[group.by]] %in% Gly.types] = "Gly"
  object$classification[object@meta.data[[group.by]] %in% Gaba.types] = "GABA"
  
  # Order based on group
  if(ac.order) object$ac.order = factor(object$seurat_clusters, 
                                        levels = arrange(Metadata(object, feature.1 = group.by, feature.2 = "classification", feature.3 = "lit_type"), classification, lit_type)[[group.by]])
  
  return(object)
}
