# Dario's functions!


compare_frequency = function(object){
  
  # 1. Pre-calculate totals to avoid repeated subsetting
  total_sc    <- sum(object$species == "scRNA-seq")
  total_xen   <- sum(object$species == "Xenium")
  
  # 2. Build the dataframe
  count.df <- data.frame(
    type = as.numeric(ExtractString(names(table(object$harmony.AC)), before = '_')), 
    
    scRNAseq = as.vector(
      table(object$harmony.AC[object$species == "scRNA-seq"]) / total_sc
    ), 
    
    Xenium_harmony = as.vector(
      table(object$harmony.AC[object$species == "Xenium"]) / total_xen
    )
  )
  
  count.df = arrange(count.df, type)
  ScatterPlot(count.df, 'scRNAseq', 'Xenium_harmony', labels = 'type', y_equals_x = T, r = T, logX = T, logY = T, max.overlaps = 10) | 
    ScatterPlot(count.df, 'scRNAseq', 'Xenium_harmony', labels = 'type', y_equals_x = T, r = T, max.overlaps = 10)

}

compute_retinal_regularity <- function(coords, concavity = 2, sample_size = 10000) {
  
  library(spatstat)
  library(concaveman)
  library(sf)
  library(ggplot2)
  
  # coords: data.frame with 'x' and 'y' columns
  
  # 1. Generate the Concave Hull
  # We convert to sf to calculate area accurately within the "patches"
  pts_sf <- st_as_sf(as.data.frame(coords), coords = c("x", "y"))
  hull <- concaveman(pts_sf, concavity = concavity)
  hull_area <- as.numeric(st_area(hull))
  num_points <- nrow(coords)
  
  # 2. NN Math for Retinal Regularity Index (Mean / SD)
  # nndist is efficient for large datasets
  dist_nn <- nndist(coords[,1], coords[,2])
  
  mu <- mean(dist_nn)
  sigma <- sd(dist_nn)
  
  # This is the "Regularity Index" (RI) used in retinal anatomy
  # Random ~ 1.7-2.0 | Regular/Mosaic ~ 4.0+
  RI <- mu / sigma
  
  # 3. Diagnostic Plotting
  # We sample points to prevent ggplot from hanging with 800k points
  plot_pts <- if(nrow(coords) > sample_size) {
    coords[sample(1:nrow(coords), sample_size), ]
  } else {
    coords
  }
  
  diag_plot <- ggplot() +
    geom_sf(data = hull, fill = "lightblue", alpha = 0.2, color = "blue", size = 0.5) +
    geom_point(data = as.data.frame(plot_pts), aes(x, y), size = 0.2, alpha = 0.3, color = "black") +
    labs(
      title = paste("Nearest neighbor regularity index (RI) =", round(RI, 3)),
      subtitle = paste("Mean NN Dist:", round(mu, 2), "µm | SD:", round(sigma, 2)),
      caption = paste("Concavity:", concavity, "| Area:", round(hull_area, 0), "µm²")
    ) +
    theme_minimal()
  
  # 4. Return results
  return(list(
    stats = data.frame(
      regularity_index = RI,
      mean_nn_dist = mu,
      sd_nn_dist = sigma,
      n_cells = num_points,
      area = hull_area
    ),
    plot = diag_plot
  ))
}

ProperCase <- function(x) {
  # 1. Convert everything to lowercase first
  x <- tolower(x)
  
  # 2. Extract the first letter, uppercase it, and paste it back to the rest
  # We use substring() to get the string from the 2nd character to the end
  paste0(toupper(substr(x, 1, 1)), substring(x, 2))
}


# Helper function for safe indexing
safe_subset <- function(m, rows, cols, impute = NA) {
  m_sub <- m[match(rows, rownames(m)), match(cols, colnames(m)), drop = FALSE]
  rownames(m_sub) <- rows
  colnames(m_sub) <- cols
  
  
  m_sub[is.na(m_sub)] = impute
  as.data.frame.matrix(m_sub)
}


barplot_blast = function(df, n = 10){
  
  # vector = df$bitscore %>% setNames(df$species2)
  # QuickBarplot(head(vector, n))
  
  df.sorted = df %>%
    group_by(species1) %>%
    slice_max(bitscore, n = n)
  
  ggbarplot(df.sorted, x = 'species2', y = 'bitscore', facet.by = 'species1', scales = "free_x") %>%
    facet(facet.by = "species1", nrow = 1, scales = 'free_x')+
    RotatedAxis()
  
}

run_phylo_pipeline <- function(prefix, target_genes, sequences, species_index, redo = FALSE) {
  
  # Define file paths based on the prefix
  raw_path      <- paste0("protein_trees/", prefix, "_raw.fa")
  aligned_path  <- paste0("protein_trees/", prefix, "_aligned.fa")
  trimmed_path  <- paste0("protein_trees/", prefix, "_trimmed.fa")
  treefile_path <- paste0(trimmed_path, '.treefile')
  
  files_exist <- file.exists(raw_path, aligned_path, trimmed_path, treefile_path)
  
  if (all(files_exist) && !redo) {
    message("Files already exist for prefix '", prefix, "'. Plotting existing tree. Use redo=TRUE to rerun.")
  } else {
    
    # Process and filter sequences
    family.sequences = AAStringSet(unlist(lapply(seq_along(sequences), function(i) {
      
      seq = sequences[[i]]
      
      # Filter by the target gene list
      seq = seq[toupper(names(seq)) %in% toupper(target_genes)]
      
      # Check if any sequences remained after filtering
      if(length(seq) == 0) return(NULL)
      
      # Name sequences using the species index
      names(seq) = paste0(rownames(species_index)[[i]], '-', names(seq))
      
      # 1. Sort by gene name, then by width (descending)
      df = data.frame(width = width(seq),
                      gene = names(seq),
                      index = 1:length(seq))
      df_sorted <- df[order(df$gene, -df$width), ]
      
      # 2. Keep only the first occurrence of each gene name (longest)
      df_longest <- df_sorted[!duplicated(df_sorted$gene), ]
      
      # 3. Subset sequences using the index
      seqs_filtered <- seq[df_longest$index]
      
      as.character(seqs_filtered)
    })))
    
    # Write raw sequences
    writeXStringSet(family.sequences, filepath = raw_path, format = "fasta")
    
    # Align with Muscle
    system(paste0('muscle -align ', raw_path, ' -output ', aligned_path))
    
    # Trim with trimAl
    system(paste0('trimal -in ', aligned_path, ' -out ', trimmed_path, ' -automated1'))
    
    # Run IQTree2 (pass -redo flag if requested)
    redo_flag <- if (redo) ' -redo' else ''
    system(paste0('iqtree2 -s ', trimmed_path, ' -m MFP -bb 1000 -alrt 1000 -nt AUTO', redo_flag))
  }
  
  # Plot the resulting treefile
  my_tree <- read.tree(treefile_path)
  plot(my_tree)
  nodelabels(my_tree$node.label, adj = c(1.2, -0.5), frame = "none", cex = 0.7)
  
  return(readAAStringSet(trimmed_path))
}

BatchVarExp = function(object, 
                       group.by, 
                       reduction = c('pca', 'harmony'), 
                       var.type = 'auto',
                       n_cells = 100000, 
                       seed = 12345){
  
  # Extract embedding
  embedding <- Embeddings(object, reduction = reduction)
  
  # Extract metadata
  meta <- object@meta.data
  
  # Subsample to N cells
  set.seed(seed)
  if(n_cells < ncol(object)){
    idx <- sample(ncol(object), n_cells)
    embedding_sub <- embedding[idx, ]
    meta_sub <- meta[idx, ]
  } else {
    embedding_sub = embedding
    meta_sub = meta
  }
  
  covariate <- meta_sub[[group.by]]
  
  # Auto-detect variable type if not specified
  if(var.type == 'auto'){
    var.type <- ifelse(is.numeric(covariate) & length(unique(covariate)) > 20,
                       'continuous', 'categorical')
    message("Detected var.type: ", var.type)
  }
  
  # Coerce based on type
  if(var.type == 'categorical'){
    covariate <- as.factor(covariate)
  } else if(var.type == 'continuous'){
    covariate <- as.numeric(covariate)
  }
  
  # Compute multivariate R²
  grand_mean <- colMeans(embedding_sub)
  SS_total <- sum(rowSums((embedding_sub - grand_mean)^2))
  
  fit <- lm(embedding_sub ~ covariate)
  SS_model <- sum(rowSums((fitted(fit) - grand_mean)^2))
  
  r.squared <- SS_model / SS_total
  r.squared
}

amari_score <- function(mat) {
  
  N_col <- ncol(mat)
  N_row <- ncol(mat)
  
  col_max <- sum(apply(mat, 2, max))  # for each column, sum best value
  row_max <- sum(apply(mat, 1, max))  # for each row, sum best value
  
  best_score = N_col + N_row
  diss <- (col_max + row_max) / best_score
  diss
}

format_p <- function(p) {
  if (p < 0.001) {
    # Convert to scientific notation string: "3.2e-04" -> "3.2 %*% 10^-4"
    formatted <- format(p, scientific = TRUE, digits = 3)
    formatted <- gsub("e", " %*% 10^", formatted)
    # Remove leading plus signs in exponent if they exist
    formatted <- gsub("\\+", "", formatted)
    return(paste0(" == ", formatted))
  } else {
    return(paste0(" == ", round(p, 3)))
  }
}

convert_ids = function(ids, dataset = "mfascicularis_gene_ensembl"){
  # 2. Connect to the Ensembl BioMart
  # We use the Macaque dataset as our starting point
  mart <- useMart("ensembl", dataset = dataset)
  
  # 4. Run the query
  # 'ensembl_gene_id' is the Macaque ID (the filter)
  # 'hsapiens_homolog_associated_gene_name' is the Human Gene Symbol (the attribute)
  ortholog_map <- getBM(
    attributes = c(
      "ensembl_gene_id", 
      "hsapiens_homolog_associated_gene_name",
      "hsapiens_homolog_orthology_type"
    ),
    filters = "ensembl_gene_id",
    values = ids,
    mart = mart
  )
  
  # 5. View results
  ortholog_map = subset(ortholog_map, hsapiens_homolog_orthology_type == 'ortholog_one2one')
  ortholog_map$hsapiens_homolog_associated_gene_name
}

process_ortho_heatmap <- function(csv_path, 
                                  title_label = "Bridge Heatmap",
                                  cluster_col = "leiden_clusters_4.5",
                                  index_ref = index, 
                                  ortho_ref = ac.ortho,
                                  ortho_levels = ORTHOTYPES, 
                                  stagger.threshold) {
  
  ac.umap = as.data.frame(fread(csv_path)) 
  ac.umap$leiden_clusters = ac.umap[[cluster_col]]
  ac.umap$species_full = convert_values(ac.umap$species, index$species %>% setNames(index$ident))
  ac.umap$species_full = factor(factor(ac.umap$species_full, levels = index$species))
  
  # Just the lizard, because it has '-1' suffix
  ac.umap$orthotypes_li = ac.ortho$type[paste0(ac.umap$species_full, '_', apply(str_split_fixed(ac.umap$V1, '-', 3)[,1:2], 1, paste, collapse = '-'))]
  
  # Chicken, opossum, mouse, and macaque ('x' suffix)
  ac.umap$orthotypes_ch = ac.ortho$type[paste0(ac.umap$species_full, '_', ExtractString(ac.umap$V1, after = '-'))]
  
  # Combine orthotype labels
  ac.umap$orthotypes = ifelse(is.na(ac.umap$orthotypes_li), as.character(ac.umap$orthotypes_ch), as.character(ac.umap$orthotypes_li))
  ac.umap$orthotypes = orthotype_labels(ac.umap$orthotypes, as.factor = FALSE)
  
  message('Using species: ', unique(ac.umap$species[!is.na(ac.umap$orthotypes)]), ' for cross-tabulation!')
  
  # Gemini was messing up here, which is why it looked different
  # # 3. Orthotype Extraction
  # # Extracting prefix (V1 parts 1 & 2) and suffix (after '-')
  # df = ac.umap
  # v1_parts <- stringr::str_split_fixed(df$V1, '-', 3)
  # v1_prefix <- apply(v1_parts[, 1:2], 1, paste, collapse = '-')
  # 
  # # Using your custom ExtractString or stringr alternative
  # v1_suffix <- stringr::str_extract(df$V1, "(?<=-).*") 
  # 
  # # Mapping from ac.ortho reference
  # df$orthotypes_li <- ortho_ref$type[paste0(df$species_full, '_', v1_prefix)]
  # df$orthotypes_ch <- ortho_ref$type[paste0(df$species_full, '_', v1_suffix)]
  # df$orthotypes <- ifelse(is.na(df$orthotypes_li), 
  #                         as.character(df$orthotypes_ch), 
  #                         as.character(df$orthotypes_li))
  
  # Tabulation and counting
  js_mat = JSMatrix(table(factor(ac.umap$orthotypes, ORTHOTYPES), ac.umap$leiden_clusters))
  
  ari = round(adj.rand.index(factor(ac.umap$orthotypes, ORTHOTYPES), ac.umap$leiden_clusters), 2)
  # message('ARI: ')

  # 5. Plotting
  p <- JSHeatmap2(
    factor(ac.umap$orthotypes, ortho_levels), ac.umap$leiden_clusters,
    stagger.threshold = stagger.threshold,
    row.order = ortho_levels,
    border.col = NA,
    title = paste0(title_label, '\n ARI: ', ari)
  ) +
    ArialFont() +
    theme(
      axis.text.y = element_text(size = 6, color = 'black'),
      axis.ticks.length.y = unit(1, "pt"),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_y_discrete(position = "right") + 
    NoLegend()
  
  # Final Arrangement
  # final_plot <- ggpubr::ggarrange(
  #   p,
  #   common.legend = TRUE,
  #   legend = 'none'
  # )
  
  return(list(data = ac.umap, matrix = js_mat, plot = p))
}


get_numbers <- function(string_vector) {
  # str_extract (singular) returns a vector, not a list
  matches <- str_extract(string_vector, "-?\\d+\\.?\\d*")
  
  # Convert to numeric (non-matches stay as NA)
  return(as.numeric(matches))
}


FindMode = function(vec){
  names(tail(sort(table(vec)), 1))
  
}

FindSilhouette = function(object, 
                          group.by, 
                          reduction = 'harmony', 
                          dist = NULL, 
                          average = TRUE, 
                          method = 'cosine'){
  
  # Load library
  library(cluster)
  
  # Dist can be precomputed
  if(is.null(dist)) dist = proxy::dist(object@reductions[[reduction]]@cell.embeddings, method = method)
  
  # Compute cell-wise scores
  if(inherits(object, 'Seurat')){
    cell.scores = as.data.frame(silhouette(as.integer(object@meta.data[[group.by]]), 
                                           dist))
  } else {
    cell.scores = as.data.frame(silhouette(as.integer(object[[group.by]]), 
                                           dist))
  }
  
  # Per cell score? 
  if(!average) return(cell.scores)
  
  # Average score per cluster
  avg.sil = as.data.frame(aggregate(sil_width ~ cluster, data = cell.scores, FUN = mean))
  avg.sil
}

ReadFile <- function(path) {
  ext <- tools::file_ext(path)
  if(ext == "rds") {
    readRDS(path)
  } else if(ext == "qs2") {
    qs2::qs_read(path)
  } else {
    stop("Unsupported file extension: ", ext, ". Expected .rds or .qs2")
  }
}

QuickBarplot = function(vector, ...){
  
  summary = data.frame(group = names(vector), 
                       value = vector)
  
  PrettyBarplot(summary, x = 'group', y = 'value', ...)
}


DotPlot4 <- function(object, 
                     coord.flip = FALSE, 
                     max.pct = 100,
                     binarization = NULL, 
                     col.high = '#584B9FFF', 
                     col.low = 'lightgrey', 
                     show = NULL,
                     max.size = 6,
                     str_width = 20,
                     gene.groups = NULL, # named vector: names = genes, values = group labels
                     facet.ncol = NULL,  # passed to facet_wrap
                     ...){
  
  args = list(...)
  args$object = object
  args$col.min = -1
  args$col.max = 2
  
  if(!is.null(gene.groups)) {
    args$features = unlist(gene.groups)
    
    # Str wrapping
    names(gene.groups) = stringr::str_wrap(names(gene.groups), width = str_width)
  }
  
  if(coord.flip) {
    message('reversing order of genes for y-axis...')
    args$features = rev(args$features)
  } 
  # Remove features that are not present
  args$features = intersect(args$features, rownames(object))
  stopifnot(length(args$features) > 0)
  # Generate DotPlot data
  dot_data <- do.call(DotPlot, args)$data
  dot_data$key = paste0(dot_data$features.plot, '-', dot_data$id)
  # Scale the 'pct.exp' column
  dot_data$pct.exp <- pmin(dot_data$pct.exp, max.pct)
  if(!coord.flip){
    dot_data$id = factor(dot_data$id, levels = rev(levels(dot_data$id)))
  }
  if(!is.null(binarization)){
    bin.melt = reshape2::melt(binarization)
    bin.melt$key = paste0(bin.melt$Var1, '-', bin.melt$Var2)
    dot_data = dot_data[match(bin.melt$key, dot_data$key),]
    dot_data$heat_value = factor(bin.melt$value, levels = c(0,1))
  }
  # Add gene group annotation if provided
  if(!is.null(gene.groups)){
    if(inherits(gene.groups, 'list')) {
      gene.names = unlist(gene.groups)
      gene.groups = rep(names(gene.groups), sapply(gene.groups, length))
      names(gene.groups) = gene.names
    }
    dot_data$gene.group = gene.groups[as.character(dot_data$features.plot)]
    
    # Preserve group order as provided
    dot_data$gene.group = factor(dot_data$gene.group, levels = unique(gene.groups))
  }
  if(is.null(show)) show = unique(dot_data$id)
  
  
  # Remove NAs
  # dot_data = na.omit(dot_data)
  
  
  # Plot
  ggplot(subset(dot_data, id %in% show), aes(x = features.plot, y = id)) +
    {if(!is.null(binarization)) geom_tile(aes(fill = heat_value))} +
    {if(!is.null(binarization)) scale_fill_manual(values = c("0" = "white", "1" = "lightgrey"))} +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    theme_cowplot() +
    {if(coord.flip) coord_flip()} +
    {if(coord.flip) theme(axis.title = element_blank(), axis.text.y = element_text(face = 'italic'))} + 
    {if(!coord.flip) theme(axis.title = element_blank(), axis.text.x = element_text(face = 'italic'))} +
    RotatedAxis() + 
    scale_color_gradient(name = 'Scaled\nexpression', 
                         low = col.low, 
                         high = col.high, 
                         limits = c(-1, 2), 
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    scale_radius(name = 'Percent\nexpressed', limits = c(0, max.pct), range = c(0, max.size)) +
    ArialFont() +
    {if(!is.null(gene.groups) &&  coord.flip) facet_grid(gene.group ~ ., scales = "free_y", space = "free_y")} +
    {if(!is.null(gene.groups) && !coord.flip) facet_grid(. ~ gene.group, scales = "free_x", space = "free_x")}
}

BestResCurve = function(object, from = 0.04, to = 0.2, by = 0.01){
  res = seq(from = from, to = to, by = by)
  values = ResolutionSweep(object, from = from, to = to, by = by, FUN = function(x) 
    mean(silhouette(as.integer(x), dist(object@reductions$harmony@cell.embeddings))[,'sil_width']))
  
  plot(res, values)
  
  # FindClusters(object, resolution = res[which.max(values)])
}

overlap_p = function(list1, list2, population = 17000){
  intersection = intersect(list1, list2)
  phyper(length(intersection)-1,
         length(list1),
         population-length(list1),
         length(list2),
         lower.tail=FALSE, 
         log.p=FALSE)
}

DegOverlap = function(obj1, obj2, group.by = 'seurat_clusters', n = 50, 
                      return.genes = FALSE, plot = FALSE, DEG_FUN = TopNDEGs, 
                      use.hyper = FALSE, p.adjust = TRUE, population = 17000, ...){
  
  if(inherits(obj1, 'Seurat')){
    de1 = DEG_FUN(obj1, group.by = group.by, n = n)
    genes1 <- split(de1$gene, de1$cluster) 
  } else {
    genes1 = obj1
  }
  
  if(inherits(obj2, 'Seurat')){
    de2 = DEG_FUN(obj2, group.by = group.by, n = n)
    genes2 <- split(de2$gene, de2$cluster)
  } else {
    genes2 = obj2
  }
  
  overlap_matrix <- sapply(genes2, function(a) {
    sapply(genes1, function(b) {
        length(intersect(a, b))
      })
  })
  
  stat_matrix <- sapply(genes2, function(a) {
    sapply(genes1, function(b) {
      if(use.hyper){
        overlap_p(a, b, population = population)
      } else {
        jaccard(a, b)
      }
      
    })
  })
  
  if(use.hyper & p.adjust){
    stat_matrix <- matrix(p.adjust(stat_matrix, method = "BH"),
              nrow = nrow(stat_matrix),
              ncol = ncol(stat_matrix),
              dimnames = dimnames(stat_matrix))
  }
  
  if(use.hyper) stat_matrix = -log10(stat_matrix)
  
  if(plot) return(Heatmap2(stat_matrix, 
                           parenthesis = overlap_matrix, 
                           label = TRUE, 
                           legend.name = ifelse(use.hyper, '-log10 pval', 'Jaccard\nindex'), 
                           ...))
  
  return(list(stat_matrix = stat_matrix, 
              genes = list(genes1 = genes1, genes2 = genes2)))
  
}

ResolutionSweep = function(object, from = 0.1, to = 2, by = 0.1, mc.cores = 1, FUN = NULL, algorithm = 4, verbose = FALSE, seed = 12345){
  
  resolutions = seq(from, to, by = by)
  
  library(parallel)
  library(pbapply)
  clusterings = pblapply(resolutions, function(res){
    obj <- FindClusters(object, resolution = res, verbose = verbose, algorithm = algorithm, seed = seed)
    return(obj$seurat_clusters)
  }) # mc.cores = mc.cores
  
  names(clusterings) = paste0('RNA_snn_res.', resolutions)
  
  # n_k = sapply(clusterings, function(x) length(unique(x)))
  
  object@meta.data = cbind(object@meta.data, clusterings)
  
  if(!is.null(FUN)){
    results = dapply(names(clusterings), function(clustering){
      FUN(object@meta.data[[clustering]])
    })
    
    return(results)
  } 

  object
}


ShannonEntropy2 = function(vector, threshold = 0, normalize = FALSE){
  
  # Early return if all zero
  if(all(vector == 0)) return(0)
  
  # Scale the vector by removing background
  vector <- vector / sum(vector)
  vector[vector<threshold] <- 0
  # Renormalize
  vector <- vector / sum(vector)
  # Calculate Shannon entropy
  Hx <- -sum(vector[vector>0] * log(vector[vector>0]))
  
  if(normalize){
    Hx = Hx / log(length(vector))
  }

  return(Hx)
}

MaxARI = function(object, group.by, from = 0.1, to = 2, by = 0.1, mc.cores = 10){
  
  resolutions = seq(from, to, by = by)
  
  library(parallel)
  clusterings = mclapply(resolutions, function(res){
    obj <- FindClusters(object, resolution = res)
    return(obj$seurat_clusters)
  }, mc.cores = mc.cores)
  
  names(clusterings) = resolutions
  
  scores = sapply(clusterings, function(x) adj.rand.index(x, object@meta.data[[group.by]]))
  print(scores)
  
  clusterings[which.max(scores)]
}


BrowseIntegration = function(object, size_type = 3, size_class = 4, clust.cols = NULL, small = FALSE, ...){
  
  if(small){
    (PrettyUmap2(object, group.by = 'species', label = FALSE, show.legend = TRUE, title = 'By species', size = size_class, ...)  + LegendLowerRight() | 
       PrettyUmap2(object, group.by = 'common_class', geom.label = geom_text_repel, max.overlaps = 100, size = size_class, cols = clust.cols, min.segment.length = 0, title = 'By class', ...))
  } else {
    (PrettyUmap2(object, group.by = 'species', label = FALSE, show.legend = TRUE, title = 'By species', size = size_class)  + LegendLowerRight() | 
       PrettyUmap2(object, group.by = 'common_class', geom.label = geom_text_repel, max.overlaps = 100, size = size_class, cols = clust.cols, min.segment.length = 0, title = 'By class')) /
    (PrettyUmap2(subset(object, species == 'Human'), group.by = 'cell_type', geom.label = geom_text_repel, max.overlaps = 100, size = size_type, min.segment.length = 0, title = 'Human (type)') | 
       PrettyUmap2(subset(object, species == 'Mouse'), group.by = 'sn_annotation', geom.label = geom_text_repel, max.overlaps = 100, size = size_type, min.segment.length = 0, title = 'Mouse (type)'))
  }
  
}


dtColors = function(n = 435, omit = c('white', 'ivory', 'grey60', 'darkgrey', 'lightyellow', 'lightcyan', 'floralwhite', 'lightcyan1'), lighten.by = 0){
  colors = standardColors()[which(!standardColors() %in% omit)]
  if(length(colors) < n) n = length(colors)
  return(colorspace::lighten(colors[1:n], lighten.by))
}

labels2colors_dt = function(vector, ...){
  return(labels2colors(vector, colorSeq = dtColors(length(vector), ...)))
}

myPalette = function(n_colors, darken = 0.2){
  # Amacrine palette
  return(colorspace::darken(as.character(paletteer_c("grDevices::Spectral", n_colors), darken)))
}

AcPalette2 = function(n = 100, lighten.by = 0.4){
  colorspace::lighten(dtColors(n), lighten.by)
}

AcPalette3 = function(n, colors = c('deeppink', 'orange', 'chartreuse', 'cyan', '#584B9FFF', 'violet'), vary.by = 0.4){
  
  n.high = ceiling(n / length(colors))
  # n.low = floor(n / length(colors))
  col.list = unlist(lapply(colors, function(color) {
    colorRampPalette(c(color, colorspace::lighten(color, vary.by)))(n.high)
  }))
  
  col.list[1:n]
}

ClusterPalette2 = function(groups, colors = c('deeppink', 'orange', 'chartreuse', 'cyan', '#584B9FFF', 'violet', 'antiquewhite2'), 
                           debug = FALSE){
  
  levels = table(groups)[unique(groups)]
  col.list = unlist(lapply(seq_along(levels), function(i) {
    if(levels[[i]] > 1) {
      if(debug) browser()
      interweave(make_hcl_ramp(colors[[i]], levels[[i]]))
    } else {
      colors[[i]]
    }
  }))
  
  col.list
}


ClusterPalette = function(groups, colors = c('deeppink', 'orange', 'chartreuse', 'cyan', '#584B9FFF', 'violet', 'antiquewhite2'), 
                          vary.by = 0.2, dark.first = TRUE, debug = FALSE){
  
  levels = table(groups)[unique(groups)]
  col.list = unlist(lapply(seq_along(levels), function(i) {
    if(levels[[i]] > 1) {
      if(dark.first) {
        # if(levels[[i]] == 10 & debug) browser()
        # interweave(make_hcl_ramp(colors[[i]], levels[[i]]))
        interweave(colorRampPalette(c(colorspace::darken(colors[[i]], vary.by),
                                      # colors[[i]],
                                      colorspace::lighten(colors[[i]], vary.by)))(levels[[i]]))
      } else {
        interweave(colorRampPalette(c(colorspace::lighten(colors[[i]], vary.by), 
                                      # colors[[i]],
                                      colorspace::darken(colors[[i]], vary.by)))(levels[[i]]))
      }
    } else {
      colors[[i]]
    }
  }))
  
  col.list
}

DcPalette = function(n_colors, colors = c("violet", "deepskyblue", "chartreuse3", "orange",  "orangered","brown", 'grey')){
  # Double cone palette
  return(colorRampPalette(colors = colors)(n_colors))
}

SeuratV5toV4 = function(object, assay = 'RNA', ...){
  
  counts = attr(object@assays[[assay]], 'layers')$counts
  features = attr(object@assays[[assay]], "meta.data")$gene_names
  if(length(features) < 1) features = rownames(object[[assay]]@features)
  meta = object@meta.data
  
  rownames(counts) = features
  colnames(counts) = rownames(meta)
  object2 = CreateSeuratObject(counts, ...)
  object2@meta.data = meta
  object2
}

kNNClassifier_fast = function(object, group.by, graph = 'integrated_snn') {
  
  g <- object@graphs[[graph]]
  labels <- object@meta.data[[group.by]]
  
  # Identify unlabeled cells
  unlabeled <- which(is.na(labels))
  if (length(unlabeled) == 0) return(labels)
  
  labeled <- which(!is.na(labels))
  
  # Extract submatrix for unlabeled cells and their neighbors
  # This avoids processing all cells in the loop
  g_sub <- g[unlabeled, labeled, drop = FALSE]
  
  # Get labels for labeled cells
  labeled_labels <- labels[labeled]
  
  # For each unlabeled cell, find the most common label among neighbors
  # Using apply for vectorization
  new_labels <- apply(g_sub, 1, function(x) {
    # Find which labeled cells are neighbors (weight > 0)
    neighbor_idx <- which(x > 0)
    
    if (length(neighbor_idx) == 0) return(NA)
    
    # Get weights and labels for neighbors
    neighbor_weights <- x[neighbor_idx]
    neighbor_labels <- labeled_labels[neighbor_idx]
    
    # Aggregate weights by label and find max
    label_scores <- tapply(neighbor_weights, neighbor_labels, sum)
    names(which.max(label_scores))
  })
  
  # Update labels
  result <- labels
  result[unlabeled] <- new_labels
  
  return(result)
}

kNNClassifier2 = function(object, group.by, graph = 'integrated_snn'){
  
  g <- object@graphs[[graph]]
  labels = object@meta.data[[group.by]] # cell type from scRNA-seq
  # labels[object$orig.ident == 'patch-seq'] = NA # set patch-seq labels to NA
  
  # kNN classification of patch-seq cells
  new_labels <- labels
  unlabeled <- which(is.na(labels))
  stopifnot(length(unlabeled) > 0)
  
  for (i in unlabeled) {
    neigh <- which(g[i, ] > 0)
    lab   <- labels[neigh]
    w     <- g[i, neigh]
    
    keep <- !is.na(lab)
    if (any(keep)) {
      lab <- lab[keep]
      w   <- w[keep]
      
      score <- tapply(w, lab, sum)
      new_labels[i] <- names(which.max(score))
    }
  }
  
  return(new_labels)
}

kNNClassifier = function(object){
  
  g <- object@graphs$integrated_snn
  labels = object$goetz.cluster # cell type from scRNA-seq
  labels[object$orig.ident == 'patch-seq'] = NA # set patch-seq labels to NA
  
  # kNN classification of patch-seq cells
  new_labels <- labels
  unlabeled <- which(object$orig.ident == 'patch-seq')
  
  for (i in unlabeled) {
    neigh <- which(g[i, ] > 0)
    lab   <- labels[neigh]
    w     <- g[i, neigh]
    
    keep <- !is.na(lab)
    if (any(keep)) {
      lab <- lab[keep]
      w   <- w[keep]
      
      score <- tapply(w, lab, sum)
      new_labels[i] <- names(which.max(score))
    }
  }
  
  return(new_labels)
}

ScoreIntegration = function(object){
  pch = subset(object, orig.ident == 'patch-seq')
  missclassified.cells = names(which(pch$goetz.cluster != pch$knn.res))
  print(DimPlot(object, cells.highlight = missclassified.cells))
  length(which(pch$goetz.cluster == pch$knn.res)) / length(which(!is.na(pch$goetz.cluster)))
}

PlotIntegration = function(object){
  TitlePlot(DimPlot(object, pt.size = .1, group.by = 'orig.ident', shuffle = TRUE) + NoLegend(), 'After Seurat') | 
    TitlePlot(DimPlotLabeled(object, pt.size = .1, group.by = 'goetz.cluster') + NoLegend(), 'Type annotations')
}

TestIntegration = function(conf.obj, patch.size = 5, scrna.size = 10, return.obj = FALSE, ...){
  pchseq = subset(conf.obj, orig.ident == 'patch-seq')
  scrna = subset(conf.obj, orig.ident == 'scRNA-seq')
  int.obj = merge(DownsampleSeurat(pchseq, 'type', size = patch.size), 
                  DownsampleSeurat(scrna, 'type', size = scrna.size))
  int.obj = quietly_run2(ClusterSeurat, int.obj, ...)
  print(PlotIntegration(int.obj))
  int.obj$knn.res = kNNClassifier(int.obj)
  
  raw.accuracy = ScoreIntegration(int.obj)
  ds = DownsampleSeurat(subset(int.obj, orig.ident == 'patch-seq'), 'type', size = 1)
  plate1.acc = mean(subset(ds, plate == 'RGCp1')$knn.res == subset(ds, plate == 'RGCp1')$goetz.cluster, na.rm = TRUE)
  plate2.acc = mean(subset(ds, plate == 'RGCp2')$knn.res == subset(ds, plate == 'RGCp2')$goetz.cluster, na.rm = TRUE)
  numNAs = length(which(is.na(ds$knn.res)))
  
  message('Accuracy: ', raw.accuracy)
  message('mean plate 1: ', plate1.acc)
  message('mean plate 2: ', plate2.acc)
  message('number of NAs: ', numNAs)
  
  if(return.obj) return(int.obj)
  
  data.frame(raw.accuracy = raw.accuracy, 
             plate1.acc = plate1.acc, 
             plate2.acc = plate2.acc, 
             numNAs = numNAs)
}


ConservationHeatmap = function(features, 
                               col.low = 'cyan', 
                               col.high = 'deeppink', 
                               min.z.score = -3,
                               max.z.score = 3, 
                               font.family = 'ArialMT', 
                               font.size = 10, 
                               column_title_rot = 45,
                               border.width = 1,
                               ...){
  
  col_fun = circlize::colorRamp2(c(min.z.score, 0, max.z.score), c(col.low, "white", col.high))
  matrix.list = readRDS("../../Ortho_Objects/scaled_OT_expression_log.rds")
  gene.matrix = do.call(rbind, matrix.list[features]) %>% as.data.frame() %>% setNames(levels(ac.ortho$type_unordered))
  gene.matrix = gene.matrix[,levels(ac.ortho@meta.data[['type']])]
  colnames(gene.matrix) = orthotype_labels(colnames(gene.matrix))
  genes = ExtractString(rownames(gene.matrix), after = '\\.')
  genes = factor(genes, levels = unique(genes))
  
  ha_bottom <- HeatmapAnnotation(
    split = anno_block(
      gp = gpar(fill = NA),
      labels = levels(factor(genes)),
      labels_gp = gpar(
        fontfamily = font.family,
        fontsize = font.size,
        fontface = "italic"
      ),
      labels_rot = column_title_rot
    ),
    which = "column"   # 🔑 THIS fixes the error
  )
  
  ht = Heatmap(t(as.matrix(gene.matrix)), 
          name = 'Scaled\nexpression',
          column_split = genes, 
          column_gap = unit(0, "mm"), 
          cluster_rows = FALSE, 
          cluster_columns = FALSE,
          row_names_gp = gpar(fontsize = font.size, fontfamily = font.family),
          column_names_gp = gpar(fontsize = font.size, fontfamily = font.family),
          column_title_gp = gpar(fontfamily = font.family, fontsize = font.size, fontface = "italic"),
          row_title_gp = gpar(fontfamily = font.family, fontsize = font.size),
          column_title_rot = column_title_rot, 
          col = col_fun, 
          border_gp = gpar(col = "black", lty = 1, lwd = border.width), 
          heatmap_legend_param = list(
            at = c(min.z.score, 0, max.z.score),       # exact tick positions
            labels = c(min.z.score, "0", max.z.score),  # exact tick labels
            border = "black",      # Black border around legend
            labels_gp = gpar(col = "black", fontfamily = font.family), # Black tick labels
            ticks_gp = gpar(col = "black"), 
            title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
          ...)
  
  ht 
}

SaturationCurve = function(object, prop.values = c(0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1), resolution = 1.5){
  lapply(prop.values, function(proportion){
    message('Trying ', proportion * 100, '% of cells...')
    if(proportion == 100){
      Harmonize(object, batch = 'animal', nPCs = 30)
    } else {
      Harmonize(SubsampleSeurat(object, proportion * ncol(object)), batch = 'animal', 
                # nPCs = 30, 
                cluster_resolution = resolution)
    }
  })
}

scale_ticks_first_last <- function(axis = "y", 
                                   expand = expansion(mult = c(0,0.01)), 
                                   do = identity, 
                                   FUN = NULL,
                                   ...) {
  
  if(is.null(FUN) & axis == 'y') FUN = scale_y_continuous
  if(is.null(FUN) & axis == 'x') FUN = scale_x_continuous
  
  if (axis == "y") {
    FUN(
      ...,
      labels = function(x) {
        x = do(x)
        x_clean <- x[!is.na(x)]
        ifelse(x == min(x_clean) | x == max(x_clean), as.character(x), "")
      }, 
      expand = expand
    )
  } else if (axis == "x") {
    FUN(
      ...,
      labels = function(x) {
        x = do(x)
        x_clean <- x[!is.na(x)]
        ifelse(x == min(x_clean) | x == max(x_clean), as.character(x), "")
      }, 
      expand = expand
    )
  }
}

#' Write gene lists to a long-format Excel table
#'
#' Converts multiple named gene lists (e.g., shared, primate-specific,
#' rodent-specific, laurasiatherian-specific) into a single long-format
#' table with one row per gene and writes it to an Excel file.
#'
#' Empty gene lists (length 0) are safely ignored.
#'
#' @param gene_lists A named list of gene lists. Each element corresponds
#'   to an enrichment category (e.g., \code{shared}, \code{primate},
#'   \code{rodent}, \code{laurasiatherian}) and must itself be a named list
#'   where names are cell types and values are character vectors of gene names.
#' @param file Character string specifying the output Excel file path.
#' @param sheet Character string specifying the name of the Excel sheet.
#'
#' @return Invisibly returns the long-format \code{data.frame}.
#'
#' @export
write_gene_lists_long <- function(
    gene_lists,
    file = "DEG_lists_long.xlsx",
    sheet = "DEGs"
) {
  stopifnot(is.list(gene_lists))
  
  library(writexl)
  
  rows <- list()
  
  for (enriched_type in names(gene_lists)) {
    x <- gene_lists[[enriched_type]]
    stopifnot(is.list(x), !is.null(names(x)))
    
    for (ct in names(x)) {
      genes <- x[[ct]]
      
      # skip empty gene vectors
      if (length(genes) == 0) next
      
      rows[[length(rows) + 1]] <- data.frame(
        cell_type = ct,
        gene = genes,
        enriched_type = enriched_type,
        stringsAsFactors = FALSE
      )
    }
  }
  
  df <- do.call(rbind, rows)
  
  write_xlsx(setNames(list(df), sheet), file)
  
  invisible(df)
}

#' Make nice species names
species_labels = function(names){
  new_names = gsub('MouseLemur', 'Mouse lemur', 
                   gsub('CatShark', 'Cat shark', 
                        gsub('TreeShrew', 'Tree shrew', names)))
  new_names
}

#' Cluster until k groups
ClusterUntil = function(object, k = 2, resolution = 0.5, verbose = FALSE){
  fraction = 0.1*resolution
  object = FindClusters(object, resolution = resolution, verbose = verbose)
  if(nClusters(object) == k){
    message('Finished!')
    return(object)
  } else if(nClusters(object) > k) {
    message('Trying resolution ', resolution - fraction, '...')
    ClusterUntil(object, k = k, resolution = resolution - fraction)
  } else if(nClusters(object) < k) {
    message('Trying resolution', resolution - fraction, '...')
    ClusterUntil(object, k = k, resolution = resolution + fraction)
  } else {
    stop('unable to resolve number of clusters')
  }
}


#' This function computes the co-clustering matrix between clusters in x given another clustering y
#' @param x clustering of interest
#' @param y clustering on which to project x
#' @param row.norm row normalize the second matrix? This will make the resulting 
#'  co-clustering frequency also row normalized. Without row.norm, the values represent 
#'  the degree to which cluster i clusters with j (0 to 1), but will not be row-normalized. 
CoclusteringMatrix = function(x, y, row.norm = TRUE){
  # Row norm matrix of proportion of x going to y
  mat1 = as.matrix(ConfusionMatrix(x, y, plot = FALSE))
  
  # Row norm matrix of proportion of y going to x
  mat2 = t(ConfusionMatrix(x, y, plot = FALSE, row.norm = row.norm))
  
  # Take product to get co-clustering frequency
  ccf = mat1 %*% mat2
  ccf
}

MergeClusters2 <- function(object, idents, meta.data = "seurat_clusters", refactor = FALSE) {
  
  # Get the metadata column
  col_data <- object@meta.data[[meta.data]]
  
  # Find cells that match any of the idents to merge
  cells_to_merge <- col_data %in% idents
  
  # Get the new label (first ident in the list)
  new_label <- idents[1]
  
  # If factor, ensure the new label is a valid level
  if (is.factor(col_data)) {
    if (!(new_label %in% levels(col_data))) {
      levels(col_data) <- c(levels(col_data), new_label)
    }
  }
  
  # Assign new label to all cells in the clusters being merged
  col_data[cells_to_merge] <- new_label
  
  # Refactor if requested (drop unused levels and renumber)
  if (refactor) {
    col_data <- droplevels(col_data)
    if (is.factor(col_data)) {
      levels(col_data) <- 0:(length(levels(col_data)) - 1)
    }
  }
  
  # Update the metadata column
  object@meta.data[[meta.data]] <- col_data
  
  # Set active identity
  Idents(object) <- meta.data
  
  return(object)
}


make_hcl_ramp <- function(base_col, n, l_pad = 50, c_scale = 1) {
  
  rgb <- grDevices::col2rgb(base_col) / 255
  luv <- as(colorspace::sRGB(rgb[1], rgb[2], rgb[3]), "polarLUV")
  coords <- colorspace::coords(luv)
  
  # Extract coordinates correctly - coords is a matrix with one row
  h <- coords[1, "H"]
  c <- coords[1, "C"] * c_scale
  l <- coords[1, "L"]
  
  # Increase the range - go from darker to lighter
  l_vals <- seq(max(10, l - l_pad), min(95, l + l_pad), length.out = n)
  
  # Create colors and handle out-of-gamut issues
  colors <- colorspace::polarLUV(
    H = rep(h, n),
    C = rep(c, n),
    L = l_vals
  )
  
  # Use fixup=TRUE to handle out-of-gamut colors
  colorspace::hex(colors, fixup = TRUE)
}

AppendGabaGly = function(object){
  object$annotated[object$cell_class2 == 'gabaAC'] = paste0('gaba', object$annotated[object$cell_class2 == 'gabaAC'])
  object$annotated[object$cell_class2 == 'glyAC'] = paste0('gly', object$annotated[object$cell_class2 == 'glyAC'])
  object
}

#' Returns clusters renumbers from 1 to length(clusters)
RenumberClustering = function(clustering, zero.base = TRUE){
  
  f = factor(clustering)
  
  if(zero.base){
    clustering = factor(f, 
                        levels = levels(f), 
                        labels = seq_along(levels(f)) - 1)
  } else {
    clustering = factor(f, 
                        levels = levels(f), 
                        labels = seq_along(levels(f)))
  }
  
  clustering
}

#' Returns object with renumbers clusters
RenumberClusters = function(object, group.by = 'seurat_clusters', ...){
  
  object@meta.data[[group.by]] = RenumberClustering(object@meta.data[[group.by]], ...)
  object
}


# Return only tip descendants (vector of tip indices)
get_tip_descendants <- function(tree, node) {
  children <- tree$edge[tree$edge[,1] == node, 2]
  out <- integer(0)
  for (child in children) {
    if (child <= length(tree$tip.label)) {
      out <- c(out, child)
    } else {
      out <- c(out, get_tip_descendants(tree, child))
    }
  }
  out
}

are_adjacent_tips <- function(tree, tips) {
  tip_nodes <- match(tips, tree$tip.label)
  
  # MRCA
  mrca <- getMRCA(tree, tip_nodes)
  
  # Nodes along paths tip -> MRCA
  paths <- lapply(tip_nodes, function(x) nodepath(tree, x, mrca))
  connecting_nodes <- unique(unlist(paths))
  
  # Check adjacency at each internal node
  for (node in connecting_nodes) {
    children <- tree$edge[tree$edge[,1] == node, 2]
    
    # Force build a logical vector manually
    child_hits <- rep(FALSE, length(children))
    for (i in seq_along(children)) {
      desc <- get_tip_descendants(tree, children[i])
      child_hits[i] <- any(desc %in% tip_nodes)
    }
    
    # Now child_hits is GUARANTEED logical vector
    if (sum(child_hits) > 1)
      return(FALSE)
  }
  
  return(TRUE)
}



plot_tanglegram = function(t1, t2){
  library(dendextend)
  dnd1 <- as.dendrogram(t1)
  dnd2 <- as.dendrogram(t2)
  ## rearrange in ladderized fashion
  # dnd1 <- ladder(dnd1)
  # dnd2 <- ladder(dnd2)
  ## plot the tanglegram
  dndlist <- dendextend::dendlist(dnd1, dnd2)
  dndlist %>% entanglement
  dndlist %>% untangle(method = "step1side") %>% entanglement
  dendextend::tanglegram(dndlist %>% untangle(method = "step1side"), 
                         fast = TRUE, 
                         margin_inner = 5, 
                         main_left = "1", 
                         main_right = "2")
}

find_fraction_nodes = function(tree, clade, verbose = FALSE){
  tree$tip.label = ExtractString(tree$tip.label, after = '_')
  this.mrca = getMRCA(tree, as.character(clade))
  # browser()
  desc <- Descendants(tree, node = this.mrca, type = "tips")[[1]]
  if(verbose) print(desc)
  fraction.found = length(clade) / length(desc)
  fraction.found
  # is.monophyletic(tree, tips = as.character(clade))
}

shuffle_values = function(lst){
  
  vals <- unlist(tf.switches.clean, recursive = TRUE, use.names = FALSE)
  vals_shuffled <- sample(vals)
  refill <- function(x, values) {
    i <- 1
    
    fill_rec <- function(item) {
      if (is.list(item)) {
        lapply(item, fill_rec)
      } else if (is.atomic(item)) {
        n <- length(item)
        out <- values[i:(i+n-1)]
        i <<- i + n
        out
      } else {
        stop("Unsupported type")
      }
    }
    
    fill_rec(x)
  }
    
  return(refill(lst, vals_shuffled))
}

#' This should work as intended
quietly_run2 <- function(FUN, ...) {
  null_out <- file(nullfile(), open = "wt")
  null_err <- file(nullfile(), open = "wt")
  
  sink(null_out, type = "output")
  sink(null_err, type = "message")
  
  # Capture result before closing sinks
  result <- tryCatch({
    suppressWarnings(suppressMessages(FUN(...)))
  }, finally = {
    sink(type = "message")
    sink(type = "output")
    close(null_err)
    close(null_out)
  })
  
  result
}

#' Doesn't suppress everything
quietly_run <- function(FUN, ...) {
  out <- NULL
  capture.output({
    out <- suppressWarnings(suppressMessages(FUN(...)))
  }, file = NULL)
  out
}

MixOrthotypes = function(object, ref.name = 'oAC1 [A2]', target.name = 'oAC42* [SAC]', 
                         ALPHAS = seq(0.1, 0.9, by = 0.2)){
  
  ints = lapply(ALPHAS, function(alpha){
    
    message('Working on alpha ', alpha, '!')
    
    SPECIES = unique(object$species)
    objects = dapply(SPECIES, function(this.species){
      message('Working on species ', this.species, '!')
      object = subset(primateAC, species == this.species)
      target = subset(object, orthotype == target.name)
      reference = subset(object, orthotype == ref.name)
      dat = reference@assays$RNA@data
      n_cells = ncol(dat)
      target.expr = target@assays$RNA@data[,sample(1:ncol(target), n_cells, replace = TRUE)]
      new_dat = (1-alpha) * dat + alpha * target.expr
      mat <- as.matrix(object@assays$RNA@data)             # direct pointer, no Seurat copy
      cols <- which(object$orthotype == ref.name)   # FAST base R logical test
      mat[, cols] <- as.matrix(new_dat)                    # dgCMatrix ← dgCMatrix, runs in C
      object@assays$RNA@data <- Matrix(mat, sparse = TRUE)             # write back once
      object
    })
    
    primate_reint = merge(objects[[1]], objects[2:length(objects)])
    primate_reint = ClusterSeurat(primate_reint, 
                                  normalize = FALSE, # turn off normalization
                                  integrate.by = 'species', 
                                  cluster_resolution = 1.2
    )
    primate_reint
    
  })
  
  ints
}


factor_unique = function(vector){
  
  factor(vector, levels = unique(vector))
}

# A matrix with row and column names
SnakePlot = function(mat, 
                     lwd1 = 6, 
                     lwd2 = 1, 
                     lty2 = 'solid',
                     pt.size = 5, 
                     grid.pt.size = 1, 
                     cols = type_cols4){
  
  # Convert to long
  df <- as.data.frame(mat) %>%
    tibble::rownames_to_column("row") %>%
    pivot_longer(-row, names_to = "col", values_to = "value")
  
  # Define numeric coordinates
  col_levels <- colnames(mat)
  row_levels <- rownames(mat)
  
  df <- df %>%
    mutate(
      col_num = as.integer(factor(col, levels = col_levels)),
      row_num = as.integer(factor(row, levels = row_levels))
    )
  
  # ---- Horizontal segments (same value within a row)
  segments_h <- df %>%
    arrange(row_num, col_num) %>%
    group_by(row) %>%
    mutate(next_value = lead(value),
           next_col = lead(col_num)) %>%
    filter(value == next_value) %>%
    transmute(
      value,
      x = col_num,
      xend = next_col,
      y = row_num,
      yend = row_num
    )
  
  segments_h_skip <- df %>%
    na.omit() %>% 
    arrange(row_num, col_num) %>%
    group_by(row) %>%
    mutate(next_value = lead(value),
           next_col = lead(col_num)) %>%
    filter(value == next_value) %>%
    transmute(
      value,
      x = col_num,
      xend = next_col,
      y = row_num,
      yend = row_num
    )
  
  # ---- Vertical segments (same value within a column)
  segments_v <- df %>%
    arrange(col_num, row_num) %>%
    group_by(col) %>%
    mutate(next_value = lead(value),
           next_row = lead(row_num)) %>%
    filter(value == next_value) %>%
    transmute(
      value,
      x = col_num,
      xend = col_num,
      y = row_num,
      yend = next_row
    )
  
  # ---- Vertical segments (not adjacent)
  segments_v_skip <- df %>%
    na.omit() %>% 
    arrange(col_num, row_num) %>%
    group_by(col) %>%
    mutate(next_value = lead(value),
           next_row = lead(row_num)) %>%
    filter(value == next_value) %>%
    transmute(
      value,
      x = col_num,
      xend = col_num,
      y = row_num,
      yend = next_row
    )
  
  # Combine both
  segments_all <- bind_rows(segments_h, segments_v)
  
  # ---- Plot
  ggplot() +
    geom_point(
      data = df,
      aes(x = col_num, y = row_num),
      color = 'grey',
      size = grid.pt.size,
      show.legend = FALSE
    ) +
    # segments first so dots draw on top
    geom_segment(
      data = segments_h,
      aes(x = x, xend = xend, y = y, yend = yend, color = value),
      linewidth = lwd1,
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_segment(
      data = segments_v,
      aes(x = x, xend = xend, y = y, yend = yend, color = value),
      linewidth = lwd1,
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_segment(
      data = segments_v_skip,
      aes(x = x, xend = xend, y = y, yend = yend, color = value),
      linewidth = lwd2,
      linetype = lty2,
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_segment(
      data = segments_h_skip,
      aes(x = x, xend = xend, y = y, yend = yend, color = value),
      linewidth = lwd2,
      linetype = lty2,
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_point(
      data = na.omit(df),
      aes(x = col_num, y = row_num, color = value),
      size = pt.size,
      show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(col_levels),
      labels = col_levels,
      expand = expansion(0.01)
    ) +
    scale_y_reverse(
      breaks = seq_along(row_levels),
      labels = row_levels,
      expand = expansion(0.03)
    ) +
    scale_color_manual(values = cols)+
    theme_minimal(base_size = 14) +
    theme(axis.text = element_text(color = 'black'))+
    RotatedAxis()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank()
    )
}

row_norm = function(table){
  
  tab = as.data.frame.matrix(table)
  
  if(length(dim(tab)) == 1) {
    return(tab / sum(tab))
  } else {
    tab / rowSums(tab)
  }
}

col_norm = function(table){
  matrix = t(t(table)/colSums(table))
}

orthotype_no_labels = function(vector, as.factor = TRUE, asterisk = FALSE, ...){
  
  old.values <- c(
    "22_A2", "27_VG3", "18_TRHDE", "41_TRHDE", "34", "26_A8", "24_SEG", "2_NNgly",
    "15", "17_ROBO3", "29_SLC35D3", "35", "30_A17", "32_A17", "9", "1_MAF*", "3*",
    "4_nNOS", "14_CA2", "6_NPY", "10_MAF*", "11*", "25*", "13", "40*", "42", "38",
    "7", "20", "21_PDGFRA", "28", "5", "36_VIP", "23_RXRG", "39", "16_NNgaba*",
    "19_CRH*", "12_nNOS*", "33_CA1*", "31_nNOS", "37", "8_SAC*"
  )
  
  if(!asterisk) old.values = gsub('\\*', '', old.values)
  
  new.levels <- c(
    "oAC1 [A2]", "oAC2 [VG3]", "oAC3", "oAC4", "oAC5", "oAC6 [A8]",
    "oAC7 [SEG]", "oAC8 [NNgly]", "oAC9", "oAC10", "oAC11",
    "oAC12", "oAC13 [A17]", "oAC14 [A17]", "oAC15", "oAC16*", "oAC17*",
    "oAC18 [nNOS]", "oAC19 [CA2]", "oAC20 [NPY]", "oAC21*", "oAC22*",
    "oAC23*", "oAC24", "oAC25", "oAC26", "oAC27", "oAC28", "oAC29",
    "oAC30 [PDGFRA]", "oAC31", "oAC32", "oAC33 [VIP]", "oAC34", "oAC35",
    "oAC36* [NNgaba]", "oAC37* [CRH]", "oAC38* [nNOS]", "oAC39* [CA1]",
    "oAC40 [nNOS]", "oAC41", "oAC42* [SAC]"
  )
  
  old.values = ExtractString(old.values, after = '_')
  new.levels = gsub('\\*', '', gsub('oAC', '', ExtractString(new.levels, after = ' ')))
  
  if(as.factor){
    factor(convert_values2(vector, key_table = data.frame(old.names = old.values,
                                                          new.names = new.levels), ...), 
           levels = new.levels)
  } else{
    convert_values2(vector, key_table = data.frame(old.names = old.values,
                                                   new.names = new.levels), ...)
  }
  
}

orthotype_labels = function(vector, as.factor = TRUE, asterisk = TRUE, ...){
  
  old.values <- c(
    "22_A2", "27_VG3", "18_TRHDE", "41_TRHDE", "34", "26_A8", "24_SEG", "2_NNgly",
    "15", "17_ROBO3", "29_SLC35D3", "35", "30_A17", "32_A17", "9", "1_MAF*", "3*",
    "4_nNOS", "14_CA2", "6_NPY", "10_MAF*", "11*", "25*", "13", "40*", "42", "38",
    "7", "20", "21_PDGFRA", "28", "5", "36_VIP", "23_RXRG", "39", "16_NNgaba*",
    "19_CRH*", "12_nNOS*", "33_CA1*", "31_nNOS", "37", "8_SAC*"
  )
  if(!asterisk) old.values = gsub('\\*', '', old.values)
  
  new.levels <- ORTHOTYPES
  
  if(as.factor){
    factor(convert_values2(vector, key_table = data.frame(old.names = old.values,
                                                                new.names = new.levels), ...), 
           levels = new.levels)
  } else{
    convert_values2(vector, key_table = data.frame(old.names = old.values,
                                                         new.names = new.levels), ...)
  }
  
}

assign_class = function(ortho.annotated){
  ortho.annotated$class_species = factor(paste0(ortho.annotated$cell_class, '-', ortho.annotated$species), 
                                         levels = outer(unique(ortho.annotated$species), c('PR', 'BC', 'HC', 'AC', 'RGC', 'MG'), function(y,x) paste0(x, '-', y)))
  
  ortho.annotated$class2_species = factor(paste0(ortho.annotated$cell_class2, '-', ortho.annotated$species), 
                                          levels = outer(unique(ortho.annotated$species), c('PR', 'onBC', 'offBC', 'HC', 'glyAC', 'gabaAC', 'RGC', 'MG'), function(y,x) paste0(x, '-', y)))
  
  ortho.annotated
}

ControlMultiplicity = function(df, multiplicity = 2, sort = FALSE){

  df.trimmed = do.call(rbind, lapply(unique(df$row_annotation), function(anno) {
    # Each gene can appear no more than n times
    # Pickes the top n instances of each gene
    x = subset(df, row_annotation == anno)
    if(!is.na(x[1,'Chicken'])){
      x = x %>%
        mutate(.orig_order = row_number()) %>%
        group_by(Chicken) %>%
        slice_head(n = multiplicity) %>%
        ungroup() %>%
        arrange(.orig_order) %>%
        select(-.orig_order)
      # x[sample(seq_len(nrow(x))),]
    } else if(!is.na(x[1,'Lizard'])){
      x = x %>%
        mutate(.orig_order = row_number()) %>%
        group_by(Lizard) %>%
        slice_head(n = multiplicity) %>%
        ungroup() %>%
        arrange(.orig_order) %>%
        select(-.orig_order)
      # x[sample(seq_len(nrow(x))),]
    } else if(!is.na(x[1,'Zebrafish'])){
      x = x %>%
        group_by(Zebrafish) %>%
        slice_head(n = multiplicity) %>%
        ungroup() %>% as.data.frame()
      # x[sample(seq_len(nrow(x))),]
    }
    if(sort){
      x = x %>% 
        arrange(desc(fraction))
      # x %>%
      #   mutate(num_nas = rowSums(is.na(.))) %>%
      #   arrange(num_nas) %>%
      #   select(-num_nas)
    }
    x
  }))
  
  df.trimmed
}


difference = function(x){
  which(x[-1] != x[-length(x)]) + 1
}

smart_sort <- function(x, na.last = TRUE, decreasing = FALSE) {
  x <- as.character(x)
  
  # capture "prefix" (everything before the final run of digits) and the final digits
  parts <- regexec("^(.*?)([0-9]+)$", x)
  matches <- regmatches(x, parts)
  
  prefix <- vapply(matches, function(m) if (length(m) == 3) m[2] else NA_character_, character(1))
  num    <- vapply(matches, function(m) if (length(m) == 3) as.numeric(m[3]) else NA_real_, numeric(1))
  
  # clean prefix: remove trailing non-alphanumeric separators like '-' or '_'
  prefix <- ifelse(is.na(prefix), NA_character_, gsub("[^A-Za-z0-9]+$", "", prefix))
  
  # create ordering: prefix then numeric; NA prefixes/numbers will respect na.last
  ord <- order(prefix, num, na.last = na.last, decreasing = decreasing, method = "radix")
  x[ord]
}

#' Get reciprocal best hits from a protein blast matrix
#' 
#' @param gnnm blast graph
#' @param species1 reference species, will be used as rows
#' @param species2 target species, will be used as columns
#' 
#' @returns a matrix with species 1 gene names as rows and species 2 gene names as columns
#' 
GetReciprocalBestHits = function(gnnm, species1, species2){
  
  # Subset to mouse vs shark
  li2ch = gnnm[startsWith(rownames(gnnm), species1), startsWith(colnames(gnnm), species2)]
  ch2li = gnnm[startsWith(rownames(gnnm), species2), startsWith(colnames(gnnm), species1)]
  
  # Row max of li2ch
  q <- summary(li2ch) %>% group_by(i) %>% slice_max(order_by = x, n = 1, with_ties = FALSE)
  
  # Row max of ch2li
  r <- summary(ch2li) %>% group_by(i) %>% slice_max(order_by = x, n = 1, with_ties = FALSE)
  
  m1 = sparseMatrix(i = q$i, j = q$j, x = q$x, dims = dim(li2ch), dimnames = dimnames(li2ch))
  m2 = sparseMatrix(i = r$i, j = r$j, x = r$x, dims = dim(ch2li), dimnames = dimnames(ch2li))
  
  # Bidirectional best hits
  m = sqrt(m1 * t(m2))
  df = GraphToDF(m)
  
  rownames(m) = gsub(paste0(species1, '_'), '', rownames(m))
  message(paste0('Rows are species ', species1))
  
  colnames(m) = gsub(paste0(species2, '_'), '', colnames(m))
  message(paste0('Columns are species ', species2))
  
  return(m)
}

MatchACs = function(object){
  # if(!exists('ac.ortho')) ac.ortho = LoadACOrtho()
  if(!exists('vertebrateAC')) vertebrateAC = readRDS('../../Ortho_Objects/vertebrateAC_unintegrated.rds')
  vertebrateAC$barcode = ExtractString(Cells(vertebrateAC), before = '_')
  
  # Exclude ACs that are not in the original tetrapod integration
  object[,!(object$cell_class == 'AC' & !Cells(object) %in% vertebrateAC$barcode)]
}

LoadSAMapGraph = function(gnnm_path = 'samap/graphs/maps_pm_sh_kf_ze_am_ch_li_op_mm_mf_v1.gnnm.mtx', 
                          gns_path = 'samap/graphs/maps_pm_sh_kf_ze_am_ch_li_op_mm_mf_v1.gns.csv'){
  
  library(Matrix)
  
  # Read in homology graph files from SAMap
  gnnm = readMM(gnnm_path)
  gns = read.csv(gns_path, row.names = 1)
  
  # Set row and column names
  rownames(gnnm) = gns[,1]
  colnames(gnnm) = gns[,1]
  
  gnnm
}


which_diff = function(x){
  which(x[-1] != x[-length(x)]) + 1
}

species_names_pretty = function(species_names){
  species_names_pretty = gsub('CatShark|Shark|Cat shark', 'Shark', 
                              gsub('MouseLemur', 'Mouse lemur', 
                                   gsub('TreeShrew', 'Tree shrew', 
                                        gsub('ZebrafishLyu', 'Zebrafish', 
                                             species_names))))
  species_names_pretty
}


pretty_dendrogram = function(dendro, space_bw_text = 0.2, line_start_buffer = 1, linetype = '22', label_size = 4){
  
  # Get data
  dend_data <- dendro_data(dendro)
  
  # Get leaf and segment data
  seg <- segment(dend_data)
  lab <- label(dend_data)
  
  # Compute the maximum height for where labels will go
  x_max <- max(seg$yend)
  
  ggplot() +
    # dendrogram branches
    geom_segment(
      data = seg,
      aes(x = -y, y = -x, xend = -yend, yend = -xend),
      lineend = "round"
    ) +
    
    # dashed lines from leaves to labels
    geom_segment(
      data = lab,
      aes(x = y+line_start_buffer, xend = x_max * space_bw_text - 1, # make it end before the end of the word
          y = -x, yend = -x),
      linetype = linetype, color = "black"
    ) +
    
    # labels to the right of dashed lines
    geom_label(
      data = lab,
      aes(x = x_max * space_bw_text + 0.02, # a bit of buffer
          y = -x, label = label),
      hjust = 1, size = label_size, 
      label.size = 0,         # removes border around label box
      label.padding = unit(0, "lines"), # make tighter box
      fill = "white"          # optional: keeps white background behind text
    ) +
    
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}

renumber = function(fact){
  
  stopifnot(inherits(fact, 'factor'))
  levels = levels(factor(fact))
  match(fact, levels)
  
}

bootstrap_operator = function(var1, var2, FUN = `+`, n_boot = 1000){
  sapply(seq_len(n_boot), function(x){
    mean1 = mean(sample(var1, size = length(var1), replace = TRUE))
    mean2 = mean(sample(var2, size = length(var2), replace = TRUE))
    FUN(mean1, mean2)
  })
}

elementwise_mean = function(mat_list){
  mean_mat <- Reduce("+", mat_list) / length(mat_list)
}

ViewNestedList = function(list){
  print(names(list))
  if(length(list[[1]]) > 1) ViewNestedList(list[[1]])
}

dapply = function(list, FUN){
  lapply(list, FUN) %>% setNames(list)
}

SummarizeMappings = function(object){
  
  library(formattable)
  
  # GLM
  object$first.glm = glue('{object$glm.anno.ac.1} ({round(object$glm.prob.ac.1, 2)})')
  object$second.glm = glue('{object$glm.anno.ac.2} ({round(object$glm.prob.ac.2, 2)})')
  object$third.glm = glue('{object$glm.anno.ac.3} ({round(object$glm.prob.ac.3, 2)})')
  
  # XGB
  object$first.xgb = glue('{object$xgb.anno.ac.1} ({round(object$xgb.prob.ac.1, 2)})')
  object$second.xgb = glue('{object$xgb.anno.ac.2} ({round(object$xgb.prob.ac.2, 2)})')
  object$third.xgb = glue('{object$xgb.anno.ac.3} ({round(object$xgb.prob.ac.3, 2)})')
  
  # TabPFN
  object$first.tab = glue('{object$tab.anno.ac.1} ({round(object$tab.prob.ac.1, 2)})')
  object$second.tab = glue('{object$tab.anno.ac.2} ({round(object$tab.prob.ac.2, 2)})')
  object$third.tab = glue('{object$tab.anno.ac.3} ({round(object$tab.prob.ac.3, 2)})')
  
  df = object[,startsWith(object$plate, 'ACp')]@meta.data[,c('nFeature_RNA', 
                                                            'glm.anno.class',
                                                            'xgb.anno.class',
                                                            'first.glm',
                                                            'first.xgb',
                                                            'first.tab',
                                                            'second.glm',
                                                            'second.xgb',
                                                            'second.tab',
                                                            'third.glm',
                                                            'third.xgb',
                                                            'third.tab',
                                                            # 'agreement',
                                                            'module.score.glm',
                                                            'module.score.xgb',
                                                            'module.score.tab',
                                                            'lit_type')]
    
  
  all_three_agree = ExtractString(df$first.glm, after = ' ') == ExtractString(df$first.xgb, after = ' ') & 
    ExtractString(df$first.glm, after = ' ') == ExtractString(df$first.tab, after = ' ')
  two_agree = ExtractString(df$first.glm, after = ' ') == ExtractString(df$first.xgb, after = ' ') | 
    ExtractString(df$first.glm, after = ' ') == ExtractString(df$first.tab, after = ' ') | 
    ExtractString(df$first.xgb, after = ' ') == ExtractString(df$first.tab, after = ' ')
  
  formattable(df, align = c("l", rep("c", ncol(df) - 1)),
              list(
                glm.anno.class = formatter("span",
                                           style = x ~ style(
                                             color = ifelse(df$glm.anno.class == 'AC', "green", "red"),
                                             font.weight = "bold"
                                           )
                ),
                xgb.anno.class = formatter("span",
                                           style = x ~ style(
                                             color = ifelse(df$xgb.anno.class == 'AC', "green", "red"),
                                             font.weight = "bold"
                                           )
                ),
                first.glm = formatter("span",
                                      style = x ~ style(
                                        color = ifelse(all_three_agree, "green", ifelse(two_agree, "gold", "red")),
                                        font.weight = "bold"
                                      )
                ),
                first.xgb = formatter("span",
                                      style = x ~ style(
                                        color = ifelse(all_three_agree, "green", ifelse(two_agree, "gold", "red")),
                                        font.weight = "bold"
                                      )
                ),
                first.tab = formatter("span",
                                      style = x ~ style(
                                        color = ifelse(all_three_agree, "green", ifelse(two_agree, "gold", "red")),
                                        font.weight = "bold"
                                      )
                ),
                module.score.glm = color_tile("white", "lightgreen"),
                module.score.xgb = color_tile("white", "lightgreen"),
                module.score.tab = color_tile("white", "lightgreen"),
                nFeature_RNA = color_tile("white", "lightgreen"),
                lit_type = formatter("span",
                                     style = x ~ style(color = ifelse(replace_na(startsWith(df$lit_type, 'AC'), FALSE), "green", "red"), 
                                                       font.weight = "bold"))
              )
  )
  
  # Kable
  # cat("\\begin{landscape}")
  # knitr::kable(subset(object, startsWith(plate, 'ACp'))@meta.data[,c('nFeature_RNA', 
  #                                                                                          'glm.anno.class', 
  #                                                                                          'xgb.anno.class', 
  #                                                                                          'first.glm', 
  #                                                                                          'first.xgb', 
  #                                                                                          'second.glm', 
  #                                                                                          'second.xgb', 
  #                                                                                          'third.glm', 
  #                                                                                          'third.xgb', 
  #                                                                                          # 'agreement', 
  #                                                                                          'module.score.glm', 
  #                                                                                          'module.score.xgb', 
  #                                                                                          'lit_type')], 
  #              align = 'c') %>% 
  #   kableExtra::kable_styling(font_size = 8, latex_options = "scale_down")
  # cat("\\end{landscape}")
}

vec2df = function(vector){
  data.frame(names = names(vector), 
             values = vector)
}

CollectSmartMatrix = function(smatrix, group.by){
  
  matrix = smatrix@matrix
  group = smatrix@col.data[,group.by]
  rownames(matrix) = group
  colnames(matrix) = group
  diag(matrix) = NA
  
  sim_df <- as.data.frame(as.table(matrix)) %>%
    dplyr::rename(sample1 = Var1, sample2 = Var2, similarity = Freq) %>%
    na.omit() %>% 
    group_by(sample1, sample2) 
    # summarise(mean_similarity = mean(similarity), .groups = "drop")
  
  sim_df
}
AverageSmartMatrix = function(smatrix, group.by){
  
  matrix = smatrix@matrix
  group = smatrix@col.data[,group.by]
  rownames(matrix) = group
  colnames(matrix) = group
  diag(matrix) = NA
  
  sim_df <- as.data.frame(as.table(matrix)) %>%
    dplyr::rename(sample1 = Var1, sample2 = Var2, similarity = Freq) %>%
    na.omit() %>% 
    group_by(sample1, sample2) %>%
    summarise(mean_similarity = mean(similarity), .groups = "drop")
  
  sim_df
}

get_root <- function(tree) {
  parents <- tree$edge[,1]
  children <- tree$edge[,2]
  # root = a parent that never appears as a child
  setdiff(parents, children)
}

getMRCA2 = function(tree, clade){
  if(length(clade) == 1){
    node = which(tree$tip.label == clade)
    node
    # This gives you the parent, not the MRCA of a single node!
    # tree$edge[tree$edge[,2] == node, 1]
  } else {
    getMRCA(tree, clade)
  }
}

GetSisterCladeDebug = function(tree, clade){
  
  # 1. Get MRCA (most recent common ancestor) node of the clade
  node <- getMRCA2(tree, clade)
  
  # 2. Find the parent node of this MRCA
  parent <- tree$edge[tree$edge[,2] == node, 1]
  
  if(length(parent) == 0) {
    message("browser")
    browser()
  }
  
  if(length(parent) < 1){
    tree = root(tree, node = sample.int(tree$Nnode, 1)+length(tree$tip.label))
    GetSisterClade(tree, clade)
  }
  
  if(parent == getRoot(tree)){
    
    message("Parent is root, rooting elsewhere...")
    # browser()

    # Re-root randomly at some internal node
    tree = root(tree, node = sample.int(tree$Nnode, 1)+length(tree$tip.label))
    # node <- getMRCA2(tree, clade)
    # parent <- tree$edge[tree$edge[,2] == node, 1]
    
    # Recursion step
    GetSisterClade(tree, clade)
    
  } else {
    
    # 3. Get the sister clade (the other child of the parent)
    # browser()
    sisters <- tree$edge[tree$edge[,1] == parent, 2]
    sister_node <- sisters[sisters != node]
    
    # 4. Extract the tips descending from that sister node
    # sister_tips <- extract.clade(tree, sister_node)$tip.label
    if(length(sister_node) > 1) {
      message("browser")
      browser()
    }
    if(sister_node <= length(tree$tip.label)){
      message("single sister")
      sister_tips = tree$tip.label[sister_node]
    } else {
      message("multiple sisters")
      sister_tips <- extract.clade(tree, sister_node)$tip.label
    }
  }
  
  sister_tips
  
}

GetSisterClade2 = function(tree, clade){
  
  # 1. Get MRCA (most recent common ancestor) node of the clade
  node <- getMRCA2(tree, clade)
    
  # 2. Find the parent node of this MRCA
  parent <- tree$edge[tree$edge[,2] == node, 1]
  
  if(length(parent) == 0) {
    message("browser")
    browser()
  }
  
  if(parent == getRoot(tree) | node == getRoot(tree)){
    
    message("Parent is root, rooting elsewhere...")
    
    # Re-root randomly at some internal node
    tree = root(tree, node = sample.int(tree$Nnode, 1)+length(tree$tip.label))
    node <- getMRCA2(tree, clade)
    parent <- tree$edge[tree$edge[,2] == node, 1]
    
    # Recursion step
    # GetSisterClade(tree, clade) 
    
  } 
  
  # 3. Get the sister clade (the other child of the parent)
  # browser()
  sisters <- tree$edge[tree$edge[,1] == parent, 2]
  sister_node <- sisters[sisters != node]
  
  # 4. Extract the tips descending from that sister node
  # sister_tips <- extract.clade(tree, sister_node)$tip.label
  if(length(sister_node) > 1) {
    message("browser")
    browser()
  }
  if(sister_node <= length(tree$tip.label)){
    message("single sister")
    sister_tips = tree$tip.label[sister_node]
  } else {
    message("multiple sisters")
    sister_tips <- extract.clade(tree, sister_node)$tip.label
  }
  
  sister_tips
  
}

GetSisterClade = function(tree, clade){
  
  # 1. Get MRCA (most recent common ancestor) node of the clade
  node <- getMRCA2(tree, clade)
  
  if(node == getRoot(tree)){
    
    # Reroot randomly at some internal node
    tree = root(tree, node = sample.int(tree$Nnode, 1)+length(tree$tip.label))
    
    # Recursion step
    GetSisterClade(tree, clade) 
  }
  else {
      
    # 2. Find the parent node of this MRCA
    parent <- tree$edge[tree$edge[,2] == node, 1]
    
    if(length(parent) == 0) browser()
    
    if(parent == getRoot(tree)){
      
      # Re-root randomly at some internal node
      tree = root(tree, node = sample.int(tree$Nnode, 1)+length(tree$tip.label))
      
      # Recursion step
      GetSisterClade(tree, clade) 
      
    } else {
      # 3. Get the sister clade (the other child of the parent)
      sisters <- tree$edge[tree$edge[,1] == parent, 2]
      sister_node <- sisters[sisters != node]
      
      # 4. Extract the tips descending from that sister node
      # sister_tips <- extract.clade(tree, sister_node)$tip.label
      if(length(sister_node) > 1) browser()
      if(sister_node <= length(tree$tip.label)){
        sister_tips = tree$tip.label[sister_node]
      } else {
        sister_tips <- extract.clade(tree, sister_node)$tip.label
      }
    }
    
    sister_tips
  }
  
}

find_k_nearest_neighbors <- function(graph, k) {
  
  library(igraph)
  #' Find the k-nearest neighbors for each node in a graph.
  #'
  #' This function uses the shortest path distance to determine the "closeness"
  #' of nodes. It assumes an unweighted graph where the distance is the
  #' number of edges in the shortest path.
  #'
  #' @param graph An igraph object representing the graph.
  #' @param k An integer specifying the number of nearest neighbors to find.
  #' @return A list where each element is a vector of k-nearest neighbors
  #'         for a corresponding node.
  
  # Ensure k is a positive integer
  if (k <= 0) {
    warning("k must be a positive integer. Returning an empty list.")
    return(lapply(V(graph)$name, function(x) character(0)))
  }
  
  # Calculate the shortest path distances between all pairs of nodes
  # The 'distances' function in igraph is efficient for this task.
  # We use 'out' mode for directed graphs, but it works for undirected too.
  dist_matrix <- distances(graph, mode = "out")
  
  # Get the names of all nodes in the graph
  nodes <- V(graph)$name
  num_nodes <- length(nodes)
  
  # Initialize a list to store the results
  all_neighbors <- vector("list", num_nodes)
  names(all_neighbors) <- nodes
  
  # Iterate through each node to find its k nearest neighbors
  for (i in 1:num_nodes) {
    current_node <- nodes[i]
    
    # Get the distances from the current node to all other nodes
    dists <- dist_matrix[i, ]
    
    # Remove the distance to itself
    dists <- dists[-i]
    
    # Sort the distances in ascending order and get the node names
    sorted_nodes <- names(sort(dists))
    
    # Select the top k nodes
    k_neighbors <- head(sorted_nodes, k)
    
    # Store the result
    all_neighbors[[current_node]] <- k_neighbors
  }
  
  return(all_neighbors)
}

are.neighbors = function(dist.mat, nodes, binary = FALSE){

  stopifnot(all(nodes %in% rownames(dist.mat)))
  k = length(nodes)
  
  # Grab the knn for each node; these compose the pool
  knn.list = t(apply(dist.mat, 1, function(row) names(head(sort(row), k))))
  pool = unique(as.vector(knn.list[nodes,]))
  
  # Jaccard overlap between pool and nodes
  if(binary){
    all(nodes %in% pool)
  } else {
    # jaccard(pool, nodes)
    length(intersect(pool, nodes))/length(pool)
  }
  
}

#' Testing for phylogenetic signal
#' 
#' @param tree an object of class phylo
#' @param variable a vector to test for, whose names are the species correpsonding to tree tips
#' 
#' @return A data.frame with LRT p-value and statistic
PhylogeneticSignalTest = function(tree, variable){
  
  library(phylolm)
  
  trait_data = data.frame(species = names(variable), 
                          trait_A = variable)
  
  # If binary, then do logistic; else do linear regression
  if(all(variable %in% c(0,1))){
    
    # Fit logistic regression assuming a tree-dependent covariance structure
    model_phyloglm <- phyloglm(trait_A ~ 1, data = trait_data, phy = tree)
    
    # Fit a standard generalized linear model (GLM) assuming independence
    model_glm <- glm(trait_A ~ 1, data = trait_data, family = "binomial")
    
    # --- Likelihood Ratio Test (LRT) ---
    LRT_statistic <- 2 * (logLik(model_phyloglm) - logLik(model_glm))
    LRT_p_value <- pchisq(LRT_statistic, df = 1, lower.tail = FALSE)
    
  } else {
    
    # Linear regression with tree
    model_phylolm <- phylolm(trait_A ~ 1, data = trait_data, phy = tree, model = "BM")
    
    # Linear regression without tree
    model_lm <- lm(trait_A ~ 1, data = trait_data)
    
    # --- Likelihood Ratio Test (LRT) ---
    LRT_statistic <- 2 * (logLik(model_phylolm) - logLik(model_lm))
    LRT_p_value <- pchisq(LRT_statistic, df = 1, lower.tail = FALSE)
    
  }
  
  cat(paste("LRT p-value:", signif(LRT_p_value, 3)))
  data.frame(stat = as.numeric(LRT_statistic), pval = as.numeric(LRT_p_value))
  
}


invert_list <- function(my_list) {
  # Stack the list to create a two-column data frame (values, names)
  stacked_list <- stack(my_list)
  
  # Use tapply to group the original names (ind) by the new elements (values)
  inverted <- tapply(stacked_list$ind, stacked_list$values, c)
  
  # Convert the result back to a named list
  as.list(inverted)
}

calculate_silhouette_coefficient <- function(dist_matrix, cluster_nodes, all_clusters) {
  
  # --- Input Validation ---
  if (!is.matrix(dist_matrix) || nrow(dist_matrix) != ncol(dist_matrix)) {
    stop("Input 'dist_matrix' must be a square matrix.")
  }
  
  if (!is.character(cluster_nodes) && !is.factor(cluster_nodes)) {
    stop("Input 'cluster_nodes' must be a vector of character strings or factors.")
  }
  
  if (!is.list(all_clusters) || !all(sapply(all_clusters, is.character))) {
    stop("Input 'all_clusters' must be a list of character vectors.")
  }
  
  all_nodes <- unlist(all_clusters)
  if (length(unique(all_nodes)) != length(all_nodes)) {
    stop("Nodes must be unique across all clusters.")
  }
  
  if (!all(cluster_nodes %in% all_nodes)) {
    missing_nodes <- cluster_nodes[!cluster_nodes %in% all_nodes]
    stop("The following nodes from 'cluster_nodes' were not found in 'all_clusters': ",
         paste(missing_nodes, collapse = ", "))
  }
  
  # --- Preparation for Calculation ---
  # Create a vector of cluster assignments for each node in the distance matrix
  cluster_assignments <- rep(NA, nrow(dist_matrix))
  names(cluster_assignments) <- rownames(dist_matrix)
  
  for (i in 1:length(all_clusters)) {
    cluster_name <- all_clusters[[i]][1] # Use first node name as cluster identifier
    cluster_assignments[all_clusters[[i]]] <- paste0("Cluster-", i)
  }
  
  # Filter out nodes not in any cluster to match the distance matrix
  valid_nodes <- names(cluster_assignments[!is.na(cluster_assignments)])
  dist_matrix <- dist_matrix[valid_nodes, valid_nodes]
  cluster_assignments <- cluster_assignments[valid_nodes]
  
  # Need to convert the matrix to a dist object for silhouette() function
  dist_object <- as.dist(dist_matrix)
  
  # --- Calculation ---
  library(cluster)
  sil_result <- silhouette(x = as.integer(as.factor(cluster_assignments)), dist = dist_object)
  
  # Get the average silhouette score for the specific cluster
  target_cluster_id <- unique(cluster_assignments[cluster_nodes])
  cluster_sil_score <- mean(sil_result[cluster_assignments == target_cluster_id, "sil_width"])
  
  return(cluster_sil_score)
}

# Not very useful or interpretable
calculate_normalized_distance <- function(dist_matrix, cluster_nodes, normalize = FALSE) {
  
  # --- Input Validation ---
  if (!is.matrix(dist_matrix) || nrow(dist_matrix) != ncol(dist_matrix)) {
    stop("Input 'dist_matrix' must be a square matrix.")
  }
  
  if (!is.character(cluster_nodes) && !is.factor(cluster_nodes)) {
    stop("Input 'cluster_nodes' must be a vector of character strings or factors.")
  }
  
  # Check if all cluster_nodes exist in the distance matrix
  if (!all(cluster_nodes %in% rownames(dist_matrix))) {
    missing_nodes <- cluster_nodes[!cluster_nodes %in% rownames(dist_matrix)]
    stop("The following nodes from 'cluster_nodes' were not found in the distance matrix: ",
         paste(missing_nodes, collapse = ", "))
  }
  
  # --- Calculation ---
  
  # 1. Calculate the average pairwise distance for the entire dataset.
  # We use the mean of the upper triangle to avoid double counting and the diagonal (which is 0).
  total_avg_dist <- mean(dist_matrix[upper.tri(dist_matrix, diag = FALSE)])
  
  # 2. Extract the sub-matrix for the specified cluster.
  sub_matrix <- dist_matrix[cluster_nodes, cluster_nodes]
  
  # 3. Calculate the average pairwise distance within the cluster.
  # We again use the upper triangle to get the average pairwise distance.
  cluster_avg_dist <- mean(sub_matrix[upper.tri(sub_matrix, diag = FALSE)])
  
  # 4. Compute the normalized distance.
  if(normalize) cluster_avg_dist / total_avg_dist else cluster_avg_dist
  
}

RFMatrix = function(list){
  rf.mat = Iterate(c(
    upgma = lapply(list, function(x) x$UPGMA),
    nj = lapply(list, function(x) x$NJ),
    mp = lapply(list, function(x) x$MP),
    ml = lapply(list, function(x) x$ML)
  ), function(x, y) {
    tryCatch({
      tips.remove = union(setdiff(x$tip.label, y$tip.label), setdiff(y$tip.label, x$tip.label))
      print(tips.remove)
      x = drop.tip(x, tip = tips.remove)
      y = drop.tip(y, tip = tips.remove)
      RF.dist(x, y, normalize = TRUE)
    }, error = function(e) {
      # Print the error message to the console
      message("An error occurred: ", conditionMessage(e))
      browser()
    })
  })
  
  colnames(rf.mat) = sub("(\\d)", "\\.\\1", colnames(rf.mat))
  rownames(rf.mat) = sub("(\\d)", "\\.\\1", rownames(rf.mat))
  rf.mat = MakeSmartMatrix(rf.mat, delimiter.row = '\\.', delimiter.col = '\\.')
  
  Heatmap2(rf.mat@matrix, 
           top_anno = HeatmapAnnotation(`Marginal\nRFD` = anno_barplot(min_max_norm(matrixStats::colMedians(rf.mat@matrix, na.rm = T)))), 
           legend.name = 'Robinson-\nFoulds\ndistance')
}

CompareTrees = function(tree1, tree2, title1 = 'tree1', title2 = 'tree2'){
  
  library(dendextend)
  library(phylogram)
  
  message('RF distance: ', RF.dist(tree1, tree2, normalize = TRUE))
  phy1 <- midpoint(tree1) #root(all.trees$all$NJ, node = 72)
  phy2 <- midpoint(tree2) #root(all.trees$all$MP, node = 49)
  ## convert phylo objects to dendrograms
  dnd1 <- as.dendrogram.phylo(phy1)
  dnd2 <- as.dendrogram.phylo(phy2)
  ## rearrange in ladderized fashion
  dnd1 <- ladder(dnd1)
  dnd2 <- ladder(dnd2)
  ## plot the tanglegram
  dndlist <- dendextend::dendlist(dnd1, dnd2)
  dndlist %>% entanglement
  dndlist %>% untangle(method = "step1side") %>% entanglement
  dendextend::tanglegram(dndlist %>% untangle(method = "step1side"), 
                         fast = TRUE, 
                         margin_inner = 5, 
                         main_left = title1, 
                         main_right = title2)
}

Col2Vec = function(x){
  x[,1] %>% setNames(rownames(x))
}

BidirectionalBestHits = function(matrix, verbose = T){
  
  row.max = as.data.frame(apply(matrix, 1, function(x) names(which.max(x)))) %>% rownames_to_column() %>% setNames(c('row', 'col'))
  row.max$value = apply(matrix, 1, max)
  col.max = as.data.frame(apply(matrix, 2, function(x) names(which.max(x)))) %>% rownames_to_column() %>% setNames(c('col', 'row'))
  col.max$value = apply(matrix, 2, max)
  df = rbind(row.max, col.max[,c(2,1,3)])
  bidirectional = df[duplicated(df), ]
  
  # No clear match!
  if(verbose){
    print('No clear match!')
    print(as.data.frame(anti_join(df, bidirectional)))
  }
  
  bidirectional
}

OnOffBcClassification = function(object, gene = 'ISL1', cutoff = 50, plot = FALSE, 
                                 plot.features = c('ISL1', 'GRM6', 'GRM5', 'GRIK1', 'PRKCA'), 
                                 remove.on.bc = NULL){
  
  on.bcs = names(which(PercentageExpressed2(subset(object, cell_class == 'BC'), gene, group.by = 'annotated')[1,,drop = TRUE] > cutoff))
  if(!is.null(remove.on.bc)) {
    # Remove if present, add if not present
    on.bcs = setdiff(on.bcs, remove.on.bc)
  }
  
  object$bc.subclass = NA
  object$bc.subclass[object$cell_class == 'BC'] = 'OFF'
  object$bc.subclass[object$annotated %in% on.bcs] = 'ON'
  # object$annotated = as.character(object$annotated)
  # object$annotated[which(object$bc.subclass == 'ON')] = paste0('on', object$annotated[which(object$bc.subclass == 'ON')])
  # object$annotated[which(object$bc.subclass == 'OFF')] = paste0('off', object$annotated[which(object$bc.subclass == 'OFF')])
  
  # Assign as cell_class2, then rename types
  object$cell_class2[which(object$bc.subclass == 'ON')] = paste0('on', object$cell_class2[which(object$bc.subclass == 'ON')])
  object$cell_class2[which(object$bc.subclass == 'OFF')] = paste0('off', object$cell_class2[which(object$bc.subclass == 'OFF')])
  
  # Rename types
  # object$cell_class2 = ExtractString(object$annotated, after = '-')
  object$annotated = paste0(object$cell_class2, '-', ExtractString(object$annotated, before = '-'))
  
  if(plot) print(DotPlot3(subset(object, cell_class == 'BC'), features = plot.features, group.by = 'annotated'))
  object
}

AddGabaGly = function(object){
  object$annotated = as.character(object$annotated)
  object$annotated[object$cell_class2 == 'gabaAC'] = paste0('gaba', object$annotated[object$cell_class2 == 'gabaAC'])
  object$annotated[object$cell_class2 == 'glyAC'] = paste0('gly', object$annotated[object$cell_class2 == 'glyAC'])
  object
}

PlotWave = function(dat.norm, roi, coef = NULL, shifted = TRUE){
  dat = data.frame(Time = as.vector(times), 
                   Roi = dat.norm[roi,], 
                   Stim = stim$V2)
  # dat.long = reshape2::melt(as.matrix(dat), ) %>% setNames(c('Time', 'Identity', 'Intensity'))
  dat.long <- dat %>% tidyr::pivot_longer(cols = c(Roi, Stim), names_to = "Identity", values_to = "Intensity")
  dat.long$Identity = factor(dat.long$Identity, levels = c('Stim', 'Roi'))
  TitlePlot(
    ggline(dat.long, 'Time', 'Intensity', size = 1, shape = NA, color = 'Identity', 
           legend = 'right', palette = c(Stim = 'grey', Roi = 'blue')) + 
    theme(axis.text.x = element_blank()), 
    paste0('ROI: ', roi, ', coef: ', round(coef, 2)))
}

CompareBinMats = function(binmat_v12, binmat_v13){
  
  sets = VennDiagram(rownames(binmat_v12), rownames(binmat_v13))
  genes.use = sets[[3]]
  
  # Intersection
  print(
    (PrettyHistogram(rowMeans(binmat_v12[genes.use,]), title = 'intersection in 1') |
     PrettyHistogram(rowMeans(binmat_v13[genes.use,]), title = 'intersection in 2')) /
    
  # Specific genes
  (PrettyHistogram(rowMeans(binmat_v12[sets[[1]],]), title = '1 only genes') |
     PrettyHistogram(rowMeans(binmat_v13[sets[[2]],]), title = '2 only genes'))
  )
  
  bin.mat = binmat_v13[genes.use,]
  hc_rows <- hclust(dist(bin.mat))
  row_order <- hc_rows$order
  ordered_row_names <- rownames(bin.mat)[row_order]
  col.order = colnames(bin.mat)
  
  print(
    (ggHeatmap(t(binmat_v12[ordered_row_names,col.order]), title = '1') + theme(axis.text.x = element_blank()))  /
    (ggHeatmap(t(binmat_v13[ordered_row_names,col.order]), title = '2') + theme(axis.text.x = element_blank()))
  )
  
  print(
    (ggHeatmap(t(binmat_v12[sets[[1]],col.order]), title = '1 only genes') + theme(axis.text.x = element_blank())) / 
    (ggHeatmap(t(binmat_v13[sets[[2]],col.order]), title = '2 only genes') + theme(axis.text.x = element_blank()))
  )
}

replace_factor <- function(factor_var, old_values, new_value) {
  
  for(old_value in old_values){
    order = levels(factor_var)
    order[order == old_value] = new_value
    order = unique(order)
    
    # Convert to character for replacement
    char_var <- as.character(factor_var)
    char_var[char_var == old_value] <- new_value
    
    # Convert back to factor (will include new level if needed)
    factor_var = factor(char_var, levels = order)
  }
  
  return(factor_var)
  
}

# Function to collapse groups in a single tree
collapse_by_centroid <- function(tree, groups) {
  
  group_names <- names(groups)
  for(i in 1:length(groups)) {
    clade_members <- groups[[i]]
    
    # Check if all members are still in the tree
    current_tips <- tree$tip.label
    available_members <- intersect(clade_members, current_tips)
    
    if(length(available_members) == 0) {
      next  # Skip if no members found
    }
    
    if(length(available_members) == 1) {
      # Single member - just rename it
      tree$tip.label[tree$tip.label == available_members[1]] <- group_names[i]
    } else {
      # Multiple members - find centroid
      dist_matrix <- cophenetic(tree)
      group_distances <- dist_matrix[available_members, available_members]
      
      # Find the centroid (tip with minimum average distance to others)
      avg_distances <- rowMeans(group_distances)
      centroid_tip <- names(which.min(avg_distances))
      
      # Keep centroid, drop others
      to_drop <- setdiff(available_members, centroid_tip)
      
      if(length(to_drop) > 0) {
        tree <- drop.tip(tree, to_drop, rooted = FALSE)
      }
      
      # Rename centroid to group name
      tree$tip.label[tree$tip.label == centroid_tip] <- group_names[i]
    }
  }
  
  # FORCE UNROOTED STATUS
  # tree$root.edge <- NULL
  # attr(tree, "order") <- NULL
  tree = unroot(tree)
  
  return(tree)
}

# Function to collapse groups in a single tree
collapse_by_MRCA <- function(tree, groups) {
  
  # Rooting is required for MRCA
  if('MG' %in% tree$tip.label) {
    tree = root(tree, 'MG')
    message('rooting tree at MG')
  }
  if('MG-1' %in% tree$tip.label) tree = root(tree, 'MG-1')
  
  group_names = names(groups)
  for(i in 1:length(groups)) {
    # Find MRCA of the group
    mrca_node <- getMRCA(tree, groups[[i]])
    if(!is.null(mrca_node)) {
      # Extract subtree and replace with single tip
      tree <- bind.tip(tree, group_names[i], where = mrca_node, edge.length = 0)
      tree <- drop.tip(tree, groups[[i]])
    }
  }
  return(tree)
}

map_support_to_edges <- function(tree) {
  edge_support <- rep(NA, nrow(tree$edge))
  node_offset <- length(tree$tip.label)
  internal_nodes = seq_along(tree$node.label) + node_offset
  
  for (i in seq_along(tree$node.label)) {
    node_num <- node_offset + i
    edge_idx <- which(tree$edge[, 2] == node_num)
    if (length(edge_idx) == 1) {
      edge_support[edge_idx] <- tree$node.label[i]
    }
  }
  
  return(edge_support)
}

LOOTrees = function(matrix, mc.cores = 1, ...){
  
  mclapply(seq_len(ncol(matrix)), function(i) {
    
    message('Removing site ', i)
    
    # Resample sites (columns)
    bs_data <- matrix[, -i]
    
    ConstructTrees(bs_data, ...)
    
  }, mc.cores = mc.cores)
  
}


LOOParsimonyTrees = function(matrix, mc.cores = 1){
  
  as.multiPhylo(mclapply(seq_len(ncol(matrix)), function(i) {
    
    # Resample sites (columns)
    bs_data <- matrix[, -i]
    
    # Convert to phyDat
    phydat_bs <- phyDat(bs_data, type = "USER", levels = c("0", "1"), compress = FALSE)
    
    # Parsimony tree search
    tree = pratchet(phydat_bs, trace = 0, minit = 100)
    
    # Node labels
    tree$node.label = (length(tree$tip.label)+1):((length(tree$tip.label)+tree$Nnode))
    
    tree
    
  }, mc.cores = mc.cores))
  
}

BootstrapParsimonyTrees = function(matrix, nIter = 100, mc.cores = 1){
  
  as.multiPhylo(mclapply(seq_len(nIter), function(iter) {
    
    # Resample sites (columns)
    bs_data <- matrix[, sample(seq_len(ncol(matrix)), replace = TRUE)]
    
    # Convert to phyDat
    phydat_bs <- phyDat(bs_data, type = "USER", levels = c("0", "1"), compress = FALSE)
    
    # Parsimony tree search
    tree = pratchet(phydat_bs, trace = 0, minit = 100)
    
    # Node labels
    tree$node.label = (length(tree$tip.label)+1):((length(tree$tip.label)+tree$Nnode))
    
    tree
    
  }, mc.cores = mc.cores))
  
}

RunIQTree2 = function(phydat, nthreads = 1, iqtree.model = 'MFP+ASC', verbose = FALSE, asr = FALSE, 
                      constraint = NULL, prefix = NULL, options = '-wsr', tool = 'iqtree2'){
  
  if(is.null(prefix)) prefix = 'temp_files/temp'
  phyfile = paste0(prefix, '.phy')
  treefile = paste0(prefix, '.treefile')
  
  write.phyDat(phydat, phyfile, format = "phylip")
  
  if(asr) {
    args = c('-s', phyfile, '-redo', '-pre', prefix, '-nt', nthreads, '-m', iqtree.model, '--ancestral', options)
  } else {
    args = c('-s', phyfile, '-redo', '-pre', prefix, '-nt', nthreads, '-m', iqtree.model, options)
  }
  
  if(!is.null(constraint)){
    args = c(args, '-g', constraint)
  }
  
  call = paste0(tool, ' ', paste0(args, collapse = ' '))
  message('Running call: ', call)
  # system2('iqtree2', 
  #         args, 
  #         stdout = ifelse(verbose, "", FALSE))
  
  system(call, show.output.on.console = ifelse(verbose, "", FALSE))
  
  tree = read.tree(treefile)
  tree
}

MultivariateR2 = function(data, cluster){
  
  # Refactor 
  cluster = factor(cluster)
  
  # Mean of entire data
  grand_mean <- colMeans(data)
  
  # Sum of squared deviances from grand mean
  total_ss <- sum(rowSums((data - rep(grand_mean, each = nrow(data)))^2))
  
  # Group means of clusters
  group_means <- aggregate(data, by = list(cluster), FUN = mean)[, -1]
  
  # Cluster sizes
  group_sizes <- table(factor(cluster))
  
  # Sum of squared deviances of group means from grand mean
  between_ss <- sum(group_sizes * rowSums((group_means - rep(grand_mean, each = length(unique(cluster))))^2))
  
  # Multivariate R^2
  r.squared = between_ss / total_ss
  
  r.squared
  
}

na_propagate <- function(data, fun, ...) {
  # Handle vectors
  if (is.vector(data) || is.factor(data)) {
    complete_idx <- !is.na(data)
    complete_data <- data[complete_idx]
    
    if (length(complete_data) == 0) {
      return(rep(NA, length(data)))
    }
    
    result <- fun(complete_data, ...)
    
    # Create full result vector with NAs
    if (is.vector(result)) {
      full_result <- rep(NA, length(data))
      full_result[complete_idx] <- result
    } else {
      # For complex objects, just return the result on complete data
      # and add an attribute showing which indices were complete
      attr(result, "complete_indices") <- which(complete_idx)
      attr(result, "original_length") <- length(data)
      full_result <- result
    }
    
    return(full_result)
  }
  
  # Handle data frames/matrices
  if (is.data.frame(data) || is.matrix(data)) {
    complete_idx <- complete.cases(data)
    complete_data <- data[complete_idx, , drop = FALSE]
    
    if (nrow(complete_data) == 0) {
      return(list(result = NULL, all_na = TRUE))
    }
    
    result <- fun(complete_data, ...)
    
    # Handle different result types
    if (is.vector(result)) {
      full_result <- rep(NA, nrow(data))
      full_result[complete_idx] <- result
    } else if ("cluster" %in% names(result)) {
      # Handle clustering objects specifically
      full_result <- result
      full_clusters <- rep(NA, nrow(data))
      full_clusters[complete_idx] <- result$cluster
      full_result$cluster <- full_clusters
    } else {
      # For other complex objects, add attributes
      attr(result, "complete_indices") <- which(complete_idx)
      attr(result, "original_nrow") <- nrow(data)
      full_result <- result
    }
    
    return(full_result)
  }
  
  # Fallback for other data types
  stop("Data must be a vector, data frame, or matrix")
}

create_dynamic_annotation <- function(col_data, 
                                      colors = NULL, 
                                      show_legend = TRUE,
                                      position = "top",
                                      ...) {
  # Start with the annotation data
  anno_list <- as.list(col_data)
  
  # Add color mappings if provided
  if (!is.null(colors)) {
    for (col_name in names(colors)) {
      if (col_name %in% names(anno_list)) {
        # Add color mapping for this column
        anno_list[[col_name]] <- anno_simple(
          anno_list[[col_name]], 
          col = colors[[col_name]]
        )
      }
    }
  }
  
  # Add other parameters
  anno_list$show_legend <- show_legend
  
  # Add any additional arguments
  extra_args <- list(...)
  anno_list <- c(anno_list, extra_args)
  
  # Create annotation based on position
  if (position == "top") {
    return(do.call(HeatmapAnnotation, anno_list))
  } else if (position == "left" || position == "right") {
    return(do.call(rowAnnotation, anno_list))
  } else {
    stop("Position must be 'top', 'left', or 'right'")
  }
}

hamming_distance <- function(x, y) {
  sum(x != y)
}

# Create Hamming distance matrix
hamming_matrix <- function(data) {
  n <- nrow(data)
  dist_matrix <- matrix(0, n, n)
  
  for(i in 1:n) {
    for(j in 1:n) {
      dist_matrix[i, j] <- hamming_distance(data[i, ], data[j, ])
    }
  }
  
  rownames(dist_matrix) <- rownames(data)
  colnames(dist_matrix) <- rownames(data)
  return(dist_matrix)
}

FlipOrder = function(list, by){
  mat = matrix(list, ncol = by, byrow = TRUE)
  as.vector(mat)
}

LoadRefFiles = function(path = '../../reference_files'){
  neuropeptides = read.table(paste0(path, "/complete_neuropeptide_symbols.txt"))$V1
  hk_genes = read.table(paste0(path, "/Eisenberg_Cell_2013_hk_genes.txt"))$V1
  tf_genes = read.csv(paste0(path, "/cisbp_20220920.csv"))$Name
}

MatchDimensions = function(x, y){
  x[match(rownames(y), rownames(x)), match(colnames(y), colnames(x))]
}
  
CheckDimensions = function(x, y){
  # Check dimensions
  dim_match <- identical(dim(x), dim(y))
  
  # Check row names
  rownames_match <- identical(rownames(x), rownames(y))
  
  # Check column names
  colnames_match <- identical(colnames(x), colnames(y))
  
  # All must be TRUE to fully match structure
  all(dim_match, rownames_match, colnames_match)
}

DotPlotAndBinarization = function(object, bin.mat, features){
  (DotPlot3(object, features = features, coord.flip = FALSE) + 
          # theme_cowplot(font_size = 10)+
          NoLegend() + theme(axis.text.x = element_blank())) /
  (ggHeatmap(t(bin.mat[features,unique(object$annotated)]), 
                  high_color = '#584B9FFF', 
                  legend.name = 'Scaled\nexpression') + 
          theme_cowplot(font_size = 10)+
          RotatedAxis() + 
          NoLegend())
  
  # return(invisible(NULL))
}

gmm_loglikelihood <- function(x, means, covariances, weights, log = TRUE) {
  # x: n x d matrix
  # means: list of k mean vectors (each of length d)
  # covariances: list of k d x d covariance matrices
  # weights: numeric vector of length k, should sum to 1
  # log: whether to return log-likelihoods (default = TRUE)
  
  n <- nrow(x)
  k <- length(means)
  d <- ncol(x)
  
  log_lik_mat <- matrix(NA, n, k)
  
  for (i in seq_len(k)) {
    mu <- means[[i]]
    sigma <- covariances[[i]]
    sigma_inv <- solve(sigma)
    log_det <- determinant(sigma, logarithm = TRUE)$modulus
    
    xc <- sweep(x, 2, mu)
    mahal <- rowSums((xc %*% sigma_inv) * xc)
    
    log_lik_mat[, i] <- log(weights[i]) - 0.5 * (d * log(2 * pi) + log_det + mahal)
  }
  
  if (log) {
    # Log-sum-exp trick for numerical stability
    max_log <- apply(log_lik_mat, 1, max)
    log_lik <- max_log + log(rowSums(exp(log_lik_mat - max_log)))
    return(log_lik)
  } else {
    return(rowSums(exp(log_lik_mat)))
  }
}

safe_scale <- function(x) {
  apply(x, 2, function(col) {
    if (all(is.na(col))) return(rep(NA, length(col)))  # skip all-NA columns
    mu <- mean(col, na.rm = TRUE)
    sd_val <- sd(col, na.rm = TRUE)
    if (sd_val > 0) {
      (col - mu) / sd_val
    } else {
      col - mu  # center only
    }
  })
}

Bootstrap = function(data, n = 500, expansion = 10, noise.alpha = 0.1){
  # noise.sd = apply(data, 2, sd)
  data = safe_scale(data)
  noise.vector = do.call(cbind, lapply(seq_len(ncol(data)), function(sd) rnorm(nrow(data)*expansion, sd = noise.alpha)))
  # data.boot = data[sample(seq_len(nrow(data)), n, replace = TRUE),] + noise.vector
  data.boot = data[rep(1:nrow(data), each = expansion),] + noise.vector
  data.boot
}

CompareBinarization = function(object, gene, snf){
  order = unique(object$annotated)
  (ggHeatmap(rbind(pct = binmat_v12[gene,order], 
                   expr = binmat_v13[gene,order], 
                   score = binmat_v14[gene,order], 
                   mv = snf[order]), 
                   high_color = '#584B9FFF', title = NULL) + NoLegend() + theme(axis.text.x = element_blank())) /
    (DotPlot3(objectList$Shark, features = c(gene)))
}

FindBackboneLength = function(tree, winsorize = FALSE){
  tip_nodes <- 1:Ntip(tree)
  edges <- tree$edge
  is_internal_edge <- !(edges[, 2] %in% tip_nodes)
  if(winsorize) {
    edge.lengths = DescTools::Winsorize(tree$edge.length)
  } else {
    edge.lengths = tree$edge.length
  }
  backbone.length <- sum(edge.lengths[is_internal_edge])
  backbone.length
}

list_objects <- function(env = .GlobalEnv, n = 10) {
  objs <- ls(env)
  obj_sizes <- sapply(objs, function(x) object.size(get(x, envir = env)))
  obj_sizes <- sort(obj_sizes, decreasing = TRUE)
  data.frame(
    Name = names(obj_sizes),
    Size = format(obj_sizes, units = "auto")
  )[1:n, ]
}

RenameAnnotations = function(object){
  object$annotated = paste0(object$cell_class2, '-', 
                            ifelse(ExtractString(object$annotated, before = '-') == '', 1,
                                   ExtractString(object$annotated, before = '-')))
  object
}

MeanDistance = function(dist_matrix, summarize = TRUE, plot = FALSE){
  
  sdist = MakeSmartMatrix(as.matrix(dist_matrix), delimiter.row = '-')
  sdist
  
  if(plot) print(SmartHeatmap(sdist))
  
  temp = sdist@matrix
  group = sdist@col.data$x
  rownames(temp) = group
  colnames(temp) = group
  diag(temp) = NA
  dist_df <- as.data.frame(as.table(temp)) %>%
    dplyr::rename(type1 = Var1, type2 = Var2, dist = Freq) %>%
    na.omit() %>% 
    group_by(type1, type2)
  
  if(summarize){
    dist_df %>% summarise(mean_dist = mean(dist), 
                          sd_dist = sd(dist),
                          .groups = "drop")
  } else {
    dist_df
  }
    
}

SquarePlot = function(myplot, legend.ratio = 0.3){
  legend <- get_legend(myplot)
  square_plot <- myplot + theme(legend.position = "none")
  
  plt = plot_grid(NULL, NULL, NULL, 
                  NULL, square_plot, legend, 
                  NULL, NULL, NULL, 
                  rel_widths = c(legend.ratio, 1, legend.ratio), 
                  rel_heights = c(legend.ratio, 1, legend.ratio), 
                  nrow = 3, ncol = 3)
  
  return(plt)
}

RemoveBatches = function(object, group.by, remove.batches = 60){
  object@meta.data[,group.by] = as.character(object@meta.data[,group.by])
  batches.remove = names(which(table(object@meta.data[,group.by]) < remove.batches))
  message("Removing batches: ", paste0(batches.remove, collapse = ", "))
  object[,!object@meta.data[,group.by] %in% batches.remove]
}

DrawBipartiteGraph = function(edges, top = 20){
  
  library(igraph)
  
  edges = edges %>% arrange(desc(bit_score)) %>% head(top)
  
  # Combine gene names into one vector of nodes
  all_nodes <- unique(c(edges$gene1, edges$gene2))
  
  # Create a logical vector indicating which set each node belongs to (TRUE = gene1 set)
  type <- all_nodes %in% edges$gene1  # required for bipartite graphs in igraph
  
  # Create the graph
  g <- graph_from_data_frame(edges, vertices = data.frame(name = all_nodes, type = type), directed = FALSE)
  
  # Add bit score as an edge attribute
  E(g)$weight <- edges$bit_score
  
  print(
    plot(
      g,
      layout = layout_as_bipartite,
      vertex.color = ifelse(V(g)$type, "lightblue", "salmon"),
      vertex.label.cex = 0.8,
      edge.width = E(g)$weight / max(E(g)$weight) * 5
    )
  )
  
  return(edges)
}

Read10X_custom = function(data.dir, gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE) {
  full.data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, "barcodes.tsv")
    gene.loc <- file.path(run, "genes.tsv")
    features.loc <- file.path(run, "features.tsv.gz")
    matrix.loc <- file.path(run, "matrix.mtx")
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing. Expecting ", 
           basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", 
           basename(path = matrix.loc))
    }
    data <- t(readMM(file = matrix.loc))
    if (has_dt) {
      cell.barcodes <- as.data.frame(data.table::fread(barcode.loc, 
                                                       header = FALSE))
    }
    else {
      cell.barcodes <- read.table(file = barcode.loc, header = FALSE, 
                                  sep = "\t", row.names = NULL)
    }
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    }
    else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    if (has_dt) {
      feature.names <- as.data.frame(data.table::fread(ifelse(test = pre_ver_3, 
                                                              yes = gene.loc, no = features.loc), header = FALSE))
    }
    else {
      feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                                yes = gene.loc, no = features.loc), header = FALSE, 
                                  stringsAsFactors = FALSE)
    }
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as.sparse(x = list_of_data[[j]])
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}

GraphToList = function(graph){
  apply(graph, 1, function(x) names(which(x > 0)))
}

GraphToDF = function(graph){
  
  list = as.data.frame(summary(graph))
  list$i = rownames(graph)[list$i]
  list$j = colnames(graph)[list$j]
  
  list %>% setNames(c('species1', 'species2', 'bitscore'))
  
}

DropSCT = function(object){
  object[['SCT']] = NULL
  object
}

DropIntegrated = function(object){
  object[['integrated']] = NULL
  object
}

ComputeCopheneticDistance = function(tree, metadata, plot = FALSE){
  dist_matrix <- cophenetic.phylo(tree)
  sdist = MakeSmartMatrix(as.matrix(dist_matrix), delimiter.row = '-')
  sdist@col.data$z = metadata[match(rownames(dist_matrix), metadata$annotated), 'cell_class2']
  if(plot) print(SmartHeatmap(sdist))
  
  temp = sdist@matrix
  group = sdist@col.data$z
  rownames(temp) = group
  colnames(temp) = group
  diag(temp) = NA
  dif.df <- as.data.frame(reshape2::melt(temp)) %>%
    setNames(c('sample1', 'sample2', 'distance')) %>%
    na.omit() %>% 
    group_by(sample1, sample2) %>%
    summarise(distance = mean(distance), .groups = "drop")
  
  acast(dif.df, sample1 ~ sample2, value.var = "distance")
}

plot_tree <- function(tree, annotation_df = NULL, layout = 'equal_angle', 
                      category_regex = "-", use.ape = TRUE, background.color = 'grey', 
                      title = NULL, max.bl = NULL, return.df = FALSE, 
                      show.tip.label = FALSE) {
  
  library(ape)
  library(ggtree)
  library(treeio)
  
  color.key = data.frame(old.names = c('gabaAC', 'glyAC', 'achAC', 'PR', 'Rod', 'Cone', 'MG', 'HC', 'BC', 'RGC'), 
                         new.names = c('cyan', 'cyan4', 'cyan', 'red', 'red', 'red', 'magenta', 'gold', 'chartreuse2', 'blue'))
  
  # Read and convert tree
  # tree <- read.tree(newick_file)
  td <- as.treedata(tree)
  tib = as.tibble(td)
  
  # Extract labels and infer category from label
  if(is.null(annotation_df)){
    tib <- tib %>%
      mutate(
        cell_class2 = ExtractString(label, after = '-')
      )
  } else {
    tib <- tib %>%
      left_join(annotation_df, by = c("label" = "annotated"))
  }
  
  if(!is.null(max.bl)){
    tree$edge.length[tree$edge.length > max.bl] = max.bl
  }
  
  # Plot tree with daylight layout and colored branches
  if(use.ape){
    # edge_colors
    tib$node_id = paste0(tib$parent, '-', tib$node)
    ape.data = as.data.frame(tree$edge)
    ape.data$node_id = paste0(ape.data[,1], '-', ape.data[,2])
    tib.sorted = tib[match(ape.data$node_id, tib$node_id),]
    tib.sorted$color = convert_values(tib.sorted$cell_class2, 
                                      key_table = color.key)
    tib.sorted$color[is.na(tib.sorted$color)] = background.color
    if(return.df) return(tib.sorted)
    plot.phylo(tree, type = 'unrooted', show.tip.label = show.tip.label, edge.color = tib.sorted$color, main = title)
  } else {
    ggtree(as.treedata(tib), layout = layout, aes(color = cell_class2)) +
      geom_tiplab(size = 3) +
      coord_cartesian(clip = 'off')+
      theme(legend.position = "right")
  }
  
}

ComputeSimilarity = function(object, group.by = 'annotated', batch = 'animal'){
  
  # expressed.genes = names(which(PercentageExpressed2(object, rownames(object)) > 0))
  expressed.genes = names(which(apply(object@assays$RNA@data, 1, function(x) var(x) != 0)))
  res = do.call(rbind, lapply(unique(object@meta.data[[batch]]), function(this.batch){
    cor.mat = CorrelationHeatmap(subset(object, animal == this.batch), 
                                 group.by = group.by, 
                                 annotate.by = 'cell_class2', 
                                 return.matrix = TRUE, 
                                 features = expressed.genes, 
                                 label = TRUE)
    smatrix = MakeSmartMatrix(cor.mat, delimiter.row = '-', delimiter.col = '-')
    
    # Average similarity of gaba AC to RGC
    scores.achAC.rgc = as.vector(smatrix[smatrix@row.data$x == 'achAC', smatrix@col.data$x == 'RGC']@matrix)
    scores.gabaAC.rgc = as.vector(smatrix[smatrix@row.data$x == 'gabaAC', smatrix@col.data$x == 'RGC']@matrix)
    scores.glyAC.rgc = as.vector(smatrix[smatrix@row.data$x == 'glyAC', smatrix@col.data$x == 'RGC']@matrix)
    scores.bc.rgc = as.vector(smatrix[smatrix@row.data$x == 'BC', smatrix@col.data$x == 'RGC']@matrix)
    
    print((PrettyHistogram2(scores.gabaAC.rgc, scores.bc.rgc, vline = c(mean(scores.gabaAC.rgc), mean(scores.bc.rgc))) + scale_x_continuous(limits = c(min(smatrix@matrix),1))) /
            (PrettyHistogram2(scores.glyAC.rgc, scores.bc.rgc, vline = c(mean(scores.glyAC.rgc), mean(scores.bc.rgc))) + scale_x_continuous(limits = c(min(smatrix@matrix),1))))
    
    ach.score = mean(fisher_z(scores.achAC.rgc)) - mean(fisher_z(scores.bc.rgc))
    gaba.score = mean(fisher_z(scores.gabaAC.rgc)) - mean(fisher_z(scores.bc.rgc))
    gly.score = mean(fisher_z(scores.glyAC.rgc)) - mean(fisher_z(scores.bc.rgc))
    # message('gaba AC score:', gaba.score)
    # message('gly AC score:', gly.score)
    data.frame(batch = this.batch, 
               ach.score = ach.score,
               gaba.score = gaba.score, 
               gly.score = gly.score)
  }))
  
  res
}

fisher_z <- function(r) {
  # if (any(abs(r) >= 1)) stop("Correlation coefficients must be between -1 and 1 (exclusive).")
  0.5 * log((1 + r) / (1 - r))
}

HighlightCluster = function(object, group.by, cluster){
  DimPlot(object, cells.highlight = WhichCells(object, expr = annotated == cluster))
}

SubsetSeuratGenes2 = function(object, features){
  object@assays$RNA@counts = object@assays$RNA@counts[features, ]
  object@assays$RNA@data = object@assays$RNA@data[features, ]
  object
}

GetR2 = function(object){
  type = adonis2(t(GetAssayData(object)) ~ type, 
          data = object@meta.data, 
          method = "euclidean", 
          by = "margin")$R2[[1]]
  species = adonis2(t(GetAssayData(object)) ~ species, 
          data = object@meta.data, 
          method = "euclidean", 
          by = "margin")$R2[[1]]
  technology = adonis2(t(GetAssayData(object)) ~ technology, 
          data = object@meta.data, 
          method = "euclidean", 
          by = "margin")$R2[[1]]
  data.frame(type = type, species = species, technology = technology)
}

AddModuleScoreDebug = function(object, features, pool = NULL, nbin = 24, ctrl = 1, 
          k = FALSE, assay = NULL, name = "Cluster", seed = 1, search = FALSE, 
          ...) 
{
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
                                         i))
    }
    cluster.length <- length(x = features)
  }
  else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning("The following features are not present in the object: ", 
                paste(missing.features, collapse = ", "), 
                ifelse(test = search, yes = ", attempting to find updated synonyms", 
                       no = ", not searching for symbol synonyms"), 
                call. = FALSE, immediate. = TRUE)
        if (search) {
          tryCatch(expr = {
            updated.features <- UpdateSymbolList(symbols = missing.features, 
                                                 ...)
            names(x = updated.features) <- missing.features
            for (miss in names(x = updated.features)) {
              index <- which(x == miss)
              x[index] <- updated.features[miss]
            }
          }, error = function(...) {
            warning("Could not reach HGNC's gene names database", 
                    call. = FALSE, immediate. = TRUE)
          })
          missing.features <- setdiff(x = x, y = rownames(x = object))
          if (length(x = missing.features) > 0) {
            warning("The following features are still not present in the object: ", 
                    paste(missing.features, collapse = ", "), 
                    call. = FALSE, immediate. = TRUE)
          }
        }
      }
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }
  # if (!all(LengthCheck(values = features))) {
  #   warning(paste("Could not find enough features in the object from the following feature lists:", 
  #                 paste(names(x = which(x = !LengthCheck(values = features)))), 
  #                 "Attempting to match case..."))
  #   features <- lapply(X = features.old, FUN = CaseMatch, 
  #                      match = rownames(x = object))
  # }
  # if (!all(LengthCheck(values = features))) {
  #   stop(paste("The following feature lists do not have enough features present in the object:", 
  #              paste(names(x = which(x = !LengthCheck(values = features)))), 
  #              "exiting..."))
  # }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
    ])
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

FindControlGenes = function(object, features, pool = NULL, nbin = 24, ctrl = 1){
  
  pool <- pool %||% rownames(x = object)
  assay.data = GetAssayData(object = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  features = lapply(features, function(x) intersect(x, rownames(object)))
  cluster.length = length(features)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  # ctrl.use = vector()
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.use
}

SmartHeatmap2 = function(smatrix, top_anno = NULL, left_anno = NULL, colors = NULL, show_annotation_legend = TRUE, ...){
  
  if(!is.null(top_anno)){
    top_annotation = create_dynamic_annotation(smatrix@col.data[, top_anno, drop = FALSE], 
                                         colors = colors, 
                                         show_legend = show_annotation_legend)
  } else{
    top_annotation = NULL
  }
  
  if(!is.null(left_anno)){
    left_annotation = create_dynamic_annotation(smatrix@row.data[, left_anno, drop = FALSE], 
                                         colors = colors, 
                                         position = 'left',
                                         show_legend = show_annotation_legend)
  } else {
    left_annotation = NULL
  }
  
  Heatmap2(smatrix@matrix, 
           top_annotation = top_annotation, 
           left_annotation = left_annotation,
           ...)
}

SmartHeatmap = function(smatrix, annotate.by = NULL, colors = NULL, show_annotation_legend = TRUE, ...){
  
  # if(is.null(annotation_names)) names = paste0('anno', 1:4)
  
  # if(ncol(smatrix@col.data) == 3){
  #   top_anno = HeatmapAnnotation(anno1 = smatrix@col.data[,1],
  #                                anno2 = smatrix@col.data[,2],
  #                                anno3 = smatrix@col.data[,3])
  # } else if(ncol(smatrix@col.data) == 2){
  #   top_anno = HeatmapAnnotation(anno1 = smatrix@col.data[,1],
  #                                anno2 = smatrix@col.data[,2])
  # }
  
  # selected_cols <- c("anno1", "anno2", "anno3")  # or any vector of column names
  # anno_list <- as.list(smatrix@col.data[, annotate.by, drop = FALSE])
  # top_anno <- do.call(HeatmapAnnotation, anno_list, show_legend = show_annotation_legend)
  
  # annotations_list <- lapply(1:ncol(smatrix@col.data), function(i) {
  #   HeatmapAnnotation(name = colnames(smatrix@col.data)[i], 
  #                     annotation = smatrix@col.data[, i])
  # })
  
  top_anno = create_dynamic_annotation(smatrix@col.data[, annotate.by, drop = FALSE], 
                                       colors = colors, 
                                       show_legend = show_annotation_legend)
  
  Heatmap2(smatrix@matrix, 
           top_annotation = top_anno, 
           ...)
}


CombineAndRefactor = function(list1, list2, reverse = FALSE){
  if(reverse) {
    factor(paste0(list1, '-', list2), 
           levels = outer(levels(list2), levels(list1), function(x, y) paste0(y, '-', x)) %>% as.character)
  } else {
    factor(paste0(list1, '-', list2), 
           levels = outer(levels(list1), levels(list2), function(x, y) paste0(x, '-', y)) %>% as.character)
  }
}

population_variance <- function(x) {
  mean((x - mean(x))^2)
}

plot_tree_with_bars <- function(tree, data_vector, bar_color = "steelblue") {
  
  library(ggtree)
  if (is.null(names(data_vector))) {
    stop("data_vector must be a named numeric vector with names matching tip labels of the tree.")
  }
  
  # Create data frame
  bar_data <- data.frame(label = names(data_vector), value = as.numeric(data_vector))
  
  # Base tree
  p <- ggtree(tree) +
    geom_tree() +
    geom_tiplab(align = TRUE, linetype = "dotted") +
    geom_facet(
      data = bar_data,
      mapping = aes(x = value),
      orientation = 'y', 
      geom = geom_col,
      panel = "Barplot",
      width = 1,
      fill = bar_color
    ) +
    theme_tree2()
  
  return(p)
}

PrettyBarplot = function(df, x, y, fill = 'lightgrey', remove.x = FALSE, ...){
  ggbarplot(df, x, y, fill = fill, ...) + 
    RotatedAxis()+
    # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))+
    scale_y_continuous(expand = expansion(mult = c(0,0.01)))+
    theme(plot.title = element_text(hjust = 0.5))+
    {if(remove.x) theme(axis.title.x = element_blank())}+
    coord_cartesian(clip = 'off')+
    ArialFont()
}


ArialFont <- function(...) {
  theme(
    text = element_text(family = "ArialMT", ...)
  )
}

# ArialFont <- function() {
#   base <- theme_get()
#   
#   maybe_text <- function(el) {
#     if (is.null(el) || inherits(el, "element_blank")) {
#       el
#     } else {
#       element_text(
#         family = "ArialMT",
#         colour = el$colour %||% "black",
#         face   = el$face,
#         size   = el$size,
#         hjust  = el$hjust,
#         vjust  = el$vjust,
#         angle  = el$angle,
#         lineheight = el$lineheight
#       )
#     }
#   }
#   
#   theme(
#     text          = maybe_text(base$text),
#     axis.text     = maybe_text(base$axis.text),
#     axis.title    = maybe_text(base$axis.title),
#     axis.text.x   = maybe_text(base$axis.text.x),
#     axis.title.x = maybe_text(base$axis.title.x),
#     axis.text.y   = maybe_text(base$axis.text.y),
#     axis.title.y = maybe_text(base$axis.title.y)
#   )
# }


# ArialFont = function(){
#   text_settings = element_text(family = "ArialMT", colour = 'black')
#   theme(text = text_settings, 
#         axis.text = text_settings, 
#         axis.title = text_settings,
#         axis.text.x = text_settings,
#         axis.title.x = text_settings,
#         axis.text.y = text_settings, 
#         axis.title.y = text_settings)
# }
  
interweave <- function(v) {
  n <- ceiling(length(v) / 2)
  # first_half <- v[1:n]
  # second_half <- v[(n + 1):length(v)]
  first_half = 1:n
  second_half = (n+1):(2*n)
  
  # Interleave using mapply, recycling shorter half if needed
  result <- c(rbind(first_half, second_half))
  
  # Trim
  result = result[1:length(v)]
  v[result]
  # result[!is.na(result)]  # Remove any NAs that appear due to unequal length
}

flip_var = function(df, var){
  df[[var]] = factor(df[[var]], levels = rev(levels(df[[var]])))
  df
}

coord_flip_discrete <- function(p, discrete_axis) {
  
  # Reverse factor levels of the discrete axis
  if(is.factor(p$data[[discrete_axis]])){
    message('detecting factor')
    p$data[[discrete_axis]] <- factor(as.character(p$data[[discrete_axis]]),
                                      levels = rev(levels(p$data[[discrete_axis]])))
  } else {
    message('detecting non-factor')
    p$data[[discrete_axis]] <- factor(p$data[[discrete_axis]],
                                      levels = rev(unique(p$data[[discrete_axis]])))
  }
  
  p + coord_flip()
}

merge_by_rownames_intersect <- function(df_list) {
  common_rownames <- Reduce(intersect, lapply(df_list, rownames))
  trimmed_list <- lapply(df_list, function(df) df[common_rownames, , drop = FALSE])
  merged <- do.call(cbind, trimmed_list)
  rownames(merged) <- common_rownames
  return(merged)
}

merge_by_rownames <- function(df_list) {
  merged <- Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE), df_list)
  rownames(merged) <- merged$Row.names
  merged <- merged[, -1, drop = FALSE]
  return(merged)
}

violin_plot_matrix_gg <- function(mat, fill = "lightblue") {
  if (!is.matrix(mat)) stop("Input must be a matrix")
  
  # Convert matrix to long format
  df <- melt(mat, varnames = c("row", "column"), value.name = "value")
  
  # Remove NA values
  df <- df[!is.na(df$value), ]
  
  # Plot
  ggplot(df, aes(x = column, y = value)) +
    geom_violin(fill = fill, color = "black", trim = FALSE) +
    theme_minimal() +
    xlab("") +
    ylab("Value")
}

violin_plot_matrix <- function(mat, col = "lightblue", border = "black", add_axis = TRUE, ...) {
  # Check if input is a matrix
  if (!is.matrix(mat)) stop("Input must be a matrix")
  
  # Set up empty plot
  x_positions <- seq_len(ncol(mat))
  y_range <- range(mat, na.rm = TRUE)
  plot(NA, xlim = range(x_positions) + c(-0.5, 0.5), ylim = y_range,
       xaxt = "n", xlab = "", ylab = "Value", ...)
  
  # Loop over columns to draw violins
  for (i in seq_len(ncol(mat))) {
    data <- mat[, i]
    data <- data[!is.na(data)]
    if (length(data) < 2) next  # skip if not enough data
    
    d <- density(data, na.rm = TRUE)
    d$y <- d$y / max(d$y) * 0.4  # scale width
    
    polygon(c(i - d$y, rev(i + d$y)),
            c(d$x, rev(d$x)),
            col = col, border = border)
  }
  
  # Add x-axis labels
  if (add_axis) {
    axis(1, at = x_positions, labels = colnames(mat))
  }
}

min_max_norm <- function(mat) {
  apply(mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) {
      return(rep(0, length(x)))  # Avoid divide-by-zero if constant column
    } else {
      (x - rng[1]) / diff(rng)
    }
  })
}

pad_and_bind <- function(vec1, vec2, colnames = c("vec1", "vec2")) {
  max_len <- max(length(vec1), length(vec2))
  vec1_padded <- c(vec1, rep(NA, max_len - length(vec1)))
  vec2_padded <- c(vec2, rep(NA, max_len - length(vec2)))
  mat <- cbind(vec1_padded, vec2_padded)
  colnames(mat) <- colnames
  return(mat)
}

quantile_map <- function(target, source) {
  # Get empirical quantiles from source
  probs <- ecdf(source)(source)  # Empirical CDF of source
  probs_target <- ecdf(source)(target)  # Probabilities for target values
  
  # Get sorted values and corresponding probs from source
  source_sorted <- sort(source)
  source_probs <- ecdf(source_sorted)(source_sorted)
  
  # Remove duplicate probabilities (needed for interpolation)
  unique_probs <- !duplicated(source_probs)
  source_probs <- source_probs[unique_probs]
  source_sorted <- source_sorted[unique_probs]
  
  # Interpolate to match quantiles
  target_mapped <- approx(x = source_probs, y = source_sorted, xout = probs_target, rule = 2)$y
  
  return(target_mapped)
}


MakeSSObject = function(table){
  object = CreateSeuratObject(table)
  object = NormalizeData(object)
  
  # Replace data with log1p TPM
  object@assays$RNA@data = as.sparse(log1p(object@counts))
  object = ScaleData(object)
  object
}

LoadACOrtho = function(){
  readRDS('../../Ortho_Objects/vertebrateAC_BIGSEURAT.rds')
}

SaveACOrtho = function(object){
  saveRDS(object, '../../Ortho_Objects/vertebrateAC_BIGSEURAT.rds')
}

HistogramAnnotation <- function(plot, label, x_value, color = 'red', linewidth = 1, ...) {
  data = plot$data$variable
  max_y = max(table(cut(data, breaks = 20)))
  plot + 
    # geom_vline(xintercept = x_value, color = "red", linetype = "dashed", size = 1.2) +  # Add the vertical line
    geom_segment(aes(x = x_value, y = 0, xend = x_value, yend = max_y), color = color, linewidth = linewidth, ...) +  # Add the vertical line segment
    annotate("text", x = x_value, y = max_y, label = label, color = color, hjust = 0.5, vjust = -1) # Add the label
}

ggHeatmap = function(mat, title = NULL, low_color = "white", high_color = "red", 
                     legend.name = NULL, xlab = NULL, ylab = NULL, 
                     label = FALSE, border = FALSE, min.value = NULL, 
                     max.value = NULL) {
  
  # Convert matrix to data frame for ggplot
  df <- reshape2::melt(mat)
  colnames(df) <- c("Y", "X", "Value")
  
  # Maintain order
  df$X = factor(df$X, levels = unique(sort(df$X)))
  df$Y = factor(df$Y, levels = rev(unique(sort(df$Y))))
  
  if(is.null(min.value)) min.value = min(df$Value) 
  if(is.null(max.value)) max.value = max(df$Value)
  
  # print(min.value)
  # print(max.value)

  # Generate heatmap
  p <- TitlePlot(ggplot(df, aes(x = X, y = Y, fill = Value)) +
    {if(border) geom_tile(color = 'black')}+
    {if(!border) geom_tile()}+
    {if(label) geom_text(aes(label = round(Value, 2)), color = "black", size = 4)}+  # Labels in tiles
    scale_y_discrete(expand = expansion(mult = c(0,0)))+
    scale_x_discrete(expand = expansion(mult = c(0,0)))+
    scale_fill_gradient2(name = legend.name, 
                         low = low_color, high = high_color, 
                         limits = c(min.value-0.01, max.value+0.01), 
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'black'),
      axis.text.y = element_text(color = 'black'),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.key.height = unit(0.8, "cm")
    ) +
    labs(x = xlab, y = ylab, fill = "Value"), title = title)
  
  return(p)
}

colorRampMinMax = function(data, cols){
  colorRamp2(c(min(data), max(data)), cols)
}

add_rectangle_annotation <- function(p, index, axis = "row", alpha = 0.3, color = 'red') {
  if (!(axis %in% c("column", "row"))) {
    stop("axis must be either 'column' or 'row'")
  }
  
  # Extract plot data
  plot_data <- ggplot_build(p)$data[[1]]  # Extract first layer (geom_tile)
  
  # Get x and y ranges
  x_range <- range(plot_data$x)
  y_range <- range(plot_data$y)
  
  if (axis == "column") {
    rect_data <- data.frame(xmin = index - 0.5, 
                            xmax = index + 0.5, 
                            ymin = y_range[1], 
                            ymax = y_range[2])
  } else {  # axis == "row"
    rect_data <- data.frame(xmin = x_range[1] - 0.5, 
                            xmax = x_range[2] + 0.5, 
                            ymin = index - 0.5, 
                            ymax = index + 0.5)
  }
  
  # Add annotation layer
  p + # geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                # fill = NA, color = "red", size = 1)
    annotate("rect",
               ymin = rect_data$ymin, 
               ymax = rect_data$ymax,
               xmin = rect_data$xmin, 
               xmax = rect_data$xmax,
               alpha = alpha, 
               fill = color)
}

ReadGTF = function(path, as.df = TRUE, ...){
  gtf = rtracklayer::import(path, format = "gtf", ...)
  if(as.df) as.data.frame(gtf) else gtf
}

BalanceClasses = function(object, group.by){
  
  # refactor to avoid removed factors
  object@meta.data[[group.by]] = factor(object@meta.data[[group.by]])
  tabulation = table(object@meta.data[[group.by]])
  min.class = names(tabulation)[which.min(tabulation)]
  message("Downsampling to ", min(tabulation), ' cells in each ', group.by, ' class')
  
  DownsampleSeurat(object, group.by = group.by, size = min(tabulation))
}

FindKneePoint <- function(vector, plot = FALSE, smooth = FALSE) {
  
  data = as.matrix(data.frame(seq_along(vector), 
                              sort(as.matrix(vector))) %>% setNames(c('x', 'y')))
  n <- nrow(data)
  
  if(smooth){
    
    # Apply smoothing
    # spline_fit <- smooth.spline(data[,1], data[,2])
    # data[,2] = predict(spline_fit, data[,1])$y
    smoothed <- ksmooth(data[,1], data[,2], kernel = "normal", x.points = data[,1], bandwidth = 2)$y
    cbind(data[,2], smoothed)
    data[,2] = smoothed
  }
  
  # Calculate distances from each point to the line formed by the first and last points
  distances <- sapply(2:(n - 1), function(i) {
    # Vector from first to last point
    line_vector <- c(data[n, 1] - data[1, 1], data[n, 2] - data[1, 2])
    # Vector from first point to current point
    point_vector <- c(data[i, 1] - data[1, 1], data[i, 2] - data[1, 2])
    
    # Cross product magnitude (2D equivalent)
    cross_product_magnitude <- abs(det(matrix(c(point_vector, line_vector), nrow = 2)))
    
    # Distance calculation
    distance <- cross_product_magnitude / sqrt(sum(line_vector^2))
    return(distance)
  })
  
  # Find index of maximum distance
  max_index <- which.max(distances) + 1 # correct for shift in distances
  
  # The cutoff point could be the knee point or knee point + 1
  d1 = diff(data[,2])
  if(d1[max_index] > d1[max_index-1]) knee.point = max_index+1 else knee.point = max_index
  
  if(plot){
    plot(data)
    points(knee.point, data[knee.point,2], col = "red", pch = 19, cex = 1.5)
  }
  
  # Return the knee point
  knee.point
}

# Uses derivative method
find_elbow = function(vector, plot = FALSE){
  
  index = seq_along(vector)
  
  # add two padding since each derivate removes one
  vector = c(vector[1]-diff(vector)[1], vector)
  vector = c(vector[1]-diff(vector)[1], vector)
  
  d1 = diff(vector)
  
  # second derivative
  d2 = diff(d1)
  
  # Find the max of the second derivative
  elbow.index = which.max(d2)
  
  if(plot){
    plot(c(-1, 0, index), vector)
    abline(v = elbow.index, col = "blue", lwd = 2, lty = 2)
  }
  
  elbow.index
}

CollapseBinarization = function(binarization, expansion, data){
  pct.vector = as.numeric(tapply(binarization, rep(1:nrow(data), each = expansion), mean) > 0.5) %>% 
    setNames(rownames(data))
  cluster_means = tapply(data[,1], pct.vector, mean)
  if(cluster_means[['1']] < cluster_means[['0']]) {
    pct.vector = -1*pct.vector + 1
  }
  pct.vector %>% setNames(rownames(data))
}

BinarizeExpressionMV2 = function(gene, 
                                 avg.expr, 
                                 pct.expr, 
                                 noise.alpha = 0.25, 
                                 plot = FALSE, 
                                 verbose = FALSE, 
                                 seed = 12345, 
                                 expansion = 10, 
                                 bootstrap = TRUE, 
                                 components = 3, 
                                 fast = FALSE, 
                                 compute.r.squared = TRUE, 
                                 nstart = 1){
  
  library(mclust)
  library(MASS)
  
  set.seed(seed)
  avg.expr = as.data.frame(avg.expr)
  pct.expr = as.data.frame(pct.expr)
  data = cbind(t(avg.expr[gene,]), t(pct.expr[gene,])) %>% as.data.frame %>% setNames(c('expr', 'pct'))
  
  # Check correlation between expr and pct
  r.squared.pct.expr = (cor(na.omit(data))[1,2])^2
  
  # Bootstrap data with gaussian noise if desired
  data.bootstrap = if(bootstrap) Bootstrap(data, noise.alpha = noise.alpha, expansion = expansion) else data
  
  tryCatch({
    
    # K-means seems to work best
    k.vector <- na_propagate(data, function(x) {
      CollapseBinarization(kmeans(x, centers = 2, nstart = nstart, iter.max = 10)$cluster-1, expansion = 1, na.omit(data))
    })
    # k.vector = na_propagate(kmeans_result, function(x) CollapseBinarization(x, expansion = 1, data))
    if(compute.r.squared) {
      r.squared.k = vegan::adonis2(na.omit(data) ~ na.omit(k.vector), method = "euclidean", by = "margin")$R2[[1]]
    } else {
      # k.res = kmeans(na.omit(data), centers = 2, nstart = nstart, iter.max = 10)
      # r.squared.k = k.res$betweenss/k.res$totss
      
      r.squared.k = MultivariateR2(na.omit(data), na.omit(k.vector))
    }
    
    if (plot) plot(data, col = convert_values(k.vector, key_table = c('0' = 'red', '1' = 'blue')))
    
    if(!fast){
    
      # 3 gaussian model that gets merged
      gmm_result3 = Mclust(data.bootstrap, G = 3, modelNames = "VVV", verbose = verbose)
      merged <- clustCombi(gmm_result3)
      binarization.bootstrap3 <- merged$classification[[2]]-1
      mv3.vector = CollapseBinarization(binarization.bootstrap3, expansion, data)
      # cluster_means3 = tapply(data.bootstrap[,1], binarization.bootstrap3, mean)
      # if(cluster_means3['1'] < cluster_means3['0']) binarization.bootstrap3 = -1*binarization.bootstrap3 + 1
    
      # 2 gaussian model 
      gmm_result = Mclust(data.bootstrap, G = 2, modelNames = "VVV", verbose = verbose)
      merged <- clustCombi(gmm_result)
      binarization.bootstrap <- merged$classification[[2]]-1
      mv.vector = CollapseBinarization(binarization.bootstrap, expansion, data)
      # cluster_means = tapply(data.bootstrap[,1], binarization.bootstrap, mean)
      # if(cluster_means['1'] < cluster_means['0']) binarization.bootstrap = -1*binarization.bootstrap + 1
      
      # Expr model
      expr.model = mclust::Mclust(data.bootstrap[,1], G = 2, verbose = verbose)
      # expr.vector = as.numeric(tapply(as.numeric(expr.model$z[, 2] > 0.5), rep(1:nrow(data), each = expansion), mean) > 0.5) %>% setNames(colnames(avg.expr))
      expr.vector = CollapseBinarization(as.numeric(expr.model$z[, 2] > 0.5), expansion, data)
      
      # Pct model
      pct.model = mclust::Mclust(data.bootstrap[,2], G = 2, verbose = verbose)
      # pct.vector = as.numeric(tapply(as.numeric(pct.model$z[, 2] > 0.5), rep(1:nrow(data), each = expansion), mean) > 0.5) %>% setNames(colnames(avg.expr))
      pct.vector = CollapseBinarization(as.numeric(pct.model$z[, 2] > 0.5), expansion, data)
      
      # Optional plot
      if (plot) plot(data.bootstrap, col = convert_values(binarization.bootstrap, key_table = c('0' = 'red', '1' = 'blue')))
      
      # Fit single Gaussian to same data
      single_gaussian <- MASS::fitdistr(data.bootstrap, "normal")
      log_likelihood_single <- as.numeric(logLik(single_gaussian))
      
      # Extract log-likelihood from Mclust
      log_likelihood_mixture <- gmm_result$loglik
      
      # Perform the Likelihood Ratio Test
      LRT_statistic <- -2 * (log_likelihood_single - log_likelihood_mixture)
      p_value <- pchisq(LRT_statistic, df = 6, lower.tail = FALSE)
      if(verbose) message('LRT p-value: ', p_value)
      
      # BIC
      bic_single = -2*(log_likelihood_single)+5*log(nrow(data.bootstrap))
      bic_mixture = -2*(log_likelihood_mixture)+11*log(nrow(data.bootstrap))
      
      # Compute r.squared (fraction of variance explained)
      r.squared.mv = vegan::adonis2(data ~ mv.vector, method = "euclidean", by = "margin")$R2[[1]]
      r.squared.mv3 = vegan::adonis2(data ~ mv3.vector, method = "euclidean", by = "margin")$R2[[1]]
      r.squared.expr = vegan::adonis2(data ~ expr.vector, method = "euclidean", by = "margin")$R2[[1]]
      r.squared.pct = vegan::adonis2(data ~ pct.vector, method = "euclidean", by = "margin")$R2[[1]]
      
      # Return binarized vectors and and metrics
      list(mv.vector = mv.vector, 
           mv3.vector = mv3.vector, 
           expr.vector = expr.vector,
           pct.vector = pct.vector,
           k.vector = k.vector,
           r.squared.mv = r.squared.mv, 
           r.squared.mv3 = r.squared.mv3,
           r.squared.expr = r.squared.expr, 
           r.squared.pct = r.squared.pct, 
           r.squared.k = r.squared.k, 
           p_value = p_value, 
           delta.bic = bic_mixture-bic_single, 
           r.squared.pct.expr = r.squared.pct.expr)
      
    } else {
      list(mv.vector = rep(0,nrow(data)), 
           mv3.vector = rep(0,nrow(data)), 
           expr.vector = rep(0,nrow(data)), 
           pct.vector = rep(0,nrow(data)), 
           k.vector = k.vector,
           r.squared.mv = 0, 
           r.squared.mv3 = 0,
           r.squared.expr = 0, 
           r.squared.pct = 0, 
           r.squared.k = r.squared.k, 
           p_value = 0, 
           delta.bic = -1, 
           r.squared.pct.expr = r.squared.pct.expr)
    }
    
  }, error = function(e) {
    
    print(paste0('Binarization failed for gene: ', gene, '!'))
    
    list(rep('error', 11))
  })
}

BinarizeExpressionMV = function(gene, avg.expr, pct.expr, noise.alpha = 0.25, plot = FALSE, verbose = FALSE, seed = 12345){
  
  # set.seed(seed)
  data = cbind(t(avg.expr[gene,]), t(pct.expr[gene,]))
  data.bootstrap = Bootstrap(data, noise.alpha = noise.alpha)
  if(plot) {
    plot(data.bootstrap, 
                main = paste0('(', bimodality_coefficient(data.bootstrap[,1]), 
                              ', ', bimodality_coefficient(data.bootstrap[,1]), ')'))
  }
  
  # Find value of noise.alpha that maximizes BC? 
  gmm_result <- tryCatch({
    # mvnormalmixEM(data.bootstrap, k = 2)
    ClusterR::GMM(data.bootstrap, gaussian_comps = 2, full_covariance_matrices = FALSE, seed = seed)
  }, error = function(e) {
    list('error', 'error', 'error') # print(e)
  })
  
  # pr = predict(gmm_result, newdata = data.bootstrap)
  posteriors <- ClusterR::predict_GMM(data.bootstrap, gmm_result$centroids, 
                                      gmm_result$covariance_matrices, 
                                      gmm_result$weights)
  
  # Use that index to threshold the posterior
  high_idx <- which.max(rowMeans(gmm_result$centroids))#which.max(sapply(gmm_result$mu, mean))
  binarization.bootstrap <- as.numeric(posteriors$cluster_proba[, high_idx] > 0.5)
  
  if(plot) plot(dat, col = convert_values(binarization.bootstrap, key_table = c('0' = 'red', '1' = 'blue')))
  vector_encoding = as.numeric(tapply(binarization.bootstrap, rep(1:nrow(data), each = 10), mean) > 0.5) %>% 
    setNames(colnames(avg.expr))
  
  # Fit single gaussian
  single_gaussian <- MASS::fitdistr(data.bootstrap, "normal")
  log_likelihood_single <- as.numeric(logLik(single_gaussian))
  # log_likelihood_mixture <- gmm_result$loglik
  # log_likelihood_mixture = sum(log(rowSums(exp(gmm_result$Log_likelihood))))
  log_likelihood_mixture = sum(gmm_loglikelihood(data.bootstrap, 
                                             means = list(gmm_result$centroids[1,], 
                                                          gmm_result$centroids[2,]), 
                                             covariances = list(gmm_result$covariance_matrices[1,,],
                                                                gmm_result$covariance_matrices[2,,]),
                                             weights = gmm_result$weights))
  
  # gmm_loglikelihood gets same loglik as mixtools 
  # log_likelihood_mixture = gmm_loglikelihood(data.bootstrap,
  #                                            means = test$mu,
  #                                            covariances = test$sigma,
  #                                            weights = test$lambda)
  
  # Perform the Likelihood Ratio Test
  LRT_statistic <- -2 * (log_likelihood_single - log_likelihood_mixture)
  p_value <- pchisq(LRT_statistic, df = 6, lower.tail = FALSE)
  if(verbose) message('LRT p-value: ', p_value)
  
  # Compute r.squared (fraction of variance explained)
  # summary(lm(as.numeric(data) ~ as.numeric(vector_encoding)))$r.squared
  r.squared = vegan::adonis2(data ~ vector_encoding, method = "euclidean", by = "margin")$R2[[1]]
  
  list(vector_encoding, p_value, r.squared)
}

BinarizeExpression3 = function(object, 
                               features, 
                               group.by = 'annotated', 
                               use.pct = TRUE,
                               method = 'gmm',
                               pct.cutoff = 0, 
                               fc.cutoff = 2, 
                               verbose = FALSE, 
                               gmm.posterior = 0.5, 
                               pseudocount = 0.05, 
                               r2 = TRUE,
                               ...){

  if(use.pct){
    
    # Percentage expression
    avg.expr = PercentageExpressed2(object, group.by = group.by, features = features)
    
    patterns.list = apply(avg.expr, 1, function(row){
      row.sorted = sort(row)
      if(verbose) print(row.sorted)
      
      if(method == 'curvature'){
        message('using method curvature')
        cut.off = FindKneePoint(row.sorted, ...)
      } else if(method == 'gmm'){
        message('using method gmm')
        vector = BimodalityTest2(row, alpha = 1, ...)[[1]]
        if(length(vector) <= 1) return(rep(NA, length(row))) else return(vector)
        cut.off = which(vector == 1)[[1]] #which(gmm$posterior[,which.max(gmm$mu)] > gmm.posterior)[1] # first point above 50% probability
      } else if(method == 'kmeans'){
        
      } else {
        stop('method must be one of curvature, gmm, or kmeans')
      }
      
      diff = diff(row.sorted)
      # diff = sapply(seq(1, length(row.sorted)-1), function(i){
      #   row.sorted[i+1] - row.sorted[i]
      # }, USE.NAMES = TRUE)
      
      # if(verbose) print(diff)
      
      if(diff[cut.off-1]>pct.cutoff){
        pattern = names(row.sorted)[(cut.off):length(row.sorted)]
        names(row) %in% pattern
      } else {
        rep(NA, length(row.sorted))
      }
    })
  } else {
    
    # Average expression (not in log space)
    avg.expr = AverageExpression(object, assay = 'RNA', group.by = group.by, features = features)$RNA
    
    patterns.list = apply(avg.expr, 1, function(row){
      row.sorted = sort(row)
      
      fc = sapply(seq(1, length(row.sorted)-1), function(i){
        (row.sorted[i+1] + pseudocount*max(row.sorted))/(row.sorted[i] + pseudocount*max(row.sorted))
      }, USE.NAMES = TRUE)
      
      if(max(fc) > fc.cutoff){
        pattern = names(row.sorted)[(which.max(fc)+1):length(row.sorted)]
        names(row) %in% pattern
      } else {
        NA
      }
    })
  }
  
  # patterns.df = na.omit(do.call(rbind, patterns.list))
  patterns.df = t(patterns.list)
  colnames(patterns.df) = colnames(avg.expr)
  
  numeric.patterns = matrix(as.numeric(patterns.df), 
                            nrow = nrow(patterns.df), 
                            ncol = ncol(patterns.df))
  rownames(numeric.patterns) = rownames(patterns.df)
  colnames(numeric.patterns) = colnames(patterns.df)
  
  no.var = apply(numeric.patterns, 1, function(x) var(x) == 0)
  message('Removing ', length(which(no.var)), ' binarized genes with zero variance!')
  numeric.patterns = numeric.patterns[which(!no.var),]
  avg.expr = avg.expr[which(!no.var),]
  
  if(r2){
    sapply(seq_len(nrow(avg.expr)), function(i) print(paste0('R^2 for gene ', rownames(avg.expr)[[i]], ': ', summary(lm(as.numeric(avg.expr[i,]) ~ numeric.patterns[i,]))$r.squared)))
  }
  
  numeric.patterns
}

capitalize_first <- function(s) {
  if (is.na(s) | s == "") return(s)  # Handle NA or empty strings
  paste0(toupper(substr(s, 1, 1)), tolower(substr(s, 2, nchar(s))))
}

PositiveCell = function(object, features, return.object = FALSE, plot.func = FeaturePlot, ...){
  gene_data = object@assays$RNA@counts[features, , drop = FALSE]
  object$positive = apply(gene_data, 2, function(x) all(x > 0))
  
  if(return.object){
    print(FeaturePlot(object, features = 'positive', ...))
    object
  } 
  else {
    FeaturePlot(object, features = 'positive', ...)
  }
}

DownloadEnsembl2 = function(path, 
                            speciesName = "Macaca_fascicularis", 
                            assembly = "Macaca_fascicularis", 
                            release = 112, 
                            version = '_6.0', 
                            genome = FALSE, 
                            primary = FALSE){
  
  setwd(path)
  
  if(primary) {
    fasta.command = paste0('ftp://ftp.ensembl.org/pub/release-', release, '/fasta/', tolower(speciesName), '/dna/', speciesName, '.', assembly, version, '.dna.primary_assembly.fa.gz')
  } else {
    fasta.command = paste0('ftp://ftp.ensembl.org/pub/release-', release, '/fasta/', tolower(speciesName), '/dna/', speciesName, '.', assembly, version, '.dna.toplevel.fa.gz')
  }
  
  gtf.command = paste0('ftp://ftp.ensembl.org/pub/release-', release, '/gtf/', tolower(speciesName), '/', speciesName, '.', assembly, version, '.', release, '.gtf.gz')
  pep.command = paste0('ftp://ftp.ensembl.org/pub/release-', release, '/fasta/', tolower(speciesName), '/pep/', speciesName, '.', assembly, version, '.pep.all.fa.gz')
  cdna.command = paste0('ftp://ftp.ensembl.org/pub/release-', release, '/fasta/', tolower(speciesName), '/cdna/', speciesName, '.', assembly, version, '.cdna.all.fa.gz')
  
  # wget ftp://ftp.ensembl.org/pub/release-110/gtf/danio_rerio/Danio_rerio.GRCz11.110.gtf.gz
  # wget ftp://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz
  # curl -O ftp://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz
  
  if(genome) system2('curl', args = paste0('-O ', fasta.command))
  system2('curl', args = paste0('-O ', gtf.command))
  system2('curl', args = paste0('-O ', pep.command))
  system2('curl', args = paste0('-O ', cdna.command))
    
}

loop = function(list, FUN){
  lapply(seq_along(list), function(x) FUN(names(list)[[x]],list[[x]]))
}

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

CellClassScoring2 = function(object, de_list, seed = 1, verbose = TRUE, compute.null = FALSE){
  
  score.names = paste0(names(de_list), '.score')
  lapply(seq_along(de_list), function(i){
    message('Found ', length(which(de_list[[i]] %in% rownames(object)))/length(de_list[[i]])*100, '% for ', names(de_list)[[i]])
  })
  object = AddModuleScore(object, de_list, seed = seed, ctrl = 1)
  colnames(object@meta.data)[startsWith(colnames(object@meta.data), 'Cluster')] = score.names
  
  if(verbose){
    print(VlnPlot(object, features = score.names, group.by = 'annotated', stack = TRUE, flip = TRUE, same.y.lims = TRUE))
    print(VlnPlot(subset(object, cell_class == 'AC'), group.by = 'annotated', features = score.names, stack = TRUE, flip = TRUE, same.y.lims = TRUE))
  }
  
  # Null distribution computation by permuting cell type labels
  if(compute.null){
    null.dist = do.call(cbind, lapply(seq_len(100), function(i) {
      object$permuted = as.character(sample(object$annotated))
      MeanMetadata(object, 'permuted', 'RGC.score')
    }))
    
    quantiles = quantile(as.vector(null.dist), probs = c(0.025, 0.975))
    
    if(verbose) print(PrettyHistogram2(as.vector(null.dist), vline = quantiles))
    
  } else {
    quantiles = c(NA, NA)
  }
  
  means = MeanMetadata(object, 'annotated', 'RGC.score')
  message('Significant RGC scores:')
  
  df = data.frame(type = names(means), rgc.score = means)
  
  if(verbose) {
    print(ggbarplot(df, 'type', 'rgc.score')+
          {if(compute.null) geom_hline(yintercept = quantiles[2], linetype = 'dashed')}+
          RotatedAxis())
  }
  
  return(object)
}

find_nearest_neighbors <- function(similarity_matrix, n, threshold = 0) {
  if (is.null(rownames(similarity_matrix)) || is.null(colnames(similarity_matrix))) {
    stop("Distance matrix must have row and column names.")
  }
  
  apply(similarity_matrix, 1, function(d) {
    sorted <- sort(d, index.return = TRUE, decreasing = TRUE)
    scores = sorted$x[2:(n + 1)] # Exclude self (first index)
    indices = sorted$ix[2:(n + 1)]
    indices.threshold = indices[scores > threshold]
    final = scores[scores > threshold]
    names(final) = rownames(similarity_matrix)[indices.threshold]
    return(final)
  })
}

AUPRC = function(dist.null, dist.observed, values.check = NULL, n = 10, plot = TRUE){
  if(is.null(values.check)){
    values.check = seq(min(c(dist.null, dist.observed), na.rm = TRUE), max(c(dist.null, dist.observed), na.rm = TRUE), length.out = n)
  }
  
  summary = do.call(rbind, lapply(values.check, function(value){
    TP = length(which(dist.observed > value))
    FP = length(which(dist.null > value))
    FN = length(which(dist.observed <= value))
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    data.frame(threshold = value, TP = TP, FP = FP, FN = FN, precision = precision, recall = recall)
  }))
  
  if(plot){
    print(ggplot(summary, aes(x = recall, y = precision)) +
            geom_line()+
            # geom_line(aes(y = precision, color = "Precision"), size = 1) +
            # geom_line(aes(y = recall, color = "Recall"), size = 1) +
            labs(title = NULL, 
                 x = "Recall",
                 y = "Precision") +
            theme_minimal())
  }
  
  return(summary)
}

Interleave = function(list, value = NULL){
  null_positions <- seq(2, length(list) * 2 - 1, by = 2)
  
  # Insert NULLs at specified positions
  for (pos in null_positions) {
    list <- append(list, list(value), after = pos - 1)
  }
  
  list
}

Heatmap2 = function(matrix, 
                    title = NULL,
                    legend.name = 'value', 
                    col = c('white', 'red'), 
                    column_title = NULL, 
                    rect_gp = gpar(col = "grey", lwd = 0),
                    border_gp = gpar(col = "black", lty = 1), 
                    label = FALSE, 
                    cluster_columns = FALSE, 
                    cluster_rows = FALSE,
                    min.value = NULL,
                    max.value = NULL,
                    parenthesis = NULL, 
                    label.fun = round,
                    label.size = 10,
                    font.family = 'ArialMT', 
                    column.font.size = 10, 
                    row.font.size = 10,
                    font.size = 10,
                    label.matrix = NULL,
                    ...){
  
  if(is.null(label.matrix)) labels = matrix else labls = label.matrix
  # if(!is.null(max.value)){
  #   matrix[matrix > max.value] = max.value
  # }
  # if(!is.null(min.value)){
  #   matrix[matrix < min.value] = min.value
  # }
  
  # if(is.null(col)){
  # }
  
  # Assign colors
  if(is.null(min.value)) min.value = min(matrix)
  if(is.null(max.value)) max.value = max(matrix)
  cmap = colorRamp2(c(min.value, max.value), col)
  
  matrix = label.fun(matrix, 2)
  if(!is.null(parenthesis)){
    # labels.matrix = paste0(matrix + '\n(', parenthesis, ')')
    labels.matrix <- matrix(paste0(matrix, "\n(", parenthesis, ")"),
                     nrow = nrow(matrix), ncol = ncol(matrix),
                     dimnames = dimnames(matrix))
    label.fun = identity
  } else {
    labels.matrix = matrix
  }
  
  Heatmap(matrix, 
          name = legend.name, 
          col = cmap, 
          column_title = title, 
          rect_gp = rect_gp,
          border_gp = border_gp, 
          row_names_gp = gpar(fontsize = row.font.size, fontfamily = font.family),
          column_names_gp = gpar(fontsize = column.font.size, fontfamily = font.family),
          column_title_gp = gpar(fontfamily = font.family, fontsize = column.font.size),
          row_title_gp = gpar(fontfamily = font.family, fontsize = row.font.size),
          cluster_columns = cluster_columns, 
          cluster_rows = cluster_rows,
          cell_fun = if(label) function(j, i, x, y, width, height, fill) { 
            if(!is.na(labels[i, j])) {
              grid.text(labels.matrix[i, j], x, y, gp = gpar(fontsize = label.size))
            }
          } else NULL, 
          heatmap_legend_param = list(
            border = "black",      # Black border around legend
            # at = seq(-2, 2, by = 1), # Define tick positions
            labels_gp = gpar(col = "black", fontfamily = font.family), # Black tick labels
            ticks_gp = gpar(col = "black"), 
            title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
          ...)
  
}

BimodalityTest2 = function(input, 
                           return.vector = TRUE, 
                           plot = FALSE,
                          seed = 42, 
                          replace.na = 0, 
                          noise = TRUE, 
                          verbose = FALSE, 
                          alpha = 0.05, # p-value threshold
                          noise.sd = NULL, # stdev of gaussian noise
                          noise.alpha = 0.05, # noise level parameter
                          smart.noise = FALSE, # higher values get slightly more noise? 
                          smart.noise.alpha = 0.05,
                          bootstrap = TRUE,
                          n.bootstraps = 500){
  
  library(MASS)
  library(mixtools)
  set.seed(seed)
  
  res = tryCatch({
    
    orig.vector = unlist(input)
    if(is.null(noise.sd)){
      stdev = sd(orig.vector)
      noise.sd = noise.alpha*stdev # 5% of the feature's stdev
    }
    
    # Bootstrap
    if(bootstrap){
      simulated.data = sample(orig.vector, n.bootstraps, replace = TRUE)
    } else {
      simulated.data = orig.vector
    }
    
    # Add noise
    if(noise) {
      if(smart.noise){
        noise.vec = mapply(rnorm, n = 1, mean = 0, sd = (noise.sd + smart.noise.alpha*simulated.data))
        # noise.vec = sapply(seq_along(simulated.data), rnorm, n = 1, mean = 0, sd = (noise.sd + simulated.data/20))
      } else {
        noise.vec = rnorm(length(simulated.data), mean = 0, sd = noise.sd)
      }
    } else {
      noise.vec = rep(0, length(simulated.data))
    }
    
    # Fit a single Gaussian distribution
    single_gaussian <- fitdistr(simulated.data + noise.vec, "normal")
    log_likelihood_single <- as.numeric(logLik(single_gaussian))
    
    # Fit a mixture of two Gaussian distributions
    output <- capture.output(
      model <- mixtools::normalmixEM(simulated.data + noise.vec, verb = verbose), 
      file = NULL)
    
    if(plot) plot(model, density=TRUE, breaks = 30)
    log_likelihood_mixture <- model$loglik
    
    # Perform the Likelihood Ratio Test
    LRT_statistic <- -2 * (log_likelihood_single - log_likelihood_mixture)
    p_value <- pchisq(LRT_statistic, df = 3, lower.tail = FALSE)
    if(verbose) message('LRT p-value: ', p_value)
    
    # if(return.pval) return(p_value)
    
    # if(p_value < alpha){
      # as.numeric(model$posterior[,which.max(model$mu)] > 0.5)
    vector_encoding = as.numeric(orig.vector > suppressWarnings(FindCutpoint(model)))
    list(vector_encoding, p_value)
    # } else {
      # 'ns'
    # }
    
  }, error = function(e) {
    list('error', 'error') # print(e)
  })
  
  res
}

BimodalityTest = function(input, return.vector = TRUE, plot = FALSE,
                           seed = 42, replace.na = 0, noise = TRUE, verbose = FALSE, 
                           alpha = 0.05, noise.sd = 1, return.pval = TRUE, smart.noise = FALSE){
  
  library(MASS)
  library(mixtools)
  set.seed(seed)
  
  vector_encoding = tryCatch({
    
    orig.vector = unlist(input)
    
    # Bootstrap
    
    # Add noise
    if(noise) {
      if(smart.noise){
        noise.vec = mapply(rnorm, length(unlist(input)), mean = 0, sd = noise.sd+unlist(input)/smart.noise.sd)
      } else {
        noise.vec = rnorm(length(unlist(input)), mean = 0, sd = noise.sd)
      }
    } else {
      noise.vec = rep(0, length(unlist(input)))
    }
    
    # Fit a single Gaussian distribution
    single_gaussian <- fitdistr(unlist(input), "normal")
    log_likelihood_single <- as.numeric(logLik(single_gaussian))
    
    # Fit a mixture of two Gaussian distributions
    output <- capture.output(
      model <- mixtools::normalmixEM(unlist(input) + noise.vec, verb = verbose), 
      file = NULL)
    
    if(plot) plot(model, density=TRUE, breaks = 30)
    log_likelihood_mixture <- model$loglik
    
    # Perform the Likelihood Ratio Test
    LRT_statistic <- -2 * (log_likelihood_single - log_likelihood_mixture)
    p_value <- pchisq(LRT_statistic, df = 3, lower.tail = FALSE)
    if(verbose) message('LRT p-value: ', p_value)
    
    if(return.pval) return(p_value)
    
    if(p_value < alpha){
      as.numeric(model$posterior[,which.max(model$mu)] > 0.5)
      # as.numeric(orig.vector > FindCutpoint(model))
    } else {
      'ns'
    }
    
  }, error = function(e) {
    'error' # print(e)
  })
  vector_encoding
}

plot_mds <- function(diss_matrix, return.plot = TRUE) {
  
  # Perform MDS
  mds <- cmdscale(as.dist(diss_matrix), eig = TRUE)
  
  # Create a data frame
  mds_df <- data.frame(MDS1 = mds$points[, 1], MDS2 = mds$points[, 2], Label = rownames(diss_matrix))
  
  # Calculate variance explained
  # var_expvlained <- sum(attr(mds, "eig")[1:2]) / sum(attr(mds, "eig")) * 100
  var1 = mds$eig[[1]]/sum(mds$eig)
  var2 = mds$eig[[2]]/sum(mds$eig)
  
  # Plot using ggplot2
  plt = ggplot(mds_df, aes(x = MDS1, y = MDS2, label = Label)) +
    geom_point() +
    geom_text(vjust = 1.5) +
    labs(title = "MDS",
         x = paste0("Dim 1 (", round(var1*100, 2), '% explained)'), 
         y = paste0("Dim 2 (", round(var2*100, 2), '% explained)'),) +
    theme_minimal()
  
  if(return.plot) plt else mds_df
}

custom_median <- function(x) {
  x <- sort(x, na.last = NA)  # Sort the values and remove NAs
  n <- length(x)
  
  if (n %% 2 == 1) {
    # If odd, return the middle element
    return(x[(n + 1) / 2])
  } else {
    # If even, return the lower of the two middle elements
    return(x[n / 2])
  }
}

LegendTopLeft = function(x=0, y = 1, just = c(0,1), line.spacing = 1, size = 10){
  theme(
    legend.position = c(x, y),  # Place the legend inside at the lower right
    legend.justification = just,  # Adjust the anchor point of the legend
    legend.background = element_rect(fill = "transparent", colour = NA), 
    legend.key.height = unit(line.spacing, "lines"), 
    legend.text = element_text(size = size)
  )
}

LegendTopRight = function(x=1, y = 1, just = c(1,1), line.spacing = 1, size = 10){
  theme(
    legend.position = c(x, y),  # Place the legend inside at the lower right
    legend.justification = just,  # Adjust the anchor point of the legend
    legend.background = element_rect(fill = "transparent", colour = NA), 
    legend.key.height = unit(line.spacing, "lines"), 
    legend.text = element_text(size = size)
  )
}

LegendLowerLeft = function(x=0.01, y = 0.01, just = c(0,0), line.spacing = 1, size = 10, background.col = 'transparent', background.alpha = 0){
  theme(
    legend.position = c(x, y),  # Place the legend inside at the lower right
    legend.justification = just,  # Adjust the anchor point of the legend
    legend.background = element_rect(fill = alpha(background.col, background.alpha), colour = NA), 
    legend.key.height = unit(line.spacing, "lines"), 
    legend.text = element_text(size = size)
  )
}

LegendLowerRight = function(x=1, y = 0, just = c(1,0), line.spacing = 1, size = 10){
  theme(
    legend.position = c(x, y),  # Place the legend inside at the lower right
    legend.justification = just,  # Adjust the anchor point of the legend
    legend.background = element_rect(fill = "transparent", colour = NA), 
    legend.key.height = unit(line.spacing, "lines"), 
    legend.text = element_text(size = size)
    
    # legend.spacing.y = unit(5, "pt")         # reduces vertical gap
    # legend.text = element_text(margin = margin(r = -5, unit = "pt")),
    # legend.spacing.x = unit(0, "pt")          # reduces horizontal gap (for horizontal legends)
  )  
}

#' Install a Package from CRAN or Bioconductor
#'
#' This function attempts to install a given package from CRAN. If the installation
#' from CRAN fails, it will attempt to install the package from Bioconductor using `BiocManager`.
#' If the package is already installed, the function will skip the installation.
#'
#' @param pkg A character string specifying the name of the package to install.
#'
#' @return The function does not return a value. It attempts to install the specified package.
#' It will print messages to indicate the success or failure of the installation attempts.
#'
#' @details
#' The function first checks if the package is already installed using `requireNamespace()`. If the package
#' is not installed, it first tries to install it from CRAN using `install.packages()`. If the installation
#' from CRAN fails (e.g., the package is only available on Bioconductor), the function will check if
#' `BiocManager` is installed. If not, it installs `BiocManager` and then attempts to install the package
#' from Bioconductor.
#'
#' @examples
#' \dontrun{
#' # Install ggplot2 from CRAN or Bioconductor if it's not already installed:
#' install_package("ggplot2")
#'
#' # Install a Bioconductor package (e.g., "GenomicFeatures"):
#' install_package("GenomicFeatures")
#' }
#'
#' @seealso \code{\link{install.packages}}, \code{\link{BiocManager::install}}
#'
#' @export
install_package <- function(pkg) {
  # Check if the package is already installed
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Try to install from CRAN
    tryCatch({
      message(paste("Attempting to install", pkg, "from CRAN..."))
      install.packages(pkg)
    }, error = function(e) {
      message(paste("Failed to install", pkg, "from CRAN. Trying BiocManager..."))
      
      # Check if BiocManager is installed
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      
      # Try to install from Bioconductor
      tryCatch({
        BiocManager::install(pkg)
      }, error = function(e) {
        message(paste("Failed to install", pkg, "from Bioconductor."))
      })
    })
  } else {
    message(paste("Package", pkg, "is already installed."))
  }
}


DrawPhylogeneticTreeFromNewick = function(newick, key, as_ylab = TRUE, layout = 'rectangular', species.prune = NULL, rootedge = 50, interval = 50){
  
  # Convert to common name
  tree = ape::read.tree(newick)
  tree$tip.label = key$CommonName[match(gsub("_", " ", tree$tip.label), key$LatinName)]
  # dendro = phylogram::as.dendrogram.phylo(tree)
  # sorted.dendro = reorder(dendro, 1:20)
  
  # Prune species if necessary
  # pruned.dendro = dendextend::prune(sorted.dendro, species.prune)
  
  # phylo_obj = ape::as.phylo(as.hclust(pruned.dendro))
  phylo_obj = tree
  breaks = seq(0, max(phylo_obj$edge.length)+interval, by = interval)
  max_time = max(phylo_obj$edge.length)
  shift = max_time-breaks[length(breaks)-1]
  
  ggtree(phylo_obj, layout = layout) + 
    theme_tree2() + 
    geom_tiplab(as_ylab=as_ylab) + 
    geom_rootedge(rootedge = rootedge) +
    scale_x_continuous(name = "Time (millions of years ago)",
                       breaks = seq(shift-interval, max_time, by = 50),  # Breaks at 50 MYA intervals
                       labels = rev(breaks))+  # Adjust labels so 0 MYA is present
    theme(axis.text = element_text(color = 'black'))
}

rotate_coords <- function(coords, angle, radians = TRUE) {
  # Convert the angle from degrees to radians
  if(!radians) angle_radians <- angle * pi / 180 else angle_radians = angle
  
  # Create the rotation matrix
  rotation_matrix <- matrix(c(cos(angle_radians), -sin(angle_radians),
                              sin(angle_radians), cos(angle_radians)), 
                            nrow = 2, byrow = TRUE)
  
  # Original coordinates
  # coords <- cbind(x, y)
  
  # Rotate the coordinates
  rotated_coords <- t(rotation_matrix %*% t(coords))
  
  # Return the rotated coordinates
  return(data.frame(rotated_coords[, 1], rotated_coords[, 2]) %>% setNames(colnames(coords)))
}

UMIplot <- function(counts, cutoff = 200, title = NULL){
  count_data=data.frame(rank=1:length(counts), 
                        count=rev(sort(counts)))
  cutoff.rank=nrow(count_data[count_data$count>=cutoff,])
  # plot=plot(log10(count+1) ~ rank, data=count_data, log="x", type = "S", main=title)
  # abline(v=cutoff.rank, lty=2)
  ggplot(count_data, aes(x = rank, y = log10(count + 1))) +
    geom_step() +  # equivalent to type = "S"
    scale_x_log10() +  # equivalent to log = "x"
    geom_vline(xintercept = cutoff.rank, linetype = 'dashed') +
    labs(title = title, x = "Rank", y = "log10(Count + 1)") +
    theme_bw()
}

MakeMutuallyExclusive = function(list){
  all = unlist(list)
  duplicates = getDuplicates(all)
  lapply(list, function(x) x[!x %in% duplicates])
}

setRowNames <- function(object, names) {
  rownames(object) <- names
  return(object)
}

# Function to compute softmax
softmax <- function(x) {
  exp_x <- exp(x)        # Exponentiate each element
  return(exp_x / sum(exp_x))  # Normalize by the sum of exponentiated elements
}

GetTrainingSet = function(object, group.by, proportion = 0.6, max.size = 100, seed = 12345){
  
  Idents(object) = group.by
  clusters = factor(object@meta.data[[group.by]])
  training_cells = unlist(sapply(seq_along(levels(clusters)), function(index) 
    WhichCells(object, idents = levels(clusters)[[index]], 
               downsample = floor(table(clusters)[[index]])*proportion)))
  object_train = object[,training_cells]
  object_train = DownsampleSeurat(object_train, group.by = group.by, size = max.size, seed = seed)
  return(object_train)
}

PieChart = function(df, label = FALSE, label.func = geom_text_repel){
  ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
    geom_col(color = "black") +
    {if(label) label.func(aes(label = Var1), position = position_stack(vjust = 0.5))} +
    coord_polar(theta = "y") +
    # scale_fill_brewer() +
    labs(fill = "Cluster") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.background = element_rect(fill = "white")) 
}

FindGenes = function(object, prefix){
  rownames(object)[startsWith(rownames(object), prefix)]
}

DimPlotLabeled = function(object, group.by = 'seurat_clusters', ...){
  DimPlot(object, group.by = group.by, label = TRUE, ...) + NoLegend()
}

convert_values2 = function(data_vector, key_table, keep.missing = FALSE) {
  
  # Ensure the key_table has the correct columns
  if(inherits(key_table, 'data.frame')){
    if (!all(c("old.names", "new.names") %in% names(key_table))) {
      stop("key_table must have 'old.names' and 'new.names' columns")
    }
    
    # Make a full key including stuff not in key table
    full_key = rbind(data.frame(old.names = setdiff(data_vector, key_table$old.names), 
                                new.names = setdiff(data_vector, key_table$old.names)), # anything that is not there, will get renamed back to itself
                     key_table)
    
    # Perform the conversion
    if(keep.missing) {
      converted_vector <- full_key$new.names[match(data_vector, full_key$old.names)]
    } else {
      converted_vector <- key_table$new.names[match(data_vector, key_table$old.names)]
    }
    
  } else if(is.vector(key_table) && !is.null(names(key_table))){
    
    converted_vector = key_table[as.character(data_vector)]
    
  } else {
    stop('key_table must be named vector or data frame')
  }
  
  # Return the converted vector
  converted_vector
}

convert_values = function(data_vector, key_table, keep.missing = FALSE) {
  
  # Ensure the key_table has the correct columns
  if(inherits(key_table, 'data.frame')){
    if (!all(c("old.names", "new.names") %in% names(key_table))) {
      stop("key_table must have 'old.names' and 'new.names' columns")
    }
    
    # Make a full key including stuff not in key table
    full_key = rbind(data.frame(old.names = setdiff(data_vector, key_table$old.names), 
                                new.names = setdiff(data_vector, key_table$old.names)), 
                     key_table)
    
    # Perform the conversion
    if(keep.missing) {
      converted_vector <- full_key$new.names[match(data_vector, full_key$old.names)]
    } else {
      converted_vector <- key_table$new.names[match(data_vector, key_table$old.names)]
    }
    
  } else if(is.vector(key_table) && !is.null(names(key_table))){
    
    converted_vector = key_table[as.character(data_vector)]
    
  } else {
    stop('key_table must be named vector or data frame')
  }
  
  # Return the converted vector
  converted_vector
}

# Returns clusters that have at least X cells across each (or one) condition
MoreThanXCells = function(object, annotate.by = 'seurat_clusters', group.by = NULL, more.than = 30, allow.zero = FALSE){
  tabulation = table(object@meta.data[[annotate.by]], object@meta.data[[group.by]])
  print(tabulation)
  
  if(allow.zero) {
    names(which(apply(as.data.frame.matrix(tabulation), 1, function(x) all(x > more.than | x == 0))))
  } else {
    names(which(apply(as.data.frame.matrix(tabulation), 1, function(x) all(x > more.than))))
  }
}

exchange_factor_level <- function(factor_var, old_level, new_level) {
  # Check if the input is a factor
  if (!is.factor(factor_var)) {
    stop("The input variable must be a factor.")
  }
  
  # Check if the old level exists in the factor
  if (!(old_level %in% levels(factor_var))) {
    stop("The old level does not exist in the factor.")
  }
  
  # Check if the new level already exists in the factor
  if (new_level %in% levels(factor_var)) {
    stop("The new level already exists in the factor.")
  }
  
  # Replace the old level with the new level
  levels(factor_var)[levels(factor_var) == old_level] <- new_level
  
  return(factor_var)
}

geneDotPlotStack = function(object, genes){
  StackedPlots(lapply(genes, geneDotPlotFast))
}

PlainYAxis = function(){
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())
}

PlainXAxis = function(){
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())
}

get_legend <- function(my_plot) {
  plot_grob <- ggplotGrob(my_plot)
  legend <- plot_grob$grobs[[which(sapply(plot_grob$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

add_nulls_between <- function(lst) {
  # Use lapply to insert NULLs
  interleaved <- lapply(lst, function(x) list(x, NULL))
  # Flatten the list using do.call and c
  flattened <- do.call(c, interleaved)
  # Remove the last NULL added by the function
  return(flattened[-length(flattened)])
}

smartReadRDS = function(filepath){
  objectname = basename(gsub('.rds', '', filepath))
  if(!exists(objectname)){
    readRDS(filepath)
  } else {
    message(objectname, ' already loaded')
  }
}

VlnPlot2 = function(object, features = NULL, split.by = NULL, group.by = NULL, 
                    fill.by = 'feature', cols = NULL, 
                    idents = NULL, stack = FALSE, pt.size = 0, combine = TRUE, ...){
  meta = Metadata(object, split.by, group.by)
  VlnPlot(object, features = features, split.by = split.by, group.by = group.by, stack = stack, 
          idents = idents, fill.by = fill.by, cols = cols[match(meta[[group.by]], names(cols))], 
          pt.size = pt.size, combine = combine, ...)
}

RunDESeq2 = function(object, ident.1, ident.2){
  
  colData = object@meta.data
  colData$condition = factor(colData$condition, levels = c(ident.2, ident.1))
  dds <- DESeqDataSetFromMatrix(countData = object@assays$RNA@counts,
                                colData = colData,
                                design= ~ condition)
  dds <- DESeq(dds)
  result_name = resultsNames(dds)[2]
  res <- results(dds, name=result_name)
  result = data.frame(gene = rownames(res), 
                      cluster = str_split_fixed(result_name, '_', 4)[,2], 
                      avgExpr = res$baseMean, 
                      avg_log2FC = res$log2FoldChange, 
                      stat = res$stat,
                      p_val = res$pvalue,
                      p_val_adj = res$padj) %>% arrange(p_val_adj, avg_log2FC)
  
  # sorted.res = res %>% 
  # as.data.frame %>% 
  # setNames(c('gene', 'cluster', 'avgExpr', 'avg_log2FC', 'stat', 'auc', 'p_val', 'p_val_adj', 'pct_in', 'pct_out')) %>% 
  # arrange(padj, log2FoldChange)
  
  result
}

DendroOrder2 = function(object, group.by = 'seurat_clusters', ...){
  dendro = CorrelationHeatmap(object, group.by = group.by, return.dendrogram = TRUE, ...)
  object$dendro.order = factor(object@meta.data[[group.by]], levels = labels(dendro))
  object@tools$BuildClusterTree = dendro # edit 3/19/25: used to be ape::as.phylo(dendro)
  # PlotClusterTree(object)
  return(object)
}

ConvertGeneSymbolsFromMatrix = function(matrix, list){
  subset_list = list[list %in% colnames(matrix)]
  subset_matrix = matrix[,match(subset_list, colnames(matrix))]
  # column_vector = as.numeric(colnames(subset_matrix) %in% list)
  new_genes = unlist(apply(subset_matrix, 2, function(x) names(which(x == 1))))
  return(new_genes)
}

CelltypeProportionBarplotTwoGroup = function(object, celltype = 'lit_type', sample = 'sample_core', group = 'group', cols = NULL){
  
  library(rstatix)
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  tabulation = table(data[[sample]], data[[celltype]])
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  # Add metadata
  metadata = Metadata(object, sample, group)
  prop_data = melted %>% mutate(group = metadata[match(melted$Sample, metadata[[sample]]),group])

  stat.test <- prop_data %>%
    group_by(Type) %>%
    t_test(Proportion ~ group) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  
  bp = ggbarplot(prop_data, x = "Type", y = "Proportion", 
                 position = position_dodge(),
                 # fill = "tissue", 
                 color = group,
                 add = c("mean_sd", "jitter")) + 
    # NoLegend() + 
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    # {if(!is.null(cols)) scale_color_manual(values = cols)}+
    RotatedAxis()
  
  # Add p-values onto the bar plots
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = 'Type', dodge = 0.8)
  
  print(bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01
  ))
  
  return(stat.test)
}

CelltypeProportionBarplotGenotype = function(object, celltype = 'lit_type', sample = 'sample_core', group = 'genotype'){
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  tabulation = table(data[[sample]], data[[celltype]])
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  # Add metadata
  metadata = Metadata(seurat, 'sample_core', 'tissue', 'treatment', 'genotype')
  prop_data = cbind(melted, metadata[match(melted$Sample, metadata$sample_core),-1])
  
  # x = eval(parse(text = group))
  # fm <- expr(Proportion ~ tissue)#!!sym(group))
  stat.test <- prop_data %>%
    group_by(Type) %>%
    t_test(Proportion ~ genotype) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  
  bp = ggbarplot(prop_data, x = "Type", y = "Proportion", 
                 position = position_dodge(),
                 # fill = "tissue", 
                 color = group,
                 add = c("mean_sd", "jitter")) + 
    # NoLegend() + 
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    scale_color_manual(values = conditions_palette6)+
    RotatedAxis()
  
  # Add p-values onto the bar plots
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = 'Type', dodge = 0.8) 
  
  print(bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01
  ))
  
  return(stat.test)
}

CelltypeProportionBarplotTreatment = function(object, celltype = 'lit_type', sample = 'sample_core', group = 'treatment'){
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  tabulation = table(data[[sample]], data[[celltype]])
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  # Add metadata
  metadata = Metadata(seurat, 'sample_core', 'tissue', 'treatment', 'genotype')
  prop_data = cbind(melted, metadata[match(melted$Sample, metadata$sample_core),-1])
  
  # x = eval(parse(text = group))
  # fm <- expr(Proportion ~ tissue)#!!sym(group))
  stat.test <- prop_data %>%
    group_by(Type) %>%
    t_test(Proportion ~ treatment) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  
  bp = ggbarplot(prop_data, x = "Type", y = "Proportion", 
                 position = position_dodge(),
                 # fill = "tissue", 
                 color = group,
                 add = c("mean_sd", "jitter")) + 
    # NoLegend() + 
    scale_y_continuous(expand = expansion(mult = c(0, .2))) + 
    scale_color_manual(values = conditions_palette6)+
    RotatedAxis()
  
  # Add p-values onto the bar plots
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = 'Type', dodge = 0.8) 
  
  print(bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01
  ))
  
  return(stat.test)
}

CelltypeProportionBarplotTissue = function(object, celltype = 'lit_type', sample = 'sample_core', group = 'tissue'){
  
  library(rstatix)
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  tabulation = table(data[[sample]], data[[celltype]])
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  # Add metadata
  metadata = Metadata(object, sample, 'tissue')
  prop_data = melted %>% mutate(group = metadata[match(melted$Sample, metadata[[sample]]),group])
  # prop_data = cbind(melted, metadata[match(melted$Sample, metadata[[sample]]),-1]) %>% setNames(c('Sample', 'Type', 'Proportion'))
  
  # x = eval(parse(text = group))
  # fm <- expr(Proportion ~ tissue)#!!sym(group))
  stat.test <- prop_data %>%
    group_by(Type) %>%
    t_test(Proportion ~ group) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  
  bp = ggbarplot(prop_data, x = "Type", y = "Proportion", 
                 position = position_dodge(),
                 # fill = "tissue", 
                 color = group,
                 add = c("mean_sd", "jitter")) + 
    # NoLegend() + 
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    scale_color_manual(values = conditions_palette6)+
    RotatedAxis()
  
  # Add p-values onto the bar plots
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_sd", x = 'Type', dodge = 0.8) 
  
  print(bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01
  ))
  
  return(stat.test)
}

getDuplicates = function(vector){
  unique(vector[duplicated(vector)])
}

ShuffleObject = function(object, group.by = 'annotated'){
  
  levels = levels(object@meta.data[[group.by]])
  
  # Find index for each annotation
  idx = lapply(levels, function(this.level){
    which(object@meta.data[[group.by]] == this.level)
  })
  
  shuffled = lapply(idx, sample)
  
  return(Cells(object)[unlist(shuffled)])
  # shuffled.object = object[, Cells(object)[unlist(shuffled)] ]
  # shuffled.object
  
}

SeuratHeatmap = function(object, features, 
                         group.by = 'seurat_clusters', 
                         color = c("white", "white", 'red'), 
                         max.z.score = 2, 
                         assay = 'RNA', 
                         rotate = FALSE, 
                         annotate.by = NULL, 
                         annotation_cols = NULL, 
                         cluster_rows = FALSE, 
                         cluster_columns = FALSE,
                         scaled.expression = TRUE,
                         label = FALSE,
                         legend.name = 'scaled\nexpr',
                         border = TRUE, 
                         ha = NULL,
                         show.anno.legend = TRUE,
                         shuffle = FALSE,
                         add.rectangles = FALSE,
                         row.breaks = NULL, 
                         col.breaks = NULL, 
                         lwd = 0.6, 
                         lty = 1, 
                         ...){
  
  # Using scale.data slot because this is at the cell level, not averaged
  if(scaled.expression) {
    # Check that scaled data exists
    stopifnot(nrow(object@assays[[assay]]@scale.data) > 0)
    scaled.expr = as.data.frame(object@assays[[assay]]@scale.data)
    color.map = colorRamp2(c(-max.z.score, 0, max.z.score), color)
  } else {
    scaled.expr = as.data.frame(object@assays[[assay]]@data)
    color.map = colorRamp2(c(0, max.z.score), color)
  }
  colData = object@meta.data
  genes = rownames(scaled.expr)
  colData = colData[order(colData[[group.by]]),]
  scaled.expr = scaled.expr[match(features, genes),rownames(colData)]
  
  if(!is.null(annotate.by)){
    
    if(length(annotate.by) > 1){
      
      # First annotation
      metadata = Metadata(object, feature.1 = group.by, feature.2 = annotate.by[[1]])
      annotation1 = metadata[match(colData[[group.by]], metadata[[group.by]]), annotate.by[[1]] ]
      
      # Second annotation
      metadata = Metadata(object, feature.1 = group.by, feature.2 = annotate.by[[2]])
      annotation2 = metadata[match(colData[[group.by]], metadata[[group.by]]), annotate.by[[2]] ]
      
      if(is.null(annotation_cols)){
        ha = HeatmapAnnotation(anno1 = annotation1, anno2 = annotation2, 
                               annotation_name_side = "left", show_legend = show.anno.legend)
      } else{
        ha = HeatmapAnnotation(anno1 = annotation1, anno2 = annotation2, annotation_name_side = "left", 
                               col = list(anno1 = annotation_cols[[1]], 
                                          anno2 = annotation_cols[[2]]), 
                               show_legend = show.anno.legend)
      }
    } else {
      
      # Only one annotation
      metadata = Metadata(object, feature.1 = group.by, feature.2 = annotate.by)
      annotation = metadata[match(colData[[group.by]], metadata[[group.by]]), annotate.by]
      if(is.null(annotation_cols)){
        ha = HeatmapAnnotation(anno = annotation, annotation_name_side = "left", 
                               show_legend = show.anno.legend)
      } else{
        ha = HeatmapAnnotation(anno = annotation, annotation_name_side = "left", col = list(anno = annotation_cols), 
                               show_legend = show.anno.legend)
      }
    }
    
  }
  
  if(shuffle){
    idx = ShuffleObject(object, group.by = group.by)
    scaled.expr = scaled.expr[,idx]
  }
  
  if(rotate) scaled.expr = t(scaled.expr)
  
  draw(
    Heatmap(as.matrix(scaled.expr),
          name = "mat",
          col = color.map,
          cluster_rows = cluster_rows,
          cluster_columns = cluster_columns, 
          cell_fun = if(label) function(j, i, x, y, width, height, fill) { 
            if(!is.na(scaled.expr[i, j])) {
              grid.text(round(scaled.expr[i, j], 2), x, y, gp = gpar(fontsize = 10))
            }},
          border = border, 
          top_annotation = ha, 
          heatmap_legend_param = list(
            title = legend.name,
            border = "black",      # Black border around legend
            # at = seq(-2, 2, by = 1), # Define tick positions
            labels_gp = gpar(col = "black"), # Black tick labels
            ticks_gp = gpar(col = "black")), 
          ...), 
    padding = unit(c(0,0,0,0), "mm")
  )
  
  if(add.rectangles){
    
    if(length(unique(row.breaks)) > length(unique(col.breaks))) labels = unique(row.breaks) else labels = unique(col.breaks)
    widths.x = table(factor(col.breaks, levels = labels))
    widths.y = table(factor(row.breaks, levels = labels))
    x.breaks = c(0,cumsum(widths.x)) # start at zero
    y.breaks = c(0,cumsum(widths.y)) # start at zero
    # y.breaks = c(sapply(unique(row.breaks), function(x) which(row.breaks == x)[1])-1, length(row.breaks))
    # x.breaks = c(sapply(unique(col.breaks), function(x) which(col.breaks == x)[1])-1, length(col.breaks))
    
    for(i in seq(1, length(x.breaks)-1)) {
      y1 = y.breaks[[i]]/length(row.breaks)
      y2 = y.breaks[[i+1]]/length(row.breaks)
      x1 = x.breaks[[i]]/length(col.breaks)
      x2 = x.breaks[[i+1]]/length(col.breaks)
      
      # Make a box 
      # lwd = 0.5
      
      # left
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(x1, x1), c(1-y1, 1-y2), gp = gpar(lty = lty, lwd = lwd))
      })
      
      # top
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(x1, x2), c(1-y1, 1-y1), gp = gpar(lty = lty, lwd = lwd))
      })
      
      # right
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(x2, x2), c(1-y1, 1-y2), gp = gpar(lty = lty, lwd = lwd))
      })
      
      # bottom
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(x1, x2), c(1-y2, 1-y2), gp = gpar(lty = lty, lwd = lwd))
      })
      
    }
  }
  
}

volcanoPlot <- function(table, fdr_cutoff = 0.01, fc_cutoff = 0.25, max_fdr = 1e-100, max_fc = 2, 
                        labels = FALSE, name1 = "", name2 = "", anno.y = 0.1, color.signif = 'red', ...){
  # max_fc = max(abs(table$avg_log2FC))
  table$p_val_adj[table$p_val_adj < max_fdr] = max_fdr
  table$avg_log2FC[table$avg_log2FC < -max_fc] = -max_fc
  table$avg_log2FC[table$avg_log2FC > max_fc] = max_fc
  data_signif = subset(table, p_val_adj < fdr_cutoff & abs(avg_log2FC) > fc_cutoff)
  data_ns = subset(table, p_val_adj >= fdr_cutoff | abs(avg_log2FC) <= fc_cutoff)
  plot=ggplot(table, aes(y=(-log10(p_val_adj)), x=avg_log2FC))+
    geom_point()+
    ylab("-log10 FDR")+
    xlab("log2 fold change")+
    geom_point(data = data_signif, color = color.signif)+
    geom_point(data = data_ns, color = "black")+
    {if(labels) geom_text_repel(data = data_signif, label = paste0("italic('", data_signif$gene, "')"), parse = TRUE,...)}+
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 'dashed')+
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 'dashed')+
    annotation_custom(
      grob = textGrob(paste0(nrow(subset(data_signif, avg_log2FC < 0)), '\n', name1), x = 0.1, y = anno.y, hjust = 0.5, vjust = 1, gp = gpar(col = "grey")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    annotation_custom(
      grob = textGrob(paste0(nrow(subset(data_signif, avg_log2FC > 0)), '\n', name2), x = 0.9, y = anno.y, hjust = 0.5, vjust = 1, gp = gpar(col = "grey")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    # annotate('text', x = 0.2, y = 0.8, label = nrow(subset(data_signif, avg_log2FC < 0)))+
    scale_x_continuous(limits = c(-max_fc, max_fc))+
    theme_cowplot() +
    ArialFont()
  return(plot)
}

FindMarkersFast = function(object, group.by = 'seurat_clusters', ident.1 = NULL, ident.2 = NULL, p_val_adj_cutoff = 0.05, avg_log2FC_cutoff = 0.25){
  output = wilcoxauc(object, group_by = group.by, groups_use = c(ident.1, ident.2)) %>% setNames(c('gene', 'cluster', 'avgExpr', 'avg_log2FC', 'stat', 'auc', 'p_val', 'p_val_adj', 'pct_in', 'pct_out'))
  if(is.factor(object@meta.data[[group.by]])) output$cluster = factor(output$cluster, levels = levels(object@meta.data[[group.by]]))
  output %>% filter(p_val_adj <= p_val_adj_cutoff & abs(avg_log2FC) > avg_log2FC_cutoff & cluster == ident.1) %>% arrange(cluster, p_val_adj, desc(abs(avg_log2FC)))
}

FindAllMarkersFast = function(object, group.by = 'seurat_clusters', p_val_adj_cutoff = 0.05, avg_log2FC_cutoff = 0.25, assay = 'RNA'){
  output = wilcoxauc(object, group_by = group.by, seurat_assay = assay) %>% setNames(c('gene', 'cluster', 'avgExpr', 'avg_log2FC', 'stat', 'auc', 'p_val', 'p_val_adj', 'pct_in', 'pct_out'))
  if(is.factor(object@meta.data[[group.by]])) output$cluster = factor(output$cluster, levels = levels(object@meta.data[[group.by]]))
  output %>% filter(p_val_adj <= p_val_adj_cutoff & abs(avg_log2FC) > avg_log2FC_cutoff) %>% arrange(cluster, p_val_adj, -avg_log2FC)
}

JSHeatmap2 = function(list1, list2, xlab = NULL, ylab = NULL, ari = FALSE, union.threshold = 0, ...){
  args = list(...)
  JSHeatmap(JSMatrix(table((list1), as.character(list2)), union.threshold = union.threshold), ...)+
    {if(is.null(args$title) & ari) ggtitle(paste0('ARI = ', round(adj.rand.index(list1, list2), 2)))}+ 
    {if(!is.null(args$title) & ari) ggtitle(paste0(args$title, ' (ARI = ', round(adj.rand.index(list1, list2), 2), ')'))}+ 
    # {if(is.null(args$title) & !ari) }+  # do nothing
    {if(!is.null(args$title) & !ari) ggtitle(paste0(args$title))}+ 
    xlab(xlab)+
    ylab(ylab)+
    ArialFont()+
    theme(axis.title = element_text(color = 'black'), 
          axis.text = element_text(color = 'black'))
}

ConfusionMatrix = function(true, predicted, row.norm = TRUE, plot = TRUE, max.value = 1, 
                           col.high = "#584B9FFF", col.low = 'white', 
                           legend_name = 'Percentage\nmapping', fill.na = NULL){
  raw_table = table(true, predicted)
  if(row.norm){
    matrix = raw_table/rowSums(raw_table)
  } else {
    matrix = t(t(raw_table)/colSums(raw_table))
  }
  
  melted = reshape2::melt(matrix) %>% setNames(c('True', 'Predicted', 'Percentage'))
  plt = ggplot(melted, aes(x = Predicted, y = True))+
    geom_point(aes(colour = Percentage,  size=Percentage))+
    scale_color_gradient(legend_name, low=col.low, high = col.high, limits=c(0, max.value), na.value = 'grey') +
    scale_radius(legend_name, limits=c(0, max.value)) +
    # scale_size(range = c(1, max.size), limits = c(0, max.perc))+
    theme_bw() +
    RotatedAxis() +
    ylab('True')+
    xlab('Predicted')+
    theme(axis.title = element_text(color = 'black'), 
          axis.text = element_text(color = 'black'),
          plot.title = element_text(hjust = 0.5))
  
  if(!plot){
    matrix = as.data.frame.matrix(matrix)
    if(!is.null(fill.na)) matrix[is.na(matrix)] = fill.na
    return(matrix)
  }
    
  return(plt)
}

AddSpacing = function(plt.list){
  
}

ConvertGeneSymbols2 = function(object, orthology_graph){
  
  # Convert matrix into a key
  key = as.data.frame(ConvertGeneSymbolsFromMatrix(orthology_graph, rownames(object))) %>% 
    rownames_to_column() %>% 
    setNames(c('old.names', 'new.names'))
  
  # Save original feature names
  original.features = rownames(object)
  object@misc$orig.features = original.features
  
  # Rename features
  converted.features = convert_values(rownames(object), key)
  converted.features[is.na(converted.features)] = original.features[is.na(converted.features)]
  object = RenameFeatures(object, new.names = converted.features)
  
  return(object)
}

PlotBinaryMask = function(full_mask, mask){
  grid.grabExpr(draw(Heatmap(rbind(full_mask, mask), 
                             cluster_rows = F, 
                             cluster_columns = F, 
                             col = c("grey", "blue"), 
                             column_names_side = "bottom", 
                             row_names_side = "left", 
                             rect_gp = gpar(col = "black", lwd = 1),
                             name = "ProbExp                  "
  )))
}

GrabSamples = function(dir, prefix, column){
  matrix.list = lapply(list.files(path = dir, pattern = prefix, full.names = TRUE), function(x) as.matrix(read.csv(x, row.names = 1)))
  return(do.call(cbind, lapply(matrix.list, function(matrix) matrix[,column])))
}

StdTable2 = function(dir, prefix){
  matrix.list = lapply(list.files(path = dir, pattern = prefix, full.names = TRUE), function(x) as.matrix(read.csv(x, row.names = 1)))
  sd = apply(simplify2array(matrix.list), 1:2, sd)
  return(sd)
}

StdTable = function(filename){
  matrix.list = lapply(paste0(filename, "_", 0:9, ".csv"), function(x) as.matrix(read.csv(x, row.names = 1)))
  sd = apply(simplify2array(matrix.list), 1:2, sd)
  return(sd)
}

get_quartiles <- function(x, na.rm = TRUE) {
  q1 <- quantile(x, probs = 0.25, na.rm = na.rm)
  q3 <- quantile(x, probs = 0.75, na.rm = na.rm)
  return(c(Q1 = q1, Q3 = q3))
}

MedianTable2 = function(dir, prefix){
  matrix.list = lapply(list.files(path = dir, pattern = prefix, full.names = TRUE), function(x) as.matrix(read.csv(x, row.names = 1)))
  median = apply(simplify2array(matrix.list), 1:2, median)
  return(median)
}

IQRTable2 = function(dir, prefix){
  matrix.list = lapply(list.files(path = dir, pattern = prefix, full.names = TRUE), function(x) as.matrix(read.csv(x, row.names = 1)))
  median = apply(simplify2array(matrix.list), 1:2, get_quartiles)
  return(median)
}

IQRTable = function(dir, prefix){
  matrix.list = lapply(paste0(filename, "_", 0:9, ".csv"), function(x) as.matrix(read.csv(x, row.names = 1)))
  median = apply(simplify2array(matrix.list), 1:2, get_quartiles)
  return(median)
}

MedianTable = function(filename){
  matrix.list = lapply(paste0(filename, "_", 0:9, ".csv"), function(x) as.matrix(read.csv(x, row.names = 1)))
  median = apply(simplify2array(matrix.list), 1:2, median)
  return(median)
}

MeanTable2 = function(dir, prefix){
  matrix.list = lapply(list.files(path = dir, pattern = prefix, full.names = TRUE), function(x) as.matrix(read.csv(x, row.names = 1)))
  message('Averaging ', length(matrix.list), ' matrices!')
  mean = apply(simplify2array(matrix.list), 1:2, mean)
  return(mean)
}

MeanTable = function(filename){
  matrix.list = lapply(paste0(filename, "_", 0:9, ".csv"), function(x) as.matrix(read.csv(x, row.names = 1)))
  mean = apply(simplify2array(matrix.list), 1:2, mean)
  return(mean)
}

DEGTree = function(object, group.by = "annotated", title = NULL, mc.cores = 1, ...){
  types.use = unique(object@meta.data[[group.by]])
  ndegs = Iterate(types.use, function(i,j) nDEGs(wilcoxauc(object, group_by = group.by, groups_use = c(i,j)), presto = TRUE, ...))
  
  dist = matrix(ndegs, nrow = length(types.use), ncol = length(types.use), dimnames = list(types.use, types.use))
  # dist[lower.tri(dist)] = unlist(ndegs)
  # dist[upper.tri(dist)] = t(dist)[upper.tri(t(dist))]
  # dist = ndegs
  dist[is.na(dist)] = 0
  
  print(Heatmap(dist, 
          col = colorRamp2(c(0, max(dist, na.rm = TRUE)), c("white", "red")),
          cluster_rows = hclust(as.dist(dist)), 
          cluster_columns = hclust(as.dist(dist)), 
          clustering_method_columns = "manhattan",
          clustering_method_rows = "manhattan",
          show_row_names = TRUE, 
          show_column_names = TRUE, 
          row_dend_side = "left", 
          show_column_dend = TRUE, 
          column_title = title,
          border_gp = gpar(col = "black", lty = 1),
          # column_title = title,
          cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.0f", dist[i, j]), x, y, gp = gpar(fontsize = 10))},
          heatmap_legend_param = list(title = "# DEGs"),
          top_annotation = NULL))
  
  return(dist)
}

Iterate = function(list, FUN, list2 = NULL, return.matrix = TRUE, ...){
  
  if(is.null(list2)){
    outlist = lapply(seq_along(list[1:(length(list)-1)]), function(i) {
      lapply((i+1):length(list), function(j){
        FUN(list[[i]], list[[j]], ...)
      })
    })
  } else {
     matrix <- sapply(list2, function(a) {
      sapply(list, function(b) FUN(a, b))
    })
    
    return(matrix)
  }
  
  if(return.matrix) {
    matrix = matrix(, nrow = length(list), ncol = length(list))
    matrix[lower.tri(matrix)] = unlist(outlist)
    matrix[upper.tri(matrix)] = t(matrix)[upper.tri(t(matrix))]
    rownames(matrix) = names(list)
    colnames(matrix) = names(list)
    return(matrix) 
  } else {
    return(outlist)
  }
}

RemoveYAxis = function(plt){
  return(plt + theme(axis.text.y = element_blank(), axis.title.y = element_blank()))
}

SeuratBootstrap = function(object, group.by = "seurat_clusters"){
  clusters = unique(object@meta.data[[group.by]])
  
  bootstrap_cols = sapply(clusters, function(cluster) {
    this_cols = which(object@meta.data[[group.by]] == cluster)
    return(sample(this_cols, length(this_cols), replace = TRUE))
  }) %>% unlist
  
  # Doesn't subset unless it's shorter than original object
  boot_object = suppressWarnings(object[,bootstrap_cols[1:length(bootstrap_cols)-1]])
  
  boot_object = RenameCells(boot_object, new.names = make.unique(Cells(boot_object)))
  
  return(boot_object)
}

leaveOneOutSupport = function(matrix, ref.tree){
  
  loo.trees.data = mclapply(seq_len(ncol(matrix)), function(col) {
    phydat = phyDat(matrix[,-col], type="USER", levels=c("0", "1"), compress = FALSE)
    tree = pratchet(phydat, trace=0, minit=100)
    list(tree = tree, parsimony = parsimony(tree, phydat))
  }, mc.cores = 16)
  
  message("Generated ", length(loo.trees.data), " LOO trees!")
  
  loo.trees = as.multiPhylo(lapply(loo.trees.data, function(x) x$tree)) %>% setNames(rownames(matrix))
  treeBS = plotBS(ref.tree, loo.trees, main="LOO support", type = "unrooted")
  ref.tree$loo.support = treeBS$node.label
  return(ref.tree)
}

FindGini = function(node){
  # TF switches at node
  tfs = gained.summary[[as.character(node)]]
  
  # Descendants
  descendants = treeRooted$tip.label[Descendants(treeRooted, as.numeric(node))[[1]] ]
  
  lapply(tfs, function(tf){
    
    data = data.table(label = orthotypes %in% descendants, expr = as.character(tf.expr.mask[tf,]))
    data.frame(tf = tf, 
               node = node, 
               gini = as.numeric(gini_impurities(data)[2,3]))
  })
}

nMonophyletic = function(...){
  length(which(FindMonophyletic(...)))
}

# For each TF, test whether it's marked OTs are monophyletic!
FindMonophyletic = function(datExpr, tree, exclude.single = TRUE){
  
  # Whether to exclude TF events specific to a single orthotype
  min = ifelse(exclude.single, 1, 0)
  max = ifelse(exclude.single, ncol(datExpr)-1, ncol(datExpr))
  
  monophyletic = sapply(rownames(datExpr), function(gene){
    orthotypes.marked = names(which(datExpr[gene,] > 0))
    ifelse(length(orthotypes.marked) > min & length(orthotypes.marked) < max, is.monophyletic(phy = tree, tips = orthotypes.marked), NA)
  })
  return(monophyletic)
}

PrettyTree = function(tree, gene, cols = c("0" = "black", "1" = "red", "0.5" = "gold")){
  # TitlePlot(
  ggtree(as.treedata(tree), layout = "equal_angle", aes(color = eval(parse(text = gene)))) + 
    geom_tippoint(aes(color=eval(parse(text = gene))), show.legend = TRUE)+
    geom_nodepoint(aes(color=eval(parse(text = gene))), show.legend = TRUE)+
    # scale_color_gradient2(low = "black", high ="red", mid = "yellow", midpoint = 0.5)+
    scale_color_manual(name = gene, values = cols)+
    # geom_tiplab(geom = "text")+
    geom_text_repel(aes(label=name), max.overlaps = Inf, show.legend = FALSE, seed = 42)+
    # geom_treescale()+
    theme(legend.position = "none")
  # gene)
}

ShuffleRows = function(matrix){
  shuffled = t(apply(matrix, 1, function(row) sample(row)))
  colnames(shuffled) = colnames(matrix)
  return(shuffled)
}

MaxParsimonyTree = function(datExpr, minit){
  phydat = phyDat(t(datExpr), type="USER", levels=c("0", "1"), compress = FALSE)
  treeBinary  <- pratchet(phydat, trace = 0, minit=minit)
  treeBinary  <- acctran(treeBinary, phydat)
  treeBinary  <- di2multi(treeBinary, tol = 1e-8)
  if(inherits(treeBinary, "multiPhylo")){
    treeBinary <- unique(treeBinary)
  }
  
  # treeBinary$node.label = (length(treeBinary$tip.label)+1):((length(treeBinary$tip.label)+treeBinary$Nnode))
  # treeBinary <- root(phy = treeBinary, node = 71)
  # plot(treeBinary)
  
  # message("Parsimony score: ", parsimony(treeBinary, phydat))
  return(list(tree = treeBinary, 
              parsimony = parsimony(treeBinary, phydat), 
              monophyletic = FindMonophyletic(datExpr, treeBinary)))
}

BinarizeExpression = function(matrix, return.vector = TRUE, plot = FALSE, seed = 42, replace.na = 0, noise = TRUE, verbose = FALSE, alpha = 0.05, noise.sd = 1){
  
  library(MASS)
  library(mixtools)
  set.seed(seed)
  
  # Remove rows with zero variance, replace other NAs with replace.na
  matrix = matrix %>% replace(is.na(.), replace.na)
  input = matrix[apply(matrix, 1, var) != 0, ]
  
  vector_encoding = tryCatch({
    
    # Fit a single Gaussian distribution
    single_gaussian <- fitdistr(unlist(input), "normal")
    log_likelihood_single <- as.numeric(logLik(single_gaussian))
    
    # Fit a mixture of two Gaussian distributions
    model = if(noise) mixtools::normalmixEM(unlist(input) + rnorm(length(unlist(input)), mean = 0, noise.sd = noise.sd)) else mixtools::normalmixEM(unlist(input), k = 2)
    if(plot) plot(model, density=TRUE, breaks = 30)
    log_likelihood_mixture <- model$loglik
    
    # Perform the Likelihood Ratio Test
    LRT_statistic <- -2 * (log_likelihood_single - log_likelihood_mixture)
    p_value <- pchisq(LRT_statistic, df = 3, lower.tail = FALSE)
    if(verbose) message('LRT p-value: ', p_value)
    
    if(p_value < alpha){
      mask = matrix(model$posterior[,which.max(model$mu)], 
                    # ifelse(model$posterior[,which.max(model$mu)] > 0.8, 1, # 1 if confidently expressed, 0.5 if not confident, and 0 if confidently unexpressed
                    # ifelse(model$posterior[,which.max(model$mu)] > 0.2, 0.5, 0)),
                    nrow = nrow(input), 
                    ncol = ncol(input), 
                    dimnames = dimnames(input))
      vector_encoding = apply(mask, 2, mean)
      if(return.vector) vector_encoding else mask
    } else {
      NULL
    }
    
  }, error = function(e) {
    print(e)
  })
  vector_encoding
}

UpperCase_genes = function(object, integration = FALSE){
  
  object@misc$orig.features = rownames(object@assays$RNA@counts)
  
  rownames(object@assays$RNA@counts) <- make.unique(toupper(rownames(object@assays$RNA@counts)))
  rownames(object@assays$RNA@data)<- make.unique(toupper(rownames(object@assays$RNA@data)))
  rownames(object@assays$RNA@scale.data)<- make.unique(toupper(rownames(object@assays$RNA@scale.data)))
  object@assays$RNA@var.features <- make.unique(toupper(object@assays$RNA@var.features))
  #rownames(object@assays$RNA@meta.features)<- toupper(rownames(object@assays$RNA@meta.features))
  
  if(integration){
    rownames(object@assays$integrated@counts)<- make.unique(toupper(rownames(object@assays$integrated@counts)))
    rownames(object@assays$integrated@data)<- make.unique(toupper(rownames(object@assays$integrated@data)))
    rownames(object@assays$integrated@scale.data)<- make.unique(toupper(rownames(object@assays$integrated@scale.data)))
    object@assays$integrated@var.features <- make.unique(toupper(object@assays$integrated@var.features))
    #rownames(object@assays$integration@meta.features)<- toupper(rownames(object@assays$integration@meta.features))
  }
  
  return(object)
}

PcaContribution = function(pca, dim){
  summary = as.data.frame(pca$rotation) %>% arrange(desc(!!sym(dim)))
  contributors = summary[,dim]
  names(contributors) = rownames(summary)
  return(contributors)
}

VennDiagram3 <- function(list1, list2, list3, population = 17000, name1 = NULL, name2 = NULL, name3 = NULL) {
  library(eulerr)
  
  # Create a summary of the unique elements in each set and their intersections
  summary = list(
    unique(list1[!list1 %in% c(list2, list3)]), # Only in list1
    unique(list2[!list2 %in% c(list1, list3)]), # Only in list2
    unique(list3[!list3 %in% c(list1, list2)]), # Only in list3
    intersect(list1, list2)[!intersect(list1, list2) %in% list3], # In list1 and list2, but not list3
    intersect(list1, list3)[!intersect(list1, list3) %in% list2], # In list1 and list3, but not list2
    intersect(list2, list3)[!intersect(list2, list3) %in% list1], # In list2 and list3, but not list1
    Reduce(intersect, list(list1, list2, list3)) # In all three lists
  )
  
  # Named numeric vector for `eulerr`
  input = c(
    "A" = length(unique(list1)) - length(intersect(list1, c(list2, list3))),
    "B" = length(unique(list2)) - length(intersect(list2, c(list1, list3))),
    "C" = length(unique(list3)) - length(intersect(list3, c(list1, list2))),
    "A&B" = length(intersect(list1, list2)) - length(Reduce(intersect, list(list1, list2, list3))),
    "A&C" = length(intersect(list1, list3)) - length(Reduce(intersect, list(list1, list2, list3))),
    "B&C" = length(intersect(list2, list3)) - length(Reduce(intersect, list(list1, list2, list3))),
    "A&B&C" = length(Reduce(intersect, list(list1, list2, list3)))
  )
  
  # Default or user-provided names for the sets
  if (is.null(name1) | is.null(name2) | is.null(name3)) {
    names(input) = c(deparse(substitute(list1)), 
                     deparse(substitute(list2)), 
                     deparse(substitute(list3)), 
                     paste0(deparse(substitute(list1)), "&", deparse(substitute(list2))),
                     paste0(deparse(substitute(list1)), "&", deparse(substitute(list3))),
                     paste0(deparse(substitute(list2)), "&", deparse(substitute(list3))),
                     paste0(deparse(substitute(list1)), "&", deparse(substitute(list2)), "&", deparse(substitute(list3))))
  } else {
    names(input) = c(name1, name2, name3, 
                     paste0(name1, "&", name2), 
                     paste0(name1, "&", name3), 
                     paste0(name2, "&", name3), 
                     paste0(name1, "&", name2, "&", name3))
  }
  names(summary) = names(input)
  
  # Fit and plot the Euler diagram
  fit <- euler(input)
  print(
    plot(fit, 
         quantities = TRUE,
         fill = c('deepskyblue', 'orangered', 'forestgreen'),
         lty = 1,
         labels = list(font = 4))
  )
  
  # Hypergeometric p-value for the overlap of all three lists
  hyper_pval = phyper(length(Reduce(intersect, list(list1, list2, list3))) - 1, 
                      length(unique(list1)), 
                      population - length(unique(list1)), 
                      length(unique(list2)) + length(unique(list3)), 
                      lower.tail = FALSE, 
                      log.p = FALSE)
  
  summary$pval = hyper_pval
  
  return(summary)
}


VennDiagram = function(list1, list2, population = 17000, name1 = NULL, name2 = NULL){
  library(eulerr)
  
  summary = list(unique(list1[!list1 %in% list2]),
                 unique(list2[!list2 %in% list1]), 
                 intersect(list1,list2))
  
  # Input in the form of a named numeric vector
  input = c(
    "A" = length(unique(list1)) - length(intersect(list1,list2)),
    "B" = length(unique(list2)) - length(intersect(list1,list2)),
    "A&B" = length(intersect(list1,list2))
  )
  
  if(is.null(name1) | is.null(name2)){
    names(input) = c(deparse(substitute(list1)), deparse(substitute(list2)), paste0(deparse(substitute(list1)), "&", deparse(substitute(list2))))
  } else{
    names(input) = c(name1, name2, paste0(name1, "&", name2))
  }
  names(summary) = c(deparse(substitute(list1)), deparse(substitute(list2)), paste0(deparse(substitute(list1)), "&", deparse(substitute(list2))))
  fit <- euler(input)
  
  print(
    plot(fit, 
         quantities = TRUE,
         # fill = "transparent",
         fill = c('deepskyblue', 'orangered'),
         lty = 1,
         labels = list(font = 4))
  )
  
  # Hypergeoemtric p-value
  hyper_pval = phyper(length(intersect(list1,list2))-1, # number of white balls drawn
                      length(unique(list1)), # number of white balls
                      population-length(unique(list1)), # number of black balls
                      length(unique(list2)), # number of balls drawn
                      lower.tail=FALSE, 
                      log.p=FALSE)
  
  summary$pval = hyper_pval
  
  return(summary)
}

SAMapAlignmentHeatmap3 = function(MappingTable, type_cols, species_cols, 
                                  cols = c("white", "#fbd9d3", "#ffb09c", "red", "darkred"), 
                                  add.breaks = TRUE, show_annotation_legend = FALSE,
                                  max.value = 1,
                                  table_order = NULL, 
                                  type.order = NULL, 
                                  species.order = NULL,
                                  species.use = NULL, 
                                  font.size = 12,
                                  font.family = 'sans',
                                  return.matrix = FALSE,
                                  label = 5,
                                  minor.breaks = gpar(col = "grey", lwd = 0.5),
                                  ...){
  
  # Cell type ordering
  species = str_split_fixed(colnames(MappingTable), "_", 2)[,1]
  types = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,1]
  idents = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,2]
  
  if(!is.null(species.use)){
    MappingTable = MappingTable[species %in% species.use, species %in% species.use]
    
    # re-parse
    species = str_split_fixed(colnames(MappingTable), "_", 2)[,1]
    types = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,1]
    idents = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,2]
  }
  
  if(is.null(type.order)) type.order = unique(types)
  if(is.null(species.order)) species.order = species.use
  sort_key = data.frame(species = species, types = types, idents = idents) %>% 
    arrange(factor(types, levels = type.order), factor(species, levels = species.order))
  table_order = paste0(sort_key$species, '_', sort_key$types, '.', sort_key$idents)
  table_order = sub("\\.$", "", table_order)
  
  if(is.null(table_order)) {
    cnames = factor(colnames(MappingTable), levels = table_order) 
  } else {
    table_order = as.vector(sapply(names(type_cols), function(type) sapply(names(species_cols), function(species) paste0(species, "_", type))))
    cnames = factor(colnames(MappingTable), levels = table_order)
  }
  MappingTable = MappingTable[order(cnames), order(cnames)]
  
  # Set same species alignment scores to NA since these are never computed
  MappingTable[outer(species[order(cnames)], species[order(cnames)], FUN = function(x,y) x == y)] = NA
  
  # Max value
  MappingTable.orig = MappingTable
  if(!is.null(max.value)) MappingTable[which(MappingTable > max.value)] = max.value
  
  # Type ordering
  row_km = factor(factor(types[order(cnames)], levels = names(type_cols)))
  
  if(return.matrix) return(MappingTable)
  
  print(Heatmap(MappingTable, 
                name = "Alignment\nscore", 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                cell_fun = if(!is.null(label)) function(j, i, x, y, width, height, fill) { 
                  if(!is.na(MappingTable.orig[i, j]))
                  grid.text(round(MappingTable.orig[i, j], 2), x, y, gp = gpar(fontsize = label))
                  } else {
                    NULL
                  },
                top_annotation = HeatmapAnnotation(species = species[order(cnames)], 
                                                   `cell type` = row_km,
                                                   col = list(`cell type` = type_cols, species = species_cols), 
                                                   border = TRUE, 
                                                   show_annotation_name = TRUE, 
                                                   show_legend = FALSE,
                                                   annotation_legend_param = list(title_gp = gpar(fontfamily = font.family, fontsize = font.size))), 
                left_annotation = rowAnnotation(species = species[order(cnames)], 
                                                `cell type` = row_km,
                                                col = list(`cell type` = type_cols, species = species_cols), 
                                                border = TRUE, 
                                                show_legend = show_annotation_legend,
                                                show_annotation_name = FALSE,
                                                annotation_legend_param = list(title_gp = gpar(fontfamily = font.family, fontsize = font.size))),
                col = cols, 
                heatmap_legend_param = list(
                    title = "Scaled\nexpression",
                    border = "black", # Black border around legend
                    legend_gp = gpar(col = "black"),
                    labels_gp = gpar(col = "black", fontfamily = font.family, fontsize = font.size-2), # Black tick labels
                    ticks_gp = gpar(col = "black"), 
                    title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
                column_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2),  
                row_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2), 
                rect_gp = minor.breaks,
                border_gp = gpar(col = "black", lty = 1), 
                ...
  ))
  
  # types = sort_key$types
  types = types[order(cnames)]
  
  # Add horizontal bars separating rows
  if(add.breaks){
    breaks = c(sapply(unique(types), function(x) which(types == x)[1])-1, length(types))
    for(index in breaks) {
      x_coord = index/length(types)
      
      # Vertical
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
      })
      
      # Horizontal
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(0, 1), c(1-x_coord, 1-x_coord), gp = gpar(lty = 1, lwd = 1))
      })
    }
  }
}

SAMapAlignmentHeatmap2 = function(MappingTable, type_palette, species_palette, 
                                  cols = c("white", "#fbd9d3", "#ffb09c", "red", "darkred"), 
                                  add.breaks = TRUE, show_annotation_legend = FALSE, 
                                  table_order = NULL, species.use = NULL, 
                                  type.order = NULL, species.order = NULL,
                                  rect_gp = gpar(col = "grey", lwd = 0.5),
                                  max.value = 1, font.family = 'sans', font.size = 12,
                                  ...){
  
  # Cell type ordering
  species = str_split_fixed(colnames(MappingTable), "_", 2)[,1]
  types = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,1]
  idents = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,2]
  
  if(!is.null(species.use)){
    MappingTable = MappingTable[species %in% species.use, species %in% species.use]
    
    # re-parse
    species = str_split_fixed(colnames(MappingTable), "_", 2)[,1]
    types = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,1]
    idents = str_split_fixed(str_split_fixed(colnames(MappingTable), "_", 2)[,2], '\\.', 2)[,2]
  }
  
  if(is.null(type.order)) type.order = unique(types)
  if(is.null(species.order)) species.order = species.use
  sort_key = data.frame(species = species, types = types, idents = idents) %>% 
    arrange(factor(types, levels = type.order), factor(species, levels = species.order))
  table_order = paste0(sort_key$species, '_', sort_key$types, '.', sort_key$idents)
  table_order = sub("\\.$", "", table_order)
  
  # table_order = as.vector(sapply(names(pr_palette), function(type) sapply(names(species_palette2), function(species) paste0(species, "_", type))))
  if(is.null(table_order)) cnames = factor(colnames(MappingTable), levels = sort(unique(colnames(MappingTable)))) else cnames = factor(colnames(MappingTable), levels = table_order)
  MappingTable = MappingTable[order(cnames), order(cnames)]
  
  # Set same species alignment scores to NA since these are never computed
  MappingTable[outer(species[order(cnames)], species[order(cnames)], FUN = function(x,y) x == y)] = NA
  
  # Type ordering
  row_km = factor(factor(types[order(cnames)], levels = names(type_palette)))
  
  MappingTable[MappingTable > max.value] = max.value
  
  print(Heatmap(MappingTable, 
                name = "Alignment\nscore", 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = HeatmapAnnotation(species = factor(species[order(cnames)], species.order), 
                                                   `cell type` = row_km,
                                                   col = list(`cell type` = type_palette, 
                                                              species = species_palette), 
                                                   border = TRUE, 
                                                   show_annotation_name = TRUE, 
                                                   show_legend = FALSE), 
                left_annotation = rowAnnotation(species = factor(species[order(cnames)], species.order), 
                                                `cell type` = row_km,
                                                col = list(`cell type` = type_palette, 
                                                           species = species_palette), 
                                                border = TRUE, 
                                                show_legend = show_annotation_legend,
                                                show_annotation_name = FALSE),
                col = cols, #c("white", "lightblue", "darkblue"), 
                # show_row_names = FALSE, 
                # show_column_names = FALSE, 
                rect_gp = rect_gp,
                border_gp = gpar(col = "black", lty = 1), 
                heatmap_legend_param = list(
                  title = "Alignment\nscore",
                  border = "black", # Black border around legend
                  legend_gp = gpar(col = "black"),
                  labels_gp = gpar(col = "black", fontfamily = font.family, fontsize = font.size-2), # Black tick labels
                  ticks_gp = gpar(col = "black"), 
                  title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
                ...
  ))
  
  types = sort_key$types
  
  # Add horizontal bars separating rows
  if(add.breaks){
    breaks = c(sapply(unique(types), function(x) which(types == x)[1])-1, length(types))
    for(index in breaks) {
      x_coord = index/length(types)
      
      # Vertical
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
      })
      
      # Horizontal
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(0, 1), c(1-x_coord, 1-x_coord), gp = gpar(lty = 1, lwd = 1))
      })
    }
  }
  
  df = data.frame(types = types, 
                  species = species[order(cnames)], 
                  idents = idents[order(cnames)])
  
  return(SmartMatrix(MappingTable, df, df))
}

SAMapAlignmentHeatmap = function(MappingTable, cols = c("white", "#fbd9d3", "#ffb09c", "red", "darkred"), add.breaks = TRUE, show_annotation_legend = FALSE, ...){
  
  # Cell type ordering
  table_order = as.vector(sapply(names(pr_palette), function(type) sapply(names(species_palette2), function(species) paste0(species, "_", type))))
  cnames = factor(colnames(MappingTable), levels = table_order)
  MappingTable = MappingTable[order(cnames), order(cnames)]
  types = str_split_fixed(colnames(MappingTable), "_", 2)[,2]
  species = str_split_fixed(colnames(MappingTable), "_", 2)[,1]
  
  # Set same species alignment scores to NA since these are never computed
  MappingTable[outer(species, species, FUN = function(x,y) x == y)] = NA
  
  # Type ordering
  row_km = factor(factor(types, levels = names(pr_palette)))
  
  rotate_and_save_image <- function(image_path, angle, output_dir) {
    image <- image_read(image_path)
    image <- image_rotate(image, angle)
    output_path <- file.path(output_dir, paste0("rotated_", basename(image_path)))
    image_write(image, output_path)
    return(output_path)
  }

  rotated_image_paths <- sapply(image_paths, rotate_and_save_image, angle = -90, output_dir = "../../figures/animals/black/")

  print(Heatmap(MappingTable, name = "Alignment\nscore", cluster_rows = FALSE, cluster_columns = FALSE, 
                width = unit(4, "in"), height = unit(4, "in"), 
                top_annotation = HeatmapAnnotation(species = anno_image(rotated_image_paths[match(str_split_fixed(colnames(MappingTable), "_", 2)[,1], names(image.files2))], 
                                                                        border = FALSE, 
                                                                        space = unit(0.1, "mm")),
                                                   `cell type` = row_km,
                                                   col = list(`cell type` = pr_palette), 
                                                   border = TRUE, 
                                                   show_annotation_name = TRUE, 
                                                   show_legend = FALSE), 
                left_annotation = rowAnnotation(species = anno_image(image.files[match(str_split_fixed(colnames(MappingTable), "_", 2)[,1], names(image.files2))], 
                                                                     border = FALSE, 
                                                                     space = unit(0.1, "mm")),
                                                `cell type` = row_km,
                                                col = list(`cell type` = pr_palette), 
                                                border = TRUE, 
                                                show_legend = show_annotation_legend,
                                                show_annotation_name = FALSE),
                col = cols, #c("white", "lightblue", "darkblue"), 
                show_row_names = FALSE, show_column_names = FALSE, 
                border_gp = gpar(col = "black", lty = 1), 
                ...
  ))
  
  # Add horizontal bars separating rows
  if(add.breaks){
    breaks = c(sapply(unique(types), function(x) which(types == x)[1])-1, length(types))
    for(index in breaks) {
      x_coord = index/length(types)
      
      # Vertical
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
      })
      
      # Horizontal
      decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
        grid.lines(c(0, 1), c(1-x_coord, 1-x_coord), gp = gpar(lty = 1, lwd = 1))
      })
    }
  }
}

RMSE = function(x, y){
  sqrt(mean((x-y)^2))
}

# Inspired by https://github.com/satijalab/seurat/issues/2520
PrettyUmap2 = function(object, group.by = "seurat_clusters", 
                       color.by = NULL, 
                       remove.space.x = TRUE, remove.space.y = TRUE, 
                       nbreaks = 20, geom.label = geom_text_repel, 
                       angle = 0, pad = 1.5, 
                       cols = NULL, label = TRUE, show.legend = FALSE, 
                       shuffle = TRUE, title = NULL,
                       flip.x = FALSE, flip.y = FALSE, 
                       density = FALSE, pt.alpha = 1, pt.size = 1, 
                       rasterise = FALSE, raster.dpi = 500,
                       centroid_fun = median, legend.point.size = 4, 
                       legend.font.size = 12, remove.times = 1,
                       return.object = FALSE, 
                       ...){
  
  if(inherits(object, 'Seurat')){
    data <- as.data.frame(Embeddings(object, reduction = 'umap')[(colnames(object)), c(1, 2)])
    colnames(data) = c('UMAP_1', 'UMAP_2')
    
    # Add grouping information
    data$group.by <- object@meta.data[[group.by]]

    # Add colors
    if(is.null(color.by)) {
      data$color.by <- object@meta.data[[group.by]]
    } else {
      data$color.by <- object@meta.data[[color.by]]
    }
  } else {
    data = object
    
    data$group.by = data[[group.by]]
    # data$color.by = data[[group.by]]
    
    # Add colors
    if(is.null(color.by)) {
      data$color.by <- object[[group.by]]
    } else {
      data$color.by <- object[[color.by]]
    }
  }
  
  
  # Remove rows where groups are NA
  # data = na.omit(data)
  
  remove_space_x = function(data){
    # largest x gap
    hist_x = hist(data$UMAP_1, plot = FALSE, breaks=nbreaks)
    rle = rle(hist_x$counts)
    zero_indices = which(rle$values == 0)
    if(length(zero_indices) > 0){
      counts = rle$lengths[zero_indices]
      longest_zero_index = zero_indices[which.max(counts)]
      firstZero = sum(rle$lengths[1:(longest_zero_index-1)])+1
      secondZero = firstZero + counts[which.max(counts)] - 1
      # secondZero = sum(rle$lengths[1:(longest_zero_index-2 + counts[which.max(counts)])])
      # Check
      stopifnot(hist_x$counts[firstZero] == 0)
      stopifnot(hist_x$counts[secondZero] == 0)
      x_break1 = hist_x$mids[firstZero]
      x_break2 = hist_x$mids[secondZero]
      data[data$UMAP_1 < x_break1,"UMAP_1"] = data[data$UMAP_1 < x_break1,"UMAP_1"] + (x_break2 - x_break1)
    }
    
    data
  }
  
  # Remove empty space by binning and removing chunks with no counts
  if(remove.space.x){
    for(i in seq_len(remove.times)){
      data = remove_space_x(data)
    }
  }
  
  remove_space_y = function(data){
    # largest y gap
    hist_y = hist(data$UMAP_2, plot = FALSE, breaks=nbreaks)
    rle = rle(hist_y$counts)
    zero_indices = which(rle$values == 0)
    if(length(zero_indices) > 0){
      counts = rle$lengths[zero_indices]
      longest_zero_index = zero_indices[which.max(counts)]
      firstZero = sum(rle$lengths[1:(longest_zero_index-1)])+1
      secondZero = firstZero + counts[which.max(counts)] - 1
      # Check 
      stopifnot(hist_y$counts[firstZero] == 0)
      stopifnot(hist_y$counts[secondZero] == 0)
      y_break1 = hist_y$mids[firstZero]
      y_break2 = hist_y$mids[secondZero]
      data[data$UMAP_2 < y_break1,"UMAP_2"] = data[data$UMAP_2 < y_break1,"UMAP_2"] + (y_break2 - y_break1)
    }
    
    data
  }
  
  if(remove.space.y){
    for(i in seq_len(remove.times)){
      data = remove_space_y(data)
    }
  }
  
  # Rotate by angle
  if(angle != 0){
    data[,c("UMAP_1", "UMAP_2")] = rotate_coords(data[,c('UMAP_1', 'UMAP_2')], angle, radians = TRUE) #spdep::Rotation(data[,c("UMAP_1", "UMAP_2")], angle) %>% 
      # as.data.frame %>% 
      # setNames(c("UMAP_1", "UMAP_2"))
  }
  
  # Flip across x-axis
  if(flip.y) data[,"UMAP_2"] = -1*data[,"UMAP_2"]
  
  # Flip across y-axis
  if(flip.x) data[,"UMAP_1"] = -1*data[,"UMAP_1"]
  
  if(return.object) {
    object@reductions$umap@cell.embeddings[, c(1, 2)] <- as.matrix(data[colnames(object), c('UMAP_1', 'UMAP_2')])
    return(object)
  }
  
  # Compute centroids for labeling
  raw_centroids = data.frame(x = sapply(unique(data$group.by), function(x) centroid_fun(data$UMAP_1[data$group.by == x])),
                             y = sapply(unique(data$group.by), function(x) centroid_fun(data$UMAP_2[data$group.by == x])),
                             NAME = unique(data$group.by))
  
  centroids <- do.call(rbind, lapply(unique(data$group.by), FUN = function(x) {
    group.data <- data[data$group.by == x, c('UMAP_1', 'UMAP_2')]
    nearest.point <- RANN::nn2(data = group.data[, 1:2], query = as.matrix(raw_centroids[raw_centroids$NAME == x, 1:2]), k = 1)$nn.idx
    data.frame(x = group.data[nearest.point, 1], 
               y = group.data[nearest.point, 2], 
               NAME = x)
  }))
  
  
  if(shuffle){
    data = data[sample(1:nrow(data), nrow(data)),]
  } else {
    data = data %>% arrange(desc(color.by))
  }
  
  # Plot
  plt = theme_umap(
    ggplot(data, aes(x=UMAP_1, y=UMAP_2))+
      {if(rasterise) ggrastr::rasterise(geom_point(show.legend = show.legend, aes(color = color.by), alpha = pt.alpha, size = pt.size, shape = 16, stroke = 0), dpi = raster.dpi)}+
      {if(!rasterise) geom_point(show.legend = show.legend, aes(color = color.by), alpha = pt.alpha, size = pt.size, shape = 16, stroke = 0)}+
      {if(density) geom_density2d(show.legend = FALSE, alpha = 1, color = 'black', linewidth = 0.1, linetype = 'dashed', bins = 5)}+
      # {if(remove.space.x) scale_x_break(c(x_break1, x_break2))}+
      # {if(remove.space.y) scale_y_break(c(y_break1, y_break2))}+
      {if(label) geom.label(data = centroids, aes(label = NAME, x = x, y = y), ...)}+
      # {if(repel) geom_label_repel(data = centroids, aes(label = NAME, x = x, y = y), ...)}+
      {if(!is.null(cols)) scale_color_manual(values = cols)}+
      # guides(colour = guide_legend(override.aes = list(alpha = 1)))+
      scale_y_continuous(limits = c(min(data$UMAP_2)-pad, max(data$UMAP_2)+pad))+
      scale_x_continuous(limits = c(min(data$UMAP_1)-pad, max(data$UMAP_1)+pad))+
      theme(legend.text = element_text(size = legend.font.size))+# text = element_text(family = 'Arial'))+
      ArialFont()+
      guides(color = guide_legend(override.aes = list(size = legend.point.size, alpha = 1)))
  )
  
  TitlePlot(plt + theme(legend.title=element_blank()), title = title)
}

# Inspired by https://github.com/satijalab/seurat/issues/2520
PrettyUmap = function(object, group.by = "seurat_clusters", 
                      remove.space.x = TRUE, remove.space.y = TRUE, 
                      nbreaks = 20, geom.label = geom_text_repel, angle = 0, pad = 0.5, 
                      cols = NULL, label = TRUE, show.legend = FALSE, shuffle = TRUE, 
                      density = FALSE, alpha = 0.03, pt.size = 3, raster.dpi = 500,
                      ...){
  
  if(inherits(object, 'Seurat')){
    data <- as.data.frame(Embeddings(object = object[["umap"]])[(colnames(object)), c(1, 2)])
    
    # Add colors
    data$color.by <- object@meta.data[[group.by]]
  } else {
    data = object
    data$color.by = data[[group.by]]
  }
  
  # Remove empty space by binning and removing chunks with no counts
  if(remove.space.x){
    # largest x gap
    hist_x = hist(data$UMAP_1, plot = FALSE, breaks=nbreaks)
    rle = rle(hist_x$counts)
    zero_indices = which(rle$values == 0)
    if(length(zero_indices) > 0){
      counts = rle$lengths[zero_indices]
      longest_zero_stretch = 
        firstZero = sum(rle$lengths[1:(which.max(rle$lengths)-1)])+1
      secondZero = sum(rle$lengths[1:which.max(rle$lengths)])
      x_break1 = hist_x$mids[firstZero]
      x_break2 = hist_x$mids[secondZero]
      data[data$UMAP_1 < x_break1,"UMAP_1"] = data[data$UMAP_1 < x_break1,"UMAP_1"] + (x_break2 - x_break1)
    }
  }
  
  if(remove.space.y){
    # largest y gap
    hist_y = hist(data$UMAP_2, plot = FALSE, breaks=nbreaks)
    rle = rle(hist_y$counts)
    firstZero = sum(rle$lengths[1:(which.max(rle$lengths)-1)])
    secondZero = sum(rle$lengths[1:which.max(rle$lengths)])
    y_break1 = hist_y$mids[firstZero]
    y_break2 = hist_y$mids[secondZero]
    data[data$UMAP_2 < y_break1,"UMAP_2"] = data[data$UMAP_2 < y_break1,"UMAP_2"] + (y_break2 - y_break1)
  }
  
  # Rotate by angle
  data[,c("UMAP_1", "UMAP_2")] = spdep::Rotation(data[,c("UMAP_1", "UMAP_2")], angle) %>% 
    as.data.frame %>% 
    setNames(c("UMAP_1", "UMAP_2"))
  
  # Compute centroids for labeling
  centroids = data.frame(x = sapply(unique(data$color.by), function(x) mean(data$UMAP_1[data$color.by == x])),
                         y = sapply(unique(data$color.by), function(x) mean(data$UMAP_2[data$color.by == x])), 
                         NAME = unique(data$color.by))
  
  if(shuffle){
    data = data[sample(1:nrow(data), nrow(data)),]
  } else {
    data = data %>% arrange(desc(color.by))
  }
  
  # Plot
  plt = theme_umap(
    ggplot(data, aes(x=UMAP_1, y=UMAP_2))+
      ggrastr::rasterise(geom_point(show.legend = show.legend, aes(color = color.by), alpha = alpha, size = pt.size), dpi = raster.dpi)+
      {if(density) geom_density2d(show.legend = FALSE, alpha = 1, color = 'black', linewidth = 0.1, linetype = 'dashed', bins = 5)}+
      # {if(remove.space.x) scale_x_break(c(x_break1, x_break2))}+
      # {if(remove.space.y) scale_y_break(c(y_break1, y_break2))}+
      {if(label) geom.label(data = centroids, aes(label = NAME, x = x, y = y), point.size = NA, ...)}+
      # {if(repel) geom_label_repel(data = centroids, aes(label = NAME, x = x, y = y), ...)}+
      {if(!is.null(cols)) scale_color_manual(values = cols)}+
      guides(colour = guide_legend (override.aes = list(alpha = 1)))+
      scale_y_continuous(limits = c(min(data$UMAP_2)-pad, max(data$UMAP_2)+pad))+
      scale_x_continuous(limits = c(min(data$UMAP_1)-pad, max(data$UMAP_1)+pad))
  ) 
  
  plt + theme(legend.title=element_blank())
}

DrawSankey4 = function(filepath, min.value = 0.1, species = c('ze', 'ch', 'li', 'op', 'rn', 'hs'), axis.width = 0.1){
  library(ggforce)
  
  species1 = 'ze'
  species2 = 'ch'
  species3 = 'li'
  species4 = 'op'
  species5 = 'rn'
  species6 = 'hs'
  
  color.mapping = pr_palette
  MappingTable = reshape2::melt(as.matrix(read.csv(filepath, row.names = 1)))
  # print(head(MappingTable))
  MappingTable = rbind(MappingTable[startsWith(as.character(MappingTable$Var1), species1) & startsWith(as.character(MappingTable$Var2), species2),], 
                       MappingTable[startsWith(as.character(MappingTable$Var1), species2) & startsWith(as.character(MappingTable$Var2), species3),], 
                       MappingTable[startsWith(as.character(MappingTable$Var1), species3) & startsWith(as.character(MappingTable$Var2), species4),], 
                       MappingTable[startsWith(as.character(MappingTable$Var1), species4) & startsWith(as.character(MappingTable$Var2), species5),], 
                       MappingTable[startsWith(as.character(MappingTable$Var1), species5) & startsWith(as.character(MappingTable$Var2), species6),]
  )
  
  # print(head(MappingTable))
  # MappingTable = MappingTable %>% arrange(factor(Var1, paste0(species1, "_", names(color.mapping))),
  #                                         factor(Var2, paste0(species2, "_", names(color.mapping))))
  
  # Ordering of strata
  # MappingTable$Var1 = factor(MappingTable$Var1, levels = paste0(species1, '_', names(pr_palette)))
  # MappingTable$Var2 = factor(MappingTable$Var2, levels = paste0(species2, '_', names(pr_palette)))
  # MappingTable$Var3 = factor(MappingTable$Var3, levels = paste0(species3, '_', names(pr_palette)))
  # MappingTable$Var4 = factor(MappingTable$Var4, levels = paste0(species4, '_', names(pr_palette)))
  # MappingTable$Var5 = factor(MappingTable$Var5, levels = paste0(species5, '_', names(pr_palette)))
  # MappingTable$Var6 = factor(MappingTable$Var6, levels = paste0(species6, '_', names(pr_palette)))
  MappingTable = gather_set_data(MappingTable, 1:2)
  
  print((MappingTable))
  MappingTable$names = factor(str_split_fixed(MappingTable$y, '_', 2)[,2], levels = names(pr_palette))
  
  # Color nodes
  colors <- rep(pr_palette, 6)
  names(colors) = c(paste0(species1, '_', names(pr_palette)), 
                    paste0(species2, '_', names(pr_palette)), 
                    paste0(species3, '_', names(pr_palette)), 
                    paste0(species4, '_', names(pr_palette)), 
                    paste0(species5, '_', names(pr_palette)), 
                    paste0(species6, '_', names(pr_palette)))
  print(colors)
  MappingTable = MappingTable[MappingTable$value > min.value,]
  # return(MappingTable)
  
  # Plot
  sn = ggplot(MappingTable, aes(x, id = id, split = names, value = value)) +
    geom_parallel_sets(axis.width = axis.width, fill = 'grey', alpha = 0.4) +
    geom_parallel_sets_axes(aes(fill = y), axis.width = axis.width, color = "black") +
    geom_parallel_sets_labels(colour = 'black',
                              angle = 0
                              # hjust = c(rep(1, length(unique(MappingTable$Var1))), rep(0, length(unique(MappingTable$Var2)))),
                              # nudge_x = c(rep(-0.08, length(unique(MappingTable$Var1))), rep(0.08, length(unique(MappingTable$Var2))))
    ) +
    scale_fill_manual(values = colors) +
    theme_void() +
    NoLegend() +
    scale_x_discrete(expand = expansion(add = 0.5))
  
  return(sn)
}

DrawSankey3 = function(input, min.value = 0.1, species1 = 'ch', species2 = 'li', axis.width = 0.07, label = TRUE, text_nudge = 0.07){
  library(ggforce)
  
  color.mapping = pr_palette
  MappingTable = reshape2::melt(input) #as.matrix(read.csv(filepath, row.names = 1)))
  MappingTable = MappingTable[startsWith(as.character(MappingTable$Var1), species1) & startsWith(as.character(MappingTable$Var2), species2), ]
  MappingTable = MappingTable %>% arrange(factor(Var1, paste0(species1, "_", names(color.mapping))),
                                          factor(Var2, paste0(species2, "_", names(color.mapping))))
  
  # Ordering of strata
  MappingTable$Var1 = factor(MappingTable$Var1, levels = paste0(species1, '_', names(pr_palette)))
  MappingTable$Var2 = factor(MappingTable$Var2, levels = paste0(species2, '_', names(pr_palette)))
  MappingTable = gather_set_data(MappingTable, 1:2)
  
  # MappingTable$Var1 = str_split_fixed(MappingTable$Var1, '_', 2)[,2]
  # MappingTable$Var2 = str_split_fixed(MappingTable$Var2, '_', 2)[,2]
  MappingTable$names = factor(str_split_fixed(MappingTable$y, '_', 2)[,2], levels = names(pr_palette))
  
  # Color nodes
  colors <- rep(pr_palette, 2)
  names(colors) = c(paste0(species1, '_', names(pr_palette)), 
                    paste0(species2, '_', names(pr_palette)))
  
  MappingTable = MappingTable[MappingTable$value > min.value,]
  
  # Plot
  sn = ggplot(MappingTable, aes(x, id = id, split = names, value = value)) +
    geom_parallel_sets(axis.width = axis.width, fill = 'grey', alpha = 0.4) +
    geom_parallel_sets_axes(aes(fill = y), axis.width = axis.width, color = "black") +
    {if(label) geom_parallel_sets_labels(colour = 'black', 
                                         angle = 0,
                                         hjust = c(rep(1, length(unique(MappingTable$Var1))), rep(0, length(unique(MappingTable$Var2)))),
                                         nudge_x = c(rep(-text_nudge, length(unique(MappingTable$Var1))), rep(text_nudge, length(unique(MappingTable$Var2)))))} +
    scale_fill_manual(values = colors) +
    theme_void() +
    NoLegend() +
    scale_x_discrete(expand = expansion(add = 0.5))
  
  return(sn)
}

GetTypeList = function(gene_pairs, PHOTORECEPTORS, species_ids){
  # gene.summary = list()
  type.list = list()
  
  # This recovers all 335 entries
  for(i in 1:(length(species_ids)-1)){
    type.list[[ species_ids[i] ]] = list()
    for(k in 1:(length(PHOTORECEPTORS))){
      type.list[[ species_ids[i] ]][[ PHOTORECEPTORS[k] ]] = list()
      for(j in (i+1):length(species_ids)){
        type.list[[ species_ids[i] ]][[ PHOTORECEPTORS[k] ]][[ species_ids[j] ]] = list()
        for(l in (1):length(PHOTORECEPTORS)){
          # element = 1
          column.name = paste0(species_ids[i], "_", PHOTORECEPTORS[k], ".", species_ids[j], "_", PHOTORECEPTORS[l])
          if(column.name %in% names(gene_pairs)){
            genes.1 = ExtractString(gene_pairs[[column.name]], paste0(species_ids[i], "_"), ";")
            genes.1 = genes.1[genes.1 != ""]
            genes.2 = ExtractString(gene_pairs[[column.name]], paste0(species_ids[j], "_"))
            genes.2 = genes.2[genes.2 != ""]
            type.list[[ species_ids[i] ]][[ PHOTORECEPTORS[k] ]][[ species_ids[j] ]][[ PHOTORECEPTORS[l] ]] = setNames(data.frame(genes.1, genes.2), names(species_ids)[c(i,j)])
            # element = element + 1
          }
        }
      }
    }
  }
  
  return(type.list)
}

DrawSankey2 = function(filepath, min.value = 0.1, species1, species2){
  
  color.mapping = pr_palette
  # color.mapping = c("grey", "violet", "blue", "green", "red", "magenta", "cyan", 'black', 'brown') %>% setNames(c('rod', 'UV', 'blue', 'green', 'red', 'principle', 'accessory', 'AC', 'HC'))
  MappingTable = reshape2::melt(as.matrix(read.csv(filepath, row.names = 1)))
  MappingTable = MappingTable[startsWith(as.character(MappingTable$Var1), species1) & startsWith(as.character(MappingTable$Var2), species2), ]
  MappingTable = MappingTable %>% arrange(factor(Var1, paste0(species1, "_", names(color.mapping))),
                                          factor(Var2, paste0(species2, "_", names(color.mapping))))
  
  alignment = list()
  alignment$nodes = data.frame(name = unique(c(MappingTable$Var1, MappingTable$Var2)))
  nodeIDs = alignment$nodes$name
  names(nodeIDs) = as.numeric(rownames(alignment$nodes))-1 # Zero index
  
  # Filtration
  alignment$links = MappingTable
  
  # Convert to index
  alignment$links$Var1 = as.numeric(names(nodeIDs[match(alignment$links$Var1, nodeIDs)]))
  alignment$links$Var2 = as.numeric(names(nodeIDs[match(alignment$links$Var2, nodeIDs)]))
  alignment$links = alignment$links[!is.na(alignment$links$value),]
  alignment$links = alignment$links[alignment$links$value >= min.value,]
  # alignment$links$value[is.na(alignment$links$value)] = 0
  rownames(alignment$links) = 1:nrow(alignment$links)
  
  # Color nodes
  alignment$nodes$node_color = str_split_fixed(alignment$nodes$name, "_", 2)[,2]
  colors <- paste(unique(color.mapping[match(alignment$nodes$node_color, names(color.mapping))]), collapse = '", "')
  colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
  alignment$nodes$display = ""
  
  # Plot
  sn = sankeyNetwork(Links = alignment$links, 
                     Nodes = alignment$nodes, Source = 'Var1',
                     Target = 'Var2', Value = 'value', NodeID = 'node_color',
                     fontSize = 12, nodeWidth = 30, NodeGroup = "node_color", 
                     colourScale = colorJS, fontFamily = "arial")
  
  saveNetwork(sn, "sn.html")
  webshot::webshot("sn.html", gsub(".csv", ".png", filepath), vwidth = 300, vheight = 400, zoom = 2)
  
  return(sn)
}

DrawSankey = function(filepath, min.value = 0.1, species1, species2){
  
  color.mapping = c("purple", "blue", "green", "red", "grey", "cyan", "magenta", 'black', 'brown') %>% setNames(c('UV', 'blue', 'green', 'red', 'rod', 'accessory', 'principle', 'AC', 'HC'))
  MappingTable = reshape2::melt(as.matrix(read.csv(filepath, row.names = 1)))
  alignment = list()
  alignment$nodes = data.frame(name = unique(MappingTable$Var2))
  nodeIDs = alignment$nodes$name
  names(nodeIDs) = as.numeric(rownames(alignment$nodes))-1 # Zero index
  
  # Filtration
  MappingTable = MappingTable[startsWith(as.character(MappingTable$Var1), species1) & startsWith(as.character(MappingTable$Var2), species2), ]
  alignment$links = MappingTable %>% arrange(factor(Var1, paste0(species1, "_", names(color.mapping))),
                                             factor(Var2, paste0(species2, "_", names(color.mapping))))
  
  # Convert to index
  alignment$links$Var1 = as.numeric(names(nodeIDs[match(alignment$links$Var1, nodeIDs)]))
  alignment$links$Var2 = as.numeric(names(nodeIDs[match(alignment$links$Var2, nodeIDs)]))
  alignment$links = alignment$links[!is.na(alignment$links$value),]
  alignment$links = alignment$links[alignment$links$value >= min.value,]
  # alignment$links$value[is.na(alignment$links$value)] = 0
  rownames(alignment$links) = 1:nrow(alignment$links)
  
  # Color nodes
  alignment$nodes$node_color = str_split_fixed(alignment$nodes$name, "_", 2)[,2]
  colors <- paste(unique(color.mapping[match(alignment$nodes$node_color, names(color.mapping))]), collapse = '", "')
  colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
  alignment$nodes$display = ""
  
  # Plot
  sn = sankeyNetwork(Links = alignment$links, 
                     Nodes = alignment$nodes, Source = 'Var1',
                     Target = 'Var2', Value = 'value', NodeID = 'node_color',
                     fontSize = 12, nodeWidth = 30, NodeGroup = "node_color", 
                     colourScale = colorJS, fontFamily = "arial")
  
  saveNetwork(sn, "sn.html")
  webshot::webshot("sn.html", gsub(".csv", ".png", filepath), vwidth = 300, vheight = 400, zoom = 2)
  
  return(sn)
}

RenameFeatures = function(object, new.names){
  
  if(length(unique(new.names)) != length(new.names)) {
    message("Found ", length(new.names)-length(unique(new.names))," duplicates...making unique!")
    new.names = make.unique(new.names)
  }
  
  # Before features
  before.renamed = rownames(object@assays$RNA@counts)
  
  rownames(object@assays$RNA@counts) <- (new.names)
  rownames(object@assays$RNA@data) <- (new.names)
  rownames(object@assays$RNA@meta.features) <- (new.names)
  rownames(object@assays$RNA@scale.data) <- new.names[match(rownames(object@assays$RNA@scale.data), before.renamed)]
  
  # Would need to match order for this
  if('integrated' %in% names(object@assays)){
    integrated.features = rownames(object@assays$integrated@data)
    
    # No counts in integrated!
    # rownames(object@assays$integrated@counts) <- new.names[match(integrated.features, before.renamed)]
    rownames(object@assays$integrated@data) <- new.names[match(integrated.features, before.renamed)]
    rownames(object@assays$integrated@meta.features) <- new.names[match(integrated.features, before.renamed)]
  }
  
  return(object)
}

SelectFeatures = function(object, features = NULL, nfeatures = 20000, min.pct.expressed = 10, group.by = "seurat_clusters", assay = "RNA", return.pct = FALSE){
  
  object = FindVariableFeatures(object, selection.method = "vst", nfeatures = nrow(object@assays[[assay]]@data), verbose = FALSE)
  
  # Filter out lowly expressed genes
  if(min.pct.expressed > 0){
    if(is.null(features)) features = object@assays[[assay]]@var.features
    pct.expressed = PercentageExpressed2(object, features = features, group.by = group.by)
    pct.expressed$max.exp = apply(pct.expressed, 1, max)
    if(return.pct) return(pct.expressed)
    features.use = features[pct.expressed$max.exp >= min.pct.expressed]
    features.use = head(features.use, nfeatures)
    message("Using top ", length(features.use)," variable features expressed in at least ", min.pct.expressed, "% cells of at least one cluster")
  } else {
    message('min.pct.expr is zero, giving back var.features...')
    features.use = object@assays[[assay]]@var.features
  }
  
  return(features.use)
}

PseudoBulkPCA = function(object, 
                         assay = "RNA", 
                         group.by = "seurat_clusters", 
                         features = NULL, 
                         title = NULL, 
                         remove.genes = NULL,
                         return.plot = TRUE, 
                         min.pct.expressed = 0, 
                         nfeatures = 2000, 
                         axes = c('PC1', 'PC2'),
                         scale = FALSE, 
                         binary = NULL, 
                         size = 1, 
                         shape = 1, 
                         cols = NULL, 
                         label = TRUE, 
                         annotate.by = NULL, 
                         return.dist = TRUE, 
                         nPCs = 5, 
                         verbose = FALSE, 
                         font.size = 4, 
                         dist.method = 'euclidean'
                         ){
  
  if(inherits(object, "Seurat")){
    
    if(!is.null(features)){
      avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = features)[[assay]]))
    } else if(min.pct.expressed == 0){
      object = FindVariableFeatures(object, nfeatures = nfeatures)
      features.use = object@assays[[assay]]@var.features
      avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = features.use)[[assay]]))
    } else {
      # Feature selection to avoid noisy lowly expressed genes that get amplified after scaling
      features.use = SelectFeatures(object, nfeatures = nfeatures, min.pct.expressed = min.pct.expressed, group.by = group.by)
      
      # Get average expression
      avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = features.use)[[assay]]))
    }
  } else {
    avg_exp = object
  }
  
  # Remove any genes? 
  if(!is.null(remove.genes)){
    avg_exp = avg_exp[!rownames(avg_exp) %in% remove.genes,]
  }
  
  message("Feature matrix of size: ", paste0(dim(avg_exp), collapse = 'x'))
  
  if(!is.null(binary)){
    avg_exp_old = avg_exp
    avg_exp[avg_exp_old >= binary] = 1
    avg_exp[avg_exp_old < binary] = 0
  }
  
  # Compute principle components
  res.pca <- prcomp(t(avg_exp), scale = scale)
  var_exp = (res.pca$sdev^2)/sum(res.pca$sdev^2)*100
  res.pca$var_exp = var_exp
  
  if(verbose) print(barplot(var_exp))
  
  # Plot data projected onto PC1 and PC2
  df = as.data.frame(res.pca$x)
  df$label = rownames(df)
  if(!is.null(annotate.by)) {
    if(inherits(object, 'Seurat')){
      df$group = Metadata(object, group.by, annotate.by)[,annotate.by]
    } else {
      df$group = annotate.by
    }
  } else {
    df$group = df$label
  }
  
  # Return top nPCs
  if(return.dist){
    if(length(nPCs) > 1){
      res.pca$dist = lapply(nPCs, function(nPC) as.matrix(dist(df[,1:nPC], method = dist.method)))
    } else {
      res.pca$dist = as.matrix(dist(df[,1:nPCs], method = dist.method))
    }
  }
  
  if(!return.plot) return(res.pca)
  
  axis1 = as.numeric(gsub('PC', '', axes[[1]]))
  axis2 = as.numeric(gsub('PC', '', axes[[2]]))
  
  TitlePlot(
    ggplot(df, aes(x = !!sym(axes[[1]]), y = !!sym(axes[[2]]), label = label, color = group))+
    geom_point(size = size, shape = shape)+
    {if(label) geom_text_repel(show.legend = FALSE, size = font.size)}+
    {if(!is.null(cols)) scale_color_manual(values = cols)}+
    labs(x = paste0(axes[[1]], " (", signif(var_exp[axis1], 2), "%)"), y = paste0(axes[[2]], " (", signif(var_exp[axis2], 2), "%)"))+
    theme_dario()+
    theme(legend.key=element_rect(fill="white")),
    title = title
  )
  
}

SAMapHeatmap4 = function(objectList, 
                         gene.key, 
                         type_palette, 
                         species_palette,
                         species.use = c("Chicken", "Zebrafish"), 
                         types.use = c("red", "green", "blue", "UV", "rod"), 
                         min.z.score = -2,
                         max.z.score = 2, 
                         col.low = 'grey', 
                         col.high = "#584B9FFF",
                         rotate = FALSE,
                         dotplot = FALSE, 
                         dot.scale.factor = 0.05,
                         color.dot.by.species = FALSE,
                         pseudocount = 1,
                         max.pct = 100,
                         col_fun = NULL,
                         names.show = NULL,
                         font.family = 'ArialMT',
                         font.size = 12,
                         show_annotation_legend = FALSE,
                         types.order = NULL,
                         add.breaks = FALSE,
                         mybreaks = NULL,
                         add.rectangles = FALSE, 
                         row.breaks = NULL,
                         col.breaks = NULL,
                         species = NULL, 
                         lwd = 0.6,
                         lty = 1, 
                         ...){
  
  args = list(...)
  
  # Set col_fun
  if(is.null(col_fun)){
    col_fun = circlize::colorRamp2(c(min.z.score, 0, max.z.score), c("blue", "white", "red"))
  }
  
  # Scale data
  scaled.expr = lapply(species.use, function(species) {
    message('scaling for species: ', species)
    # if(species == 'Zebrafish') browser()
    # Edit: do scaling with all the clusters, then subset to the clusters plotted. That way if we have only one cluster, the scaled value isn't NA
    # norm.expr = AverageExpression(subset(objectList[[species]], annotated %in% types.use), group.by = "annotated", slot = "data", assay = 'RNA')$RNA
    norm.expr = AverageExpression(objectList[[species]], group.by = "annotated", slot = "data", assay = 'RNA')$RNA
    scaled.expr = t(scale(t(norm.expr)))
    scaled.expr = scaled.expr[match(gene.key[[species]], rownames(scaled.expr)),]
    
    # Now we subset to types of interest
    scaled.expr = scaled.expr[,intersect(types.use, colnames(scaled.expr)), drop = FALSE]
    colnames(scaled.expr) = colnames(scaled.expr)#paste0(species, ' ', colnames(scaled.expr))
    
    # Truncate values to max.z.score
    scaled.expr[scaled.expr < min.z.score] = min.z.score
    scaled.expr[scaled.expr > max.z.score] = max.z.score
    scaled.expr
  })
  
  # Percentage expression
  pct.expr = lapply(species.use, function(species) {
    message('percentage expressed for species: ', species)
    # if(species == 'Zebrafish') browser()
    pct.expr = PercentageExpressed2(subset(objectList[[species]], annotated %in% types.use), 
                                    features = setdiff(na.omit(unique(gene.key[[species]])), ''), 
                                    group.by = 'annotated')
    pct.expr = pct.expr[match(gene.key[[species]], rownames(pct.expr)),,drop = FALSE]
    colnames(pct.expr) = colnames(pct.expr)#paste0(species, ' ', colnames(pct.expr))
    pct.expr
  })
  
  # Combine 
  # gene.key = cbind(gene.key[,species.use], row_annotation = gene.key$row_annotation) # Sort by species so that gene names appear in same order
  if(is.null(names.show)) {
    gene.names = apply(gene.key, 1, function(x) paste0(x, collapse = "/"))
  } else {
    gene.names = make.unique(gene.key[[names.show]])
  }
  
  combined.expr = do.call(cbind, scaled.expr)
  rownames(combined.expr) = gene.names
  
  combined.pct = do.call(cbind, lapply(pct.expr, function(x) as.matrix(x)))
  rownames(combined.pct) = gene.names
  
  # Ordering matrix
  if(is.null(types.order)) types.order = as.vector(apply(outer(types.use, species.use, FUN = function(x, y) paste0(y, " ", x)), 1, function(x) x))
  
  # not_found = setdiff(types.order, colnames(combined.expr))
  not_found = setdiff(types.order, colnames(combined.pct))
  if(length(not_found) > 0) stop('Did not find ', paste0(not_found, collapse = ', '))
  combined.expr = combined.expr[, types.order] #order(factor(colnames(combined.expr), levels = types.order))]
  combined.pct = combined.pct[, types.order] #order(factor(colnames(combined.pct), levels = types.order))]
  print(dim(combined.expr))
  
  # Ordering for annotations
  types = str_split_fixed(colnames(combined.expr), " ", 2)[,2]
  classes = ExtractString(types, before = '_', after = '-')
  if(is.null(species)) species = str_split_fixed(colnames(combined.expr), " ", 2)[,1]
  row_split = paste0(gene.key$row_annotation, ' markers')
  
  if("one2one" %in% colnames(gene.key)){
    asterisks = ifelse(gene.key$one2one == 1, "*", "")
    rownames(combined.expr) = paste0(rownames(combined.expr), asterisks)
  }
  
  if(rotate){
    # Transpose
    combined.expr = t(combined.expr)
    combined.pct = t(combined.pct)
  } 
  
  if(dotplot){
    min.z.score = -1
    if(color.dot.by.species){
      col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c('white', 'grey25'))
    } else{
      col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c(col.low, col.high))
    }
    dot_col = function(species, expression) circlize::colorRamp2(c(min.z.score, max.z.score), c('white', species_palette[[species]]))(expression)
    cell_fun = function(j, i, x, y, w, h, fill){
      if(color.dot.by.species){
        grid.circle(x=x,
                    y=y,
                    r=unit(combined.pct[i, j]/100 * dot.scale.factor, "cm"),
                    gp=gpar(fill = dot_col(species[[j]], combined.expr[i, j]), col = NA))
      } else {
        grid.circle(x=x,
                    y=y,
                    r=unit(combined.pct[i, j]/100 * dot.scale.factor, "cm"),
                    gp=gpar(fill = col_fun(combined.expr[i, j]), col = NA))
      } 
    }
    
    rect_gp = gpar(type = "none")
    
    # Set this for color legend; edit: setting above
    # col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c('white', 'grey25'))
    
    # Create a dot legend
    lgd.values = c(0+pseudocount, 25+pseudocount, 50+pseudocount, 75+pseudocount, 100+pseudocount)
    lgd.values[lgd.values > (max.pct + pseudocount)] = (max.pct + pseudocount)
    lgd = Legend(labels = seq(0,100, by = 25), title = "Percent\nexpressed",
                 graphics = list(
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[1]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[2]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[3]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[4]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[5]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA))
                 ), title_gp = gpar(fontfamily = "sans", fontface = "plain"))
  } else {
    cell_fun = NULL
    rect_gp = gpar(col = NA)
    lgd = NULL
  }
  
  # Plot
  if(rotate){
    draw(Heatmap(combined.expr,
                 name = "mat",
                 col = col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE, 
                 # row_split = factor(types, levels = unique(types)),
                 # column_split = factor(row_split, levels = unique(row_split)),
                 column_title = NULL,
                 row_title = NULL, 
                 border = TRUE, 
                 cell_fun = cell_fun,
                 rect_gp = rect_gp,
                 # left_annotation = rowAnnotation(
                 #   # `Cell type` = anno_simple(gp = gpar(fill = type_palette[match(unique(types), names(type_palette))], 
                 #                                                                    # col = 'white')),
                 #                                                          #labels = unique(types), 
                 #                                                          #labels_gp = gpar(col = "black", fontsize = 10)),
                 #                                 Species = factor(species, levels = species.use),
                 #                                 # annotation_name_side = "bottom", 
                 #                                 show_annotation_name = FALSE,
                 #                                 show_legend = show_annotation_legend, 
                 #                                 border = TRUE, 
                 #                                 col = list(`Cell type` = type_palette, Species = species_palette), 
                 #                                 annotation_legend_param = list(title_gp = gpar(fontfamily = font.family, fontsize = font.size), border = 'black')), 
                 heatmap_legend_param = list(
                   title = "Scaled\nexpression",
                   border = "black", # Black border around legend
                   legend_gp = gpar(col = "black"),
                   labels_gp = gpar(col = "black", fontfamily = font.family, fontsize = font.size-2), # Black tick labels
                   ticks_gp = gpar(col = "black"), 
                   title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
                 row_names_gp = gpar(fontsize = 10, family = font.family),
                 column_names_gp = gpar(fontsize = 10, fontface = "italic", family = font.family),
                 column_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2),  
                 row_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2), 
                 ...), 
         # annotation_legend_list = lgd
    )
  } else {
    draw(Heatmap(combined.expr,
                 name = "mat",
                 col = col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE, 
                 row_split = factor(row_split, levels = unique(row_split)),
                 # column_split = factor(types, levels = unique(types)),
                 column_title = NULL,
                 row_title = NULL,   
                 border = TRUE, 
                 cell_fun = cell_fun,
                 # cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", combined.pct[i, j]), x, y, gp = gpar(fontsize = 10))},
                 rect_gp = rect_gp,
                 # top_annotation = HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = type_palette[match(unique(types), names(type_palette))], 
                 #                                                                       col = 'white'),
                 #                                                             labels = unique(types), 
                 #                                                             labels_gp = gpar(col = "black", fontsize = 10)),
                 #                                    Species = factor(species, levels = species.use),
                 #                                    show_annotation_name = FALSE,
                 #                                    # annotation_name_side = 'none', 
                 #                                    show_legend = show_annotation_legend, 
                 #                                    border = TRUE, 
                 #                                    col = list(`Cell type` = type_palette, Species = species_palette), 
                 #                                    annotation_legend_param = list(title_gp = gpar(fontfamily = font.family, fontsize = font.size), border = 'black')), 
                 heatmap_legend_param = list(
                   title = "Scaled\nexpression",
                   border = "black",      # Black border around legend
                   legend_gp = gpar(col = "black"),
                   labels_gp = gpar(col = "black", fontfamily = font.family, fontsize = font.size-2), # Black tick labels
                   ticks_gp = gpar(col = "black"), 
                   title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
                 column_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2),  
                 row_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2), 
                 row_names_gp = gpar(fontsize = 10, fontface = "italic", family = font.family),
                 column_names_gp = gpar(fontsize = 10, family = font.family),
                 ...), 
         annotation_legend_list = lgd)
  }
  
  # Add horizontal bars separating rows
  if(add.breaks){
    breaks = c(sapply(unique(mybreaks), function(x) which(mybreaks == x)[1])-1, length(mybreaks))
    for(index in breaks) {
      x_coord = index/length(types)
      
      # # Vertical
      # decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
      #   grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
      # })
      
      # Horizontal
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(0, 1), c(1-x_coord, 1-x_coord), gp = gpar(lty = 1, lwd = 0.75))
      })
    }
  }
  
  if(add.rectangles){
    y.breaks = c(sapply(unique(row.breaks), function(x) which(row.breaks == x)[1])-1, length(row.breaks))
    x.breaks = c(sapply(unique(col.breaks), function(x) which(col.breaks == x)[1])-1, length(col.breaks))
    
    
    for(i in seq(1, length(x.breaks)-1)) {
      y1 = y.breaks[[i]]/length(row.breaks)
      y2 = y.breaks[[i+1]]/length(row.breaks)
      x1 = x.breaks[[i]]/length(col.breaks)
      x2 = x.breaks[[i+1]]/length(col.breaks)
      
      # Make a box 
      # lwd = 0.5
      
      # left
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(x1, x1), c(1-y1, 1-y2), gp = gpar(lty = lty, lwd = lwd))
      })
      
      # top
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(x1, x2), c(1-y1, 1-y1), gp = gpar(lty = lty, lwd = lwd))
      })
      
      # right
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(x2, x2), c(1-y1, 1-y2), gp = gpar(lty = lty, lwd = lwd))
      })
      
      # bottom
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(x1, x2), c(1-y2, 1-y2), gp = gpar(lty = lty, lwd = lwd))
      })
      
      # Horizontal
      # decorate_heatmap_body("mat", row_slice = 1, {
      #   grid.lines(c(0, 1), c(1-y1, 1-y1), gp = gpar(lty = 1, lwd = 0.2))
      # })
    }
  
    # for(i in seq(1, length(x.breaks)-1)) {
    #   y1 = y.breaks[[i]]/length(row.breaks)
    #   y2 = y.breaks[[i+1]]/length(row.breaks)
    #   x1 = x.breaks[[i]]/length(col.breaks)
    #   x2 = x.breaks[[i+1]]/length(col.breaks)
    #   
    #   # Draw rectangle
    #   grid.rect(
    #     x = (x1 + x2) / 2,
    #     y = 1-((y1 + y2) / 2), # flip across y-axis?
    #     width  = abs(x2 - x1),
    #     height = abs(y2 - y1),
    #     gp = gpar(col = "black", lwd = 0.5, fill = NA)
    #   )
    #   
    # }
  }
}

SAMapHeatmap3 = function(objectList, 
                         gene.key, 
                         type_palette, 
                         species_palette,
                         species.use = c("Chicken", "Zebrafish"), 
                         types.use = c("red", "green", "blue", "UV", "rod"), 
                         min.z.score = -2,
                         max.z.score = 2, 
                         col.low = 'grey', 
                         col.high = "#584B9FFF",
                         rotate = FALSE,
                         dotplot = FALSE, 
                         dot.scale.factor = 0.05,
                         color.dot.by.species = FALSE,
                         pseudocount = 1,
                         max.pct = 100,
                         col_fun = NULL,
                         names.show = NULL,
                         font.family = 'ArialMT',
                         font.size = 12,
                         show_annotation_legend = FALSE,
                         types.order = NULL,
                         add.breaks = TRUE,
                         ...){
  
  args = list(...)
  
  # Set col_fun
  if(is.null(col_fun)){
    col_fun = circlize::colorRamp2(c(min.z.score, 0, max.z.score), c("blue", "white", "red"))
  }
  
  # Scale data
  scaled.expr = lapply(species.use, function(species) {
    message('scaling for species: ', species)
    # if(species == 'Zebrafish') browser()
    # Edit: do scaling with all the clusters, then subset to the clusters plotted. That way if we have only one cluster, the scaled value isn't NA
    # norm.expr = AverageExpression(subset(objectList[[species]], annotated %in% types.use), group.by = "annotated", slot = "data", assay = 'RNA')$RNA
    norm.expr = AverageExpression(objectList[[species]], group.by = "annotated", slot = "data", assay = 'RNA')$RNA
    scaled.expr = t(scale(t(norm.expr)))
    scaled.expr = scaled.expr[match(gene.key[[species]], rownames(scaled.expr)),]
    
    # Now we subset to types of interest
    scaled.expr = scaled.expr[,intersect(types.use, colnames(scaled.expr)), drop = FALSE]
    colnames(scaled.expr) = paste0(species, ' ', colnames(scaled.expr))
    
    # Truncate values to max.z.score
    scaled.expr[scaled.expr < min.z.score] = min.z.score
    scaled.expr[scaled.expr > max.z.score] = max.z.score
    scaled.expr
  })
  
  # Percentage expression
  pct.expr = lapply(species.use, function(species) {
    message('percentage expressed for species: ', species)
    # if(species == 'Zebrafish') browser()
    pct.expr = PercentageExpressed2(subset(objectList[[species]], annotated %in% types.use), 
                                    features = setdiff(na.omit(unique(gene.key[[species]])), ''), 
                                    group.by = 'annotated')
    pct.expr = pct.expr[match(gene.key[[species]], rownames(pct.expr)),,drop = FALSE]
    colnames(pct.expr) = paste0(species, ' ', colnames(pct.expr))
    pct.expr
  })
  
  # Combine 
  if(is.null(names.show)) {
    gene.names = apply(gene.key, 1, function(x) paste0(x, collapse = "/"))
  } else {
    gene.names = make.unique(gene.key[[names.show]])
  }
  
  combined.expr = do.call(cbind, scaled.expr)
  rownames(combined.expr) = gene.names
  
  combined.pct = do.call(cbind, pct.expr)
  rownames(combined.pct) = gene.names
  
  # Ordering matrix
  if(is.null(types.order)) types.order = as.vector(apply(outer(types.use, species.use, FUN = function(x, y) paste0(y, " ", x)), 1, function(x) x))
  combined.expr = combined.expr[, order(factor(colnames(combined.expr), levels = types.order))]
  combined.pct = combined.pct[, order(factor(colnames(combined.pct), levels = types.order))]
  
  # Ordering for annotations
  types = str_split_fixed(colnames(combined.expr), " ", 2)[,2]
  classes = ExtractString(types, before = '_', after = '-')
  species = str_split_fixed(colnames(combined.expr), " ", 2)[,1]
  row_split = paste0(gene.key$row_annotation, ' markers')
  
  if("one2one" %in% colnames(gene.key)){
    asterisks = ifelse(gene.key$one2one == 1, "*", "")
    rownames(combined.expr) = paste0(rownames(combined.expr), asterisks)
  }
  
  if(rotate){
    # Transpose
    combined.expr = t(combined.expr)
    combined.pct = t(combined.pct)
  } 
  
  if(dotplot){
    min.z.score = -1
    col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c(col.low, col.high))
    dot_col = function(species, expression) circlize::colorRamp2(c(min.z.score, max.z.score), c('white', species_palette[[species]]))(expression)
    cell_fun = function(j, i, x, y, w, h, fill){
      if(color.dot.by.species){
        grid.circle(x=x,
                    y=y,
                    r=unit(combined.pct[i, j]/100 * dot.scale.factor, "cm"),
                    gp=gpar(fill = dot_col(species[[j]], combined.expr[i, j]), col = NA))
      } else {
        grid.circle(x=x,
                    y=y,
                    r=unit(combined.pct[i, j]/100 * dot.scale.factor, "cm"),
                    gp=gpar(fill = col_fun(combined.expr[i, j]), col = NA))
      } 
    }
    
    rect_gp = gpar(type = "none")
    
    # Set this for color legend; edit: setting above
    # col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c('white', 'grey25'))
    
    # Create a dot legend
    lgd.values = c(0+pseudocount, 25+pseudocount, 50+pseudocount, 75+pseudocount, 100+pseudocount)
    lgd.values[lgd.values > (max.pct + pseudocount)] = (max.pct + pseudocount)
    lgd = Legend(labels = seq(0,100, by = 25), title = "Percent\nexpressed",
                 graphics = list(
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[1]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[2]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[3]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[4]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[5]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA))
                 ), title_gp = gpar(fontfamily = "sans", fontface = "plain"))
  } else {
    cell_fun = NULL
    rect_gp = gpar(col = NA)
    lgd = NULL
  }
  
  # Plot
  if(rotate){
    draw(Heatmap(combined.expr,
                 name = "mat",
                 col = col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE, 
                 # row_split = factor(types, levels = unique(types)),
                 # column_split = factor(row_split, levels = unique(row_split)),
                 column_title = NULL,
                 row_title = NULL, 
                 border = TRUE, 
                 cell_fun = cell_fun,
                 rect_gp = rect_gp,
                 # left_annotation = rowAnnotation(
                 #   # `Cell type` = anno_simple(gp = gpar(fill = type_palette[match(unique(types), names(type_palette))], 
                 #                                                                    # col = 'white')),
                 #                                                          #labels = unique(types), 
                 #                                                          #labels_gp = gpar(col = "black", fontsize = 10)),
                 #                                 Species = factor(species, levels = species.use),
                 #                                 # annotation_name_side = "bottom", 
                 #                                 show_annotation_name = FALSE,
                 #                                 show_legend = show_annotation_legend, 
                 #                                 border = TRUE, 
                 #                                 col = list(`Cell type` = type_palette, Species = species_palette), 
                 #                                 annotation_legend_param = list(title_gp = gpar(fontfamily = font.family, fontsize = font.size), border = 'black')), 
                 heatmap_legend_param = list(
                   title = "Scaled\nexpression",
                   border = "black", # Black border around legend
                   legend_gp = gpar(col = "black"),
                   labels_gp = gpar(col = "black", fontfamily = font.family, fontsize = font.size-2), # Black tick labels
                   ticks_gp = gpar(col = "black"), 
                   title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
                 column_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2),  
                 row_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2), 
                 ...), 
         annotation_legend_list = lgd
    )
  } else {
    draw(Heatmap(combined.expr,
                 name = "mat",
                 col = col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE, 
                 row_split = factor(row_split, levels = unique(row_split)),
                 column_split = factor(types, levels = unique(types)),
                 column_title = NULL,
                 row_title = NULL,   
                 border = TRUE, 
                 cell_fun = cell_fun,
                 # cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", combined.pct[i, j]), x, y, gp = gpar(fontsize = 10))},
                 rect_gp = rect_gp,
                 top_annotation = HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = type_palette[match(unique(types), names(type_palette))], 
                                                                                       col = 'white'),
                                                                             labels = unique(types), 
                                                                             labels_gp = gpar(col = "black", fontsize = 10)),
                                                    Species = factor(species, levels = species.use),
                                                    show_annotation_name = FALSE,
                                                    # annotation_name_side = 'none', 
                                                    show_legend = show_annotation_legend, 
                                                    border = TRUE, 
                                                    col = list(`Cell type` = type_palette, Species = species_palette), 
                                                    annotation_legend_param = list(title_gp = gpar(fontfamily = font.family, fontsize = font.size), border = 'black')), 
                 heatmap_legend_param = list(
                   title = "Scaled\nexpression",
                   border = "black",      # Black border around legend
                   legend_gp = gpar(col = "black"),
                   labels_gp = gpar(col = "black", fontfamily = font.family, fontsize = font.size-2), # Black tick labels
                   ticks_gp = gpar(col = "black"), 
                   title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
                 column_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2),  
                 row_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2), 
                 ...), 
         annotation_legend_list = lgd)
  }
  
  # Add horizontal bars separating rows
  if(add.breaks){
    breaks = c(sapply(unique(classes), function(x) which(classes == x)[1])-1, length(classes))
    for(index in breaks) {
      x_coord = index/length(types)
      
      # # Vertical
      # decorate_heatmap_body("Alignment\nscore", row_slice = 1, {
      #   grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
      # })
      
      # Horizontal
      decorate_heatmap_body("mat", row_slice = 1, {
        grid.lines(c(0, 1), c(1-x_coord, 1-x_coord), gp = gpar(lty = 1, lwd = 1))
      })
    }
  }
}

# SAMap genes heatmap!
SAMapHeatmap2 = function(objectList, 
                        gene.key, 
                        type_palette, 
                        species_palette,
                        species.use = c("Chicken", "Zebrafish"), 
                        types.use = c("red", "green", "blue", "UV", "rod"), 
                        min.z.score = -2,
                        max.z.score = 2, 
                        col.low = 'grey', 
                        col.high = "#584B9FFF",
                        rotate = FALSE,
                        dotplot = FALSE, 
                        dot.scale.factor = 0.05,
                        color.dot.by.species = FALSE,
                        pseudocount = 1,
                        max.pct = 100,
                        col_fun = NULL,
                        names.show = NULL,
                        font.family = 'ArialMT',
                        font.size = 12,
                        show_annotation_legend = FALSE,
                        types.order = NULL,
                        ...){
  
  args = list(...)
  
  # Set col_fun
  if(is.null(col_fun)){
    col_fun = circlize::colorRamp2(c(min.z.score, 0, max.z.score), c("blue", "white", "red"))
  }
  
  # Scale data
  scaled.expr = lapply(species.use, function(species) {
    message('scaling for species: ', species)
    # if(species == 'Zebrafish') browser()
    # Edit: do scaling with all the clusters, then subset to the clusters plotted. That way if we have only one cluster, the scaled value isn't NA
    # norm.expr = AverageExpression(subset(objectList[[species]], annotated %in% types.use), group.by = "annotated", slot = "data", assay = 'RNA')$RNA
    norm.expr = AverageExpression(objectList[[species]], group.by = "annotated", slot = "data", assay = 'RNA')$RNA
    scaled.expr = t(scale(t(norm.expr)))
    scaled.expr = scaled.expr[match(gene.key[[species]], rownames(scaled.expr)),]
    
    # Now we subset to types of interest
    scaled.expr = scaled.expr[,intersect(types.use, colnames(scaled.expr)), drop = FALSE]
    colnames(scaled.expr) = paste0(species, ' ', colnames(scaled.expr))
    
    # Truncate values to max.z.score
    scaled.expr[scaled.expr < min.z.score] = min.z.score
    scaled.expr[scaled.expr > max.z.score] = max.z.score
    scaled.expr
  })
  
  # Percentage expression
  pct.expr = lapply(species.use, function(species) {
    message('percentage expressed for species: ', species)
    # if(species == 'Zebrafish') browser()
    pct.expr = PercentageExpressed2(subset(objectList[[species]], annotated %in% types.use), features = setdiff(unique(gene.key[[species]]), ''), group.by = 'annotated')
    pct.expr = pct.expr[match(gene.key[[species]], rownames(pct.expr)),,drop = FALSE]
    colnames(pct.expr) = paste0(species, ' ', colnames(pct.expr))
    pct.expr
  })
  
  # Combine 
  if(is.null(names.show)) {
    gene.names = apply(gene.key, 1, function(x) paste0(x, collapse = "/"))
  } else {
    gene.names = make.unique(gene.key[[names.show]])
  }
  
  combined.expr = do.call(cbind, scaled.expr)
  rownames(combined.expr) = gene.names
  
  combined.pct = do.call(cbind, pct.expr)
  rownames(combined.pct) = gene.names
  
  # Ordering matrix
  if(is.null(types.order)) types.order = as.vector(apply(outer(types.use, species.use, FUN = function(x, y) paste0(y, " ", x)), 1, function(x) x))
  combined.expr = combined.expr[, order(factor(colnames(combined.expr), levels = types.order))]
  combined.pct = combined.pct[, order(factor(colnames(combined.pct), levels = types.order))]
  
  # Ordering for annotations
  types = str_split_fixed(colnames(combined.expr), " ", 2)[,2]
  species = str_split_fixed(colnames(combined.expr), " ", 2)[,1]
  row_split = paste0(gene.key$row_annotation, ' markers')
  
  if("one2one" %in% colnames(gene.key)){
    asterisks = ifelse(gene.key$one2one == 1, "*", "")
    rownames(combined.expr) = paste0(rownames(combined.expr), asterisks)
  }
  
  if(rotate){
    # Transpose
    combined.expr = t(combined.expr)
    combined.pct = t(combined.pct)
  } 
  
  if(dotplot){
    min.z.score = -1
    col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c(col.low, col.high))
    dot_col = function(species, expression) circlize::colorRamp2(c(min.z.score, max.z.score), c('white', species_palette[[species]]))(expression)
    cell_fun = function(j, i, x, y, w, h, fill){
      if(color.dot.by.species){
        grid.circle(x=x,
                  y=y,
                  r=unit(combined.pct[i, j]/100 * dot.scale.factor, "cm"),
                  gp=gpar(fill = dot_col(species[[j]], combined.expr[i, j]), col = NA))
      } else {
        grid.circle(x=x,
                    y=y,
                    r=unit(combined.pct[i, j]/100 * dot.scale.factor, "cm"),
                    gp=gpar(fill = col_fun(combined.expr[i, j]), col = NA))
        } 
    }
    
    rect_gp = gpar(type = "none")
    
    # Set this for color legend; edit: setting above
    # col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c('white', 'grey25'))
    
    # Create a dot legend
    lgd.values = c(0+pseudocount, 25+pseudocount, 50+pseudocount, 75+pseudocount, 100+pseudocount)
    lgd.values[lgd.values > (max.pct + pseudocount)] = (max.pct + pseudocount)
    lgd = Legend(labels = seq(0,100, by = 25), title = "Percent\nexpressed",
                 graphics = list(
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[1]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[2]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[3]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[4]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                   function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[5]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA))
                 ), title_gp = gpar(fontfamily = "sans", fontface = "plain"))
  } else {
    cell_fun = NULL
    rect_gp = gpar(col = NA)
    lgd = NULL
  }
  
  # Plot
  if(rotate){
    draw(Heatmap(combined.expr,
                  # name = "mat",
                  col = col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE, 
                  column_split = factor(row_split, levels = unique(row_split)), 
                  row_split = factor(types, levels = unique(types)),
                  column_title = NULL,
                  row_title = NULL, 
                  border = TRUE, 
                  cell_fun = cell_fun,
                  rect_gp = rect_gp,
                  left_annotation = rowAnnotation(`Cell type` = anno_block(gp = gpar(fill = type_palette[match(unique(types), names(type_palette))], col = 'white'),
                                                                           labels = unique(types), 
                                                                           labels_gp = gpar(col = "black", fontsize = 10)),
                                                  Species = factor(species, levels = species.use),
                                                  # annotation_name_side = "bottom", 
                                                  show_annotation_name = FALSE,
                                                  show_legend = show_annotation_legend, 
                                                  border = TRUE, 
                                                  col = list(`Cell type` = type_palette, Species = species_palette), 
                                                  annotation_legend_param = list(title_gp = gpar(fontfamily = font.family, fontsize = font.size), border = 'black')), 
                 heatmap_legend_param = list(
                   title = "Scaled\nexpression",
                   border = "black",      # Black border around legend
                   legend_gp = gpar(col = "black"),
                   labels_gp = gpar(col = "black", fontfamily = font.family, fontsize = font.size-2), # Black tick labels
                   ticks_gp = gpar(col = "black"), 
                   title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
                 # column_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2),  
                 # row_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2), 
                  ...), 
         annotation_legend_list = lgd
         )
  } else {
    draw(Heatmap(combined.expr,
                  # name = "mat",
                  col = col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE, 
                  row_split = factor(row_split, levels = unique(row_split)),
                  column_split = factor(types, levels = unique(types)),
                  column_title = NULL,
                  row_title = NULL,   
                  border = TRUE, 
                  cell_fun = cell_fun,
                 # cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", combined.pct[i, j]), x, y, gp = gpar(fontsize = 10))},
                  rect_gp = rect_gp,
                  top_annotation = HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = type_palette[match(unique(types), names(type_palette))], col = 'white'),
                                                                              labels = unique(types), 
                                                                              labels_gp = gpar(col = "black", fontsize = 10)),
                                                     Species = factor(species, levels = species.use),
                                                     show_annotation_name = FALSE,
                                                     # annotation_name_side = 'none', 
                                                     show_legend = show_annotation_legend, 
                                                     border = TRUE, 
                                                     col = list(`Cell type` = type_palette, Species = species_palette), 
                                                     annotation_legend_param = list(title_gp = gpar(fontfamily = font.family, fontsize = font.size), border = 'black')), 
                 heatmap_legend_param = list(
                   title = "Scaled\nexpression",
                   border = "black",      # Black border around legend
                   legend_gp = gpar(col = "black"),
                   labels_gp = gpar(col = "black", fontfamily = font.family, fontsize = font.size-2), # Black tick labels
                   ticks_gp = gpar(col = "black"), 
                   title_gp = gpar(fontfamily = font.family, fontsize = font.size)),
                 # column_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2),  
                 # row_title_gp = gpar(fontfamily = font.family, fontsize = font.size+2), 
                  ...), 
         annotation_legend_list = lgd)
  }
}

# SAMap genes heatmap!
SAMapHeatmap = function(objectList, 
                        gene.key, 
                        species.use = c("Chicken", "Zebrafish"), 
                        types.use = c("red", "green", "blue", "UV", "rod"), 
                        max.z.score = 2, 
                        rotate = FALSE,
                        return.data = FALSE,
                        ...){
  
  # Scale data
  scaled.expr = lapply(species.use, function(species) {
    norm.expr = AverageExpression(subset(objectList[[species]], annotated %in% types.use), group.by = "annotated", slot = "data")$RNA
    scaled.expr = t(scale(t(norm.expr)))
    scaled.expr = scaled.expr[match(gene.key[[species]], rownames(scaled.expr)),]
    colnames(scaled.expr) = paste0(colnames(scaled.expr), "_", species)
    
    # Truncate values to max.z.score
    scaled.expr[scaled.expr < -max.z.score] = -max.z.score
    scaled.expr[scaled.expr > max.z.score] = max.z.score
    scaled.expr
  })
  
  # Combine 
  gene.names = do.call(cbind, lapply(scaled.expr, function(x) rownames(x)))
  gene.names = apply(gene.names, 1, function(x) paste0(x, collapse = "/"))
  combined.expr = do.call(cbind, scaled.expr)
  rownames(combined.expr) = gene.names
  
  # Ordering matrix
  types.order = as.vector(apply(outer(types.use, species.use, FUN = function(x, y) paste0(x, "_", y)), 1, function(x) x))
  combined.expr = combined.expr[, order(factor(colnames(combined.expr), levels = types.order))]
  
  # Ordering for annotations
  types = (str_split_fixed(colnames(combined.expr), "_", 2)[,1])
  species = str_split_fixed(colnames(combined.expr), "_", 2)[,2]
  row_split = gene.key$row_annotation
  
  if("one2one" %in% colnames(gene.key)){
    asterisks = ifelse(gene.key$one2one == 1, "*", "")
    rownames(combined.expr) = paste0(rownames(combined.expr), asterisks)
  }
  
  if(return.data) return(combined.expr)
  
  if(rotate){
    print(Heatmap(t(combined.expr),
                  name = "mat",
                  col = colorRamp2(c(-max.z.score, 0, max.z.score), c("blue", "white", "red")),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE, 
                  heatmap_legend_param = list(title = ""),
                  column_split = factor(row_split, levels = unique(row_split)), 
                  border = TRUE, 
                  left_annotation = rowAnnotation(species = anno_image(image.files[match(species, names(image.files))], 
                                                                       border = FALSE, 
                                                                       space = unit(0, "mm")),
                                                  `cell type` = str_split_fixed(colnames(combined.expr), "_", 2)[,1],
                                                  annotation_name_side = "bottom", 
                                                  show_legend = FALSE, 
                                                  border = TRUE, 
                                                  col = list(`cell type` = type_palette
                                                  )), 
                  ...))
    
    # Add horizontal bars separating rows
    types = rev(types)
    breaks = (sapply(unique(types), function(x) which(types == x)[1])-1)[-1]
    for(index in breaks) {
      for(slice in seq_along(unique(row_split))){
        x_coord = index/length(types)
        
        # Horizontal
        decorate_heatmap_body("mat", column_slice = slice, {
          grid.lines(c(0, 1), c(x_coord, x_coord), gp = gpar(lty = 1, lwd = 1))
        })
      }
    }
  } else {
    print(Heatmap(combined.expr,
                  name = "mat",
                  col = colorRamp2(c(-max.z.score, 0, max.z.score), c("blue", "white", "red")),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE, 
                  row_split = factor(row_split, levels = unique(row_split)), 
                  border = TRUE, 
                  heatmap_legend_param = list(title = ""),
                  # rect_gp = gpar(col = "black", lwd = 1),
                  # left_annotation = rowAnnotation(#row = anno_textbox(gene.key$row_annotation, gene.key$row_annotation), 
                  #                     row = anno_block(labels = unique(gene.key$row_annotation), 
                  #                                labels_gp = gpar(col = "white", fontsize = 10))),
                  top_annotation = HeatmapAnnotation(species = anno_image(image.files[match(species, names(image.files))], 
                                                                          border = FALSE, 
                                                                          space = unit(0, "mm")),
                                                     `cell type` = str_split_fixed(colnames(combined.expr), "_", 2)[,1],
                                                     annotation_name_side = "left", 
                                                     show_legend = FALSE, 
                                                     border = TRUE, 
                                                     # gp = gpar(col = "black"), 
                                                     col = list(`cell type` = pr_palette
                                                     )), 
                  ...))
    
    # Add horizontal bars separating rows
    breaks = (sapply(unique(types), function(x) which(types == x)[1])-1)[-1]
    for(index in breaks) {
      for(slice in seq_along(unique(row_split))){
        x_coord = index/length(types)
        
        # Vertical
        decorate_heatmap_body("mat", row_slice = slice, {
          grid.lines(c(x_coord, x_coord), c(0, 1), gp = gpar(lty = 1, lwd = 1))
        })
      }
    }
  }
}

# From Seurat DotPlot function
DarPlot = function(data, dot.scale = 6, scale.min = NA, scale.max = NA){
  plot <- ggplot(data = data, mapping = aes_string(x = "features.plot", y = "id")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = "avg.exp.scaled")) + 
    scale_color_gradient(low = "lightgrey", high = "blue") + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    theme_cowplot() + coord_flip() + RotatedAxis() + theme(axis.title = element_blank())
  
  plot
}

TopN = function(table, group.by, sort.by, n = 10){
  table %>%
    group_by(!!sym(group.by)) %>%
    arrange(desc(!!sym(sort.by))) %>% 
    slice_head(n = n) %>%
    ungroup() -> top
}

TopNDEGs = function(object, group.by = "annotated", n = 10, avg_log2FC_cutoff = 0.25, 
                    p_val_adj_cutoff = 0.05, only.positive = TRUE, 
                    sort.by = 'avg_log2FC', assay = 'RNA', gene.list = FALSE){
  
  if(inherits(object, 'Seurat')){
    # de_table = wilcoxauc(object, group_by = group.by) %>% 
    #   setNames(c('gene', 'cluster', 'avgExpr', 'avg_log2FC', 'stat', 'auc', 'p_val', 'p_val_adj', 'pct_in', 'pct_out'))
    de_table = FindAllMarkersFast(object, group.by = group.by, avg_log2FC_cutoff = avg_log2FC_cutoff, p_val_adj_cutoff = p_val_adj_cutoff, assay = assay)
  } else {
    de_table = object
  }
  
  if(only.positive){
    de_table %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff) %>%
      arrange(desc(!!sym(sort.by))) %>% 
      slice_head(n = n) %>%
      ungroup() -> top
  } else {
    de_table %>%
      group_by(cluster) %>%
      dplyr::filter(abs(avg_log2FC) > avg_log2FC_cutoff & p_val_adj < p_val_adj_cutoff) %>%
      arrange(desc(abs(avg_log2FC))) %>% 
      slice_head(n = n) %>%
      ungroup() -> top
  }
  
  if(is.factor(object@meta.data[[group.by]])){
    # top = top %>% arrange(factor(rlang::sym(c("annotated")), levels = levels(object@meta.data[[group.by]])))
    top$cluster = factor(top$cluster, levels = levels(object@meta.data[[group.by]]))
    top = top[order(top$cluster),]
  }
  
  if(gene.list) return(split(top$gene, top$cluster))
  
  top %>% as.data.frame 
}

CelltypeProportionBarplot = function(object, 
                                     x = "annotated", 
                                     y = "animal", 
                                     group.by = NULL,
                                     return.table = FALSE, 
                                     show.all = FALSE, 
                                     normalize = TRUE, 
                                     facet = TRUE, 
                                     outline = 'grey20', 
                                     subset.to = NULL,
                                     ...){
  
  if(inherits(object, "Seurat")) {
    data = object@meta.data
  } else {
    data = object
  }
  
  # Remove unused factors
  if(is.factor(data[[x]])) data[[x]] = factor(factor(data[[x]], levels = (levels(data[[x]]))))
  
  tabulation = table(data[[y]], data[[x]])
  normalized.tabulation = tabulation / rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation) %>% setNames(c("Sample", "Type", "Proportion"))
  
  # Add group metadata
  if(!is.null(group.by)) {
    metadata = Metadata(object, y, group.by)
    melted$group = metadata[match(melted$Sample, metadata[[y]]),group.by]
  }
  
  if(!is.null(subset.to)) melted = subset(melted, Type %in% subset.to)
  
  
  if(show.all){
    if(facet){
      
      p = ggplot(melted, aes(y = Proportion, x = Sample, fill = Type))+
        geom_bar(stat = "identity", position="dodge", color = outline) + 
        theme_bw() + 
        facet_wrap(~ Type, ...)+
        scale_fill_discrete(name = "Sample")+
        ylab("Proportion") 
    } else {
      p = ggplot(melted, aes(y = Proportion, x = Type, fill = Sample))+
        geom_bar(stat = "identity", position="dodge") + 
        theme_bw() + 
        scale_fill_discrete(name = "Sample")+
        ylab("Proportion") 
    }
    
  } else if(!is.null(group.by)){
    p = ggbarplot(melted, 
                  y = "Proportion", 
                  x = 'group', 
                  fill = 'group', 
                  add = c("mean_se", "jitter"), 
                  facet.by = 'Type',
                  ...) + 
      NoLegend() +
      # facet_wrap(~ Type, ...)+
      scale_y_continuous(expand = expansion(mult = c(0, .1)))
  } else {
    p = ggbarplot(melted, 
                  y = "Proportion", 
                  x = "Type", 
                  fill = "Type", 
                  add = c("mean_se", "jitter"), 
                  ...) + 
      NoLegend() +
      scale_y_continuous(expand = expansion(mult = c(0, .1)))
    # rremove("xlab")
  }
  if(return.table) {
    if(normalize) {
      normalized.tabulation
    } else {
      tabulation
    }
  } else {
    p
  }
}

PercentageExpressed2 = function(object, features, group.by = NULL){
  
  # object = SubsetSeuratGenes(object, features = features)
  features.use = intersect(features, rownames(object))
  
  if(!is.null(group.by)){ 
    countDat = object@assays$RNA@counts[features.use, , drop = FALSE]
    groups = object@meta.data[[group.by]]
    count.list = suppressWarnings(split.default(as.data.frame(countDat), groups) |> lapply(as.matrix))
    pct.exp = as.data.frame(do.call(cbind, lapply(count.list, function(datExpr){
      pct.exp = apply(datExpr, 1, function(row) length(row[row > 0])/length(row) * 100)
    })))
    
    # Sort by levels of group.by
    if(is.factor(object@meta.data[[group.by]])) {
      # refactor in case some were subsetted out
      object@meta.data[[group.by]] = factor(object@meta.data[[group.by]])
      pct.exp = pct.exp[,levels(object@meta.data[[group.by]])]
    }
    
    pct.exp = pct.exp[match(features, rownames(pct.exp)),,drop = FALSE]
    rownames(pct.exp) = features
    
  } else {
    datExpr = object@assays$RNA@counts[features.use, , drop = FALSE]
    pct.exp = apply(datExpr, 1, function(row) length(row[row > 0])/length(row) * 100)
    
    pct.exp = pct.exp[match(features, names(pct.exp))]
    names(pct.exp) = features
  }
  
  return(pct.exp)
}

PercentageExpressed = function(object, features, group.by = NULL){
  
  # object = SubsetSeuratGenes(object, features = features)
  
  if(!is.null(group.by)){
    obj.list = SplitObject(object, split.by = group.by)
    pct.exp = as.data.frame(do.call(cbind, lapply(obj.list, function(obj){
      datExpr = obj@assays$RNA@counts[features, , drop = FALSE]
      pct.exp = apply(datExpr, 1, function(row) length(row[row > 0])/length(row) * 100)
    })))
    
    # Sort by levels of group.by
    if(is.factor(object@meta.data[[group.by]])) {
      # refactor in case some were subsetted out
      object@meta.data[[group.by]] = factor(object@meta.data[[group.by]])
      pct.exp = pct.exp[,levels(object@meta.data[[group.by]])]
    }
  } else {
    datExpr = object@assays$RNA@counts[features, , drop = FALSE]
    pct.exp = apply(datExpr, 1, function(row) length(row[row > 0])/length(row) * 100)
  }
  
  return(pct.exp)
  
  # tabulation = (object@assays$RNA@counts[feature,] > 0) %>% table
  # pct.exp = as.numeric((tabulation[["TRUE"]]/sum(tabulation)))*100
  # return(pct.exp)
}

# From https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

BlastOrthologyTable = function(file1 = "maps/zech/ch_to_ze.txt", 
                               file2 = "maps/zech/ze_to_ch.txt", 
                               keyfile1 = "Keys/Chicken_key.csv",
                               keyfile2 = "Keys/Zebrafish_key.csv", 
                               one2one = TRUE, 
                               return.both = FALSE, 
                               as.edge.list = TRUE){
  
  # Blast columns
  blast_cols <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
  )
  
  message('Reading files...')
  orthology_file1 = fread(file1) %>% as.data.frame() %>% setNames(blast_cols)
  orthology_file2 = fread(file2) %>% as.data.frame() %>% setNames(blast_cols)
  key1 = fread(keyfile1, header = TRUE) %>% as.data.frame()
  key2 = fread(keyfile2, header = TRUE) %>% as.data.frame()
  
  # Convert transcripts to genes
  message('Converting transcripts to genes...')
  orthology_file1$qseqid = key1$symbol[match(orthology_file1$qseqid, key1[, 2])]
  orthology_file1$sseqid = key2$symbol[match(orthology_file1$sseqid, key2[, 2])]
  orthology_file2$qseqid = key2$symbol[match(orthology_file2$qseqid, key2[, 2])]
  orthology_file2$sseqid = key1$symbol[match(orthology_file2$sseqid, key1[, 2])]
  
  # Only top hit per gene combination (e.g. collapse isoforms)
  if (one2one) {
    
    # Get best 1 to 2 blast
    best1 = orthology_file1 %>%
      filter(evalue < 1e-6) %>%
      dplyr::select(qseqid, sseqid, bitscore) %>%
      group_by(qseqid) %>%
      filter(bitscore == max(bitscore, na.rm = TRUE)) %>%
      unique()
    
    best2 = orthology_file2 %>%
      filter(evalue < 1e-6) %>%
      dplyr::select(qseqid, sseqid, bitscore) %>%
      group_by(qseqid) %>%
      filter(bitscore == max(bitscore, na.rm = TRUE)) %>%
      unique()
    
  } else {
    
    best1 = orthology_file1 %>%
      filter(evalue < 1e-6) %>%
      dplyr::select(qseqid, sseqid, bitscore) %>% 
    
    best1 <- orthology_file1 %>%
      # unique() %>% # remove duplicate rows, which happens when isoforms have same sequence
      rowwise() %>% # 
      mutate(pair = paste(sort(c(qseqid, sseqid)), collapse = "-")) %>% # make gene pair pair column
      ungroup() %>% # ungroup if grouped
      group_by(pair) %>% # group by gene pair
      slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>% # find the best bit score for each gene pair
      mutate(across(everything(), ~replace_na(.x, 0))) # replace nas with zero
    
    best2 <- orthology_file2 %>%
      unique() %>%
      rowwise() %>%
      mutate(pair = paste(sort(c(qseqid, sseqid)), collapse = "-")) %>%
      ungroup() %>%
      group_by(pair) %>%
      slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
      mutate(across(everything(), ~replace_na(.x, 0)))
  }
  
  # Binary masks
  matrix1 = reshape2::acast(best1, qseqid ~ sseqid, value.var = "bitscore")
  matrix2 = reshape2::acast(best2, qseqid ~ sseqid, value.var = "bitscore")
  
  if (one2one) {
    message('Making binary...')
    matrix1[!is.na(matrix1)] = 1
    matrix1[is.na(matrix1)] = 0
    matrix2[!is.na(matrix2)] = 1
    matrix2[is.na(matrix2)] = 0
  }
  
  # Subset to common genes
  message('Subsetting to common genes...')
  gene1 = intersect(rownames(matrix1), colnames(matrix2))
  gene2 = intersect(rownames(matrix2), colnames(matrix1))
  final1 = matrix1[gene1, gene2]
  final2 = t(matrix2[gene2, gene1])
  
  if (one2one) {
    
    message('Finding bidirectional blast hits...')
    
    # Find bidirectional blast hits
    A = final1 * final2
    
    # Remove genes with no ortholog
    A = A[rowSums(A) > 0, ]
    A = A[, colSums(A) > 0]
    
    # Underscores to dashes for compatibility with Seurat
    colnames(A) = gsub("_", "-", colnames(A))
    rownames(A) = gsub("_", "-", rownames(A))
    
    # Get ambiguous mappings
    dupl_rows = apply(A, 1, function(x) sum(x) > 1)
    dupl_cols = apply(A, 2, function(x) sum(x) > 1)
    
    # Remove ambiguous orthologs
    new = A[!dupl_rows, !dupl_cols]
    
    # Remove rows/columns with zeroes
    zero_rows = apply(new, 1, function(x) sum(x) == 0)
    zero_cols = apply(new, 2, function(x) sum(x) == 0)
    new2 = new[!zero_rows, !zero_cols]
    
  } else if (return.both) {
    
    message('Returning both graphs...')
    if (as.edge.list) {
      edges1 <- as.data.frame(as.table(final1)) %>%
        filter(Freq != 0) %>%
        dplyr::rename(gene1 = Var1, gene2 = Var2, bit_score = Freq) %>%
        mutate(gene1 = paste0('gene1:', gene1),
               gene2 = paste0('gene2:', gene2))
      edges2 <- as.data.frame(as.table(final2)) %>%
        filter(Freq != 0) %>%
        dplyr::rename(gene1 = Var1, gene2 = Var2, bit_score = Freq) %>%
        mutate(gene1 = paste0('gene1:', gene1),
               gene2 = paste0('gene2:', gene2))
      return(edges1, edges2)
    } else {
      return(list(final1, final2))
    }
    
  } else {
    
    message('Finding element-wise geometric mean...')
    
    new = sqrt(final1 * t(final2))
    new2 = new
    
    if (as.edge.list) {
      edges <- as.data.frame(as.table(new2)) %>%
        filter(Freq != 0) %>%
        dplyr::rename(gene1 = Var1, gene2 = Var2, bit_score = Freq) %>%
        mutate(gene1 = paste0('gene1:', gene1),
               gene2 = paste0('gene2:', gene2))
      return(edges)
    }
  }
  
  return(new2)
}

# With italicized gene symbols
DotPlot3 = function(object, 
                    coord.flip = TRUE, 
                    max.pct = 100,
                    binarization = NULL, 
                    col.high = '#584B9FFF', 
                    col.low = 'lightgrey', 
                    show = NULL, # which clusters to show (scaling is done on all the clusters, not just these)
                    max.size = 6,
                    ...){
  args = list(...)
  args$object = object
  args$col.min = -1
  args$col.max = 2
  
  if(coord.flip) {
    message('reversing order of genes for y-axis...')
    args$features = rev(args$features) # reverse gene order to make it top to bottom
  } 
  
  # Remove features that are not present
  args$features = intersect(args$features, rownames(object))
  stopifnot(length(args$features) > 0)
  
  # Generate DotPlot data
  dot_data <- do.call(DotPlot, args)$data
  dot_data$key = paste0(dot_data$features.plot, '-', dot_data$id)
  
  # Scale the 'pct.exp' column so the max is 50%
  dot_data$pct.exp <- pmin(dot_data$pct.exp, max.pct) # Cap at 50%
  
  if(!coord.flip){
    dot_data$id = factor(dot_data$id, levels = rev(levels(dot_data$id))) # reverse type order if types are on the y-axis
  }
  
  if(!is.null(binarization)){
    bin.melt = reshape2::melt(binarization)
    bin.melt$key = paste0(bin.melt$Var1, '-', bin.melt$Var2)
    dot_data = dot_data[match(bin.melt$key, dot_data$key),]
    dot_data$heat_value = factor(bin.melt$value, levels = c(0,1))
  }
  
  if(is.null(show)) show = unique(dot_data$id)
  
  # Plot with modified data
  ggplot(subset(dot_data, id %in% show), aes(x = features.plot, y = id)) +
    {if(!is.null(binarization)) geom_tile(aes(fill = heat_value))}+
    {if(!is.null(binarization)) scale_fill_manual(values = c("0" = "white", "1" = "lightgrey"))}+
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    theme_cowplot()+
    {if(coord.flip) coord_flip()}+
    {if(coord.flip) theme(axis.title = element_blank(), axis.text.y = element_text(face = 'italic'))}+ 
    {if(!coord.flip) theme(axis.title = element_blank(), axis.text.x = element_text(face = 'italic'))}+
    RotatedAxis() + 
    scale_color_gradient(name = 'Scaled\nexpression', 
                         low = col.low, 
                         high = col.high, 
                         limits = c(-1, 2), 
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    scale_radius(name = 'Percent\nexpressed', limits = c(0, max.pct), range = c(0, max.size))+
    ArialFont()
  
}

DotPlot2 = function(object, coord.flip = TRUE, max.pct = 100, ...){
  args = list(...)
  if(coord.flip) args$features = rev(args$features) # reverse gene order to make it top to bottom
  
  # Generate DotPlot data
  dot_data <- DotPlot(object, col.min = -1, col.max = 2, ...)$data
  
  # Scale the 'pct.exp' column so the max is 50%
  dot_data$pct.exp <- pmin(dot_data$pct.exp, max.pct) # Cap at 50%
  
  # Plot with modified data
  ggplot(dot_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
    geom_point() +
    theme_cowplot()+
    # scale_size_continuous(range = c(0, 6))  # Adjust dot size scale as needed
  # DotPlot(object, col.min = -1, col.max = 2, ...) + 
    theme(axis.title = element_blank(), axis.text.y = element_text(face = 'italic'))+ # Set this before coord flipping
    {if(coord.flip) coord_flip()}+ 
    RotatedAxis() + 
    guides(color = guide_colorbar(title = "Scaled\nexpression"), 
           size = guide_legend(title = "Percent\nexpressed"))+
    scale_color_gradient(name = 'Averageexpression', low = "lightgrey", high = "#584B9FFF", limits = c(-1, 2))+
    scale_radius(name = 'Percentexpressed', limits = c(0, max.pct), range = c(0, 6), )
}

SeuratToH5ad2 = function(object, filepath, counts.only = TRUE){
  
  # Turn to character, factor leads to some issues downstream
  if('annotated' %in% colnames(object@meta.data)) object$annotated = as.character(object$annotated)
  
  if(counts.only){
    new_object = CreateSeuratObject(object@assays$RNA@counts)
    new_object = TransferMetadata(object, new_object)
    object = new_object
  }
  
  tempfile = gsub("h5ad", "h5Seurat", filepath)
  SeuratDisk::SaveH5Seurat(object, filename = tempfile, overwrite = TRUE)
  SeuratDisk::Convert(tempfile, dest = "h5ad", overwrite = TRUE)
  file.remove(tempfile)
  
  return(object)
}

SeuratToH5ad = function(object, filepath, genes.remove = NULL, types.remove = NULL, downsample = 100, seed = 12345){
  
  # Unfactor as this causes issues downstream
  object$annotated = as.character(object$annotated)
  
  if(!is.null(types.remove)){
    object = subset(object, annotated %in% types.remove, invert = TRUE)
  }
  if(!is.null(downsample)) object = DownsampleSeurat(object, group.by = "annotated", size = downsample, seed = seed)
  message("Some cell ids: ", paste0(head(Cells(object)), collapse = ", "))
  annotated = object$annotated
  animal = object$animal
  object = CreateSeuratObject(object@assays$RNA@counts[!rownames(object@assays$RNA@counts) %in% genes.remove,])
  object$annotated = annotated
  object$animal = animal
  print(table(object$annotated))
  
  tempfile = gsub("h5ad", "h5Seurat", filepath)
  SeuratDisk::SaveH5Seurat(object, filename = tempfile, overwrite = TRUE)
  SeuratDisk::Convert(tempfile, dest = "h5ad", overwrite = TRUE)
  file.remove(tempfile)
  
  return(object)
}

ExtractString = function(vector, before = NULL, after = NULL){
  
  library(stringr)
  
  if(is.null(after)){
    return(str_split_fixed(vector, before, 2)[,2])
  } else if(is.null(before)){
    return(str_split_fixed(vector, after, 2)[,1])
  } else {
    return(str_split_fixed(str_split_fixed(vector, before, 2)[,2], after, 2)[,1])
  }
}

UpdateCellClass = function(object, annotation = major_annotation){
  # Convert to major cell class
  object$cell_class = as.character(object$cell_class)
  object$cell_class[object$cell_class == "GabaAC" | object$cell_class == "GlyAC"] = "AC"
  object$cell_class[object$cell_class == "BP"] = "BC"
  
  # Everything else goes to other
  object$cell_class[!object$cell_class %in% Annotation(annotation)] = "Other"
  # object$cell_class[object$cell_class == "MicroG"] = "Other"
  object$cell_class = factor(object$cell_class, unique(Annotation(annotation)))
  
  return(object)
}

ProportionBarplot = function(data, x, y, return.table = FALSE){
  tabulation = table((data[[x]]), (data[[y]]))
  normalized.tabulation = tabulation/rowSums(tabulation)
  melted = reshape2::melt(normalized.tabulation)
  p = ggplot(melted, aes(y = value, x = Var2, fill = Var1))+
    geom_bar(stat = "identity", position="dodge") + 
    rremove("xlab") + ylab("Proportion")
  if(return.table) {
    normalized.tabulation
  } else {
    p + theme_dario()
  }
}

PrettyBoxplot2 = function(data, x, y, color = NULL, alpha = 0.5, width = 0.5, size = 1, jitter = TRUE){
  ggplot(data, aes(x = eval(parse(text = x)), y = eval(parse(text = y)), color = eval(parse(text = color))))+
    geom_boxplot(linetype = "dashed", outlier.shape = NA, width = width) +
    stat_boxplot(aes(ymin = after_stat(lower), ymax = after_stat(upper)), outlier.shape = NA, width = width) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = width/2) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = width/2) +
    {if(jitter) geom_jitter(width = width/2, alpha = alpha, size = size)}+
    theme_cowplot()+
    labs(x = x, y = y)+
    theme(legend.position = "none", axis.title.x = element_blank())
}

PrettyBoxplot = function(data){
  ggplot(reshape2::melt(data), aes(x = variable, y = value, color = variable))+
    geom_boxplot(linetype = "dashed", outlier.shape = NA) +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5) +
    theme_cowplot()+
    theme(legend.position = "none", axis.title.x = element_blank())
}

ClusterDendrogram = function(object, assay = "RNA", group.by = "seurat_clusters"){
  Idents(object) = group.by
  DefaultAssay(object) = assay
  # object = FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  object = BuildClusterTree(object = object)
  return(object)
}

ComputeM1 = function(original_clusters, permuted_clusters){
  
  tableMatrix = table(original_clusters, permuted_clusters)
  
  # Proportion of type in each cluster (row-normalized)
  rnMatrix = tableMatrix/rowSums(tableMatrix)
  
  m1.scores = sapply(seq_along(unique(original_clusters)), function(index){
    proportions = rnMatrix[index,]
    score = 1-ShannonEntropy(proportions)/log2(length(proportions))
    return(score)
  })
  
  return(m1.scores)
}

ComputeM2 = function(original_clusters, permuted_clusters){
  tableMatrix = table(original_clusters, permuted_clusters)
  
  # Proportion of cluster from each type (column-normalized)
  cnMatrix = t(t(tableMatrix)/colSums(tableMatrix))
  
  m2.scores = sapply(seq_along(unique(original_clusters)), function(index){
    raw_score = sum((tableMatrix[index,] * (1-cnMatrix[index,])))/sum(tableMatrix[index,])
    score = 1-raw_score
    return(score)
  })
  
  return(m2.scores)
}

MatchClusters = function(reference, target, p.adj.threshold = 1e-10, overlap.threshold = 30, 
                         jaccard.threshold = 0, return.key = FALSE, bidirectional = FALSE, 
                         keep.missing = FALSE){
  
  # Overlap statistics
  overlap.stats = OverlapStatistics(table(reference, target))
  overlap.signif = subset(overlap.stats, pval < p.adj.threshold & overlap > overlap.threshold & jaccard > jaccard.threshold)
  
  if(bidirectional){
    overlap.signif.one2one = BidirectionalBestHits(table(reference, target))
    colnames(overlap.signif.one2one) = c('ident1', 'ident2')
  } else {
    # Throw out multi-mapping clusters
    message('Found the following multimapping clusters')
    print(subset(overlap.signif, ident1 %in% getDuplicates(overlap.signif$ident1)) %>% arrange(ident1))
    overlap.signif.one2one = subset(overlap.signif, !ident1 %in% getDuplicates(overlap.signif$ident1))
  }
  
  if(return.key) return(overlap.signif.one2one)
  
  # Assign new names
  transferred = convert_values(reference, key = overlap.signif.one2one[,c('ident1', 'ident2')] %>% 
                                 setNames(c('old.names', 'new.names')), 
                               keep.missing = keep.missing)
  return(transferred)
}

# Given an OrthoType cluster, it finds the corresponding species clusters that match it given the thresholds 
CorrespondingCluster = function(overlap.list, cluster, log.p = 10, overlap = 30, jaccard = 0){
  if(inherits(overlap.list, 'list')){
    filtered.list = lapply(overlap.list, function(table) table[table$log.p >= log.p & 
                                                                 table$overlap >= overlap & 
                                                                 table$jaccard >= jaccard,])
    results = do.call(rbind, lapply(filtered.list, function(x) subset(x, ident1 == cluster)))
    
  } else if(inherits(overlap.list, 'data.frame')) {
    table = overlap.list
    filtered.list = table[table$log.p >= log.p & table$overlap >= overlap & table$jaccard >= jaccard,]
    results = subset(filtered.list, ident1 == cluster)
  } else {
    stop('should be a data.frame or list')
  }
  
  rownames(results) = NULL
  return(results)
}

OverlapStatistics = function(table){
  pval = HyperMatrix(table, log.p = FALSE)
  log.p = HyperMatrix(table, log.p = TRUE)
  jaccard = JSMatrix(table)
  row.norm.prop = table/rowSums(table)
  col.norm.prop = t(t(table)/colSums(table))
  
  # Summary table
  data=reshape2::melt(table)
  data$pval=reshape2::melt(pval)$value
  data$padj = p.adjust(data$pval, method = 'fdr')
  data$log.p=reshape2::melt(log.p)$value
  data$log.padj = -log10(data$padj)
  data$jaccard=reshape2::melt(jaccard)$value
  data$row.norm.prop = reshape2::melt(row.norm.prop)$value
  data$col.norm.prop = reshape2::melt(col.norm.prop)$value
  colnames(data)=c("ident1", "ident2", "overlap", 'pval', 'padj', 'log.p', "log.padj", "jaccard", 'row.norm.prop', 'col.norm.prop')
  data
}

HyperMatrix = function(table, log.p = TRUE){
  pval = table
  for(column in 1:ncol(table)){
    for(row in 1:nrow(table)){
      pval[row,column]=phyper(table[row,column]-1,
                              sum(table[row,]),
                              sum(rowSums(table))-sum(table[row,]),
                              sum(table[,column]),
                              lower.tail=FALSE, 
                              log.p=FALSE)
    }
  }
  
  # Return negative log10 p-value
  if(log.p) pval = -log10(pval)
  
  return(pval)
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

TitlePlot = function(plot, title, hjust = 0.5){
  plot + ggtitle(title) + theme(plot.title = element_text(hjust = hjust, face = 'plain'))
}

TransferMetadata = function(from, to){
  nCount_RNA = to$nCount_RNA
  nFeature_RNA = to$nFeature_RNA
  to@meta.data = from@meta.data
  to$nCount_RNA = nCount_RNA
  to$nFeature_RNA = nFeature_RNA
  
  return(to)
}

UpdateEnrichment = function(object, species, enrichment.order = c("CD90+", "CD73-", "NEUN+", "NEUN-", "CHX10+", "NEUN-\nCHX10-", "NEUN-\nCHX10+", "NEUN+\nCHX10-", "CHX10-\nCD73-\nCD133-", "NONE")){
  
  # Read in species metadata
  species_metadata = read.csv("../../Species_Objects/sample_metadata2.csv")
  
  if(species == "Human") {
    object$region = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,1])
    object$enrichment = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,2])
    object$animal = gsub("A", "", gsub("B", "", gsub("C", "", toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,3]))))
    object$donor = object$animal
  } else if(species == "Goldfish"){
    object$enrichment = "NONE"
    object$animal = object$orig.ident
  } else if(species == 'Macaque'){
    object$enrichment = toupper(str_split_fixed(species_metadata$description[match(object$orig.ident, species_metadata$channel)], "_", 3)[,2])
    object$animal = substr(toupper(str_split_fixed(species_metadata$description[match(object$orig.ident, species_metadata$channel)], "_", 3)[,3]), 0, 2)
    object$region = toupper(str_split_fixed(species_metadata$description[match(object$orig.ident, species_metadata$channel)], "_", 3)[,1])
  } else if(species == "Macaque.old") {
    # Region 
    object$region = "Periphery"
    object$region[grepl("Fovea", object$orig.ident, fixed = TRUE)] = "Fovea"
    
    # Animal
    object$animal = substr(object$tag, 0, 2)
    object$animal[startsWith(object$tag, "M1Per")] = "M5"
    object$animal[startsWith(object$tag, "M2Per")] = "M6"
    object$animal[startsWith(object$tag, "M3Per")] = "M7"
    
    # Enrichment type
    object$enrichment = "NONE"
    object$enrichment[grepl("CD90|Mixed", object$tag)] = "CD90+"
    object$enrichment[grepl("CD73-", object$tag)] = "CD73-"
    
    # object$batch = paste0(object$region, "_", object$animal, "_", object$enrichment)
    object$batch = paste0(object$region, "_", object$animal)
  } else if(species == "Marmoset") {
    object$region = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,1])
    object$enrichment = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,2])
    object$animal = gsub("A", "", gsub("B", "", gsub("C", "", toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 3)[,3]))))
  } else {
    object$enrichment = toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 2)[,1])
    object$animal = gsub("A", "", gsub("B", "", gsub("C", "", toupper(str_split_fixed(species_metadata$description[match(object$orig.file, species_metadata$channel)], "_", 2)[,2]))))
  }
  
  # Make sure they are all uniform
  object$enrichment[toupper(object$enrichment) == "CD90"] = "CD90+"
  object$enrichment[toupper(object$enrichment) == "CHX10"] = "CHX10+"
  object$enrichment[toupper(object$enrichment) == "CHX10+"] = "CHX10+"
  object$enrichment[toupper(object$enrichment) == "NEUN"] = "NEUN+"
  object$enrichment[toupper(object$enrichment) == "NEUN+"] = "NEUN+"
  object$enrichment[toupper(object$enrichment) == "NEUN-"] = "NEUN-"
  object$enrichment[toupper(object$enrichment) == "CD73"] = "CD73-"
  object$enrichment[toupper(object$enrichment) == "ALL"] = "NONE"
  object$enrichment[toupper(object$enrichment) == "ISL2B"] = "ISL2B+"
  object$enrichment[toupper(object$enrichment) == "NEUN-CHX10-"] = "NEUN-\nCHX10-"
  object$enrichment[toupper(object$enrichment) == "NEUN-CHX10+"] = "NEUN-\nCHX10+"
  object$enrichment[toupper(object$enrichment) == "NEUN+CHX10-"] = "NEUN+\nCHX10-"
  
  # Order
  object$enrichment = factor(object$enrichment, levels = enrichment.order)
  
  return(object)
}

ReadAmacrineData2 = function(){
  objectList = lapply(FINALFILES, function(file) {
    print(file)
    readRDS(file)
  })
  
  objectList$Mouse$animal = objectList$Mouse$batch
  objectList
}

ReadAmacrineData = function(speciesList, mc.cores = 1, yan = FALSE){
  objectList = mclapply(speciesList, function(species) {
    message('Working on ', species)
    if(species == "Mouse") {
      if(yan){
        YanAC = readRDS("../../Species_Reference/YanAC_v4.rds")
        
        # Yan et al. cell depletion strategy
        # BCs and Müller glia - GFP, rods - CD73, cones - CD133
        # objectList$Mouse$enrichment = "NONE"
        YanAC$enrichment = "CHX10-\nCD73-\nCD133-"
        YanAC$seurat_clusters = YanAC$cluster_no
        
        # Assign literature types
        yan_otherdata = fread("../../reference_files/MouseAC_other_meta.txt", sep = " ", fill = TRUE)
        YanAC$lit_type = yan_otherdata$Notes[match(YanAC$cluster_no, yan_otherdata$Cluster)]
        YanAC$lit_type[YanAC$cluster_no == 10] = "nGnG-2_CCK"
        YanAC$lit_type[YanAC$cluster_no == 63] = "PENK_SST"
        YanAC$lit_type[YanAC$lit_type == ""] = NA
        
        # Original classification
        YanAC$orig.classification = ifelse(YanAC$classification == "GABAergic", "GABA", 
                                           ifelse(YanAC$classification == "Glycinergic", "Gly", 
                                                  ifelse(YanAC$classification == "Both", "Both", "nGnG")))
        
        return(YanAC)
      } else {
        return(readRDS("../../Species_Reference/MouseACref_v4.rds"))
      }
    } else {
      readRDS(paste0("../../Species_Objects/", species, "AC_v5.rds"))
    }
  }, mc.cores = mc.cores)
  names(objectList) = speciesList
  
  return(objectList)
}

ScaledExpression = function(object, scale.within = NULL, ...){
  avg.expr = AverageExpression(object, ...)$RNA
  
  
  scaled.expr = t(scale(t(avg.expr)))
}

LogAvgExpr = function(object, ...){
  args = list(...)
  return(as.data.frame(log1p(AverageExpression(object, ...)[[args$assay]])))
}

GetAnnotationColors = function(vector, annotation){
  return(Colors(annotation)[match(vector, Annotation(annotation))])
}

FormatTitle = function(){
  return(theme(plot.title = element_text(hjust = 0.5, face = "plain")))
}

theme_cowplot2 = function (font_size = 14, font_family = "", line_size = 0.5, 
                           rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14) 
{
  half_line <- font_size/2
  small_size <- rel_small * font_size
  theme_grey(base_size = font_size, base_family = font_family) %+replace% 
    theme(line = element_line(color = "black", size = line_size, 
                              linetype = 1, lineend = "butt"), rect = element_rect(fill = NA, 
                                                                                   color = NA, size = line_size, linetype = 1), text = element_text(family = font_family, 
                                                                                                                                                    face = "plain", color = "black", size = font_size, 
                                                                                                                                                    hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9, 
                                                                                                                                                    margin = margin(), debug = FALSE), axis.line = element_line(color = "black", 
                                                                                                                                                                                                                size = line_size, lineend = "square"), axis.line.x = NULL, 
          axis.line.y = NULL, axis.text = element_text(color = "black", 
                                                       size = small_size), axis.text.x = element_text(margin = margin(t = small_size/4), 
                                                                                                      vjust = 1), axis.text.x.top = element_text(margin = margin(b = small_size/4), 
                                                                                                                                                 vjust = 0), axis.text.y = element_text(margin = margin(r = small_size/4), 
                                                                                                                                                                                        hjust = 1), axis.text.y.right = element_text(margin = margin(l = small_size/4), 
                                                                                                                                                                                                                                     hjust = 0), axis.ticks = element_line(color = "black", 
                                                                                                                                                                                                                                                                           size = line_size), axis.ticks.length = unit(half_line/2, 
                                                                                                                                                                                                                                                                                                                       "pt"), axis.title.x = element_text(margin = margin(t = half_line/2), 
                                                                                                                                                                                                                                                                                                                                                          vjust = 1), axis.title.x.top = element_text(margin = margin(b = half_line/2), 
                                                                                                                                                                                                                                                                                                                                                                                                      vjust = 0), axis.title.y = element_text(angle = 90, 
                                                                                                                                                                                                                                                                                                                                                                                                                                              margin = margin(r = half_line/2), vjust = 1), 
          axis.title.y.right = element_text(angle = -90, margin = margin(l = half_line/2), 
                                            vjust = 0), legend.background = element_blank(), 
          legend.spacing = unit(font_size, "pt"), legend.spacing.x = NULL, 
          legend.spacing.y = NULL, legend.margin = margin(0, 
                                                          0, 0, 0), legend.key = element_blank(), legend.key.size = unit(1.1 * 
                                                                                                                           font_size, "pt"), legend.key.height = NULL, legend.key.width = NULL, 
          legend.text = element_text(size = rel(rel_small)), 
          legend.text.align = NULL, legend.title = element_text(hjust = 0), 
          legend.title.align = NULL, legend.position = "right", 
          legend.direction = NULL, legend.justification = c("left", 
                                                            "center"), legend.box = NULL, legend.box.margin = margin(0, 
                                                                                                                     0, 0, 0), legend.box.background = element_blank(), 
          legend.box.spacing = unit(font_size, "pt"), panel.background = element_blank(), 
          panel.border = element_blank(), panel.grid = element_blank(), 
          panel.grid.major = NULL, panel.grid.minor = NULL, 
          panel.grid.major.x = NULL, panel.grid.major.y = NULL, 
          panel.grid.minor.x = NULL, panel.grid.minor.y = NULL, 
          panel.spacing = unit(half_line, "pt"), panel.spacing.x = NULL, 
          panel.spacing.y = NULL, panel.ontop = FALSE, strip.background = element_rect(fill = "grey80"), 
          strip.text = element_text(size = rel(rel_small), 
                                    margin = margin(half_line/2, half_line/2, half_line/2, 
                                                    half_line/2)), strip.text.x = NULL, strip.text.y = element_text(angle = -90), 
          strip.placement = "inside", strip.placement.x = NULL, 
          strip.placement.y = NULL, strip.switch.pad.grid = unit(half_line/2, 
                                                                 "pt"), strip.switch.pad.wrap = unit(half_line/2, 
                                                                                                     "pt"), plot.background = element_blank(), plot.title = element_text(
                                                                                                       size = rel(rel_large), hjust = 0.5, vjust = 1, 
                                                                                                       margin = margin(b = half_line)), plot.subtitle = element_text(size = rel(rel_small), 
                                                                                                                                                                     hjust = 0, vjust = 1, margin = margin(b = half_line)), 
          plot.caption = element_text(size = rel(rel_tiny), 
                                      hjust = 1, vjust = 1, margin = margin(t = half_line)), 
          plot.tag = element_text(face = "bold", hjust = 0, 
                                  vjust = 0.7), plot.tag.position = c(0, 1), plot.margin = margin(half_line, 
                                                                                                  half_line, half_line, half_line), complete = TRUE)
}

SmartSaveRDS = function(object, path, overwrite = FALSE){
  if(!basename(path) %in% list.files(dirname(path))){
    message("Writing file to ", path)
    saveRDS(object, path)
  } else if(overwrite){
    message("Overwriting file ", path)
    saveRDS(object, path)
  }
  else {
    warning("File already exists, skipped writing file...")
  }
}

NoRotatedXAxis = function(angle = 0, hjust = 0.5, vjust = 1){
  return(theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust)))
}

RotatedXAxis = function(angle = 90, hjust = 1, vjust = 0.5){
  return(theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust)))
}

PrepStack = function(plt.list, nrow = NULL, ncol = NULL, remove.legends = TRUE){
  
  # Remove x axis labels for all but last row
  plt.list[1:(length(plt.list)-ncol)] = lapply(plt.list[1:(length(plt.list)-ncol)], function(plot) {
    plot = plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    return(plot)
  })
  
  # Remove y axis labels for all but first column
  if(ncol > 1){
    plt.list[1:length(plt.list) %% ncol != 1] = lapply(plt.list[1:length(plt.list) %% ncol != 1], function(plot) {
      plot = plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
      return(plot)
    })
  }
  
  # Remove all legends
  if(remove.legends){
    plt.list = lapply(plt.list, function(plot) plot + NoLegend())
  }
  
  return(plt.list)
}

StackedPlots = function(plt.list, height = 2, remove.titles = TRUE, no.legend = FALSE, combine.legend = TRUE, ...){
  
  # Modify plots
  # plt.list = lapply(plt.list, function(plot) {
  #   gene = plot$labels$title
  #   plot = plot + ylab(gene) + theme(axis.title.x = element_blank(), plot.title = element_blank())
  #   return(plot)
  # })
  
  # Remove axis labels for all but last plot
  plt.list[1:(length(plt.list)-1)] = lapply(plt.list[1:(length(plt.list)-1)], function(plot) {
    plot = plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    return(plot)
  })
  
  # Remove legend from all except middle panel
  if(combine.legend) {
    middle.element = ceiling(length(plt.list)/2)
    plt.list[-middle.element] = lapply(plt.list[-middle.element], function(plot) plot + NoLegend()) 
  }
  
  # Remove all titles
  if(remove.titles){
    plt.list = lapply(plt.list, function(plot) {
      plot = plot + theme(plot.title = element_blank()) 
      return(plot)
    })
  }
  
  # Remove legend
  if(no.legend){
    plt.list = lapply(plt.list, function(plot) {
      plot = plot + NoLegend() 
      return(plot)
    })
  }
  
  # plt = ggarrange(plotlist = plt.list,
  #                 nrow = length(plt.list),
  #                 # align = "v",
  #                 # axis = "bt",
  #                 heights = c(rep(1, length(plt.list)-1), height),
  #                 ...)
  
  plt = patchwork::wrap_plots(plotlist = plt.list, ncol = 1)
  
  return(plt)
}

ViolinPlotFlipped = function(object, height = 2, ...){
  plt.list = VlnPlot(object, ..., combine = FALSE) 
  
  # Flip coordinates
  plt.list = lapply(plt.list, function(plot) plot + scale_y_continuous(position = "right") + coord_flip())
  
  # Modify plots
  plt.list = lapply(plt.list, function(plot) {
    gene = plot$labels$title
    plot = plot + ylab(gene) + theme(axis.title.y = element_blank(), 
                                     plot.title = element_blank(), 
                                     legend.position = "none")
    return(plot)
  })
  
  # Remove axis labels for all but last plot
  plt.list[2:(length(plt.list))] = lapply(plt.list[2:(length(plt.list))], function(plot) {
    plot = plot + theme(axis.text.y = element_blank())
    return(plot)
  })
  
  # plt = ggarrange(plotlist = plt.list, 
  #           common.legend = TRUE, 
  #           nrow = length(plt.list), 
  #           legend = 'right', 
  #           align = "v", 
  #           heights = c(rep(1, length(plt.list)-1), height))
  
  plt = patchwork::wrap_plots(plotlist = plt.list, nrow = 1)
  return(plt)
}

ViolinPlot = function(object, height = 2, ...){
  plt.list = VlnPlot(object, ..., combine = FALSE)
  
  # Modify plots
  plt.list = lapply(plt.list, function(plot) {
    gene = plot$labels$title
    plot = plot + ylab(gene) + theme(axis.title.x = element_blank(), plot.title = element_blank())
    return(plot)
  })
  
  # Remove axis labels for all but last plot
  plt.list[1:(length(plt.list)-1)] = lapply(plt.list[1:(length(plt.list)-1)], function(plot) {
    plot = plot + theme(axis.text.x = element_blank())
    return(plot)
  })
  
  plt = ggarrange(plotlist = plt.list, 
                  common.legend = TRUE, 
                  nrow = length(plt.list), 
                  legend = 'right', 
                  align = "v", 
                  heights = c(rep(1, length(plt.list)-1), height))
  return(plt)
}

GetSpecies = function(object, field = "species"){
  return(unique(object@meta.data[[field]]))
}

# Annotate clusters using a dataframe with old and new names
AnnotateClusters = function(object, key, from = "seurat_clusters", to = "annotated", factor = FALSE){
  new.names = key$new.names[match(as.character(object@meta.data[,from]), as.character(key$old.names))]
  if(factor) return(factor(new.names, levels = sort(unique(key$new.names)))) else return(new.names)
}

theme_umap = function(plt, remove.axes = FALSE, title = NULL, rename.axes = TRUE){
  
  # Correct xlab ylab
  if(rename.axes) plt = plt + labs(x = 'UMAP 1', y = 'UMAP 2')
  
  if(remove.axes) {
    plt2 = plt + theme_void()
  } else{
    plt2 = plt + theme_dario() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank())
  }
  
  # Center title 
  plt2 = plt2 + 
    {if(!is.null(title)) ggtitle(title)}+ 
    theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.title=element_blank()) + 
    ArialFont()
  
  return(plt2)
}

TransferSACLabels = function(object){
  AC_model <- TrainXGBoost(PengAC, object, train.clusters = "lit_type")
  # saveRDS(AC_model, "../../data/train_Peng_test_Rat.rds")
  # train_Peng_test_Rat <- BuildConfusionMatrix(object, PengAC, model = AC_model)
  object$xgb.labels = PredictLabels(object, model = AC_model, scale.by.model = TRUE)
  object$seurat_clusters = factor(object$seurat_clusters)
  print(JSHeatmap(JSMatrix(table(object$xgb.labels, object$seurat_clusters)), heatmap = FALSE))
  # print(plotConfusionMatrix(table(as.character(object$xgb.labels), object$seurat_clusters)))
  return(object)
}

plotSACMarkers = function(object, group.by = "seurat_clusters", feature.plot = FALSE){
  if(feature.plot){
    plot_grid(
      {if(any(ON_SAC_markers %in% rownames(object))) FeaturePlot(object, features = ON_SAC_markers, ncol = length(ON_SAC_markers))}, 
      {if(any(OFF_SAC_markers %in% rownames(object))) FeaturePlot(object, features = OFF_SAC_markers, ncol = length(OFF_SAC_markers))}, 
      ncol = 1, nrow = 2)
  } else {
    plot_grid(
      {if(any(ON_SAC_markers %in% rownames(object))) VlnPlot(object, group.by = group.by, features = ON_SAC_markers, ncol = length(ON_SAC_markers))}, 
      plot_grid(
        {if(any(OFF_SAC_markers %in% rownames(object))) VlnPlot(object, group.by = group.by, features = OFF_SAC_markers, ncol = length(OFF_SAC_markers))}, 
        VlnPlot(object, features = c("ON.SAC.score", "OFF.SAC.score")),
        ncol = 2, rel_widths = c(3,2)
      ),
      ncol = 1, nrow = 2)
  }
}

plotXGBLabels = function(object){
  plot_grid(
    plot_grid(
      DimPlot(object, group.by = "seurat_clusters", label = TRUE) + NoLegend(), 
      DimPlot(object, group.by = "xgb.labels", label = TRUE) + NoLegend()
    ),
    plotSACMarkers(object, group.by = "xgb.labels"), 
    ncol = 1, nrow = 2)
}

SubclusterSACs = function(object, group.by = "animal", resolution = 0.2, harmonize = TRUE){
  if(harmonize == TRUE){
    object = Harmonize(object, batch = group.by, cluster_resolution = resolution, show.plots = FALSE)
  } else {
    object = ClusterSeurat(object, cluster_resolution = resolution)
  }
  print(BrowseSeurat(object, batch = group.by))
  object = TransferSACLabels(object)
  print(plotSACMarkers(object))
  print(plotSACMarkers(object, feature.plot = TRUE))
  return(object)
}

BarPlot = function(data, x, y, color = NULL){
  p = ggplot(data, aes(x=eval(parse(text = x)), y=eval(parse(text = y)))) + 
    geom_col(color = "black") +
    xlab(x)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    ylab(y)+
    theme_cowplot()+
    RotatedAxis()
  
  return(p)
}

stackedBarGraph2 = function(data, x, y, normalize = TRUE, position = position_stack(), 
                            label = FALSE, label.threshold = 0.05, 
                            as.factor = FALSE, border.color = 'black', percent = TRUE){
  
  if(inherits(data, "Seurat")) data = data@meta.data
  
  # data[,x] = as.character(data[,x])
  # data[,y] = as.character(data[,y])
  cross.tabulation = table(data[,x], data[,y])
  
  if(normalize) cross.tabulation = cross.tabulation/rowSums(cross.tabulation)
  melted = reshape2::melt(cross.tabulation)
  melted$Var2 = melted$Var2
  melted$label = melted$Var2
  melted$label[melted$value < label.threshold] = NA
  
  if(as.factor) melted$Var1 = factor(melted$Var1, levels = levels(data[,x]))
  
  p = ggplot(melted, aes(fill=Var2, y=value, x=Var1, label = label)) + 
    {if(normalize) geom_bar(position="stack", stat="identity", color = border.color)}+
    {if(!normalize) geom_bar(position='stack', stat="identity", color = border.color)}+
    theme_minimal()+
    # xlab(x)+
    labs(fill=y)+
    {if(label) geom_text(size = 3, position = position_stack(vjust = 0.5))}+
    # scale_x_continuous(breaks = seq(0, nClusters(seurat)-1, 1))+
    # geom_text_repel()+
    # scale_fill_manual(feature.1)+
    # scale_x_discrete()
    # scale_y_continuous(expand = c(0, 0), 
    #                    breaks = c(0, 1),  # only first and last
    #                    labels = c(0, 1),
    #                    minor_breaks = waiver()  # keeps automatic minor ticks
    #                    )+
    {if(percent) scale_ticks_first_last(do = function(x) x*100)}+
    {if(!percent) scale_ticks_first_last()}+
    ylab('Proportion')+
    theme_cowplot()+
    # theme_dario()+
    RotatedAxis() + 
    theme(axis.title.x = element_blank())
  # theme(plot.margin=unit(c(.2,.5,.2,.2),"cm"))
  return(p)
}

CorrelationHeatmap = function(object, 
                              assay = "RNA", 
                              group.by = "seurat_clusters", 
                              method = 'pearson',
                              annotate.by = NULL, 
                              annotate.numeric = NULL,
                              cluster_rows = TRUE, 
                              features = NULL, 
                              title = NULL, 
                              label = FALSE, 
                              cols = c("white", "red", 'darkred'), 
                              annotation_cols = NULL, 
                              plot = TRUE, 
                              return.dendrogram = FALSE, 
                              return.matrix = FALSE,
                              return.plot = FALSE,
                              label.size = 6,
                              ...){
  
  if(!is.null(features)){
    avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = features)[[assay]]))
  } else {
    message("Using top 2000 variable features!")
    object = FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
    avg_exp = as.data.frame(log1p(AverageExpression(object, verbose = FALSE, group.by = group.by, assay = assay, features = object@assays[[assay]]@var.features)[[assay]]))
  }
  # print(dim(avg_exp))
  type_cor = cor(avg_exp, method = method)
  
  # Categorical
  if(!is.null(annotate.by) & !is.null(annotate.numeric)){
    
    metadata = Metadata(object, feature.1 = group.by, feature.2 = annotate.by)
    annotation = metadata[match(colnames(type_cor), metadata[[group.by]]), annotate.by]
    
    mean.meta = MeanMetadata(object, feature.1 = group.by, feature.2 = annotate.numeric)
    mean.meta = mean.meta[match(colnames(type_cor), names(mean.meta))]
    
    if(is.null(annotation_cols)){
      ha = HeatmapAnnotation(anno = annotation, 
                             barplot = anno_barplot(mean.meta, 
                                                    gp = gpar(fill = "grey"), # Bar color
                                                    border = TRUE),
                             annotation_name_side = "left")
    } else{
      ha = HeatmapAnnotation(anno = annotation, 
                             barplot = anno_barplot(mean.meta, 
                                                    gp = gpar(fill = "grey"), # Bar color
                                                    border = TRUE),
                             annotation_name_side = "left", 
                             col = list(anno = annotation_cols))
    }
  } else if(!is.null(annotate.by)){
    
    metadata = Metadata(object, feature.1 = group.by, feature.2 = annotate.by)
    annotation = metadata[match(colnames(type_cor), metadata[[group.by]]), annotate.by]
    
    if(is.null(annotation_cols)){
      ha = HeatmapAnnotation(anno = annotation, 
                             annotation_name_side = "left")
    } else{
      ha = HeatmapAnnotation(anno = annotation, 
                             annotation_name_side = "left", 
                             col = list(anno = annotation_cols))
    }
  } else if(!is.null(annotate.numeric)){
    
    mean.meta = MeanMetadata(object, feature.1 = group.by, feature.2 = annotate.numeric)
    annotation = mean.meta[match(colnames(type_cor), names(mean.meta))]
    ha = HeatmapAnnotation(barplot = anno_barplot(annotation, 
                                                  gp = gpar(fill = "grey"), # Bar color
                                                  border = TRUE),
                           annotation_name_side = "left")
  } else {
    ha = NULL
  }
  
  if(inherits(cluster_rows, 'dendrogram')){
    ht = Heatmap(type_cor, 
                 name = "cor", 
                 col = cols, 
                 cluster_rows  = as.hclust(cluster_rows),
                 cluster_columns = as.hclust(cluster_rows),
                 border_gp = gpar(col = "black", lty = 1),
                 #show_row_names = TRUE, show_column_names = TRUE, row_dend_side = "left", 
                 show_column_dend = TRUE, 
                 column_title = title,
                 cell_fun = if(label) function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", type_cor[i, j]), x, y, gp = gpar(fontsize = label.size))} else NULL,
                 heatmap_legend_param = list(title = paste0(method, " R")),
                 top_annotation = ha, 
                 ...)
    
  } else if(cluster_rows == TRUE){
    ht = Heatmap(type_cor, 
                 name = "cor", 
                 col = cols, #colorRamp2(seq(0,1,0.1), rev(rainbow(11))), #colorRamp2(c(-1,0,1), c("blue", "white", "red")),
                 cluster_rows = hclust(as.dist(1 - type_cor)),
                 cluster_columns = hclust(as.dist(1 - type_cor)), 
                 border_gp = gpar(col = "black", lty = 1),
                 #show_row_names = TRUE, show_column_names = TRUE, row_dend_side = "left", 
                 show_column_dend = TRUE, 
                 column_title = title,
                 cell_fun = if(label) function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", type_cor[i, j]), x, y, gp = gpar(fontsize = label.size))} else NULL,
                 heatmap_legend_param = list(title = paste0(method, " R")),
                 top_annotation = ha, 
                 ...)
    
  } else if(cluster_rows == FALSE){
    ht = Heatmap(type_cor, 
                 name = "cor", 
                 col = cols, #colorRamp2(seq(0,1,0.1), rev(rainbow(11))),
                 cluster_rows = FALSE, 
                 # cluster_columns = FALSE, 
                 # show_row_names = TRUE, 
                 # show_column_names = TRUE, 
                 row_dend_side = "left", 
                 cell_fun = if(label) function(j, i, x, y, width, height, fill) { grid.text(sprintf("%.2f", type_cor[i, j]), x, y, gp = gpar(fontsize = label.size))} else NULL,
                 show_column_dend = TRUE, 
                 column_title = title,
                 heatmap_legend_param = list(title = paste0(method, " R")),
                 top_annotation = ha, 
                 ...)
    
  } else {
    stop('could not parse cluster_rows')
  }
  
  if(plot) ht = draw(ht)
  
  if(return.dendrogram) return(row_dend(ht))
  if(return.matrix) return(type_cor)
  if(return.plot) return(ht)
}

EvaluateModel = function(table, verbose = TRUE){
  
  TP = table[1,1]
  TN = table[2,2]
  FP = table[1,2]
  FN = table[2,1]
  accuracy = (TP+TN)/sum(rowSums(table))
  sensitivity = (TP)/(TP + FN) # true positive rate
  specificity = (TN)/(TN + FP) # true negative rate
  precision = TP/(TP+FP) # proportion of predicted positives that are truly positive
  recall = TP/(TP+FN) # same as sensitivity
  accuracy = (TP+TN)/sum(rowSums(table))
  if(verbose) {
    message(paste0(capture.output(table), collapse = '\n'))
    message("accuracy: ", accuracy, "\ntrue positive rate (AKA sensitivity or recall): ", sensitivity, "\nfalse positive rate: ", 1-specificity, '\nprecision: ', precision)
  }
  
  list(accuracy = accuracy, precision = precision, recall = recall, sensitivity = sensitivity, specificity = specificity)
}

ScatterHistogram = function(object, x, y, color = NULL, offset = -0.2){
  
  if(inherits(object, "Seurat")){
    data = object@meta.data
  } else {
    data = object
  }
  
  p1 = gghistogram(data, x = x, xlab = FALSE) + theme(axis.text.x = element_blank()) + 
    theme(plot.margin = unit(c(0,0,0,0), 'mm'))
  p2 = NULL
  p3 = ggscatter(data, x, y, color = color, title = NULL) + 
    theme(legend.position = "bottom", plot.title = element_blank()) + 
    coord_cartesian(clip = 'off') + 
    theme(plot.margin = unit(c(0,0,0,0), 'mm'))
  p4 = gghistogram(data, x = y, ylab = FALSE) + 
    coord_flip() + 
    RotatedAxis() + 
    theme(axis.text.y = element_blank()) + 
    theme(plot.margin = unit(c(0,0,0,0), 'mm'))
  
  plot_grid(
    p1, NULL, NULL, 
    NULL, NULL, NULL, 
    p3, NULL, p4,
    ncol = 3, 
    nrow = 3, 
    rel_widths = c(2, offset, 1), 
    rel_heights = c(1, offset, 2), 
    align = "hv", 
    axis = "btlr")
  
  # patchworked
  # (p1 | plot_spacer()) / (p3 | p4)
}

MixtureHistogram = function(df, variable, model, intersects = NULL, color.1 = "magenta", color.2 = "blue", legend.name = NULL, xlab = NULL, ylab = "Density", plot.points = TRUE, logY = FALSE){
  ggplot(df, aes(x=eval(parse(text = variable)))) +
    geom_histogram(bins = 20, aes(y=..density..), colour="black", fill = "white", boundary=0) + 
    # geom_density(aes(y=..density..)) +
    geom_function(fun = function(x) model$lambda[1]*dnorm(x, mean = model$mu[1], sd = model$sigma[1]), linewidth = 1, linetype = "solid", color = color.1)+
    geom_function(fun = function(x) model$lambda[2]*dnorm(x, mean = model$mu[2], sd = model$sigma[2]), linewidth = 1, linetype = "solid", color = color.2)+
    geom_vline(xintercept = c(intersects), linetype = "dashed", color = "grey")+
    scale_x_continuous(breaks = seq(0,1, by = 0.25), expand = c(0,0))+
    {if(plot.points) geom_point(aes(fill = eval(parse(text = legend.name)), y = 0), shape = 21, color = "black", size = 3)}+
    scale_fill_gradient(name = legend.name, limits = c(0,1), low = color.1, high = color.2)+
    {if(logY) scale_y_continuous(trans = "log10")}+
    ylab(ylab)+
    xlab(xlab)+
    theme_dario()
}

ComputePosteriors = function(xi, model){
  sigma1 = model$sigma[1]
  sigma2 = model$sigma[2]
  mu1 = model$mu[1]
  mu2 = model$mu[2]
  lambda1 = model$lambda[1]
  lambda2 = model$lambda[2]
  
  p_xi_1 = dnorm(xi, mean = mu1, sd = sigma1)
  p_1 = lambda1
  
  p_xi_2 = dnorm(xi, mean = mu2, sd = sigma2)
  p_2 = lambda2
  
  # Bayes Theorem
  results.df = data.frame(p_1_xi = (p_xi_1 * p_1)/(p_xi_2 * p_2 + p_xi_1 * p_1), 
                          p_2_xi = (p_xi_2 * p_2)/(p_xi_2 * p_2 + p_xi_1 * p_1))
  results.df$check = results.df$p_1_xi + results.df$p_2_xi
  
  return(results.df)
}

FindCutpoint = function(model){
  sigma1 = model$sigma[1]
  sigma2 = model$sigma[2]
  mu1 = model$mu[1]
  mu2 = model$mu[2]
  alpha = model$lambda[1]
  beta = model$lambda[2]
  
  # Intersection point of two Gaussians is where probability flips (adapted from https://stats.stackexchange.com/questions/311592/how-to-find-the-point-where-two-normal-distributions-intersect)
  A = (-1/sigma1^2+1/sigma2^2)
  B = 2*(-mu2/sigma2^2 + mu1/sigma1^2)
  C = (mu2^2/sigma2^2) - (mu1^2/sigma1^2) + log(sigma2^2/sigma1^2) + 2*log(alpha/beta)
  discriminant = B^2 - 4*A*C
  intersect1 = (-B + sqrt(discriminant))/(2*A)
  intersect2 = (-B - sqrt(discriminant))/(2*A)
  intersect = c(intersect1, intersect2)[which(c(intersect1, intersect2) > min(c(mu1,mu2)) & c(intersect1, intersect2) < max(c(mu1,mu2)))]
  return(intersect)
}

UpdateLitTypes = function(object, field = "lit_type"){
  previous_names = object@meta.data[[field]]
  object@meta.data[[field]][previous_names == "AII"] = "A2"
  object@meta.data[[field]][previous_names == "CA-I"] = "CA1"
  object@meta.data[[field]][previous_names == "CA-II"] = "CA2"
  return(object)
}

# ConservedEnrichedHeatmap3 = function(object1, object2, object3, geneList, color, ...){
#   
#   DoHeatmap(object1, features = unlist(geneList), ...) + 
#     # scale_x_discrete(expand = c(0,0))+
#     # NoLegend() +
#     # theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#     # scale_fill_gradientn(colors = c("white", "white", color))
#   # coord_flip()
#   
#   DoHeatmap(object2, features = unlist(geneList), ...) + 
#     # NoLegend() +
#     # theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#     # scale_x_discrete(expand = c(0,0))+
#     # scale_fill_gradientn(colors = c("white", "white", color))
#   # coord_flip()
#   
#   DoHeatmap(object3, features = unlist(geneList), ...)  
#     # NoLegend() +
#     # scale_x_discrete(expand = c(0,0))+
#     # theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#     # scale_fill_gradientn(colors = c("white", "white", color))
#   # coord_flip()
#   
#   panel = plot_grid(p1, NULL, p2, NULL, p3, 
#                     # rel_widths = c(3, -0.05, 3, -0.05, 3),
#                     rel_widths = c(3, 0, 3, 0, 3),
#                     ncol = 5) 
#   
#   return(panel)
# }

ConservedEnrichedHeatmap2 = function(object1, object2, object3, geneList, color, ...){
  
  p1 = DoHeatmap(object1, features = unlist(geneList), ...) + 
    scale_x_discrete(expand = c(0,0))+
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", color))
  # coord_flip()
  
  p2 = DoHeatmap(object2, features = unlist(geneList), ...) + 
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_gradientn(colors = c("white", "white", color))
  # coord_flip()
  
  p3 = DoHeatmap(object3, features = unlist(geneList), ...) + 
    NoLegend() +
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", color))
  # coord_flip()
  
  panel = plot_grid(p1, NULL, p2, NULL, p3, 
                    # rel_widths = c(3, -0.05, 3, -0.05, 3),
                    rel_widths = c(3, 0, 3, 0, 3),
                    ncol = 5) 
  
  return(panel)
}

ConservedEnrichedHeatmap = function(object, sharedList, primateList, rodentList, laurasiaList, ...){
  # Shared
  s1 = DoHeatmap(object, features = (unlist(sharedList)),  ...) + 
    scale_x_discrete(expand = c(0,0))+
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", "magenta"))
  # coord_flip()
  
  # Species-specific 
  p1 = DoHeatmap(object, features = unlist(primateList),  ...) + 
    scale_x_discrete(expand = c(0,0))+
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", "blue"))
  # coord_flip()
  
  p2 = DoHeatmap(object, features = unlist(rodentList),  ...) + 
    NoLegend() +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_gradientn(colors = c("white", "white", "turquoise"))
  # coord_flip()
  
  p3 = DoHeatmap(object, features = unlist(laurasiaList),  ...) + 
    NoLegend() +
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_gradientn(colors = c("white", "white", "orange"))
  # coord_flip()
  
  panel = plot_grid(s1,p1,p2,p3, 
                    ncol = 4) 
  
  return(panel)
}

theme_dario = function(margin = 10, linewidth = 1){
  dario.theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                      panel.border = element_rect(fill = NA, colour = "black", linewidth=linewidth),
                      panel.background = element_blank(), 
                      plot.margin = margin(t = margin, r = margin, b = margin, l = margin),
                      axis.line = element_blank())+
                      # panel.background = element_rect(fill = 'white', colour = "black", linewidth=1), 
                      # axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), 
                      # plot.title = element_text(hjust = 0.5))  +
                ArialFont()
  return(dario.theme)
}

PlotCelltypeProportions = function(objectList, proportions, celltype, color = "orange", reference.line){
  
  SAC.clusters = lapply(objectList, function(object) {
    metadata = Metadata(object, feature.1 = "seurat_clusters", feature.2 = "lit_type")
    return(metadata$seurat_clusters[which(metadata$lit_type == celltype)])
  })
  
  SAC.clusters = SAC.clusters[match(names(proportions), names(SAC.clusters))]
  
  SAC.props = data.frame(species = factor(c(names(SAC.clusters)), levels = phylogenetic_order), 
                         proportion = 100*c(unlist(sapply(seq_along(proportions), function(i) sum(proportions[[i]][SAC.clusters[[i]] ])))), 
                         group = c(rep("Test", length(objectList))))
  
  ggplot(SAC.props, aes(x = species, y = proportion, fill = group)) + 
    ylab(paste0(celltype, " proportion (%)"))+
    geom_col(color = "black")+
    geom_hline(yintercept = reference.line, linetype = "dashed")+
    scale_fill_manual(values = c(color))+
    scale_y_continuous(expand = c(0,0))+
    theme_cowplot()+
    theme(legend.position = "none")+
    RotatedAxis() 
  # coord_flip()
}

# Downsample seurat by variable identity class
DownsampleSeurat = function(object, group.by, size, seed = 12345, verbose = FALSE){
  # previous_ident = Idents(object)
  
  group_ident = paste(c(group.by), collapse = '.')
  if(length(group.by) > 1){
    object[[group_ident]] = apply(object@meta.data[, group.by], 1, paste, collapse = "-")
  }
  Idents(object) = group_ident
  object = subset(object, downsample = size, seed = seed)
  # object = subset(object, cells = WhichCells(seurat, downsample = size, seed = seed))
  # Idents(object) = previous_ident
  
  if(verbose) print(t(t(table(object@meta.data[[group_ident]]))))
  
  return(object)
}

ComputeIntegrationLISI  = function(object, group.by = "species", plot = TRUE, subsample.types = TRUE, nPCs = 20, method = 'seurat'){
  
  if(subsample.types) {
    # Downsample to minimum category size
    downsample = min(table(object[[group.by]]))
    object = DownsampleSeurat(object, group.by, downsample)
  }
  
  if(method == 'seurat'){
    iLISI_observed <- compute_lisi(object@reductions$pca@cell.embeddings[,1:nPCs], 
                                   object@meta.data, c(group.by))
  } else if(method == 'harmony') {
    iLISI_observed <- compute_lisi(object@reductions$harmony@cell.embeddings[,1:nPCs], 
                                   object@meta.data, c(group.by))
  } else {
    stop('method should be seurat or harmony')
  }
  
  iLISI_expected = SimpsonIndex(table(object@meta.data[[group.by]]), invert = TRUE)
  
  # Normalized iLISI
  object$normalized.iLISI = (iLISI_observed[[group.by]]-1)/(iLISI_expected-1)
  
  if(plot){
    print(plot_grid(ClusterBatchPlot(object, batch = group.by, shuffle = TRUE), 
                    FeaturePlot(object, features = "normalized.iLISI"), 
                    VlnPlot(object, features = "normalized.iLISI", pt.size = 0)+NoLegend(),
                    VlnPlot(object, features = "normalized.iLISI", group.by = group.by, pt.size = 0)+NoLegend(),
                    ncol = 2, nrow = 2))
  }
  
  return(object)
}

ComputeLittypeLISI = function(object, types.to.consider = c("A2", "VG3", "SAC", "A17"), 
                              subset.types = FALSE, subsample.types = TRUE, nPCs = 20, plot = TRUE){
  knownAC = object
  knownAC$species_littype[!knownAC$species_littype %in% types.to.consider] = "Other"
  
  if(subset.types) knownAC = subset(knownAC, species_littype %in% types.to.consider)
  if(subsample.types) {
    # Downsample to minimum category size
    downsample = min(table(knownAC[["species_littype"]]))
    knownAC = DownsampleSeurat(knownAC, "species_littype", downsample)
  }
  
  cLISI_observed <- compute_lisi(knownAC@reductions$pca@cell.embeddings[,1:nPCs], 
                                 knownAC@meta.data, c('species_littype'))
  
  cLISI_expected = SimpsonIndex(table(knownAC$species_littype), invert = TRUE) # length(unique(knownAC$species_littype)) 
  
  # Normalized cLISI
  knownAC$normalized.cLISI = (cLISI_observed$species_littype-1)/(cLISI_expected-1)
  if(plot){
    print(plot_grid(DimPlot(knownAC, group.by = "species_littype", label = TRUE, shuffle = TRUE) + NoLegend(), 
                  ClusterBatchPlot(knownAC, batch = "species", group.by = "species_littype", shuffle = TRUE), 
                  ClusterFeaturePlot(knownAC, features = "normalized.cLISI", group.by = "species_littype"), 
                  plot_grid(VlnPlot(subset(knownAC, species_littype != "Other"), features = "normalized.cLISI", group.by = "species_littype", pt.size = 0, log = F) + NoLegend(),
                            VlnPlot(subset(knownAC, species_littype != "Other"), features = "normalized.cLISI", group.by = "species", pt.size = 0, log = F) + NoLegend(), 
                            ncol = 2, align = "h", rel_widths = c(1,2)), 
                  ncol = 2, nrow = 2, rel_widths = c(1,1.2)))
  }
  return(knownAC)
}

# ComputeLittypeLISITest = function(object, types.to.consider = c("A2", "VG3", "SAC", "A17")){
#   knownAC = object
#   knownAC$species_littype = sample(types.to.consider, ncol(knownAC), replace = TRUE)
#   
#   cLISI_observed <- compute_lisi(knownAC@reductions$umap@cell.embeddings, 
#                                  knownAC@meta.data, c('species_littype'))
#   
#   cLISI_expected = SimpsonIndex(table(knownAC$species_littype), invert = TRUE) # length(unique(knownAC$species_littype)) 
#   
#   # Normalized cLISI (no normalization for now)
#   knownAC$normalized.cLISI = cLISI_observed$species_littype # (cLISI_observed$species_littype-1)/(species_littype-1)
#   print(plot_grid(DimPlot(knownAC, group.by = "species_littype", label = TRUE, shuffle = TRUE) + NoLegend(), 
#                   ClusterBatchPlot(knownAC, batch = "species", group.by = "species_littype", shuffle = TRUE), 
#                   ClusterFeaturePlot(knownAC, features = "normalized.cLISI", group.by = "species_littype"), 
#                   plot_grid(VlnPlot(subset(knownAC, species_littype != "Other"), features = "normalized.cLISI", group.by = "species_littype", pt.size = 0, log = F) + NoLegend(),
#                             VlnPlot(knownAC, features = "normalized.cLISI", group.by = "species", pt.size = 0, log = F) + NoLegend(), ncol = 2, align = "h"), 
#                   ncol = 2, nrow = 2))
#   return(knownAC)
# }

SubsetSeuratGenes = function(object, features){
  counts = object@assays$RNA@counts
  newObject = CreateSeuratObject(counts[rownames(counts) %in% features,])
  newObject = NormalizeData(newObject, verbose = FALSE)
  nCount_RNA = newObject$nCount_RNA
  nFeature_RNA = newObject$nFeature_RNA
  newObject@meta.data = object@meta.data
  newObject$nCount_RNA = nCount_RNA
  newObject$nFeature_RNA = nFeature_RNA
  # newObject@meta.data[, !colnames(newObject@meta.data) %in% c("nCount_RNA", "nFeature_RNA")] = object@meta.data[, !colnames(newObject@meta.data) %in% c("nCount_RNA", "nFeature_RNA")]
  return(newObject)
}

transparent = function(){
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA))
}

TwoGeneDotPlot = function(object, type, genes, max.pct = 100, legend = FALSE, font.size = 12){
  
  # dat1 = geneDotPlotFast(object, gene1)$data %>% mutate(gene = gene1)
  # dat2 = geneDotPlotFast(object, gene2)$data %>% mutate(gene = gene2)
  
  dat = do.call(rbind, lapply(genes, function(this.gene) {
    geneDotPlotFast(ac.ortho, this.gene)$data %>% mutate(gene = this.gene)
    }))
  dat = subset(dat, OT == type)
  dat$gene = factor(dat$gene, levels = genes)
  
  plt = TitlePlot(
    ggplot(dat, aes(y=factor(species, levels = rev(levels(object$species))), x = gene, alpha=avg.exp.scaled, color = gene, size=pct.exp))+
      geom_point(shape = 16)+
      xlab(NULL)+
      ylab(NULL)+
      theme_classic(base_size = font.size, base_family = 'ArialMT')+
      {if(!legend) NoLegend()}+
      coord_cartesian(clip = 'off')+
      scale_color_discrete(guide = "none")+
      scale_alpha_continuous(limits = c(-1, 2), 
                             range = c(0.1,1),
                             oob = scales::squish)+
      theme(axis.text.y=element_text(colour="black"),
            axis.text.x=element_text(colour="black", angle=30, hjust=1, face = 'italic'))+
      transparent()+
      # ArialFont()+
      scale_radius(name = 'Percent\nexpressed', limits = c(0, max.pct), range = c(0, 6)),
    orthotype_labels(type)) + theme(plot.title = element_text(hjust = 0.9))
  
  plt
}

MultigeneDotPlot = function(object, genes, remove.titles = FALSE, no.legend = FALSE, ...){
  StackedPlots(plt.list = lapply(genes, function(gene) geneDotPlotFast(object, gene, ...)), remove.titles = remove.titles, no.legend = no.legend)
}

geneDotPlotFast = function(object, gene, group.by = "orthotype", scale.within.species = TRUE, mini = FALSE, 
                           col.low = 'lightgrey', col.high = '#584B9FFF', mc.cores = 20, max.pct = 100, 
                           binarization = NULL){
  
  stopifnot(is.factor(object$species))
  
  if(scale.within.species){
    Idents(object) = object$species
    dat = mclapply(unique(object$species), function(x) DotPlot(object, ident = x, features=gene, group.by = group.by, col.min = -1, col.max = 2)$data, mc.cores = mc.cores)
    names(dat) = unique(object$species)
    dat = do.call(rbind, dat)
    dat$species = stringr::str_split_fixed(rownames(dat), "\\.", 2)[,1]
    dat$OT = dat$id
  } else{
    Idents(object) = paste0(object$species, "_", object@meta.data[[group.by]])
    dat=DotPlot(object, features=gene)$data
    # metadata = Metadata(object, feature.1 = "species", feature.2 = group.by)
    dat$species = str_split_fixed(dat$id, "_", 2)[,1]
    dat$OT = factor(str_split_fixed(dat$id, "_", 2)[,2], levels = levels(object@meta.data[[group.by]]))
  }
  
  # Plot bigger values last
  dat = dat %>% arrange(pct.exp)
  dat$key = paste0(dat$species, '-', dat$OT)
  
  if(!is.null(binarization)){
    bin.melt = na.omit(reshape2::melt(as.matrix(binarization)))
    bin.melt$key = paste0(bin.melt$Var1, '-', bin.melt$Var2)
    dat = dat[match(bin.melt$key, dat$key),]
    dat$heat_value = factor(bin.melt$value, levels = c(0,1))
  }
  
  # Modify max.pct
  dat$pct.exp[dat$pct.exp > max.pct] = max.pct
  
  # Plot
  plt = ggplot(dat, aes(y=factor(species, levels = rev(levels(object$species))), x = OT))+
    {if(!is.null(binarization)) geom_tile(aes(fill = heat_value))}+
    {if(!is.null(binarization)) scale_fill_manual(values = c("0" = "white", "1" = "lightgrey"))}+
    geom_point(aes(color=avg.exp.scaled, size=pct.exp))+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle(gene)+
    theme_classic()+ 
    coord_cartesian(clip = 'off')+
    # guides(color = guide_colorbar(title = "Average\nexpression"),
           # size = guide_legend(title = "Percent\nexpressed"))+
    # scale_color_gradient(low=col.low, high=col.high, #limits = c(-1,2),
    #                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    theme(plot.title=element_text(hjust=0.5, face = 'italic'),
          axis.text.y=element_text(colour="black"),
          axis.text.x=element_text(colour="black", angle=45, hjust=1))+
    scale_color_gradient(name = 'Scaled\nexpression', 
                         low = col.low, 
                         high = col.high,
                         limits = c(-1, 2), 
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    scale_radius(name = 'Percent\nexpressed', limits = c(0, max.pct), range = c(0, 6)) + 
    ArialFont()
  
  if(mini){
    plt = plt + 
      ggtitle(gene)+
      theme(axis.text.x = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks = element_blank(), 
            axis.title = element_blank(),
            plot.title=element_text(hjust=0.5, face = 'italic'),
            legend.position = 'none')
  }
  
  return(plt)
}

geneDotPlot <- function(object, gene, group.by = "type", scale.within.species = TRUE){
  
  if(scale.within.species){
    Idents(object) = object$species
    dat = lapply(unique(object$species), function(x) DotPlot(object, ident = x, features=gene, group.by = group.by)$data)
    names(dat) = unique(object$species)
    dat = do.call(rbind, dat)
    dat$species = stringr::str_split_fixed(rownames(dat), "\\.", 2)[,1]
    dat$OT = dat$id
  } else{
    Idents(object) = paste0(object$species, "_", object@meta.data[[group.by]])
    dat=DotPlot(object, features=gene)$data
    # metadata = Metadata(object, feature.1 = "species", feature.2 = group.by)
    dat$species = str_split_fixed(dat$id, "_", 2)[,1]
    dat$OT = factor(str_split_fixed(dat$id, "_", 2)[,2], levels = levels(object@meta.data[[group.by]]))
  }
  
  # Plot
  return(ggplot(dat, aes(y=factor(species, levels = rev(unique(species))), x = OT, color=avg.exp.scaled, size=pct.exp))+
           geom_point()+
           xlab("OT")+
           ylab("Species")+
           ggtitle(gene)+
           theme_classic()+ 
           theme(plot.title=element_text(hjust=0.5), 
                 axis.text.x=element_text(colour="black", angle=45, hjust=1))+
           scale_color_gradient(low="grey", high="blue"))
}

# We don't want scaling across all celltypes and across all species, we want scaling within species across celltypes
# Old geneSpeciesPlot:                                                      Vsx2180  13.5364298    87.676768           Vsx2                    Pig_RBC     0.367316626 
# DotPlot(object, ident = "Pig", features="Vsx2", group.by = "OT")$data: Vsx2     13.536430     87.67677            Vsx2                    RBC         2.50000000
OTDotPlot <- function(object, gene, celltype, scale.within.species = TRUE){
  
  if(scale.within.species){
    Idents(object) = object$species
    dat = lapply(unique(object$species), function(x) DotPlot(object, ident = x, features=gene, group.by = group.by)$data)
    names(dat) = unique(object$species)
    dat = do.call(rbind, dat)
    dat$species = stringr::str_split_fixed(rownames(dat), "\\.", 2)[,1]
    dat$OT = dat$id
  } else{
    Idents(object) = object$species_OT
    dat=DotPlot(object, features=gene)$data
    dat$species=stringr::str_split_fixed(dat$id, "_", 2)[,1]
    dat$OT=stringr::str_split_fixed(dat$id, "_", 2)[,2]
  }
  
  dat=dat[dat$OT == as.character(celltype),]
  
  # Plot
  return(ggplot(dat, aes(y=factor(species, levels = rev(seurat_datasets)), x=features.plot, color=avg.exp.scaled, size=pct.exp))+
           geom_point()+
           xlab("Gene")+
           ylab("Species")+
           ggtitle(celltype)+
           theme_classic()+ 
           theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(colour="black", angle=45, hjust=1))+
           scale_color_gradient(low="grey", high="blue"))
  #return(dat)
}

ComputeCelltypeLISI = function(object, batch = "animal", nPCs = 20, harmony = FALSE){
  obj.list = SplitObject(object, split.by = batch)
  names(obj.list) = NULL
  cLISI.list = lapply(obj.list, function(each) {
    
    # Using first 20 PCs
    pcs = if(harmony) object@reductions$harmony@cell.embeddings[,1:nPCs] else object@reductions$pca@cell.embeddings[,1:nPCs]
    
    # Subset embedding 
    subset_pcs = pcs[match(Cells(each), rownames(pcs)),]
    
    # Compute cLISI
    batch.cLISI = compute_lisi(subset_pcs, each@meta.data, c("batch_clusters"))
    return(batch.cLISI)
  })
  cLISI.df = do.call(rbind, cLISI.list)
  return(cLISI.df[match(Cells(object), rownames(cLISI.df)), , drop = FALSE])
}

ClusterBatches = function(object, batch = "animal", resolution = 1.5){
  obj.list = SplitObject(object, split.by = batch)
  obj.list = lapply(obj.list, function(each) ClusterSeurat(each, cluster_resolution = resolution, do.umap = FALSE))
  names(obj.list) = NULL
  metadata = do.call(rbind, lapply(obj.list, function(each) each@meta.data))
  metadata$batch_clusters = paste0(metadata[[batch]], "_", metadata[["seurat_clusters"]])
  return(metadata$batch_clusters[match(Cells(object), rownames(metadata))])
}

getProportions = function(object, enrichment.group = c("ALL", "NONE", "CD73", "CD73-")){
  if(length(intersect(object$enrichment, c("ALL", "NONE", "CD73", "CD73-"))) == 0){
    message("No cells found!")
    return(list(NULL, 0))
  }
  object = subset(object, enrichment %in% enrichment.group)
  message("Subsetting to un-enriched samples: ", length(Cells(object)), " cells left")
  props = table(object$seurat_clusters)/sum(table(object$seurat_clusters))
  return(list(props, length(Cells(object))))
}

# Assigns each cluster in object 
AssignAnnotations = function(object, annotation, use = 'seurat_clusters', as = "lit_type"){
  object@meta.data[[as]] = NA
  for(index in seq_along(annotation)){
    object@meta.data[[as]][object@meta.data[[use]] %in% annotation[[index]] ] = names(annotation)[[index]]
  }
  
  # if(factor) object@meta.data[[as]] = factor(object@meta.data[[as]])
  
  return(object)
}

PlotOrthoFeature = function(object, gene, order = FALSE){
  new_genename = paste0(gene, ".")
  object[[new_genename]] = BigMatrix[gene,]
  FeaturePlot(object, feature = new_genename, raster = FALSE, order = order)
}

ExtractBarcodes = function(object, return.barcodes = FALSE, return.aliases = FALSE){
  if(all(endsWith(Cells(object), "-1"))){
    object_barcodes = gsub("-1", "", substrRight(Cells(object), 18))
    
    if(return.barcodes) return(object_barcodes)
    if(return.aliases) {
      # aliases = sapply(seq_along(object_barcodes), function(i) str_split_fixed(Cells(object)[i], object_barcodes[i], 2)[,1])
      aliases = substr(Cells(object), 1, nchar(Cells(object))-19)
      return(aliases)
    }
    
    # Remove ambiguous barcodes
    object_duplicates = names(table(object_barcodes))[table(object_barcodes) > 1]
    object = object[,!object_barcodes %in% object_duplicates]
    object = RenameCells(object, new.names = gsub("-1", "", substrRight(Cells(object), 18)))
    
  } else if(all(endsWith(Cells(object), "x"))){
    object_barcodes = gsub("x", "", str_split_fixed(Cells(object), ":", 2)[,2])
    
    if(return.barcodes) return(object_barcodes)
    if(return.aliases) {
      # aliases = sapply(seq_along(object_barcodes), function(i) str_split_fixed(Cells(object)[i], object_barcodes[i], 2)[,1])
      aliases = str_split_fixed(Cells(object), ":", 2)[,1]
      return(aliases)
    }
    
    # Remove ambiguous barcodes
    object_duplicates = names(table(object_barcodes))[table(object_barcodes) > 1]
    object = object[,!object_barcodes %in% object_duplicates]
    object = RenameCells(object, new.names = gsub("x", "", str_split_fixed(Cells(object), ":", 2)[,2]))
    
  } else {
    stop("Cell names did not match criteria!")
  }
  
  return(object)
}

MapBack = function(object1, object2, group.by = "seurat_clusters", convert.barcodes = FALSE, match.barcodes = FALSE, mc.cores = 1, as = NULL){
  
  if(convert.barcodes){
    object1 = ExtractBarcodes(object1)
    object2 = ExtractBarcodes(object2)
  }
  
  if(match.barcodes){
    # Search over all barcode types in object1 and see if they match all barcode types in object2
    object1$alias = ExtractBarcodes(object1, return.aliases = TRUE)
    object1$barcode = ExtractBarcodes(object1, return.barcodes  = TRUE)
    object2$alias = ExtractBarcodes(object2, return.aliases = TRUE)
    object2$barcode = ExtractBarcodes(object2, return.barcodes  = TRUE)

    # Compare barcodes
    intersects = do.call(rbind, mclapply(unique(object1$alias), function(i) {
      lapply(unique(object2$alias), function(j) {
        length(intersect(subset(object1, alias == i)$barcode, subset(object2, alias == j)$barcode))
      })
    }, mc.cores = mc.cores))
    
    rownames(intersects) = unique(object1$alias)
    colnames(intersects) = unique(object2$alias)
    
    # Find correspondence
    if(length(unique(object1$alias)) > length(unique(object2$alias))){
      # column max
      key = data.frame(object1 = rownames(intersects)[apply(intersects, 2, which.max)], 
                       object2 = colnames(intersects))
    } else {
      # row max
      key = data.frame(object1 = rownames(intersects), 
                       object2 = colnames(intersects)[apply(intersects, 1, which.max)])
    }
    
    print(key)
    
    # Match the barcodes
    object2 = RenameCells(object2, new.names = paste0(convert_values(object2$alias, key %>% setNames(c('new.names', 'old.names'))), '_', object2$barcode, '-1'))
    
    # Ensure same format in first object
    object1 = RenameCells(object1, new.names = paste0(object1$alias, '_', object1$barcode, '-1'))
  }
  
  if(length(group.by) > 1){
    object1@meta.data = cbind(object1@meta.data[,!colnames(object1@meta.data) %in% group.by], object2@meta.data[match(colnames(object1), colnames(object2)), group.by])
  } else {
    if(is.null(as)) as = group.by
    object1@meta.data[[as]] = object2@meta.data[match(colnames(object1), colnames(object2)), group.by] 
  }
  
  object1
}

SmartMerge = function(matrixList, by = NULL, all = TRUE, common_keys = NULL){
  
  if(!all){
    # Find common
    common_genes = Reduce(intersect, lapply(matrixList, function(x) x[[by]]))
    
    # Subset each item to common
    matrixList = lapply(matrixList, function(item) item[item[[by]] %in% common_genes,])
    
    # Initialize merged matrix
    first = matrixList[[1]]
    merge = first[match(common_genes, first[[by]]),]
    for(i in 2:length(matrixList)){
      current = matrixList[[i]]
      new = current[match(common_genes, current[[by]]),]
      merge = cbind(merge, new) 
    }
  } else if(!is.null(common_keys)){

    # Initialize merged matrix
    merge = matrixList[[1]][match(common_keys, rownames(matrixList[[1]])), , drop = FALSE]
    rownames(merge) = common_keys
    # merge = first[match(common_keys, first[[by]]),]
    for(i in 2:length(matrixList)){
      current = matrixList[[i]]
      new = current[match(common_keys, rownames(current)), , drop = FALSE]
      merge = cbind(merge, new) 
    }
  } else {
    common_keys = Reduce(union, lapply(matrixList, function(x) rownames(x)))
    rownames(merge) = common_keys
    
    # Initialize merged matrix
    merge = matrixList[[1]][match(common_keys, rownames(matrixList[[1]])), , drop = FALSE]
    rownames(merge) = common_keys
    # merge = first[match(common_keys, first[[by]]),]
    for(i in 2:length(matrixList)){
      current = matrixList[[i]]
      new = current[match(common_keys, rownames(current)), , drop = FALSE]
      merge = cbind(merge, new) 
    }
  }
  
  merge
}

BigOrthologMatrix2 = function(species, orthologous_genes, ortho){
  
  message("Working on ", species, "!")
  
  # Read in species data
  if(species == "Mouse") {
    # speciesAC = readRDS("../../Species_Reference/MouseACref_v4.rds")
    # speciesAC = readRDS("../../../storage/HahnObjects/BCs/MouseBC_int_ann_v3.rds")
    speciesAC = readRDS("../../Species_Objects/Mouse_initial.rds")
  } else {
    speciesAC = readRDS(paste0("../../Species_Objects/", species, "_initial.rds"))
  }
  
  # Remove cells that are not in orthotype object
  cells.use = Cells(ortho)[Cells(ortho) %in% Cells(speciesAC)]
  message("Found ", length(cells.use), " cells from ortho object")
  speciesAC = speciesAC[,cells.use]
  
  # Revert gene names to original species names
  if(!is.null(speciesAC@misc$orig.features)) {
    rownames(speciesAC@assays$RNA@data) = speciesAC@misc$orig.features
    message("Reverting to original feature names for ", species)
  }
  
  # Convert gene symbols if not human
  if(!species %in% c("Human", "TreeShrew", "Rhabdomys")) {
    speciesAC = suppressMessages(ConvertGeneSymbols(speciesAC, "../../Orthology/martMergeRefHuman.rds", species))
  } else if(species == "Rhabdomys") {
    speciesAC = suppressMessages(ConvertGeneSymbols(speciesAC, "../../Orthology/martMergeRefHuman.rds", "Mouse"))
  }
  
  # Get count matrix
  count_matrix = speciesAC@assays$RNA@data
  
  subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthologous_genes),]
  message("Found ", nrow(subset_count_matrix), " genes out of ", length(orthologous_genes), " total orthologs. Efficiency = ", nrow(subset_count_matrix)/length(orthologous_genes))
  
  # Sort by orthology list order
  subset_count_matrix = suppressWarnings(as.matrix(subset_count_matrix))
  ortholog_count_matrix = subset_count_matrix[match(orthologous_genes, rownames(subset_count_matrix)),]
  rownames(ortholog_count_matrix) = orthologous_genes
  # ortholog_count_matrix[is.na(ortholog_count_matrix)] <- 0
  # ortholog_count_matrix[1:10,1:5]
  
  return(ortholog_count_matrix)
}

BigOrthologMatrix = function(species, orthologous_genes, mart_file = '../../Orthology/martMergeRefHuman_modified.rds', use.counts = FALSE){
  
  message("Working on ", species, "!")
  # Read in species data
  if(species == "Mouse") {
    speciesAC = readRDS("../../Species_Reference/MouseACref_v4.rds")
    # speciesAC = readRDS("../../Species_Reference/YanAC_v4.rds")
  } else {
    speciesAC = readRDS(paste0("../../Species_Objects/", species, "AC_v6.rds"))
  }
  
  # Revert names to original species names
  if(!is.null(speciesAC@misc$orig.features)) {
    rownames(speciesAC@assays$RNA@counts) = speciesAC@misc$orig.features # must be the counts since this is what convert ConvertGeneSymbols uses
    message("Reverting to original feature names for ", species)
  }
  
  # Convert gene symbols if not human
  if(!species %in% c("Human", "TreeShrew", "Rhabdomys")) {
    speciesAC = suppressMessages(ConvertGeneSymbols(speciesAC, mart_file, species))
  } else if(species == "Rhabdomys") {
    speciesAC = suppressMessages(ConvertGeneSymbols(speciesAC, mart_file, "Mouse"))
  }
  
  # Get count matrix
  if(use.counts) count_matrix = speciesAC@assays$RNA@counts else count_matrix = speciesAC@assays$RNA@data
  
  subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthologous_genes),]
  message("Found ", nrow(subset_count_matrix), " genes out of ", length(orthologous_genes), " total orthologs. Efficiency = ", nrow(subset_count_matrix)/length(orthologous_genes))
  
  # Sort by orthology list order
  subset_count_matrix = suppressWarnings(as.matrix(subset_count_matrix))
  ortholog_count_matrix = subset_count_matrix[match(orthologous_genes, rownames(subset_count_matrix)),]
  rownames(ortholog_count_matrix) = orthologous_genes
  # ortholog_count_matrix[is.na(ortholog_count_matrix)] <- 0
  ortholog_count_matrix[1:10,1:5]
  
  return(ortholog_count_matrix)
}

OrthologSeurat = function(speciesObject, orthology_key, common.genes = FALSE, mart_filepath = "../../Orthology/martMergeRefHuman.txt", reference_species = "Human"){
  
  species = GetSpecies(speciesObject)
  message("Working on ", species, "!")
  
  # Read in species data
  # if(species == "Mouse") {
  #   speciesObject = readRDS("../../Species_Reference/MouseACref_v4.rds")
  #   # speciesObject = subset(speciesObject, source == "Yan")
  # } else {
  #   speciesObject = readRDS(paste0(directory, species, "AC_v5.rds"))
  # }
  
  # Revert names to original species names
  if(!is.null(speciesObject@misc$orig.features)) {
    rownames(speciesObject@assays$RNA@counts) = speciesObject@misc$orig.features
    message("Reverting to original feature names for ", species)
  }
  
  # Convert gene symbols if not human
  if(species == "Rhabdomys") {
    # Using mouse as reference
    speciesObject = (ConvertGeneSymbols(speciesObject, orthology_key, "Mouse"))
  } else if(species == "Goldfish") {
    # Using zebrafish as reference
    speciesObject = (ConvertGeneSymbols(speciesObject, orthology_key, "Zebrafish"))
  } else if(!species %in% c(reference_species, "TreeShrew")) {
    speciesObject = (ConvertGeneSymbols(speciesObject, orthology_key, species))
  }
  
  # Get count matrix
  count_matrix = speciesObject@assays$RNA@counts
  
  # Now using same conversion function for consistency
  # # Convert some names to proper column name in orthology file
  # if(species == "Peromyscus") species = "Northern American deer mouse"
  # if(species == "Macaque") species = "Crab-eating macaque"
  # if(species == "Lizard") species = "Green anole"
  # if(species == "MouseLemur") species = "Mouse Lemur"
  # if(species == "Marmoset") species = "White-tufted-ear marmoset"
  # if(species == "Rhabdomys") species = "Mouse" # Rhabdomys assembly used mouse names
  # 
  # Subset to genes in orthology key
  # if(species == "Human") {
  #   subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthology_key$`Gene name`),]
  # } else {
  #   subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthology_key[[paste0(species," gene name")]]),]
  # }
  
  # Subset to genes in orthology key
  subset_count_matrix = count_matrix[toupper(rownames(count_matrix)) %in% toupper(orthology_key$`Gene name`),]
  n_genes_found = nrow(subset_count_matrix)
  n_genes = length(unique(orthology_key$`Gene name`))
  message("Found ", n_genes_found, " genes out of ", n_genes, " total orthologs. Efficiency = ", n_genes_found/n_genes)
  
  # Sort by orthology list order
  subset_count_matrix = suppressWarnings(as.matrix(subset_count_matrix))
  
  if(common.genes){
    ortholog_count_matrix = subset_count_matrix
  } else {
    # Now returning raw matrix with all genes, including genes that were not found in every species
    ortholog_count_matrix = subset_count_matrix[match(orthology_key$`Gene name`, rownames(subset_count_matrix)),]
    rownames(ortholog_count_matrix) = orthology_key$`Gene name`
    ortholog_count_matrix[is.na(ortholog_count_matrix)] <- 0
  }
  
  # Create new seurat object with metadata
  # speciesOrtho = CreateSeuratObject(subset_count_matrix)
  speciesOrtho = CreateSeuratObject(ortholog_count_matrix)
  speciesOrtho@meta.data = speciesObject@meta.data
  
  return(speciesOrtho)
}

# convertGenes <- function(genes, in_species, out_species = "Human"){
#   in_species_column=which(toupper(colnames(orthology_key))==paste0(species, " GENE NAME"))
#   orthologs = orthology_key[match(toupper(genes), toupper(orthology_key[,in_species_column])), paste0(out_species, " Gene Name")]
#   return(orthologs)
# }

ScoreCelltype = function(object, features){
  # Scale data if not already scaled for every feature
  # if(!all(rownames(object@assays$RNA@scale.data) == rownames(object))) object = ScaleData(object, features = rownames(object))
  features_in_object = intersect(features, rownames(object@assays$RNA@scale.data))
  message("Found ", length(features_in_object), " features in this object")
  raw_score = colSums(object@assays$RNA@scale.data[features_in_object,])
  normalized_score = raw_score/length(features_in_object)
  return(normalized_score)
}

AssignName = function(object, use = 'seurat_clusters', zero.index = FALSE){
  
  # Old function
  # speciesAC$type = paste0(sprintf('%02d', speciesAC$seurat_clusters), "_", speciesAC$lit_type) 
  # speciesAC$type[is.na(speciesAC$lit_type)] = sprintf('%02d', speciesAC$seurat_clusters[is.na(speciesAC$lit_type)])
  # speciesAC$type = factor(speciesAC$type, levels = 0:length(unique(speciesAC$type)))
  
  # object$type = paste0(str_pad(object$seurat_clusters, 2, side = "left"), "_", object$lit_type)
  # object$type[is.na(object$lit_type)] = str_pad(object$seurat_clusters[is.na(object$lit_type)], 2, side = "left")
  
  if(use == 'seurat_clusters' & zero.index) {
    object$clusters = as.numeric(as.character(object@meta.data[[use]]))
  } else if(use == 'seurat_clusters' & !zero.index) {
    object$clusters = as.numeric(as.character(object@meta.data[[use]]))+1
  } else {
    object$clusters = object@meta.data[[use]]
  }
  object$type = paste0(object$clusters, "_", object$lit_type)
  object$type[is.na(object$lit_type)] = as.character(object$clusters[is.na(object$lit_type)])
  meta = Metadata(object, feature.1 = "clusters", feature.2 = "lit_type")
  object$type = factor(object$type, levels = ifelse(is.na(meta$lit_type), as.character(meta$clusters), paste0(meta$clusters, "_", meta$lit_type)))
  object$type_no = object$clusters
  Idents(object) = "type"
  return(object)
}

nSpecies = function(object, group.by = "species"){
  return(length(unique(object@meta.data[[group.by]])))
}

nClusters = function(object, clusters = "seurat_clusters"){
  return(nrow(unique(object[[clusters]])))
}

GADPlot = function(df, Gaba.gene, Gly.gene = "SLC6A9", cutoff.Gly = NULL, cutoff.GABA = NULL, zero.base.clusters = TRUE){
  
  if(!zero.base.clusters){
    df$cluster = as.numeric(as.character(df$cluster))+1
  }
  
  ggplot(df, 
         aes(x = eval(parse(text = Gaba.gene)), 
             y = eval(parse(text = Gly.gene)), 
             label = cluster)) + 
    geom_point()+
    xlab(Gaba.gene)+
    ylab(Gly.gene)+
    geom_hline(yintercept = cutoff.Gly, linetype = "dashed", color = "grey")+
    geom_vline(xintercept = cutoff.GABA, linetype = "dashed", color = "grey")+
    theme_dario() +
    scale_x_continuous(expand = expansion(mult = c(0,0)))+
    scale_y_continuous(expand = expansion(mult = c(0,0)))+
    geom_text_repel(min.segment.length = 0.05)
  
}

SubsampleSeurat = function(object, size, seed = 12345){
  set.seed(seed)
  object <- object[, sample(colnames(object), size = size, replace=F)]
  return(object)
}

ExpressionScatter2 = function(data, x, y, 
                             xlab = NULL, 
                             ylab = NULL, 
                             title = NULL, 
                             genes.to.label = NULL, 
                             logX = FALSE, 
                             logY = FALSE, 
                             label.size = 4,
                             de.genes = FALSE,
                             top.n.genes = 50,
                             label.x.npc = 0,
                             label.y.npc = 0.9,
                             assay = 'RNA',
                             ...){
  
  args = list(...)
  
  if(inherits(data, "Seurat")) data = LogAvgExpr(data, assay = assay, ...)
  
  # Set all other genes to NA to avoid labeling
  data$gene = rownames(data)
  data$gene[!data$gene %in% c(genes.to.label)] = NA
  data$gene.italic = paste0("italic(", data$gene, ")")
  data$gene.italic[!data$gene %in% c(genes.to.label)] = NA
  
  if(de.genes){
    de_res = TopNDEGs(data[,data@meta.data[,args$group.by] %in% c(x, y)], group.by = args$group.by, n = top.n.genes, avg_log2FC_cutoff = 0.25)
    genes.to.label = de_res$gene
  }
  
  # Set labeled for coloring
  # data$labeled = ifelse(!is.na(data$gene), TRUE, FALSE)
  data$labeled = FALSE
  data$labeled[!is.na(data$gene)] = TRUE
  
  data$group = TRUE
  
  plt = ggplot(data, aes_(x=as.name(x), y=as.name(y)))+
    ggrastr::rasterise(geom_point(data = subset(data, !labeled), color = "grey"), dpi = 300)+
    geom_point(data = subset(data, labeled), color = "black", shape = 21, fill = 'red')+
    # scale_color_manual(values = c("grey", "red"))+
    geom_smooth(formula = y ~ x, method = "lm", se = FALSE, linetype = "dashed", color = "black") +
    {if(!is.null(genes.to.label)) geom_text_repel(data = subset(data, !is.na(gene)),
                                                  aes(label = gene.italic),
                                                  max.overlaps = Inf, 
                                                  size = label.size, 
                                                  parse = TRUE)}+
    # labs(title = title,
    #       x = xlab,
    #      {if(!is.null(ylab)) y = ylab})+
    ggtitle(title)+
    {if(!is.null(xlab)) xlab(xlab)}+
    {if(!is.null(ylab)) ylab(ylab)}+
    stat_cor(method = "pearson", aes(label = ..r.label..), label.x.npc = label.x.npc, label.y.npc = label.y.npc)+
    {if(logX) scale_x_continuous(trans='log10')} +
    {if(logY) scale_y_continuous(trans='log10')} +
    theme_cowplot2()+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
    ArialFont()
  
  return(plt)
}

ExpressionScatter = function(data, x, y, 
                             xlab = NULL, 
                             ylab = NULL, 
                             title = NULL, 
                             genes.to.label = NULL, 
                             logX = FALSE, 
                             logY = FALSE, 
                             label.size = 4,
                             de.genes = FALSE,
                             top.n.genes = 50,
                             ...){
  
  args = list(...)
  
  if(inherits(data, "Seurat")) data = LogAvgExpr(data, assay = "RNA", ...)
  
  # Set all other genes to NA to avoid labeling
  data$gene = rownames(data)
  data$gene[!data$gene %in% c(genes.to.label)] = NA
  
  if(de.genes){
    de_res = TopNDEGs(data[,data@meta.data[,args$group.by] %in% c(x, y)], group.by = args$group.by, n = top.n.genes, avg_log2FC_cutoff = 0.25)
    genes.to.label = de_res$gene
  }
  
  # Set labeled for coloring
  # data$labeled = ifelse(!is.na(data$gene), TRUE, FALSE)
  data$labeled = FALSE
  data$labeled[!is.na(data$gene)] = TRUE
  
  data$group = TRUE
  
  plt = ggplot(data, aes_(x=as.name(x), y=as.name(y), label = as.name("gene")))+#color = labeled, group = group))+
    ggrastr::rasterise(geom_point(data = subset(data, !labeled), color = "grey"), dpi = 300)+
    geom_point(data = subset(data, labeled), color = "black", shape = 21, fill = 'red')+
    # scale_color_manual(values = c("grey", "red"))+
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
    {if(!is.null(genes.to.label)) geom_text_repel(max.overlaps = Inf, size = label.size)}+
    # labs(title = title,
    #       x = xlab,
    #      {if(!is.null(ylab)) y = ylab})+
    ggtitle(title)+
    {if(!is.null(xlab)) xlab(xlab)}+
    {if(!is.null(ylab)) ylab(ylab)}+
    stat_cor(method = "pearson", aes(label = ..r.label..))+
    {if(logX) scale_x_continuous(trans='log10')} +
    {if(logY) scale_y_continuous(trans='log10')} +
    theme_cowplot2()+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  
  return(plt)
}

#' @example
#' ScatterPlot(nTypes.df, 
#' x = "nCluster.residual", 
#' y = "nRgcCluster.residual", 
#' label = "species", 
#' ylab = "# of RGC clusters (residual)", 
#' xlab = "# of AC clusters (residual)", 
#' label.x.npc = 0.7, 
#' label.y.npc = 0.1)
ScatterPlot = function(data, x, y, 
                       xlab = NULL, 
                       ylab = NULL, 
                       title = NULL, 
                       labels = NULL, 
                       logX = FALSE, 
                       logY = FALSE, 
                       # cor.position = c("top", "left"),
                       fill = 'grey',
                       lm = FALSE, 
                       y_equals_x = FALSE, 
                       r = FALSE, 
                       r_pval = FALSE, 
                       r_break_pval = FALSE,
                       label_func = geom_text_repel, 
                       min.segment.length = 0.05, 
                       label.size = 3,
                       max.overlaps = 1, 
                       lm.color = 'grey',
                       pt.size = 1, 
                       pt.shape = 21,
                       pt.col = 'black',
                       label.x.npc = 'right', 
                       label.y.npc = 0.01,
                       ...){
  
  if(!is.null(labels)) {
    data$names = data[[labels]]
  } else {
    data$names = ""
  }
  
  plt = ggplot(data, aes(x = eval(parse(text = x)), 
                         y = eval(parse(text = y)), 
                         label = names))+
    {if(lm) geom_smooth(method = "lm", se = FALSE, color = lm.color, linetype = "dashed")} +
    {if(fill %in% colnames(data)) geom_point(aes(fill = eval(parse(text = fill))), 
                                             color = pt.col, 
                                             shape = pt.shape, 
                                             size = pt.size) 
      else geom_point(fill = fill, color = 'black', shape = 21)}+
    labs(title = title, 
         x = ifelse(is.null(xlab), x, xlab), 
         y = ifelse(is.null(ylab), y, ylab))+
    {if(fill %in% colnames(data)) labs(fill = NULL)}+
    {if(r) stat_cor(method = "pearson", aes(label = ..r.label..), label.x.npc = label.x.npc, label.y.npc = label.y.npc, hjust = 1)}+
    {if(r_pval) stat_cor(method = "pearson", label.x.npc = label.x.npc, label.y.npc = label.y.npc, hjust = 1)}+
    {if(r_break_pval) stat_cor(
      method = "pearson",
      aes(label = paste("R =", ..r.., "\np =", ..p.., sep = " ")),
      # aes(label = paste(..r.label.., ..p.label.., sep = "~`\n`~")),
      # aes(label = paste("R =", ..r.label.., "\np =", after_stat(p), sep = " ")),
      # aes(label = sprintf("%s\n%s", after_stat(r.label), after_stat(p.label))),
      label.x.npc = label.x.npc,
      label.y.npc = label.y.npc,
      hjust = 1
    )}+
    {if(y_equals_x) geom_abline(slope = 1, intercept = 0)}+
    {if(logX) scale_x_continuous(trans='log10')} +
    {if(logY) scale_y_continuous(trans='log10')} +
    {if(!is.null(labels)) label_func(max.overlaps = max.overlaps, 
                                     min.segment.length = min.segment.length, 
                                     size = label.size, 
                                     ...)}+
    theme(plot.title = element_text(hjust = 0.5))+
    ArialFont()+
    theme_dario()
  
  return(plt)
}

JSMatrix = function(matrix, union.threshold = 0){
  rowSums = rowSums(matrix)
  colSums = colSums(matrix)
  
  js_matrix = matrix
  
  for(row in seq_along(rowSums)){
    for(column in seq_along(colSums)){
      intersection = matrix[row,column]
      union = rowSums[row] + colSums[column] - intersection
      if(union == 0) {
        js_matrix[row,column] = 0
      } else {
        js_matrix[row,column] = intersection / union
      }
      
      # If the union is too small, impute as zero (no trustworthy overlap)
      if(union < union.threshold) js_matrix[row,column] = 0
    }
  }
  
  return(js_matrix)
}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

JSHeatmap = function(matrix, title = NULL, stagger.threshold = 0.1, heatmap = TRUE, 
                     col.low = "white", col.high = "#584B9FFF", max.value = NULL,
                     legend_name = NULL, row.order = NULL, column.order = NULL, 
                     xlab = NULL, ylab = NULL, border.col = 'black',
                     split = FALSE, label = FALSE){
  
  melted = reshape2::melt(matrix)
  colnames(melted) = c("row", "col", "jaccard")
  
  # Diagonalization
  row.max = apply(matrix, 1, which.max)
  match.df = data.frame(row = 1:length(row.max), 
                        col = as.integer(row.max))
  match.df$value = matrix[as.matrix(match.df)]
  match.df.sorted = match.df %>% arrange(-value)
  
  orig.rownames = rownames(matrix)
  orig.colnames = colnames(matrix)
  
  # If row.order is supplied, reorder the match.df.sorted
  if(!is.null(row.order)) match.df.sorted = match.df.sorted %>% arrange(factor(row, levels = match.df$row[match(row.order, orig.rownames)]))
  
  # Sorting
  melted.sorted = melted %>% arrange(factor(row, levels = unique(orig.rownames[match.df.sorted$row])), 
                                     factor(col, levels = unique(orig.colnames[match.df.sorted$col])))
  melted.sorted.nolow = subset(melted.sorted, jaccard > stagger.threshold)
  
  if(is.null(row.order)) {
    melted$row = factor(melted$row, levels = rev(unique(c(as.character(melted.sorted.nolow$row), setdiff(orig.rownames, melted.sorted.nolow$row)))))
  } else {
    melted$row = factor(melted$row, levels = rev(row.order))
  }
  
  if(is.null(column.order)) {
    melted$col = factor(melted$col, levels = (unique(c(as.character(melted.sorted.nolow$col), setdiff(orig.colnames, melted.sorted.nolow$col)))))
  } else {
    melted$col = factor(melted$col, levels = (column.order))
  }
  
  if(split) melted$species = str_split_fixed(melted$col, "_", 2)[,1]
  
  # Handle naming 
  if(is.null(legend_name)){
    # if(heatmap) legend_name = 'Jaccard\nindex' else legend_name = 'Jaccard\nindex'
    legend_name = 'Jaccard\nindex'
  }
  
  if(is.null(max.value)){
    if(heatmap) max.value = 0.5 else max.value = 1
  }
  
  plt = ggplot(melted, aes(x = col, y = row))+
    {if(heatmap) geom_tile(aes(fill = jaccard), color=border.col) else geom_point(aes(colour = jaccard,  size=jaccard))}+
    {if(heatmap) scale_fill_gradient(legend_name, low=col.low, high = col.high, limits=c(0, max.value), na.value = col.high) else scale_color_gradient(legend_name, low=col.low, high = col.high, limits=c(0, max.value), na.value = col.high)} +
    # {{if(!heatmap) scale_size(range = c(1, max.size), limits = c(0, max.perc))+   
    guides(fill = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    theme_bw() +
    {if(label) geom_text(aes(label = round(jaccard, 2)), color = "black")}+
    {if(heatmap) scale_y_discrete(expand = c(0,0))}+
    {if(heatmap) scale_x_discrete(expand = c(0,0))}+
    {if(split) facet_wrap(~species)}+
    {if(!heatmap) labs(size = legend_name)}+
    {if(!heatmap) scale_radius(range = c(0.3,3))}+
    RotatedAxis()+ 
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(title)+
    theme(axis.title = element_text(color = 'black'), 
          axis.text = element_text(color = 'black'),
          plot.title = element_text(hjust = 0.5))
  return(plt)
}

ConvertGeneSymbols = function(object, mart = '../../Orthology/martMergeRefHuman.rds', species, make.unique = TRUE){
  
  if(inherits(mart, "character")){
    mart = readRDS(mart)
  }
  
  if(species == "Peromyscus") species = "Northern American deer mouse"
  if(species == "Macaque") species = "Crab-eating macaque"
  if(species == "Lizard") species = "Brown anole" #species = "Green anole"
  if(species == "MouseLemur") species = "Mouse Lemur"
  if(species == "Marmoset") species = "White-tufted-ear marmoset"
  
  if(species == "Lamprey"){
    lamprey_annotations = read.csv("../../../storage/data/GeneAnotation_Lamprey.xlsx - Sheet1.csv")
    lamprey_annotations$Symbol = gsub(" ", "", lamprey_annotations$Symbol)
    
    # Custom edits
    new_genes = rownames(object@assays$RNA@counts)
    # new_genes[new_genes == "LOC103091742"] = "TFAP2B_X2"
    # new_genes[new_genes == "MSTRG.3969"] = "TFAP2B_X3"
    
    # Convert from Gene first
    sorted_lamprey_annotations = lamprey_annotations[match(rownames(object@assays$RNA@counts), lamprey_annotations$Gene),]
    new_genes[new_genes %in% sorted_lamprey_annotations$Gene] = sorted_lamprey_annotations$Symbol[new_genes %in% sorted_lamprey_annotations$Gene]
    
    # Convert from Stringtie second
    sorted_lamprey_annotations = lamprey_annotations[match(rownames(object@assays$RNA@counts), lamprey_annotations$StringTie),]
    new_genes[new_genes %in% sorted_lamprey_annotations$StringTie] = sorted_lamprey_annotations$Symbol[new_genes %in% sorted_lamprey_annotations$StringTie]
    
    # Return blanks to original 
    new_genes[new_genes == ""] = rownames(object@assays$RNA@counts)[new_genes == ""]
    
    # Custom edits
    new_genes[new_genes == "SLC32A1/VIAAT1"] = "SLC32A1"
    new_genes[new_genes == "SLC6A9/GlyT1"] = "SLC6A9"
    new_genes[new_genes == "SLC6A1/GAT1"] = "SLC6A1"
    new_genes = make.unique(new_genes)
    
    # Print some examples of converted genes
    df = data.frame(before = rownames(object@assays$RNA@counts), after = new_genes)
    message("Converted ", nrow(df[df$before != df$after,]), " genes! Some examples below:")
    message(paste0(capture.output(head(df[df$before != df$after,], 200)), collapse = "\n"))
    
    # Add to object
    object@misc$orig.features = rownames(object@assays$RNA@counts)
    rownames(object@assays$RNA@counts) = new_genes
    rownames(object@assays$RNA@data) = new_genes 
  } else if(species == "TreeShrew"){
    TOGA_ts = read.table("../../Orthology/TOGA/TreeShrew.tsv", header = T)
    TOGA_ts$t_gene_name = str_split_fixed(TOGA_ts$t_transcript, "\\.", 2)[,2]
    TOGA_ts$q_gene_name = str_split_fixed(TOGA_ts$q_transcript, "\\.", 3)[,2]
    TOGA_ts_one2one = subset(TOGA_ts, orthology_class == "one2one")
    TOGA_ts_filt = unique(TOGA_ts_one2one[,c("q_gene","orthology_class","t_gene_name","q_gene_name")])
    new_genes = TOGA_ts_filt$t_gene_name[match(rownames(object), gsub("_", "-", TOGA_ts_filt$q_gene))]
    new_genes[is.na(new_genes)] = rownames(object)[is.na(new_genes)]
    rownames(object@assays$RNA@counts) = make.unique(new_genes)
    rownames(object@assays$RNA@data) = make.unique(new_genes)
  } else {
    # Pull genes to convert from object counts matrix
    before_genes = rownames(object@assays$RNA@counts)
    
    # Read in ensembl biomart orthology table
    # mart = as.data.frame(fread(mart_filepath))
    
    # Get one-to-one orthologs between species and mouse
    conversion_key = mart[mart[["Gene name"]] != "" & 
                            # mart[[paste0(species, " gene name")]] != "" & 
                            mart[[paste0(species, " homology type")]] == "ortholog_one2one",
                          c("Gene name", paste0(species, " gene stable ID"), paste0(species, " gene name"))] %>% unique
    
    # If the gene name is blank or not present in the object, try the ensembl id instead
    conversion_key[(conversion_key[[paste0(species, " gene name")]] == "" | !conversion_key[[paste0(species, " gene name")]] %in% before_genes), paste0(species, " gene name")] = 
      conversion_key[(conversion_key[[paste0(species, " gene name")]] == "" | !conversion_key[[paste0(species, " gene name")]] %in% before_genes), paste0(species, " gene stable ID")]
    
    # Sort genes by rownames of seurat object
    sorted_key = conversion_key[match(before_genes, conversion_key[[paste0(species, " gene name")]]),]
    
    # Convert seurat object names to reference (human or mouse) names if they are in the conversion key
    converted_genes = before_genes
    converted_genes[converted_genes %in% sorted_key[[paste0(species, " gene name")]] ] = toupper(sorted_key$`Gene name`[converted_genes %in% sorted_key[[paste0(species, " gene name")]] ])
    
    # In the rare case that we introduce duplicates, make unique to avoid gene merging and keep downstream analyses consistent
    if(make.unique) converted_genes = make.unique(converted_genes)
    
    # Change underscores to dashes
    converted_genes = gsub("_", "-", converted_genes)
    
    # Make conversion df
    before_after = data.frame(before = before_genes, after = converted_genes)
    
    # Print some examples of converted genes
    message("Converted ", nrow(before_after[converted_genes != before_genes,]), " genes! Some examples below:")
    message(paste0(capture.output(head(before_after[converted_genes != before_genes,], 200)), collapse = "\n"))
    
    # Save old symbols in the object
    object@misc$orig.features = before_genes
    
    # Change seurat object rownames
    rownames(object@assays$RNA@counts) <- converted_genes
    rownames(object@assays$RNA@data) <- converted_genes
    rownames(object@assays$RNA@scale.data) <- converted_genes[match(rownames(object@assays$RNA@scale.data), before_genes)]
  }
  return(object)
}

ProfileContamination = function(object, z.score.threshold = 3){
  object[["percent.rod"]] <- PercentageFeatureSet(object, features = intersect(Rod_markers, rownames(object)))
  object[["percent.cone"]] <- PercentageFeatureSet(object, features = intersect(Cone_markers, rownames(object)))
  object[["percent.hc"]] <- PercentageFeatureSet(object, features = intersect(HC_markers, rownames(object)))
  object[["percent.bc"]] <- PercentageFeatureSet(object, features = intersect(setdiff(BC_markers, "PRKCA"), rownames(object)))
  object[["percent.ac"]] <- PercentageFeatureSet(object, features = intersect(AC_markers, rownames(object)))
  object[["percent.rgc"]] <- PercentageFeatureSet(object, features = intersect(setdiff(RGC_markers, c("POU6F2", "RBFOX3")), rownames(object)))
  object[["percent.mg"]] <- PercentageFeatureSet(object, features = intersect(MG_markers, rownames(object)))
  
  # rod.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.rod")) > 2) - 1
  # cone.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.cone")) > 2) - 1
  # hc.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.hc")) > 2) - 1
  # bc.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.bc")) > 2) - 1
  # ac.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.ac")) > 2) - 1
  # rgc.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.rgc")) > 2) - 1
  # mg.contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = "percent.mg")) > 2) - 1
  
  cell_classes = c("rod", "cone", "hc", "bc", "ac", "rgc", "mg")
  
  contamin.clusters = lapply(cell_classes, function(class) {
    contamin = which(scale(MeanMetadata(object, feature.1 = "seurat_clusters", feature.2 = paste0("percent.", class))) > z.score.threshold) - 1
    return(contamin)
  })
  
  names(contamin.clusters) = cell_classes
  object@misc$contamin.clusters = contamin.clusters
  return(object)
}

TrainXGBoost = function(train, test, train.clusters = "seurat_clusters", nfeatures = 2000){
  
  Idents(train) = train.clusters
  # Idents(test) = test.clusters
  
  DefaultAssay(train) = "RNA"
  DefaultAssay(test) = "RNA"
  
  train = FindVariableFeatures(train, selection.method = "vst", nfeatures = nfeatures)
  test = FindVariableFeatures(test, selection.method = "vst", nfeatures = nfeatures)
  common_HVGs <- intersect(VariableFeatures(train), VariableFeatures(test))
  message("Using ", length(common_HVGs), " common highly variable genes...\n")
  
  AC_model <- TrainModel(train, training_genes = common_HVGs, train_ident = train.clusters, do.scale = TRUE)
  
  return(AC_model)
}

LIBRARIES = c('paletteer', 'tidyverse', 'Seurat', 'ggplot2', 'reshape2', 
              'dplyr', 'xgboost', 'cowplot', 'pdfCluster', 'ggrepel', 
              'presto', 'scales', 'harmony', 'RColorBrewer', 'WGCNA', 
              'openxlsx', 'lisi', 'ggpubr', 'doParallel', 'circlize', 
              'viridis', 'ComplexHeatmap', 'Polychrome', 
              'HGNChelper', 'openxlsx', 'Matrix', 'patchwork', 'colorspace', 
              'ape', 'dendextend', 'qs2')

LoadLibraries = function(load.lisi = TRUE, verbose = FALSE){

  res = lapply(LIBRARIES, require, character.only = TRUE)
  names(res) = LIBRARIES
  
  if(verbose) return(res)
}

SourceFiles = function(path = '../../utils/', objects = TRUE, AC = FALSE){
  filenames = c('xgboost_train.R', 'utilFxns.R', 'xgboost_train.R', 'plottingFxns.R',
                'dario_functions.R', 'wrappers.R', 
                'objects.R', 'objects_AC.R', 'xgboost_train_DT.R', 
                'SmartMatrix.R', 'SeuratV5_functions.R')
  
  # Remove files? 
  if(!objects) filenames = setdiff(filenames, 'objects.R')
  if(!AC) filenames = setdiff(filenames, 'objects_AC.R')
  
  lapply(filenames, function(filename) source(paste0(path, filename)))
  
  invisible(NULL)
}

DoubletAnalysis = function(object, group.by = "seurat_clusters"){
  plot_grid(
    plot_grid(
      VlnPlot(object, features = "nFeature_RNA", group.by = "DF.classifications", pt.size = 0), 
      DimPlot(object, group.by = "DF.classifications") + NoLegend(), 
      DimPlot(object, group.by = group.by, label = T) + NoLegend(), 
      ncol = 3, nrow = 1, rel_widths = c(1, 2, 2)), 
    stackedBarGraph(object, feature.1 = "DF.classifications", feature.2 = group.by) + RotatedAxis(), 
    nrow = 2)
}

BrowseSeurat = function(object, batch = "animal", group.by = 'seurat_clusters'){
  if(is.null(object@meta.data[["percent.mt"]])) object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  plot = plot_grid(DimPlot(object, group.by = group.by, label = TRUE, shuffle = T, repel = T) + NoLegend(), 
                   ClusterBatchPlot(object, batch = batch, group.by = group.by, shuffle = TRUE) + NoLegend(), 
                   ClusterFeaturePlot(object, features = "nFeature_RNA", group.by = group.by) + NoLegend(),
                   VlnPlot(object, "nCount_RNA", group.by = group.by, pt.size = 0) + RotatedAxis() + NoLegend(), 
                   VlnPlot(object, "nFeature_RNA", group.by = group.by, pt.size = 0) + RotatedAxis() + NoLegend(), 
                   VlnPlot(object, "percent.mt", group.by = group.by, pt.size = 0) + RotatedAxis() + NoLegend(),
                   nrow = 2, ncol = 3, labels = LETTERS)
  return(plot)
}

BrowseSeuratSmall = function(object, batch = "orig.file"){
  plot = plot_grid(theme_umap(DimPlot(object, label = TRUE, group.by = "seurat_clusters", shuffle = TRUE)) + NoLegend(), 
                   theme_umap(ClusterBatchPlot(object, batch = batch, shuffle = TRUE)) + NoLegend(), 
                   theme_umap(ClusterFeaturePlot(object, features = "nFeature_RNA")) + NoLegend(),
                   nrow = 1, ncol = 3)
  return(plot)
}

ReprocessIntegrated = function(object, 
                               nPCs = 20, 
                               cluster_resolution = 1.5, 
                               nfeatures = 2000, 
                               selection.method = "vst", 
                               method = "seurat", 
                               k.param = 20, 
                               recompute.var.genes = FALSE, 
                               run.umap = TRUE, 
                               verbose = TRUE){
  if(method == 'seurat'){
    
    # Save old clusters
    object$old_seurat_clusters = object$seurat_clusters
    
    # Re-process using integrated assay
    DefaultAssay(object) = "integrated"
    if(recompute.var.genes) object = FindVariableFeatures(object, selection.method = selection.method, nfeatures = nfeatures)
    object <- ScaleData(object) %>% # vars.to.regress = "nCount_RNA" in ScaleData worsens batch correction somehow
      RunPCA() %>% 
      FindNeighbors(dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution, verbose = verbose)
    
    if(run.umap) object = RunUMAP(object, dims = 1:nPCs)
    
    DefaultAssay(object) = "RNA"
  } else if(method == "harmony"){
    
    # Save old clusters
    object$old_seurat_clusters = object$seurat_clusters
    
    object <- object %>%
      FindNeighbors(reduction = "harmony", dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution) %>%
      identity()
    if(run.umap) object = RunUMAP(object, reduction = "harmony", dims = 1:nPCs)
  } else {
    stop("method should be either 'seurat' or 'harmony'")
  }
  
  return(object)
}

ClusterBatchPlot <- function(object, batch = "orig.file", group.by = "seurat_clusters", ...){
  plot = DimPlot(object, group.by = batch, ...)
  plot$data$seurat_clusters <- as.factor(object@meta.data[[group.by]][match(rownames(plot$data), rownames(object@meta.data))])
  return(LabelClusters(plot, id = "seurat_clusters"))
}

ClusterFeaturePlot <- function(object, group.by = "seurat_clusters", repel = FALSE, ...){
  # args = list(...)
  plot = FeaturePlot(object, ...)
  plot$data$seurat_clusters <- as.factor(object@meta.data[[group.by]])
  return(LabelClusters(plot, id = "seurat_clusters", repel = repel))
}

plotBatches2 = function(SeuratObject, batch = "orig.file", clusters = "seurat_clusters"){
  nGroups = length(unique(SeuratObject@meta.data[,batch]))
  p = plot_grid(
    ClusterBatchPlot(SeuratObject, batch = batch, group.by = clusters, shuffle = TRUE) + NoLegend(),
    plot_grid(
      DimPlot(SeuratObject, reduction = "umap", pt.size = .1, group.by = batch, split.by = batch) + NoLegend(), 
      stackedBarGraph(SeuratObject, feature.1 = batch, feature.2 = clusters) + RotatedAxis() + NoLegend(), 
      nrow = 2, 
      rel_heights = c(1,1.5)
    ),
    ncol = 2, 
    rel_widths = c(1, (1.5 + 0.1*nGroups))
  )
  return(p)
}

plotBatches = function(SeuratObject, batch = "orig.file", clusters = "seurat_clusters"){
  p = plot_grid(
    plot_grid(
      DimPlot(SeuratObject, reduction = "umap", pt.size = .1, group.by = batch) + NoLegend(),
      DimPlot(SeuratObject, reduction = "umap", pt.size = .1, group.by = batch, split.by = batch) + NoLegend(), 
      ncol = 2, 
      rel_widths = c(1.2, nrow(unique(SeuratObject[[batch]])))
    ),
    stackedBarGraph(SeuratObject, feature.1 = batch, feature.2 = clusters) + RotatedAxis(), 
    nrow = 2
  )
  return(p)
}

FindDimensions = function(object, nPCs = 50){
    
  DefaultAssay(object) = "RNA"
  
  # Pre-processing
  object <- object %>% 
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(npcs = nPCs, verbose = FALSE)
  
  ElbowPlot(object, ndims = nPCs)
}

Harmonize = function(SeuratObject, batch = "orig.file", nPCs = 20, cluster_resolution = 0.5, 
                     k.param = 20, show.plots = FALSE, run.umap = TRUE, save.clusters = FALSE, 
                     lambda = NULL, remove.batches = 60, verbose = FALSE){
  
  # Remove tiny batches to avoid singular matrices
  SeuratObject@meta.data[,batch] = as.character(SeuratObject@meta.data[,batch])
  batches.remove = names(which(table(SeuratObject@meta.data[,batch]) < remove.batches))
  message("Removing batches: ", paste0(batches.remove, collapse = ", "))
  SeuratObject = SeuratObject[,!SeuratObject@meta.data[,batch] %in% batches.remove]
  
  # Save old clusters
  if(save.clusters) SeuratObject$old_seurat_clusters = SeuratObject$seurat_clusters
  
  DefaultAssay(SeuratObject) = "RNA"
  
  # Pre-processing
  SeuratObject <- SeuratObject %>% 
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(npcs = nPCs, verbose = FALSE)
  
  # Run harmony
  harmony <- SeuratObject %>% RunHarmony(batch, 
                                         plot_convergence = {if(show.plots) TRUE else FALSE}, 
                                         lambda = lambda)
  
  # Check PCA for mitigation of batche effects
  if(show.plots){
    print(plot_grid(
      TitlePlot(DimPlot(SeuratObject, reduction = "pca", pt.size = .1, group.by = batch) + NoLegend(), 'Before Harmony'),
      TitlePlot(DimPlot(harmony, reduction = "harmony", pt.size = .1, group.by = batch) + NoLegend(), 'After Harmony')
    ))
  }
  
  # Downstream analysis
  if(!run.umap){
    harmony <- harmony %>%
      FindNeighbors(reduction = "harmony", dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution, verbose = verbose) %>%
      identity()
  } else {
    harmony <- harmony %>%
      RunUMAP(reduction = "harmony", dims = 1:nPCs) %>%
      FindNeighbors(reduction = "harmony", dims = 1:nPCs, k.param = k.param) %>%
      FindClusters(resolution = cluster_resolution, verbose = verbose) %>%
      identity()
  }
  
  return(OneBasedClusters(harmony))
}

OneBasedClusters = function(object, group.by = 'seurat_clusters'){
  object@meta.data[[group.by]] <- as.factor(as.integer(object@meta.data[[group.by]]))
  Idents(object) <- object@meta.data[[group.by]]
  object
}

RobustResolution = function(object, from = 0.5, to = 1.5, by = 0.1, mc.cores = 10){
  
  resolutions = seq(from, to, by = by)
  
  library(parallel)
  clusterings = mclapply(resolutions, function(res){
    obj <- FindClusters(object, resolution = res)
    return(obj$seurat_clusters)
  }, mc.cores = mc.cores)
  
  names(clusterings) = resolutions
  
  # Make a data frame
  # clusterings_df = do.call(cbind, clusterings)
  Iterate(clusterings, adj.rand.index)
  
}

plotClusterRes = function(SeuratObject, dims = 1:20, mc.cores = 10, from = 0.1, to = 2, by = 0.05){
  
  # Vary resolution parameter from 0.4 to 2.0
  resolutions = seq(from, to, by = by)
  library(parallel)
  nclusters = mclapply(resolutions, function(res){
    SeuratObject <- FindClusters(SeuratObject, resolution = res)
    return(length(unique(SeuratObject$seurat_clusters)))
  }, mc.cores = mc.cores) %>% unlist
  
  res.cluster.table = data.frame(resolution = resolutions,
                                 nclusters = nclusters)
  
  p = ggplot(res.cluster.table, aes(x = resolution, y = nclusters))+
    geom_line()+ 
    ylab("Number of clusters")+
    theme_bw()
  
  print(res.cluster.table)
  
  return(p)
}

proportionTable = function(SeuratObject, feature.1, feature.2){
  freq = table(SeuratObject@meta.data[[feature.1]], SeuratObject@meta.data[[feature.2]]) 
  prop = t(freq)/colSums(freq)
  return(prop)
}

#' Find doublets
#' 
#' Function that merges similar clusters on tree if they have fewer than X number of DEGs.
#' Note that the BuildClusterTree function should have been called for this object
#' 
#' @param SeuratObject Object of class Seurat
#' @param channels the identifier corresponding to the different 10X channels each cell was in
#' @param annotations cell type annotations for homotypic doublet estimation
#' @param suffix suffix in colnames(SeuratObject) that is often added automatically, e.g. ".1"
#' @param num.cores number of cores, default is 1
FindDoublets = function(SeuratObject, 
                        channels = "orig.file", 
                        annotations = "seurat_clusters", 
                        suffix = "", 
                        num.cores = 1, #length(unique(SeuratObject$orig.file)), 
                        classify.by = NULL, 
                        nPCs = 10){
  
  if("DF.classifications" %in% colnames(SeuratObject@meta.data)) SeuratObject$DF.classifications = NULL
  
  # Split object by sequencing channel since doublets can only arise between cells in the same channel
  obj.list <- SplitObject(SeuratObject, split.by = channels)
  
  # Run DoubletFinder on each sample separately
  library(DoubletFinder)
  obj.list = lapply(obj.list, function(object){
    
    # print(head(as.character(object@meta.data[[classify.by]])))
    
    # Pre-processing
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
    object <- ScaleData(object)
    object <- RunPCA(object, nPCs = nPCs)
    object <- RunUMAP(object, dims = 1:nPCs)
    
    # pK Identification (no ground-truth)
    sweep.res.list <- paramSweep_v3(object, PCs = 1:nPCs, sct = FALSE, num.cores = 1)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    # BC maximization
    max_pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    message("Found max pK of: ", max_pK, "\n")
    
    # Homotypic Doublet Proportion Estimate
    homotypic.prop <- modelHomotypic(object@meta.data[[annotations]]) ## ex: annotations <- object@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(object@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    # Run DoubletFinder with varying classification stringencies
    # object <- doubletFinder_v3(object, PCs = 1:20, pN = 0.25, pK = max_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    if(!is.null(classify.by)){
      object <- doubletFinder_v3(object, 
                                 PCs = 1:nPCs, 
                                 pN = 0.25, 
                                 pK = max_pK, 
                                 nExp = nExp_poi.adj, 
                                 reuse.pANN = FALSE, 
                                 sct = FALSE, 
                                 annotations = as.character(object@meta.data[[classify.by]]))
    }
    else {
      object <- doubletFinder_v3(object, 
                                 PCs = 1:nPCs, 
                                 pN = 0.25, 
                                 pK = max_pK, 
                                 nExp = nExp_poi.adj, 
                                 reuse.pANN = FALSE, 
                                 sct = FALSE)
    }
    
    return(object)
  })
  
  # Extract doublet annotations
  extract_dfs = lapply(obj.list, function(object){
    columns = colnames(object@meta.data)[which(startsWith(colnames(object@meta.data), "DF."))]
    df = object@meta.data[columns]
    
    # Re-name
    if(!is.null(classify.by)){
      colnames(df) = c("DF.classifications", 
                       paste0(unlist(lapply(strsplit(columns[-1], "_"), function(x) x[[5]])), "_contribution"))
    } else {
      colnames(df) = "DF.classifications"
    }
    
    return(df)
  })
  names(extract_dfs) = NULL # to prevent rbind from appending list names to rownames
  extract_df = do.call(rbind, extract_dfs)
  
  # Save annotations for troubleshooting
  saveRDS(extract_df, "DF_annotations.rds")
  
  # rownames(extract_df) = unlist(lapply(strsplit(rownames(extract_df), "\\."), function(x) x[length(x)]))
  
  # Remove the prefix added from do.call(rbind)
  # prefix = paste0(unlist(lapply(strsplit(rownames(extract_df), "\\."), function(x) x[1])), ".")
  # prefix = paste0(unlist(lapply(names(obj.list), function(x) rep(x, nrow(obj.list[[x]])))), '.')
  # rownames(extract_df) = unlist(lapply(1:nrow(extract_df), function(i) sub(prefix[i], "", rownames(extract_df)[i])))
  
  # Add annotations to object
  if(!is.null(classify.by)){
    sorted_metadata = extract_df[match(colnames(SeuratObject), paste0(rownames(extract_df), suffix)),]
    SeuratObject = AddMetaData(SeuratObject, sorted_metadata)
  } else {
    SeuratObject$DF.classifications = extract_df[match(colnames(SeuratObject), paste0(rownames(extract_df), suffix)),]
  }
  
  return(SeuratObject)
}

#' Merge clusters fast!
#' 
#' Function that merges similar clusters on tree if they have fewer than X number of DEGs.
#' Note that the BuildClusterTree function should have been called for this object
#' 
#' 
quickMergeClusters = function(object, more.than.nDEGs = 6, log2FC.threshold = 1, p_val_adj.threshold = 0.001){
  tree <- object@tools$BuildClusterTree
  
  # Node heights
  library(ape)
  node.heights = as.integer(max(node.depth.edgelength(tree))-node.depth.edgelength(tree))
  names(node.heights) = 1:length(node.heights)
  order_nodes = sort(node.heights) #order(node.heights[-c(1:length(unique(YanAC$seurat_clusters)))])
  
  # Extract nodes below a given height (but not the tips), using 15% of max height as cutoff
  # nodes_to_check = names(order_nodes[order_nodes != 0 & order_nodes < height * max(node.heights)])
  nodes_to_check = head(names(order_nodes[order_nodes != 0]), 10)
  
  message("Checking nodes: ", paste0(nodes_to_check, collapse = ", "))
  
  node_nDEGs = lapply(nodes_to_check, function(node){
    node.markers <- wilcoxauc(object, ident.1 = 'clustertree', ident.2 = node)
    nDEGs(node.markers, log2FC.threshold = log2FC.threshold, p_val_adj.threshold = p_val_adj.threshold, presto = TRUE)
    # Same as doing: FindMarkers(YanAC, ident.1 = x, ident.2 = y)
  })
  
  message("number of DEGs: ", paste0(unlist(node_nDEGs), collapse = ", "))
  
  nodes_to_merge = nodes_to_check[which(node_nDEGs <= more.than.nDEGs)]
  
  message("Merging nodes: ", paste0(nodes_to_merge, collapse = ", "))
  
  for(node_to_merge in nodes_to_merge){
    merge_idents = sort(c(Seurat:::GetLeftDescendants(tree, node_to_merge) - 1, 
                          Seurat:::GetRightDescendants(tree, node_to_merge) - 1))
    
    # Set merged nodes to lower number of the two clusters
    object$seurat_clusters[WhichCells(object, idents = merge_idents)] = merge_idents[1] # paste0(merge_idents, collapse = "_")
  }
  
  # Make the clusters be in order
  old_clusters = sort(unique(object$seurat_clusters))
  new_clusters = 0:(length(old_clusters)-1)
  
  # Change each old cluster to its new cluster 
  for(i in 1:length(old_clusters)){
    object$seurat_clusters[object$seurat_clusters == old_clusters[i]] = new_clusters[i]
  }
  
  # Set the order and Idents
  object$seurat_clusters = factor(object$seurat_clusters, levels = new_clusters)
  Idents(object) = object$seurat_clusters
  
  return(object)
}

#' Merge clusters
#' 
#' Function that merges similar clusters on tree if they have fewer than X number of DEGs.
#' Note that the BuildClusterTree function should have been called for this object
#' 
#' 
mergeClusters = function(object, more.than.nDEGs = 6, log2FC.threshold = 1, p_val_adj.threshold = 0.001, num.nodes.to.check = 10){
  tree <- object@tools$BuildClusterTree
  
  # Node heights
  library(ape)
  node.heights = as.integer(max(node.depth.edgelength(tree))-node.depth.edgelength(tree))
  names(node.heights) = 1:length(node.heights)
  order_nodes = sort(node.heights) #order(node.heights[-c(1:length(unique(YanAC$seurat_clusters)))])
  
  # cat("Using ", height, " of ", max(node.heights), ": ", height * max(node.heights), "\n")
  # cat("Checking the ten most similar nodes...\n")
  
  # Extract nodes below a given height (but not the tips), using 15% of max height as cutoff
  # nodes_to_check = names(order_nodes[order_nodes != 0 & order_nodes < height * max(node.heights)])
  nodes_to_check = head(names(order_nodes[order_nodes != 0]), num.nodes.to.check)
  
  message("Checking nodes: ", paste0(nodes_to_check, collapse = ", "))
  
  node_nDEGs = lapply(nodes_to_check, function(node){
    node.markers <- FindMarkers(object, ident.1 = 'clustertree', ident.2 = node)
    nDEGs(node.markers, log2FC.threshold = log2FC.threshold, p_val_adj.threshold = p_val_adj.threshold)
    # Same as doing: FindMarkers(YanAC, ident.1 = x, ident.2 = y)
  })
  
  message("number of DEGs: ", paste0(unlist(node_nDEGs), collapse = ", "))
  
  nodes_to_merge = nodes_to_check[which(node_nDEGs <= more.than.nDEGs)]
  
  message("Merging nodes: ", paste0(nodes_to_merge, collapse = ", "))
  
  # order = levels(object$seurat_clusters)
  # object$seurat_clusters = as.character(object$seurat_clusters)
  
  for(node_to_merge in nodes_to_merge){
    merge_idents = sort(c(Seurat:::GetLeftDescendants(tree, node_to_merge) - 1, 
                          Seurat:::GetRightDescendants(tree, node_to_merge) - 1))
    
    # Set merged nodes to lower number of the two clusters
    object$seurat_clusters[WhichCells(object, idents = merge_idents)] = merge_idents[1] # paste0(merge_idents, collapse = "_")
  }
  
  # Make the clusters be in order
  old_clusters = sort(unique(object$seurat_clusters))
  new_clusters = 0:(length(old_clusters)-1)
  
  # Change each old cluster to its new cluster 
  for(i in 1:length(old_clusters)){
    object$seurat_clusters[object$seurat_clusters == old_clusters[i]] = new_clusters[i]
  }
  
  # for(element in 1:max){
  #   if(any(object$seurat_clusters == element)){
  #     next
  #   }
  #   object$seurat_clusters[object$seurat_clusters > element] = object$seurat_clusters[object$seurat_clusters > element] - 1
  #   cat(paste0("Re-numbering from ", element, "\n"))
  # }
  
  # Set the order and Idents
  object$seurat_clusters = factor(object$seurat_clusters, levels = new_clusters)
  Idents(object) = object$seurat_clusters
  
  return(object)
}

stackedBarGraph = function(SeuratObject, feature.1, feature.2, title = NULL, percent = FALSE){
  
  if(inherits(SeuratObject@meta.data[[feature.1]], "factor")) SeuratObject@meta.data[[feature.1]] = droplevels(SeuratObject@meta.data[[feature.1]])
  if(inherits(SeuratObject@meta.data[[feature.2]], "factor")) SeuratObject@meta.data[[feature.2]] = droplevels(SeuratObject@meta.data[[feature.2]])
  
  data <- table(SeuratObject@meta.data[[feature.1]], SeuratObject@meta.data[[feature.2]])
  melted = reshape2::melt(as.data.frame(data))
  # colnames(melted) = c("batch", "type", "composition")
  
  p = ggplot(melted, aes(fill=Var1, y=value, x=Var2)) + 
    geom_bar(position="fill", stat="identity", color = "black") +
    theme_minimal()+
    ggtitle(title)+
    xlab(feature.2)+
    labs(fill=feature.1)+
    theme_cowplot()+
    {if(percent) scale_x_continuous(labels = function(x) gsub('0\\.', '', x))}+
    # scale_fill_manual(feature.1)+
    ylab('composition')+
    theme(axis.text.x = element_text(hjust = 1, angle = 45), plot.title = element_text(hjust = 0.5))
  
  return(p)
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

nDEGs = function(table, up = TRUE, down = TRUE, 
                 log2FC.threshold = 0.25, 
                 p_val_adj.threshold = 0.05, 
                 presto = FALSE, value = FALSE){
  if(presto){
    if(up & !down){
      DEGs = (unique(subset(table, logFC > log2FC.threshold & padj < p_val_adj.threshold)$feature))
    } 
    else if(!up & down){
      DEGs = (unique(subset(table, logFC < -log2FC.threshold & padj < p_val_adj.threshold)$feature))
    } 
    else{
      DEGs = (unique(subset(table, abs(logFC) > log2FC.threshold & padj < p_val_adj.threshold)$feature))
    }
  }
  else {
    if(up & !down){
      DEGs = (unique(rownames(subset(table, avg_log2FC > log2FC.threshold & p_val_adj < p_val_adj.threshold))))
    } 
    else if(!up & down){
      DEGs = (unique(rownames(subset(table, avg_log2FC < -log2FC.threshold & p_val_adj < p_val_adj.threshold))))
    }
    else{
      DEGs = (unique(rownames(subset(table, abs(avg_log2FC) > log2FC.threshold & p_val_adj < p_val_adj.threshold))))
    } 
  }
  
  if(value) return(DEGs) else return(length(DEGs))
}

#' Label transfer
#' 
#' Performs the Seurat label transfer procedure
#' 
#' @param reference reference with labels
#' @param query query to label
#' @param dims dimensions to use from PCA
#' @param threshold prediction score threshold for a label to be transferred
labelTransfer = function(reference, query, dims = 1:30, threshold = 0.7, ...){
  anchors <- FindTransferAnchors(reference = reference, 
                                 query = query, 
                                 dims = dims, 
                                 reference.reduction = "pca", 
                                 ...)
  predictions <- TransferData(anchorset = anchors, refdata = reference$transfer, dims = dims)
  query <- AddMetaData(query, metadata = predictions)
  
  # UMAP projection
  reference <- RunUMAP(reference, dims = dims, reduction = "pca", return.model = TRUE)
  query <- MapQuery(anchorset = anchors, reference = reference, 
                    query = query, refdata = list(transfer = "transfer"), 
                    reference.reduction = "pca", reduction.model = "umap")
  
  query$predicted.transfer[query$prediction.score.max < threshold] = "Unknown" # threshold based on figure 3D of Stuart et al. 2019
  return(query)
}

#' Clustered dotplot
#' 
#' Hierarchically clusters the columns of a DotPlot
#' 
#' @param dotplot a DotPlot that needs to be clustered
ClusteredDotplot2 = function(object, 
                             features,
                             group.by = 'annotated',
                             cluster_rows = FALSE, 
                             cluster_columns = FALSE, 
                             order = NULL, 
                             dot.scale.factor = 0.2, 
                             show.legend = TRUE,
                             font.family = 'ArialMT',
                             max.pct = 100,
                             title = NULL, 
                             pseudocount = 10, # for minimum dot size 0%
                             col.low = 'grey', 
                             col.high = "#584B9FFF",
                             min.z.score = -1, 
                             max.z.score = 2,
                             rotate = FALSE,
                             species.size = NULL,
                             ...){
  
  if(!is.null(order)) {
    Idents(object) = order
    cluster_columns = FALSE
  }
  
  # This prevents duplicate gene symbols
  
  norm.expr = AverageExpression(object, features = features, group.by = group.by, slot = "data", assay = 'RNA')$RNA
  norm.expr = norm.expr[match(features, rownames(norm.expr)),] # allow for duplicate genes
  scaled.expr = t(scale(t(norm.expr)))

  # Truncate values to max.z.score
  scaled.expr[scaled.expr < min.z.score] = min.z.score
  scaled.expr[scaled.expr > max.z.score] = max.z.score
  exp_mat = scaled.expr
  
  # Percentage expression
  if(!is.null(species.size)){
    pct.list = readRDS("../../Ortho_Objects/AC_OT_percent_expressed.rds") #readRDS("../../Ortho_Objects/AC_OT_percent_expressed_list.rds")
    pct.list.features = pct.list[features]
    cutoff.list = lapply(pct.list.features, function(mat) apply(mat, 2, function(col) length(which(col > species.size))))
    percent_mat = do.call(rbind, cutoff.list) %>% as.data.frame() %>% setNames(gsub('\\*', '', levels(ac.ortho$type_unordered)))
    percent_mat = percent_mat[,levels(ac.ortho@meta.data[[group.by]])]
    
  } else {
    percent_mat = PercentageExpressed2(object, features = unique(features), group.by = group.by)
    percent_mat = percent_mat[match(features, rownames(percent_mat)),] # allow for duplicate genes
    percent_mat[percent_mat > max.pct] = max.pct
  }
  
  if(rotate){
    exp_mat = t(exp_mat)
    percent_mat = t(percent_mat)
  }
  
  col_fun = circlize::colorRamp2(c(min.z.score, max.z.score), c(col.low, col.high))
  cell_fun = function(j, i, x, y, w, h, fill){
    grid.circle(x=x, 
                y=y, 
                r=unit(((percent_mat[i, j]+pseudocount)/100) * dot.scale.factor, "cm"), 
                gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
    }
  
  # Create a dot legend
  lgd.values = c(0+pseudocount, 25+pseudocount, 50+pseudocount, 75+pseudocount, 100+pseudocount)
  lgd.values[lgd.values > (max.pct + pseudocount)] = (max.pct + pseudocount)
  lgd = Legend(labels = seq(0,100, by = 25), title = "Percent\nexpressed",
               graphics = list(
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[1]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[2]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[3]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[4]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[5]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA))
               ), title_gp = gpar(fontfamily = font.family, fontface = "plain"))
  
  
  # if(rotate){
  #   draw(
  #     Heatmap(exp_mat,
  #             heatmap_legend_param = list(title = "Scaled\nexpression", 
  #                                         title_gp = gpar(fontfamily = font.family, fontface = "plain"), 
  #                                         border = "black",      # Black border around legend
  #                                         # at = seq(-2, 2, by = 1), # Define tick positions
  #                                         labels_gp = gpar(col = "black"), # Black tick labels
  #                                         ticks_gp = gpar(col = "black")),
  #             column_title = title,
  #             col=col_fun,
  #             rect_gp = gpar(type = "none"),
  #             cell_fun = cell_fun,
  #             row_names_gp = gpar(fontsize = 10, fontface = "italic", family = font.family),
  #             column_names_gp = gpar(fontsize = 10, family = font.family),
  #             cluster_rows = cluster_columns,
  #             cluster_columns = cluster_rows,
  #             cluster_row_slices = FALSE,
  #             show_heatmap_legend = show.legend,
  #             border = "black", 
  #             use_raster = FALSE, 
  #             ...),
  #     annotation_legend_list = if(show.legend) lgd else NULL)
  # } else {
  
  if(rotate){
    rownames = gpar(fontsize = 10, family = font.family)
    colnames = gpar(fontsize = 10, fontface = "italic", family = font.family)
  } else {
    rownames = gpar(fontsize = 10, fontface = "italic", family = font.family)
    colnames = gpar(fontsize = 10, family = font.family)
  }
    draw(
      Heatmap(exp_mat,
              heatmap_legend_param = list(title = "Scaled\nexpression", 
                                          title_gp = gpar(fontfamily = font.family, fontface = "plain"), 
                                          border = "black",      # Black border around legend
                                          # at = seq(-2, 2, by = 1), # Define tick positions
                                          labels_gp = gpar(col = "black"), # Black tick labels
                                          ticks_gp = gpar(col = "black")),
              column_title = title,
              col=col_fun,
              rect_gp = gpar(type = "none"),
              cell_fun = cell_fun,
              row_names_gp = rownames,
              column_names_gp = colnames,
              cluster_rows = cluster_rows,
              cluster_columns = cluster_columns,
              cluster_row_slices = FALSE,
              show_heatmap_legend = show.legend,
              border = "black", 
              use_raster = FALSE, 
              ...),
      annotation_legend_list = if(show.legend) lgd else NULL, )
  # }
  
}

#' Clustered dotplot
#' 
#' Hierarchically clusters the columns of a DotPlot
#' 
#' @param dotplot a DotPlot that needs to be clustered
ClusteredDotplot = function(object, 
                            row.dendrogram = NULL, 
                            column.dendrogram = NULL, 
                            order = NULL, 
                            cluster.columns = TRUE, 
                            # scale.factor = 0.04, 
                            dot.scale.factor = 0.2,
                            max.pct = 100,
                            title = NULL, 
                            pseudocount = 0,
                            col.low = 'grey', 
                            col.high = "#584B9FFF",
                            ...){
  
  if(!is.null(order)) {
    Idents(object) = order
    cluster.columns = FALSE
  }
  dotplot = DotPlot(object, ...) + coord_flip()
  
  df <- dotplot$data
  exp_mat<-df %>%
    dplyr::select(-pct.exp, -avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()
  row.names(exp_mat) <- exp_mat$features.plot  
  exp_mat <- exp_mat[,-1] %>% as.matrix()
  
  percent_mat<-df %>% 
    dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
    pivot_wider(names_from = id, values_from = pct.exp) %>% 
    as.data.frame() 
  row.names(percent_mat) <- percent_mat$features.plot  
  percent_mat <- percent_mat[,-1] %>% as.matrix()
  percent_mat[percent_mat > max.pct] = max.pct
  
  col_fun = circlize::colorRamp2(c(-2, 2), c(col.low, col.high))
  cell_fun = function(j, i, x, y, w, h, fill){
    
    grid.circle(x=x,y=y,r= unit(((percent_mat[i, j]+pseudocount)/100) * dot.scale.factor, "cm"), 
                gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
  
  # Create a dot legend
  lgd.values = c(pseudocount, 25+pseudocount, 50+pseudocount, 75+pseudocount, 100+pseudocount)
  lgd.values[lgd.values > max.pct] = max.pct
  lgd = Legend(labels = seq(0,100, by = 25), title = "Percentage\nexpressed",
               graphics = list(
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[1]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[2]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[3]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[4]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA)),
                 function(x, y, w, h) grid.circle(x, y, r=unit(lgd.values[5]/100 * dot.scale.factor, "cm"), gp=gpar(fill = 'black', col = NA))
               ))
  
  if(is.null(row.dendrogram)){
    draw(
      Heatmap(exp_mat,
            heatmap_legend_param = list(title = "Scaled\nexpression"),
            column_title = title,
            col=col_fun,
            rect_gp = gpar(type = "none"),
            cell_fun = cell_fun,
            row_names_gp = gpar(fontsize = 10),
            # row_km = 4,
            cluster_rows = FALSE,
            cluster_columns = cluster.columns,
            cluster_row_slices = FALSE,
            border = "black", 
            use_raster = FALSE),
      annotation_legend_list = lgd)
  } else if(!is.null(column.dendrogram)){
    draw(
      Heatmap(exp_mat,
            heatmap_legend_param = list(title = "Scaled\nexpression"),
            column_title = title,
            col=col_fun,
            rect_gp = gpar(type = "none"),
            cell_fun = cell_fun,
            row_names_gp = gpar(fontsize = 10),
            # row_km = 4,
            cluster_rows = row.dendrogram,
            cluster_columns = column.dendrogram,
            cluster_row_slices = FALSE,
            border = "black", 
            use_raster = FALSE),
      annotation_legend_list = lgd)
  } else {
    draw(
      Heatmap(exp_mat,
            heatmap_legend_param = list(title = "Scaled\nexpression"),
            column_title = title,
            col=col_fun,
            rect_gp = gpar(type = "none"),
            cell_fun = cell_fun,
            row_names_gp = gpar(fontsize = 10),
            # row_km = 4,
            cluster_columns = cluster.columns,
            cluster_rows = FALSE,
            cluster_row_slices = FALSE,
            border = "black", 
            use_raster = FALSE),
      annotation_legend_list = lgd)
  }
  
}

PrettyHistogram2 = function(vector, vector2 = NULL, vline = NULL, logX = FALSE, logY = FALSE, 
                            title = NULL, xlab = NULL, name1 = 'group1', name2 = 'group2', 
                            line.col = 'red', fill = 'white',
                            ...){
  
  if(!is.null(vector2)){
    # ggplot(vegLengths, aes(length, fill = veg)) + 
    #   geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
    df = data.frame(variable = c(vector, vector2), group = c(rep(name1, length(vector)), rep(name2, length(vector2))))
    p = ggplot(df, aes(x=variable, fill = group)) + 
      # geom_histogram(aes(y = ..density..), ...)+
      geom_histogram(aes(y = after_stat(count / sum(count))), ...) +
      xlab(xlab)+
      ylab("Proportion")+
      scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      {if(logX) scale_x_continuous(trans='log2')} +
      {if(logY) scale_y_continuous(trans='log2', expand = expansion(mult = c(0,0.1)))} +
      {if(!is.null(vline)) geom_vline(xintercept=vline, color="grey", linetype="dashed")}+
      theme_classic()+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    df = data.frame(variable = vector)
    p = ggplot(df, aes(x=variable)) + 
      geom_histogram(color="black", fill=fill, ...)+ #aes(y=..density..))+
      xlab(xlab)+
      ylab("Frequency")+
      scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      {if(logX) scale_x_continuous(trans='log2')} +
      {if(logY) scale_y_continuous(trans='log2', expand = expansion(mult = c(0,0.1)))} +
      {if(!is.null(vline)) geom_vline(xintercept=vline, color=line.col, linetype="solid")}+
      theme_classic()+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(p)
}

PrettyHistogram = function(vector, vector2 = NULL, vline = NULL, logX = FALSE, logY = FALSE, title = NULL, xlab = NULL, fill = 'white', ...){
  # name = deparse(substitute(vector))
  df = data.frame(variable = vector)
  p = ggplot(df, aes(x=variable)) + 
    geom_histogram(color="black", fill=fill, ...)+ #bins = bins, binwidth = binwidth)+ #aes(y=..density..))+
    xlab(xlab)+
    ylab("Frequency")+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    {if(logX) scale_x_continuous(trans='log2')} +
    {if(logY) scale_y_continuous(trans='log2', expand = expansion(mult = c(0,0.1)))} +
    {if(!is.null(vline)) geom_vline(xintercept=vline, color="red3", linetype="dashed")}+
    theme_classic()+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

BINplot <- function(SeuratObject, x.limits=NULL, y.limits=NULL, title=NULL, 
                    feature.low = 200, feature.high = 7000, percent.mt.low = 0, percent.mt.high = 5){
  dat=data.frame(percent.mt=SeuratObject$percent.mt, 
                 count=SeuratObject$nCount_RNA, 
                 feature=SeuratObject$nFeature_RNA)
  plot=ggplot(dat, aes(percent.mt, feature))+
    geom_bin_2d(binwidth = c(0.1, 100))+
    scale_fill_gradient(low="grey", high="black", trans = "log")+
    geom_hline(yintercept=feature.low, color="red3", linetype="dashed")+
    geom_hline(yintercept=feature.high, color="red3", linetype="dashed")+
    geom_vline(xintercept=percent.mt.high, color="red3", linetype="dashed")+
    geom_vline(xintercept=percent.mt.low, color="red3", linetype="dashed")+
    {if(!is.null(x.limits)) scale_x_continuous(limits = x.limits)}+
    {if(!is.null(y.limits)) scale_y_continuous(limits = y.limits)}+
    {if(!is.null(title)) ggtitle(title)}+
    # stat_cor(method = "pearson", label.x = 20)+
    theme_classic()
  
  return(plot)
}



#' Convert genes 
#'
#' @param genes vector of genes
#' @param from the species from which the symbols originate
#' @param to the species to convert gene symbols to
#' 
#' @return Returns a vector of converted gene symbols
convertGenes <- function(genes, from, to){
  from_column=which(toupper(colnames(orthology_all_others))==paste0(toupper(from), ".GENE.NAME"))
  to_column=which(toupper(colnames(orthology_all_others))==paste0(toupper(to), ".GENE.NAME"))
  
  # If gene is not in orthology table, revert to original
  converted_genes = orthology_all_others[match(toupper(genes), toupper(orthology_all_others[,from_column])), to_column]
  converted_genes[is.na(converted_genes)] = genes[is.na(converted_genes)]
  
  return(converted_genes)
}

#' Annotated DotPlot
#'
#' @param object A Seurat object.
#' @param ... Other parameters to DotPlot 
#' 
#' @import ggplot2
#' @import Seurat
#'
#' @return Returns a DotPlot ggplot object
AnnotatePlot = function(plt, annotation, 
                        object_x=NULL, x_annotation = NULL, color_genes_x = FALSE,
                        object_y=NULL, y_annotation = NULL, color_genes_y = FALSE, 
                        COLOR_FUN = ColorCode){
  
  y_axis_order = layer_scales(plt)$y$range$range
  x_axis_order = layer_scales(plt)$x$range$range
  
  annotate_y = (!is.null(y_annotation) | color_genes_y)
  annotate_x = (!is.null(x_annotation) | color_genes_x)
  
  if(annotate_y){
    if(!is.null(object_y)){
      y_idents = Idents(object_y) #@meta.data[,y_annotation]
      y_axis_annotation = object_y@meta.data[match(y_axis_order, y_idents),][,y_annotation]
      colors = COLOR_FUN(annotation, y_axis_annotation)
      # print(y_axis_annotation)
    } else if(color_genes_y){
      colors = Colors(annotation)[match(y_axis_order, Genes(annotation))]
    }
  } else if(annotate_x){
    if(!is.null(object_x)){
      x_idents = Idents(object_x) #@meta.data[,y_annotation]
      x_axis_annotation = object_x@meta.data[match(x_axis_order, x_idents),][,x_annotation]
      # print(x_axis_annotation)
      colors = COLOR_FUN(annotation, x_axis_annotation)
    }
    else if(color_genes_x){
      # x_axis_annotation = Annotation(annotation_object)[match(x_axis_order, Genes(annotation_object))]
      colors = Colors(annotation)[match(x_axis_order, Genes(annotation))]
    }
  } else {
    stop('Please add only an x or y axis annotation')
  }
  
  message("\nColors: \n", paste0(colors, collapse = ", "))
  
  plt2 = plt + 
    {if(annotate_y) annotate("rect",
                             ymin = seq(0.5, length(y_axis_order)-0.5, by = 1), 
                             ymax = seq(1.5, length(y_axis_order)+0.5, by = 1),
                             xmin = 0.5, xmax = length(x_axis_order)+0.5,
                             alpha = .4, fill = colors)}+
    {if(annotate_x) annotate("rect",
                             xmin = seq(0.5, length(x_axis_order)-0.5, by = 1), 
                             xmax = seq(1.5, length(x_axis_order)+0.5, by = 1),
                             ymin = 0.5, ymax = length(y_axis_order)+0.5,
                             alpha = .4, fill = colors)} +
    theme(panel.grid = element_line(color = rgb(235, 235, 235, 100, 
                                                maxColorValue = 255), 
                                    linewidth = 1, 
                                    linetype = 1))
  return(plt2)
}

MetadataNamedVector = function(object, feature.1, feature.2){
  Metadata(object, feature.1, feature.2)[[feature.2]] %>% setNames(Metadata(object, feature.1, feature.2)[[feature.1]])
}

#' Metadata from SeuratObject
#' 
Metadata2 = function(object, feature.list, sort = TRUE){
  if(sort) {
    df = unique(object@meta.data[,feature.list]) %>% arrange(!!sym(feature.list[[1]]))
  } else {
    df = unique(object@meta.data[,feature.list])
  }
  rownames(df) = NULL
  df
}

#' Metadata from SeuratObject
#' 
Metadata = function(object, feature.1, feature.2, feature.3 = NULL, feature.4 = NULL){
  df = unique(object@meta.data[,c(feature.1, feature.2, feature.3, feature.4)]) %>% arrange(!!sym(feature.1)) #arrange(eval(parse(text=feature.1)))
  rownames(df) = NULL
  return(df)
}

DEGsPerCluster = function(object, DE_table, group.by = "seurat_clusters"){
  DEGs = lapply(Clusters(object, group.by = group.by), function(cluster) {
    subset(DE_table, group == cluster & abs(logFC) > 1 & padj < 0.001)$feature
  })
  
  nDEG_data = data.frame(cluster = Clusters(object), 
                         nDEGs = unlist(lapply(DEGs, length)), 
                         size = as.numeric(table(object@meta.data[,group.by])))
  
  return(nDEG_data)
}

Clusters = function(object, group.by = "seurat_clusters"){
  clusters = sort(unique(object@meta.data[[group.by]]))
  return((clusters))
}

#' Metadata from SeuratObject
#' 
MeanMetadata = function(object, feature.1, feature.2, as.df = FALSE, FUN = mean){
  df = object@meta.data[,c(feature.1, feature.2)] %>% 
    arrange(eval(parse(text=feature.1)))
  rownames(df) = NULL
  
  clusters = Clusters(object, group.by = feature.1)
  means = lapply(clusters, function(cluster){
    FUN(df[df[,feature.1] == cluster,feature.2])
  }) %>% unlist
  
  names(means) = Clusters(object, group.by = feature.1)
  
  if(as.df) return(data.frame(names(means), means) %>% setNames(c(feature.1, feature.2)))
  return(means)
}

#' Ordered DotPlot
#'
#' @param object A Seurat object.
#' @param ... Other parameters to DotPlot 
#' 
#' @import ggplot2
#' @import Seurat
#'
#' @return Returns a DotPlot ggplot object
OrderedDotPlot = function(object, 
                          annotation,
                          group.by = NULL, 
                          color.clusters.by = NULL, 
                          order = NULL,
                          ...){
  
  args = list(...)
  
  if(is.null(group.by)){
    labels = Idents(object)
  } else {
    labels = object[[group.by]]  
  }
  
  if(is.null(order)){
    order = CelltypeOrder(annotation)
  }
  
  # Arrange in proper order
  # object@meta.data[,color.clusters.by] = as.character(object@meta.data[,color.clusters.by])
  cell_class_df = Metadata(object, group.by, color.clusters.by) %>% 
    arrange(factor(eval(parse(text = color.clusters.by)), levels = order), # by cell class
            factor(!!sym(group.by), levels = smart_sort(unique(object@meta.data[[group.by]]))))
  
  # Set idents to order of color labels so that plot is in correct order
  if('coord.flip' %in% args){
    if(!args$coord.flip){
      Idents(object) = factor(object@meta.data[,group.by], levels = rev(cell_class_df[,group.by]))
    } 
  } else {
    Idents(object) = factor(object@meta.data[,group.by], levels = (cell_class_df[,group.by]))
  }
  
  plot <- DotPlot3(object, col.min = -1, col.max = 2, ...) + 
    # scale_color_gradient(low = "lightgrey", high = "#584B9FFF", limits = c(-1,2))+
    # scale_radius(limits = c(0,100), range = c(0, 6))+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    # guides(color = guide_colorbar(title = "Average\nexpression"), size = guide_legend(title = "Percent\nexpressed"))+
    theme(panel.grid = element_line(color = rgb(235, 235, 235, 100, 
                                                maxColorValue = 255), 
                                    linewidth = 1, 
                                    linetype = 1)) 
          # Italicize gene symbols, edit: doing this in dotplot3 now
          # axis.text.x = element_text(face = "italic"))
  
  return(plot)
}

#' Annotated DotPlot
#'
#' @param object A Seurat object.
#' @param ... Other parameters to DotPlot 
#' 
#' @import ggplot2
#' @import Seurat
#'
#' @return Returns a DotPlot ggplot object
AnnotatedDotPlot = function(object, 
                            annotation,
                            color.genes = FALSE,
                            ...){
  
  args = list(...)
  
  Idents(object) = args$group.by
  
  # Convert color.clusters.by to character since it causes issues when it's a factor
  object@meta.data[,args$color.clusters.by] = as.character(object@meta.data[,args$color.clusters.by])
  # if(is.factor(object[[args$group.by]])) object[[args$group.by]] = as.character(object[[args$group.by]]); message("converting factor to character\n")
  
  if(color.genes){
    ordered.plot = DotPlot3(object, features = Genes(annotation), ...)
    annotated.plot = AnnotatePlot(ordered.plot, annotation, color_genes_x = TRUE) # y_annotation = args$color.clusters.by)
  } else{
    ordered.plot = OrderedDotPlot(object, annotation, ...)
    annotated.plot = AnnotatePlot(ordered.plot, annotation, object_y=object, y_annotation = args$color.clusters.by)
  }
  
  return(annotated.plot)
}



# Moratorium
major_class_colors = function(annotations){
  for(i in 1:nrow(major_color_code)){
    annotations[annotations == major_color_code$annotation[i]] = major_color_code$color[i]
  }
  return(annotations)       
}

peptide_colors = function(annotations){
  for(i in 1:nrow(peptide_color_code)){
    annotations[annotations == peptide_color_code$annotation[i]] = peptide_color_code$color[i]
  }
  return(annotations)       
}

#' Color code amacrines
#' 
#' Color codes a vector of amacrine cell annotations with a consistent color code
#'
#' @param annotations a vector of amacrine cell annotations
#'
#' @return Returns a color vector
amacrine_colors = function(annotations){
  for(i in 1:nrow(amacrine_color_code)){
    annotations[annotations == amacrine_color_code$annotation[i]] = amacrine_color_code$color[i]
  }
  return(annotations)       
}

#' Amacrine order
#' 
#' A function that returns the amacrine cell default order
#' 
amacrine_order = function(){
  return(amacrine_color_code$annotation)
}

peptide_order = function(){
  return(peptide_color_code$annotation)
}

major_class_order = function(){
  return(c("Rod", "Cone", "HC", "BC", "BP", "GabaAC", "GlyAC", "RGC", "MG", "Other"))
}

paramSweep_DT <- function(seu, PCs=1:10, sct = FALSE, num.cores=6) {
  require(Seurat); require(fields); require(parallel);
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands
  
  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@assays$RNA@counts[ , real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    # cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep_v3,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,
                        mc.cores=num.cores)
    # stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep_v3,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

