# Packages
library(xgboost)
#library(philentropy)

#' A function to compute the relative tightness of clusters by comparing within cluster diameter to cross-cluster distances
#' 
#' @param coord_ids character vector corresponding to the factors of the reduced dimensional space.
#' @param df data frame with columns corresponding coord_ids and also corresponding to cluster
#' @return a list with elements cross_clust - corresponding to mean cross cluster distances, and within_clust - corresponding to mean within cluster distances
WithinCrossClusterDist = function(df, coord_ids = c("UMAP1","UMAP2")){
  
  if (!("cluster" %in% colnames(df))) stop("df must have a column corresponding to cluster ids")
  
  if (sum(coord_ids %in% colnames(df)) != length(coord_ids) ) stop("df must have columns corresponding to all entries of coord_ids")
  
  if (!is.factor(df$cluster)) df$cluster = factor(df$cluster)
  
  # Compute centroids
  df_cent = df %>% dplyr::group_by(cluster) %>% dplyr::summarise_all(list(mean))
  
  # For each cluster, compute the median within cluster distance and 
  within_clust = c()
  cross_clust = c()
  for (clust in levels(df$cluster)){
    print(clust)
    within_clust_tmp = median(sapply(rownames(subset(df, cluster==clust)),
                                     function(x){
                                       cents = as.numeric(df_cent[which(rownames(df_cent) == clust),coord_ids]);
                                       vec = as.numeric(df[x, coord_ids])
                                       distx = sqrt(sum((cents-vec)^2));
                                       return(distx)
                                     } )
    )
    within_clust = c(within_clust, within_clust_tmp)
    cross_clust_tmp = median(sapply(rownames(subset(df, cluster==clust)),
                                    function(x){
                                      cents = as.matrix(df_cent[-which(rownames(df_cent) == clust),coord_ids]);
                                      vec = as.numeric(df[x, coord_ids])
                                      distx = min(sqrt(rowSums(sweep(cents,2,vec)^2)));
                                      return(distx)
                                    } )
    )
    cross_clust = c(cross_clust, cross_clust_tmp)
  }
  
  return(list(cross_clust_dist = cross_clust, within_clust_dist = within_clust))
  
}



#' Training an xgboost model for determining classification errors
#'
#' @param train_object0 the scR object used for training
#' @param var.genes character vec, set of features used for training
#' @param do.scale bool, scales the data
#' @param scale.mean numeric vec, user supplied mean vector to use for centering
#' @param scale.var numeric vec, user supplied standard deviation vector for scaling
#'
#' @return list, with elements corresponding to
#' bst_model - the trained xgboost model
#' class_errors - the class specific errors computed from the validation set
#' @export
#'
#' @examples
XGBoost_train_KS = function(train_object0, var.genes, do.scale=FALSE,scale.mean = NULL, scale.var = NULL, min.val = NULL, max.val=NULL){
  
  predictor_Data = as.matrix(train_object0@assays$RNA@data[var.genes,])
  if (do.scale){
    if (is.null(scale.mean) & is.null(scale.var)){
      predictor_Data = t(scale(t(predictor_Data)))
    } else {
      predictor_Data = t(scale(t(predictor_Data), center=scale.mean, scale=scale.var))
    }
    
  }
  
  if (!is.null(min.val) & !is.null(max.val)){
    predictor_Data = minmax(predictor_Data, min = min.val, max = max.val)
  }
  
  max.cells.per.ident = 300; train.frac = 0.6
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using mininum of ", 0.6*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training"))
  for (i in as.character(levels(Idents(train_object0)))){
    cells.in.clust = WhichCells(train_object0,idents = i);
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(as.numeric(i)-1,length(train.temp))); validation.label = c(validation.label, rep(as.numeric(i)-1, length(validation.temp)));
  }
  train_matrix <- xgb.DMatrix(data = t(predictor_Data[,training.set]), label=training.label)
  validation_matrix <- xgb.DMatrix(data = t(predictor_Data[,validation.set]), label=validation.label)
  
  numberOfClasses <- length(unique(training.label))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "eta" = 0.2,"max_depth"=6, subsample = 0.6)
  nround    <- 200 # number of XGBoost rounds
  print(1)
  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = nround)
  
  # Predict hold-out validation set
  validation_pred <- predict(bst_model, newdata = validation_matrix)
  validation_prediction <- matrix(validation_pred, nrow = numberOfClasses,
                                  ncol=length(validation_pred)/numberOfClasses)
  
  valid_predlabels=apply(validation_prediction,2,which.max)-1
  plotConfusionMatrix(table(validation.label, valid_predlabels))
  
  A = table(validation.label, valid_predlabels)
  
  class_errors = sapply(c(1:nrow(A)), function(x){vec = A[x,]; return(sum(vec[-x])/sum(vec) )})
  
  to.return = list()
  to.return$bst_model = bst_model
  to.return$class_errors = class_errors
  to.return$Conf = A
  to.return$true_label = validation.label
  to.return$assigned_label = valid_predlabels
  
  return(to.return)
}




#' Title Computes the rao diversity index
#'
#' @param X matrix with rows as categories and columns as features
#' @param distfun character, distance metric (Default: : `"manhattan"``)
#' @param p_i 
#'
#' @return the value of the Rao diversity index
#' @export
#'
#' @examples
raoDiversity = function(X, distfun = "manhattan", p_i = NULL ){
  
  if (is.null(p_i)){
    p_i = rep(1/nrow(X), nrow(X))
  }
  
  if (length(p_i) != nrow(X)) stop("p_i needs to have the same number of elements as the rows of X")
  if (sum(p_i) != 1) stop("p_i needs to sum to 1")
  
  D = as.matrix(dist(X, method = distfun))
  rao_div = sum (D * (p_i %o% p_i))
  return(rao_div)
}



#' plotConfusionMatrix
#'
#' @param X matrix, contingency table
#' @param row.scale boolean, if row.scaling should be done. Values are TRUE, FALSE, NA
#' @param col.low 
#' @param col.high 
#' @param max.size 
#' @param ylab.use 
#' @param xlab.use 
#' @param order character. Values are "Col" or "Row". Determines of column or row ordering should be performed
#' @param x.lab.rot 
#' @param plot.return 
#' @param max.perc 
#' @param alpha.use transparency
#'
#' @return
#' @export
#'
#' @examples
plotConfusionMatrix <- function(X,row.scale=TRUE, col.low="blue", col.high="red", max.size=5, ylab.use="Known", xlab.use="Predicted", order=NULL, x.lab.rot=FALSE, plot.return=TRUE, max.perc=100, alpha.use = 0.9){
  
  if (!is.na(row.scale)){
    if (row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
    if (!row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  }
  
  X[is.na(X)] = 0
  if (max(X) > 100){
    X=X/100
  }
  
  X[X > max.perc] = max.perc
  
  orig.rownames = rownames(X)
  orig.colnames = colnames(X)
  if (!is.null(order)){
    if (order == "Row"){  
      factor.levels = c()
      for (i1 in colnames(X)){
        if (max(X[,i1]) < 50) next
        ind.sort = rownames(X)[order(X[,i1], decreasing=TRUE)]
        ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
        factor.levels = c(factor.levels, ind.sort[1])
      }
      factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
      factor.levels = factor.levels[!is.na(factor.levels)]
    } 
    
    if (order == "Col") {
      factor.levels = c()
      for (i1 in rownames(X)){
        if (max(X[i1,]) < 50) next
        ind.sort = rownames(X)[order(X[i1,], decreasing=TRUE)]
        ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
        factor.levels = c(factor.levels, ind.sort[1])
      }
      factor.levels = c(factor.levels, setdiff(rownames(t(X)), factor.levels))
      factor.levels = factor.levels[!is.na(factor.levels)]
    } 
  } else {
    factor.levels = rownames(t(X))
  }
  
  factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
  X = melt(X)
  colnames(X) = c("Known", "Predicted", "Percentage")
  #X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  #X$Predicted = factor(X$Predicted, levels = rev(factor.levels))
  
  if (!is.null(order)){
    if (order == "Row"){ 
      X$Known = factor(X$Known, levels=rev(factor.levels));
      X$Predicted = factor(X$Predicted, levels = orig.colnames)
      
    }
    if (order == "Col"){
      X$Predicted = factor(X$Predicted, levels = factor.levels);
      X$Known = factor(X$Known, levels=rev(orig.rownames));
    }
  } else {
    X$Known = factor(X$Known, levels=rev(unique(X$Known)));
    X$Predicted = factor(X$Predicted, levels=unique(X$Predicted));
  }
  
  #print(sum(is.na(X$Known)))
  
  
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low =col.low,   high = col.high, limits=c(0, max.perc ))+scale_size(range = c(1, max.size), limits = c(0,max.perc))+   theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic"))  
  
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  print(p)
  
  if (plot.return) return(p)
}

# Occupancy number
OccNo = function(p){
  return(1/sum(p^2))
}

getGroupOrder = function(groups){
  vec_order = rep(1,length(unlist(groups)))
  
  for (ind in 1:length(groups)){
    vec_order[groups[[ind]]] = ind
  }
  return(vec_order)
}

getGroupNames = function(groups){
  vec_names = c()
  
  for (ind in 1:length(groups)){
    vec_names = c(vec_names, paste0(groups[ind][[1]], collapse=","))
  }
  return(vec_names)
}


#' Title
#'
#' @param X Confusion Matrix (late clusters x late clusters), assumed to be column normalized 
#' @param corr_thresh threshold for pearson correlation to collapse clusters (default:0.6)
#'
#' @return
#' @export
#'
#' @examples
FindClusterGroups = function(X=NULL, corr_thresh = 0.6, useJSD = FALSE, JSD_thresh = 0.25){
  
  groups = list()
  is_grouped = rep(FALSE, ncol(X))
  names(is_grouped) = c(1:ncol(X))
  
  for (col_id1 in c(1:ncol(X))){
    clust1 = colnames(X)[col_id1]
    if (is_grouped[clust1]) next
    # Use Pearson correlation
    if (!useJSD){
      cor_vals = as.numeric(cor(X[,clust1],X)); names(cor_vals) = colnames(X)
      clusts_to_add = as.numeric(names(cor_vals)[cor_vals >= corr_thresh])
    } else {
      X1 = as.data.frame.matrix(X)
      jsd_vals = sapply(1:ncol(X1), function(x) JSD(t(X1[,c(col_id1,x)])))
      names(jsd_vals) = colnames(X)
      clusts_to_add = as.numeric(names(jsd_vals)[jsd_vals < JSD_thresh])
    }
    
    
    # Check if part of an existing group
    if (sum(is_grouped[clusts_to_add]) > 0){
      flag=0
      # Only add to existing group if distribution overla is significant with at least 40% of existing clusters
      for (g in 1:length(groups)){
        if (length(intersect(groups[[g]],clusts_to_add )) / length(groups[[g]]) > 0.4 ){
          groups[[g]] = union(groups[[g]], clusts_to_add)
          flag=1
        } else {
          groups[[g]] = setdiff(groups[[g]], clusts_to_add)
        }
      }
      
      if (flag==0){
        # Create new group
        groups[[length(groups) + 1]] = clusts_to_add
      }
      
      is_grouped[clusts_to_add] = TRUE
    } else {
      is_grouped[clusts_to_add] = TRUE
      groups[[length(groups) + 1]] = clusts_to_add
    }
    
  }
  
  for (l in 1:length(groups)){
    groups[[l]] = sort(groups[[l]])
  }
  
  return(groups)
  
}





#' Title
#'
#' @param W0_list a list of OT weight matrices, rows correspond to earlier stage and columns correspond to later stage. Each column sums to 1
#' @param types_in_subclass a list of types to look backwards, any number of combinations 1-45
#' @param col_str color map
#' @param thresh.use minimum weight to be rendered, default 0.1
#'
#' @return Sankey plot
#' @export 
#'
#' @examples
SubclassSankey = function(W0_list = NULL, types_in_subclass = NULL, col_str = 'rgba(0,100,0,trans)',
                          thresh.use = 0.1, title_str = "Subclass"){
  W0_list_new = list()
  types_now = types_in_subclass
  for (k in c(5:1)){
    
    W0 = W0_list[[k]]
    W0 = W0[, types_now, drop=FALSE]
    types_now = rownames(W0)[apply(W0,1,max)>thresh.use]
    W0_mat = W0[apply(W0,1,max)>thresh.use, , drop=FALSE]
    W0_list_new[[k]] = W0_mat
    
  }
  
  labels.use = c(); source.use = c(); target.use = c(); value.use = c(); color.use0 = c()
  num_nodes = 0
  for (k in c(1:5)){
    
    W0 = W0_list_new[[k]]
    W0_color = W0
    W0_color[W0 > 0] = 1
    W0_color[W0 > 0.2] = 2
    W0_color[W0 > 0.5] = 3
    W0_color[W0 > 0.8] = 4
    
    
    label1 = c(rownames(W0), colnames(W0))
    labels.use = union(labels.use, label1)
    
    source1 = c(rep(1:nrow(W0) - 1 + num_nodes,ncol(W0)))
    target1 = c(rep(1:ncol(W0) - 1 + nrow(W0) + num_nodes,each = nrow(W0)))
    value1 = as.vector(W0)
    color1 = as.vector(W0_color)
    
    
    source.use = c(source.use, source1); 
    target.use = c(target.use, target1); 
    value.use = c(value.use, value1); 
    color.use0 = c(color.use0, color1)
    num_nodes = num_nodes + nrow(W0)
  }
  
  
  #value.use[value.use < 0.05] = 0
  color.use = rep(col_str,length(value.use))
  color.use[color.use0 == 1] = gsub("trans", 0.2, col_str) 
  color.use[color.use0 == 2] = gsub("trans", 0.4, col_str) 
  color.use[color.use0 == 3] =  gsub("trans", 0.6, col_str) 
  color.use[color.use0 == 4] =  gsub("trans", 0.8, col_str)
  
  
  #opacity.use = rep(0.2, length(value.use))
  #opacity.use[value.use == 0.01] = 0.001
  
  # Compute forward weights between every pair of clusters
  
  
  
  
  library(plotly)
  
  p <- plot_ly(
    type = "sankey",
    orientation = "h",
    
    node = list(
      label = labels.use,
      pad = 5,
      thickness = 5,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    
    link = list(
      source = source.use,
      target = target.use,
      value =  value.use,
      color = color.use
      
    )
  ) %>% 
    layout(
      title = title_str,
      font = list(
        size = 10
      )
    )
  
  return(p)
  
}


##################
### Bicolor scatter plot

BiMarkerScatter = function(object, feature1 = NULL, feature2 = NULL, max.cells = 8000, transparency=TRUE, do.legend=FALSE, thresh_F1 = NA, thresh_F2 = NA){
  library(paletteer)
  
  if (feature1 %in% rownames(object@data)){
    F1 = object@data[feature1,]  
  } else {
    F1 = object@data.info[,feature1]
  }
  
  if (feature2 %in% rownames(object@data)){
    F2 = object@data[feature2,]  
  } else {
    F2 = object@data.info[,feature2]
  }
  
  
  df = data.frame(F1 = F1, F2 = F2,UMAP1 = object@dr$UMAP@cell.embeddings[,1], UMAP2 = object@dr$UMAP@cell.embeddings[,2])
  
  if (is.na(thresh_F1)) thresh_F1 = mean(df$F1) + 2.6*sd(df$F1)
  if (is.na(thresh_F1)) thresh_F2 = mean(df$F2) + 2.6*sd(df$F2)
  
  df$F1[df$F1 > thresh_F1] = thresh_F1
  df$F2[df$F2 > thresh_F2] = thresh_F2
  
  df <- df %>%
    mutate(mix1 = rgb(red = F1/max(1,max(F1)), green = 0, blue = F2/max(1,max(F2))), transparency = F1/max(F1) + F2/max(F2))
  
  # Downsample
  if (nrow(df) > max.cells) df = df[sample(rownames(df), max.cells),]
  
  p = list()
  
  if (transparency){
    p[[1]] = ggplot(df, aes(x=UMAP1, y=UMAP2)) + geom_point(aes(color=mix1, alpha=transparency), size=0.7)+ scale_color_identity() + theme_classic() +  theme(legend.position = "none")
  } else {
    p[[1]] = ggplot(df, aes(x=UMAP1, y=UMAP2)) + geom_point(aes(color=mix1), size=0.7)+ scale_color_identity() + theme_classic() +  theme(legend.position = "none")
  }
  
  
  # Legend
  d <- expand.grid(x = seq(0, 1, 0.01), y = seq(0, 1, 0.01)) %>%
    mutate(fill_val = atan(y/x),
           transparency = 10*(x + y))
  d <- d %>%
    mutate(mix1 = rgb(red = x, green = 0, blue = y))
  
  if (transparency){
    p[[2]] <- ggplot(d, aes(x, y,fill = mix1, alpha=transparency)) +
      geom_tile() +
      scale_fill_identity() +
      labs(x = feature1,
           y = feature2) +
      theme_classic() +
      theme(axis.title = element_text(size = 15),
            axis.title.y = element_text(angle = 90)) +
      coord_fixed()
  } else {
    p[[2]] <- ggplot(d, aes(x, y,fill = mix1)) +
      geom_tile() +
      scale_fill_identity() +
      labs(x = feature1,
           y = feature2) +
      theme_classic() +
      theme(axis.title = element_text(size = 15),
            axis.title.y = element_text(angle = 90)) +
      coord_fixed()
  }
  
  if (do.legend){
    return(p)
  } else {
    return(p[[1]])
  }
}


getNsigmaBeta = function(x, n=3){
  library(fitdistrplus)
  fit = fitdist(x,"beta" )
  alpha1 = fit$estimate[1]
  beta1 = fit$estimate[2]
  
  mean_x = alpha1 / (alpha1 + beta1)
  sd_x = sqrt(alpha1*beta1 / ((alpha1 + beta1)^2 * (alpha1 + beta1 + 1)))
  return(mean_x + n*sd_x)
}

minmax=function(data,min,max) {
  data2=data
  data2[data2>max]=max
  data2[data2<min]=min
  return(data2)
}
