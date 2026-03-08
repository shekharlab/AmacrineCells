
#' The SmartMatrix Class
#'
#' @slot matrix
#' @slot meta.data Contains meta-information 
setClass(
  Class = 'SmartMatrix',
  slots = c(
    matrix = 'matrix',
    row.data = 'data.frame',
    col.data = 'data.frame',
    misc = 'list'
  )
)

SmartMatrix = function(matrix, row.data = NULL, col.data = NULL, misc = list(NULL)){
  
  rownames(row.data) = rownames(matrix)
  rownames(col.data) = colnames(matrix)
  
  new("SmartMatrix", 
      matrix = matrix, 
      row.data = row.data,
      col.data = col.data, 
      misc = misc)
}

# "[.Smartmatrix" <- function(x, i, j, ...) {
#   if (missing(x = i) && missing(x = j)) {
#     return(x)
#   }
#   if (missing(x = i)) {
#     i <- NULL
#   } else if (missing(x = j)) {
#     j <- colnames(x = x)
#   }
#   if (is.logical(x = i)) {
#     if (length(i) != nrow(x = x)) {
#       stop("Incorrect number of logical values provided to subset features")
#     }
#     i <- rownames(x = x)[i]
#   }
#   if (is.logical(x = j)) {
#     if (length(j) != ncol(x = x)) {
#       stop("Incorrect number of logical values provided to subset cells")
#     }
#     j <- colnames(x = x)[j]
#   }
#   if (is.numeric(x = i)) {
#     i <- rownames(x = x)[i]
#   }
#   if (is.numeric(x = j)) {
#     j <- colnames(x = x)[j]
#   }
#   return(subset.SmartMatrix(x = x, rows = i, columns = j, ...))
# }

setMethod(
  f = "[",
  signature = "SmartMatrix",
  definition = function(x, i, j, ..., drop = FALSE) {
    
      # if (missing(x = i) && missing(x = j)) {
      #   return(x)
      # }
      # if (missing(x = i)) {
      #   i <- NULL
      # } else if (missing(x = j)) {
      #   j <- colnames(x = x)
      # }
      # if (is.logical(x = i)) {
      #   if (length(i) != nrow(x = x)) {
      #     stop("Incorrect number of logical values provided to subset features")
      #   }
      #   i <- rownames(x = x)[i]
      # }
      # if (is.logical(x = j)) {
      #   if (length(j) != ncol(x = x)) {
      #     stop("Incorrect number of logical values provided to subset cells")
      #   }
      #   j <- colnames(x = x)[j]
      # }
      # if (is.numeric(x = i)) {
      #   i <- rownames(x = x)[i]
      # }
      # if (is.numeric(x = j)) {
      #   j <- colnames(x = x)[j]
      # }
    
      x@matrix = x@matrix[i, j, drop = FALSE]
      x@row.data = x@row.data[i, , drop = FALSE]
      x@col.data = x@col.data[j, , drop = FALSE]
    
      return(x)
    }
)

SubsetCols = function(smatrix, column, value){
  
  stopifnot(length(value) == 1)
  
  new.col.data = subset(smatrix@col.data, eval(parse(text = column)) == value)
  
  smatrix2 = smatrix[,rownames(new.col.data)]
  smatrix2
}

SubsetColsIn = function(smatrix, column, values){
  
  # Identify if value is in columns or rows
  
  new.col.data = subset(smatrix@col.data, eval(parse(text = column)) %in% values)
  
  smatrix2 = smatrix[,rownames(new.col.data)]
  smatrix2
}


MakeSmartMatrix = function(table, delimiter.row = c('_'), delimiter.col = c('-'), names = c('x', 'y')){
  SmartMatrix(unclass(as.matrix(table)), 
              row.data = data.frame(x = ExtractString(rownames(table), after = delimiter.row), 
                                    y = ExtractString(rownames(table), before = delimiter.row)) %>% setNames(names), 
              col.data = data.frame(x = ExtractString(colnames(table), after = delimiter.col), 
                                    y = ExtractString(colnames(table), before = delimiter.col)) %>% setNames(names))
}

setMethod(
  f = "+",
  signature = c(e1 = "SmartMatrix", e2 = "SmartMatrix"),
  definition = function(e1, e2) {
    new("SmartMatrix", 
        matrix = e1@matrix + e2@matrix, 
        row.data = e1@row.data, 
        col.data = e1@col.data, 
        misc = e1@misc)
  }
)

setMethod(
  f = "/",
  signature = c(e1 = "SmartMatrix", e2 = "numeric"),
  definition = function(e1, e2) {
    new("SmartMatrix", 
        matrix = e1@matrix / e2, 
        row.data = e1@row.data, 
        col.data = e1@col.data, 
        misc = e1@misc)
  }
)
