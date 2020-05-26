#' @useDynLib scMatchmaker
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dist p.adjust quantile sd
NULL
#> NULL

#' Log-transformed Transcripts Per Million (logTPM) 
#'
#' logTPM normalization using log base of 2.
#'
#' @keywords internal
#' @param mat Count matrix. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param normalize_factor Normalization factor used. Default is 10000.
#' @return Returns a regular or sparse matrix.
#' @importFrom  Matrix colSums t
logTPM <- function(mat, normalize_factor = 1e4){
  return(log(x = t(x = t(x = mat) * normalize_factor / colSums(x = mat)) + 1, base = 2))
}

#' Scale Data by Row 
#'
#' Center and scale the expression matrix by its rows (features).
#'
#' @keywords internal
#' @param mat A matrix.
#' @param scale_max Maximum cutoff for scaled data. Default is 10.
#' @return Returns a matrix.
#' @importFrom  matrixStats rowMeans2 rowSds
rowScale <- function(mat, scale_max = 10) {
  row.means <- rowMeans2(x = mat)
  row.sds <- rowSds(x = mat, center = row.means)
  mat <- (mat - row.means) / row.sds
  mat <- pmin(mat, scale_max, na.rm = FALSE)
  return(mat)
}

#' Convert Matrix to Sparse Matrix
#'
#' Convert Matrix to Sparse Matrix
#'
#' @keywords internal
#' @param mat Expression data. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param zero_percent Zero-entry percentage threshold. If the number of zero entries in the returned matrices is above this number, a sparse matrix will be returned. Default is 0.7 aka 70\%.
#' @return Returns a regular or sparse matrix.
#' @importFrom  Matrix nnzero Matrix 
asSparse <- function(mat, zero_percent = 0.7){
  if(!any(class(x = mat) %in% c("dgCMatrix", "dsCMatrix"))){
    if(nnzero(x = mat, na.counted = FALSE)/length(x = mat) < (1 - zero_percent)){
      mat <- Matrix(data = mat, sparse = TRUE)
    }
  }
  return(mat)
}

#' Convert to Zero and One Range
#'
#' Convert data to values between zero and one.
#'
#' @keywords internal
#' @param x Data.
#' @param min_pct Values below this percentile will be replaced with this percentile value. Default is 0.
#' @param max_pct Values above this percentile will be replaced with this percentile value. Default is 1.
#' @return Returns new values in zero and one range.
zeroOne <- function(x, min_pct = 0, max_pct = 1){
  min_cut <- quantile(x = x, probs = min_pct, na.rm = TRUE)
  max_cut <- quantile(x = x, probs = max_pct, na.rm = TRUE)
  x[x <= min_cut] <- min_cut
  x[x >= max_cut] <- max_cut
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

#' Filter Matrix
#'
#' @keywords internal
#' @param matrix A matrix to be filtered.
#' @param low.cutoff Lower cutoff bound below which will be removed. Default is negative infinity.
#' @param high.cutoff Upper cutoff bound above which will be removed. Default is positive infinity.
#' @return Returns a filtered matrix.
filterMatrix <- function(matrix, low.cutoff = -Inf, high.cutoff = Inf){
  row.id <- apply(X = matrix, MARGIN = 1, FUN = function(x) any(x > low.cutoff & x <= high.cutoff))
  col.id <- apply(X = matrix, MARGIN = 2, FUN = function(x) any(x > low.cutoff & x <= high.cutoff))
  matrix <- matrix[row.id, col.id, drop = FALSE]
  if(nrow(x = matrix) == 0) stop("No data passed the filter.", call. = FALSE)
  return(matrix)
}

#' Match the Names of Inteactions Pairs in Expression Data.
#'
#' Match and filter the expression data and interacton pairs.
#'
#' @keywords internal
#' @param data Expression data. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param pairs Interaction pair data with at least two columns indicating gene-gene interaction partners.
#' @param partner_a Column name or index for interaction partner a. Defaulte is 1, the first column in interaction pair data.
#' @param partner_b Column name or index for interaction partner b. Defaulte is 2, the second column in interaction pairs data.
#' @return Returns a list with filtered expression data and interaction pairs.
matchNames <- function(data, pairs, partner_a = 1, partner_b = 2){
  features <- rownames(x = data)
  partner.A <- as.character(x = pairs[, partner_a])
  partner.B <- as.character(x = pairs[, partner_b])
  partner.all <- unique(x = c(partner.A, partner.B))
  partner.A.in <- partner.A %in% features
  partner.B.in <- partner.B %in% features
  partner.all.in <- apply(X = cbind(partner.A.in,partner.A.in), MARGIN = 1, FUN = all)
  if(!any(partner.all.in)) stop("No interaction pair found", call. = FALSE)
  partner.detected <- partner.all[partner.all %in% features]
  partner.not.detected <- partner.all[!partner.all %in% features]
  data <- data[partner.detected, ]
  if(length(x = partner.not.detected) > 0){
    message("Partners: ",paste(partner.not.detected, collapse = ", "),
            " not detected.\nThey will be mannually added to the data.")
    partner.add <- matrix(0, nrow = length(x = partner.not.detected), ncol = ncol(x = data), 
                          dimnames = list(partner.not.detected, colnames(x = data)))
    data <- rbind(data, partner.add)
  } 
  return(list(data = data, pairs = pairs))
} 

#' Check Cell Identity in Data
#'
#' Check if cell identities provided matched to the single cell data and reorder them. Any identity name with "_: will be subsituted as ".".
#'
#' @keywords internal
#' @param data Expression data. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param idents Cell type identity in characters or factors. 
#' @return Returns the checked identity.
checkIdents <- function(data, idents){
  if(length(x = idents) != ncol(x = data)) stop("Idents provided don't match the number of cells.", call. = FALSE)
  if(any(is.na(x = idents))) stop("NA value detected in idents.", call. = FALSE)
  if(!is.null(x = names(idents))){
    if(!all(names(x = idents) %in% colnames(x = data))) stop("Ident names don't match the cell names.", call. = FALSE)
    idents <- idents[colnames(x = data)]
  }
  message("Replace names containing _ and | with .")
  idents <- gsub(pattern = "_|[|]", replacement = ".", x = idents)
  idents <- as.factor(x = idents)
  return(idents)
}

#' Merge Two List
#'
#' Merge twp lists by the name of their elements.
#'
#' @keywords internal
#' @param a,b Two lists.
#' @return Returns the merged list.
mergeLists <- function(a, b) {
  a.names <- names(a)
  b.names <- names(b)
  m.names <- unique(c(a.names, b.names))
  sapply(m.names, function(i) {
    if (i %in% b.names) b[[i]]
    else a[[i]]
  }, simplify = FALSE)
}

#' Aggregate Interaction Data
#' 
#' Aggregate interaction data by row or column names.
#' 
#' @keywords internal
#' @param data Interaction data with rows of cell-cell paris and columns of gene-gene pairs.
#' @param aggregate Type of aggregation used:
#' \itemize{
#'   \item row, using rownames (Default).
#'   \item column, using colnames.
#'   }
#' @param cutoff Cutoff for interactions. Default is 0. Only values that are above this threshold will be considered.
#' @param counted Binarize data to aggregate. Default is FALSE.
#' @param split.by interaction pairs. Default is "_" or "|".
#' @return Returns the aggregated matrix.
aggregateName <- function(data, aggregate = c("row", "column"), counted = FALSE, cutoff = 0, split.by = "_|[|]"){
  if(counted) data[data > cutoff] <- 1
  data[data <= cutoff] <- 0
  aggregate <- match.arg(arg = aggregate)
  aggregate.data <- switch(aggregate,
                           row = rowSums(x = data),
                           column = colSums(x = data))
  name.list <- strsplit(x = names(x = aggregate.data), split = split.by)
  name.unique <- sort(x = unique(x = unlist(x = name.list)))
  aggregate.matrix <- matrix(data = 0, nrow = length(x = name.unique), ncol = length(x = name.unique), dimnames = list(name.unique, name.unique))
  for(i in 1:length(name.list)){
    aggregate.matrix[name.list[[i]][1], name.list[[i]][2]] <- aggregate.data[i]
  }
  return(aggregate.matrix)
}

#' Generate Randomized Cell Identities
#'
#' @keywords internal
#' @param idents Cell type identity in characters or factors. 
#' @param n Number of random trials. Default is 100.
#' @param seed Random seed number. Default is 1.
#' @return Returns a matrix with randomized cell identities.
randomIdents <- function(idents, n = 100, seed = 1){
  set.seed(seed = seed)
  return(replicate(n = n, expr = sample(x = idents)))
}

#' Order Matrix
#'
#' @keywords internal
#' @param matrix A matrix. 
#' @param order.dist Distance measure to be used. Default is "minkowski". See details \code{\link[stats]{dist}}.
#' @param p Power of Minkowski distance Default is 1, aka Manhattan distance. See details \code{\link[stats]{dist}}.
#' @param return.order Whether to return the matrix order. Default is FALSE.
#' @return Returns a matrix or list containing ordered matrix.
#' @details See \code{\link[seriation]{seriate}}.
#' @importFrom seriation seriate permute 
orderMatrix <- function(matrix, order.dist = "minkowski", p = 1, return.order = FALSE){
  order.use <- c(
    seriate(x = dist(x = matrix, method = order.dist, p = p)),
    seriate(x = dist(x = t(x = matrix), method = order.dist, p = p))
  )  
  ordered.matrix <- permute(x = matrix, order = order.use)
  if(return.order){
    return(list(ordered_matrix = ordered.matrix, order = order.use))
  } else {
    return(ordered.matrix)
  }
}

#' Log Command
#'
#' @keywords internal
#' @param pos Argument positions to skip. Default is 1.
#' @param frame Number of parental frame to go back. Default is 1. See \code{\link[base]{sys.frame}} for details.
#' @return Returns a list of command arguments used.
logCommand <- function(pos = 1, frame = 1){
  argument.list <- mget(x = names(x = formals(sys.function(sys.parent(n = frame)))), envir = sys.frame(which = sys.nframe()-frame))[-pos]
  argument.list$Time <- Sys.time()
  return(argument.list)
}

