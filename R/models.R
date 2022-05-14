#' @include utils.R
#'
NULL

#' Calculate the Statistics for each Cell.
#'
#' Calculate the statistics, mean or median, for each cell.
#'
#' @keywords internal
#' @param data Expression data. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param idents Cell type identity in characters or factors. 
#' @param stats Statistics used to calculate strength.
#' \itemize{
#'   \item mean, using mean (Default).
#'   \item medain, using median.
#'   }
#' @param cutoff Cutoff used to calculate the statistics (mean or median). Default is NA. If set to 0, all 0 values will be removed when calculating mean or median.
#' @return Returns a matrix with calculated statistics.
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowMedians
statsCells <- function(data, idents, stats = c("mean", "median"), cutoff = NA){
  stats <- match.arg(arg = stats)
  stats.fxn <- switch(stats,
                      mean = rowMeans,
                      median = rowMedians)
  if(stats == "median") data <- as.matrix(x = data)
  if(!is.na(x = cutoff)){
    data[data <= cutoff] <- NA
    stats.mat <- sapply(X = unique(x = idents), function(id) stats.fxn(x = data[ , idents == id, drop = FALSE], na.rm = TRUE))
    stats.mat[is.na(x = stats.mat)] <- 0
  }  else {
    stats.mat <- sapply(X = unique(x = idents), function(id) stats.fxn(x = data[ , idents == id, drop = FALSE]))
  }
  rownames(x = stats.mat) <- rownames(x = data)
  colnames(x = stats.mat) <- unique(x = idents)
  return(stats.mat)
}

#' Calculate the Strength for Each Interacting Pair.
#'
#' Calculate the strength for each interacting pair.
#'
#' @keywords internal
#' @param ligand A ligand expression data.
#' @param receptor A receptor expression data. 
#' @param f Functions used to calculate pair strength.
#' \itemize{
#'   \item +, using sum of the pairs (Default).
#'   \item *, using product of the pairs.
#'   } 
#' @param threshold Threshold used for interaction products. Default is 0, therefore interaction product (ligand*receptor) equal or below 0 will be removed.
#' @return Returns calculated interacting strength.
calcStrength.single <- function(ligand, receptor, f = c("+", "*"), threshold = 0){
  add.pairs <- function(x, y){
    pairs <- (x + y)/2
    pairs[x * y <= threshold] <- 0
    return(pairs)
  }
  ligand <- ligand[,order(colnames(x = ligand)), drop=FALSE]
  receptor <- receptor[,order(colnames(x = receptor)), drop=FALSE]
  pair.names <- expand.grid(colnames(x = ligand), colnames(x = receptor), stringsAsFactors = FALSE)
  strength <- expand.grid(ligand, receptor)
  f <- match.arg(arg = f)
  strength <- switch(f,
                     "+" = add.pairs(strength[, 1], strength[, 2]),
                     "*" = strength[, 1] * strength[, 2])
  names(x = strength) <- apply(X = pair.names, MARGIN = 1, FUN = paste, collapse = "|")
  return(strength)
}

#' Calculate the Strength for Multiple Interacting Pair.
#'
#' Calculate the strength for multiple interacting pair. See details \code{\link[scMatchmaker]{calcStrength.single}}.
#'
#' @keywords internal
#' @param ligand_mat Ligand expression matrix.
#' @param receptor_mat Receptor expression matrix. 
#' @param f Functions used to calculate pair strength.
#' \itemize{
#'   \item +, using sum of the pairs (Default).
#'   \item *, using product of the pairs.
#'   } 
#' @param threshold Threshold used for interaction products. Default is 0, therefore interaction product (ligand*receptor) equal or below 0 will be removed.
#' @return Returns a matrix with calculated interacting strength.
calcStrength <- function(ligand_mat, receptor_mat, f = c("+", "*"), threshold = 0){
  n_pairs <- nrow(x = ligand_mat)
  pair.names <- paste(rownames(x = ligand_mat), rownames(x = receptor_mat), sep = "_")
  strength <- lapply(X = 1:n_pairs, FUN = function(x) calcStrength.single(ligand = ligand_mat[x, ,drop=FALSE], receptor = receptor_mat[x, ,drop=FALSE], f = f, threshold = threshold))
  strength <- do.call(what = cbind, args = strength)
  colnames(x = strength) <- pair.names
  return(strength)
}

#' Calculate the Interaction Strength for Single Interaction
#'
#' Calculate the Interaction Strength for a single interaction pairs.
#'
#' @keywords internal
#' @param data Expression data. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param idents Cell type identity in characters or factors. 
#' @param ligand Name or index for ligand.
#' @param receptor Name or index for receptor.
#' @param stats Method used to calculate strength.
#' \itemize{
#'   \item mean, using mean (Default).
#'   \item medain, using median.
#'   }
#' @param pair.fxn Functions used to calculate pair strength.
#' \itemize{
#'   \item +, using sum of the pairs (Default).
#'   \item *, using product of the pairs.
#'   } 
#' @param cutoff Cutoff used to calculate the statistics (mean or median). Default is NA. If set to 0, all 0 values will be removed when calculating mean or median.
#' @return Returns Earth Mover's Distance.
findInteractions.single <- function(data, idents, ligand, receptor, stats = c("mean", "median"), pair.fxn = c("+", "*"), 
                                    threshold = 0, cutoff = NA){
  stats.mat <- statsCells(data = data[unique(x = c(ligand, receptor)), ,drop=FALSE], idents = idents, stats = stats, cutoff = cutoff)
  strength.mat <- calcStrength(ligand_mat = stats.mat[ligand, ,drop=FALSE], receptor_mat = stats.mat[receptor, ,drop=FALSE], f = pair.fxn, threshold = threshold)
  return(strength.mat)
}

#' Calculate the Density of Feature
#'
#' Calculate the density of an input feature.
#'
#' @keywords internal
#' @param breaks A vector giving the breakpoints between histogram cells. See details \code{\link[graphics]{hist}}.
#' @param drop0 Whether to drop all the zeros. Default is FALSE.
#' @return Returns the density.
calcDensity <- function(x, breaks, drop0 = FALSE){
  if(drop0) x <- x[x>0]
  dens <- hist(x, breaks = breaks, plot = FALSE)
  return(dens$density)
}

#' Calculate the Earth Mover's Distance Similarity for Interaction Data
#'
#' Calculate the Earth Mover's Distance Similarity for interaction data. The similarity is defined as 1-normalized Earth Mover's Distance.
#'
#' @keywords internal
#' @param data Expression data. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param idents Cell type identity in characters or factors. 
#' @param ligand Name or index for ligand.
#' @param receptor Name or index for receptor.
#' @param nbins Number of bins to use for Earth Mover's Distance calculation. Default is 100.
#' @param weighted Logic, compute weighted Earth Mover's Distance which will give higher weight/penalty on lowly expressed genes. Default is TRUE.
#' @param min_pct Values below this percentile will be replaced with this percentile's value. Default is 0.
#' @param max_pct Values above this percentile will be replaced with this percentile's value. Default is 1.
#' @param drop_zero Whether to drop all the zeros in the data. Default is FALSE.
#' @return Returns Earth Mover's Distance matrix.
#' @importFrom emdist emd2d emdw
calcEMD <- function(data, idents, ligand, receptor, nbins = 100, weighted = TRUE, min_pct = 0, max_pct = 1, drop_zero = FALSE){
  idents <- as.character(x = idents)
  split.data <- split.data.frame(x = t.data.frame(x = data), f = idents)
  nbreaks <- seq(from = min(data, na.rm = TRUE), to = max(data, na.rm = TRUE), length.out = nbins + 1)
  split.data <- lapply(X = split.data, function(x){
    t(x = apply(X = x, MARGIN = 2, FUN = calcDensity, breaks = nbreaks, drop0 = drop_zero))
  })
  gene.pairs <- paste(ligand, receptor, sep = "_")
  cell.pairs <- expand.grid(sort(x = unique(x = idents)), sort(x = unique(x = idents)), stringsAsFactors = FALSE)
  if(weighted){
    weight_breaks <- matrix(data = rev(x = nbreaks[-1]))
    emd.data <- t(x = apply(X = cell.pairs, MARGIN = 1, FUN = function(x){
      data1 <- split.data[[x[1]]]
      data2 <- split.data[[x[2]]]
      unlist(x = lapply(X = 1:length(x = gene.pairs), FUN = function(x) emdw(A = t(x = data1[ligand[x], , drop = FALSE]), wA = weight_breaks,
                                                                             B = t(x = data2[receptor[x], , drop = FALSE]), wB = weight_breaks)))
    }))
  } else {
    emd.data <- t(x = apply(X = cell.pairs, MARGIN = 1, FUN = function(x){
      data1 <- split.data[[x[1]]]
      data2 <- split.data[[x[2]]]
      unlist(x = lapply(X = 1:length(x = gene.pairs), FUN = function(x) emd2d(A = t(x = data1[ligand[x], , drop = FALSE]),
                                                                              B = t(x = data2[receptor[x], , drop = FALSE]))))
    }))
  }
  dimnames(x = emd.data) <- list(paste(cell.pairs[,1], cell.pairs[,2], sep = "|"), gene.pairs)
  emd.data <- 1 - zeroOne(x = emd.data, min_pct = min_pct, max_pct = max_pct)
  emd.data[is.na(x = emd.data)] <- 0
  return(emd.data)
}

#' Calculate the Interaction Strength
#'
#' Calculate the interaction strength and statistical significance using random permutation tests.
#'
#' @keywords internal
#' @param data Expression data. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param idents Cell type identity in characters or factors. 
#' @param ligand Name or index for ligand.
#' @param receptor Name or index for receptor.
#' @param stats Method used to calculate strength.
#' \itemize{
#'   \item mean, using mean (Default).
#'   \item medain, using median.
#'   }
#' @param pair.fxn Functions used to calculate pair strength.
#' \itemize{
#'   \item +, using sum of the pairs (Default).
#'   \item *, using product of the pairs.
#'   } 
#' @param n_perm Number of random permutations. Default is 100.
#' @param emd Whethter to run Earth Mover's Distance adjustment. Default is FALSE.
#' @param nbins Number of bins to use for Earth Mover's Distance calculation. Default is 10.
#' @param weighted Logic, compute weighted Earth Mover's Distance which will give higher weight/penalty on lowly expressed genes. Default is TRUE.
#' @param p.adjust.method p value adjustment method. Default is "none". See details \code{\link[stats]{p.adjust}}.
#' @param threshold Threshold used for interaction products. See details \code{\link[scMatchmaker]{calcStrength}}.
#' @param min_pct Earth Movre's Distance values below this percentile will be replaced with this percentile value. Default is 0.
#' @param max_pct Earth Movre's Distance Values above this percentile will be replaced with this percentile value. Default is 1.
#' @param drop_zero Whether to drop all the zeros in the data when calculating EMD. Default is FALSE.
#' @param cutoff.stats Cutoff used to calculate the statistics (mean or median). Default is NA. If set to 0, all 0 values will be removed when calculating mean or median.
#' @param return.perm Whether to return the permutation result. Default is FALSE.
#' @param n_cores Number of cores used. Default is to use all cores - 1. See details \code{\link[parallel]{makeCluster}}.
#' @param seed Random seed number for permutation. Default is 1.
#' @return Returns a list contains interaction strength and pvalues, optionally the null strength and Earth Mover's Distance similarity.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom magrittr %>% divide_by add
#' @importFrom foreach foreach %dopar%
findInteractions <- function(data, idents, ligand, receptor, stats = c("mean", "median"), pair.fxn = c("+", "*"),  n_perm = 100, emd = FALSE, nbins = 100, weighted = TRUE,
                             p.adjust.method = "none", threshold = 0, min_pct = 0, max_pct = 1, drop_zero = FALSE, cutoff.stats = NA, return.perm = FALSE, n_cores = NULL, seed = 1){
  message("Calculating Interaction Strength")
  stats.null <- findInteractions.single(data = data, idents = idents, ligand = ligand, receptor = receptor,
                                        stats = stats, pair.fxn = pair.fxn, threshold = threshold, cutoff = cutoff.stats)
  if(emd){
    emd.null <- calcEMD(data = data, idents = idents, ligand = ligand, receptor = receptor, nbins = nbins, weighted = weighted, 
                        min_pct = min_pct, max_pct = max_pct, drop_zero = drop_zero)
    interaction.null <- stats.null * emd.null
  } else {
    interaction.null <- stats.null
  }
  message("Permutating")
  ident.mat <- randomIdents(idents = idents, n = n_perm)
  if(is.null(x = n_cores)) n_cores <- detectCores() - 1
  cl <- makeCluster(spec = n_cores, setup_strategy = "sequential") # R4.0.0 Compatibility
  registerDoParallel(cl = cl)
  if(emd){
    message("Running EMD")
    interaction.perm <- foreach(i = 1:n_perm, .multicombine = TRUE, .packages = "scMatchmaker") %dopar% {
      i <- get("i")
      stat.dist <- findInteractions.single(data = data, idents = ident.mat[ ,i], ligand = ligand, receptor = receptor,
                                           stats = stats, pair.fxn = pair.fxn, threshold = threshold, cutoff = cutoff.stats)
      emd.similarity <- calcEMD(data = data, idents = ident.mat[ ,i], ligand = ligand, receptor = receptor, nbins = nbins, 
                                weighted = weighted, min_pct = min_pct, max_pct = max_pct, drop_zero = drop_zero)
      stat.dist * emd.similarity
    }
  } else {
    interaction.perm <- foreach(i = 1:n_perm, .multicombine = TRUE, .packages = "scMatchmaker") %dopar% {
      i <- get("i")
      findInteractions.single(data = data, idents = ident.mat[ ,i], ligand = ligand, receptor = receptor,
                              stats = stats, pair.fxn = pair.fxn, threshold = threshold, cutoff = cutoff.stats)
    }
  }
  stopCluster(cl = cl)
  gc()
  message("Testing and adjusting")
  interaction.pval <- 
  lapply(X = interaction.perm, FUN = function(x) interaction.null <= x) %>% 
    Reduce(f = "+") %>% 
    add(1) %>% 
    divide_by(n_perm + 1) %>%
    p.adjust(method = p.adjust.method) %>%
    `dim<-`(dim(x = interaction.null)) %>%
    `dimnames<-`(dimnames(x = interaction.null))
  to.return <- list(strength = interaction.null, pvalues = interaction.pval)
  if(return.perm) {
    message("Saving permutation result")
    to.return <- c(to.return, permute_result = list(interaction.perm))
  }
  if(emd){
    return(c(to.return, stats_null = list(stats.null), emd = list(emd.null)))
  } else {
    return(to.return)
  }
}
