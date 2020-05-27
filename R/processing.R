#' @useDynLib scMatchmaker
#' @importFrom Rcpp sourceCpp
#' @include utils.R
NULL
#> NULL

#' Load CellPhoneDB Datasets
#'
#' Load the \href{https://www.cellphonedb.org/}{CellPhoneDB} dataset for ligand-receptor interaction pairs.
#'
#' @param interaction_input Data frame (csv file), directory path or URL of the interaction_input file.
#' @param gene_input Data frame (csv file), directory path or URL of the gene_input file.
#' @param complex_input Data frame (csv file), directory path or URL of the complex_input file.
#' @param gene.symbol Gene symbols to use. The possible names are listed as column names in the gene_input file.
#' Examples include: gene_name, uniprot, hgnc_symbol, ensembl. Default is gene_name.
#' @param annotation.strategy Annotation strategy used. Default is NULL, possible options include "curated" etc.
#' @param invert Whether to perfrom invert selection on annotation strategy. Default is FALSE
#' @param path Whether inputs are in directory path format.
#' @param url Whether inputs are in URl format.
#' @return Returns the interaction pairs data.frame with the first two columns correspond to interacting partners.
#' @importFrom utils read.csv
#' @export
#' @examples \dontrun{
#' # CellPhoneDB URLs 
#' pair.data <- LoadCellPhoneDB(
#'  interaction_input = interaction_input.url,
#'  gene_input = gene_input.url,
#'  complex_input = complex_input.url,
#'  gene.symbol = "ensembl", annotation.strategy = "curated", url = TRUE)
#' }
LoadCellPhoneDB <- function(interaction_input, gene_input, complex_input, gene.symbol = "gene_name", 
                            annotation.strategy = NULL, invert = FALSE, path = FALSE, url = FALSE){
  if(path){
    interaction_input <- read.csv(file = interaction_input, header = TRUE, stringsAsFactors = FALSE)
    gene_input <- read.csv(file = gene_input, header = TRUE, stringsAsFactors = FALSE)
    complex_input <- read.csv(file = complex_input, header = TRUE, stringsAsFactors = FALSE)
  } else if(url){
    interaction_input <- read.csv(file = url(description = interaction_input), header = TRUE, stringsAsFactors = FALSE)
    gene_input <- read.csv(file = url(description = gene_input), header = TRUE, stringsAsFactors = FALSE)
    complex_input <- read.csv(file = url(description = complex_input), header = TRUE, stringsAsFactors = FALSE)
  }
  # CellPhoneDB error in hgnc_symbol P01344
  gene_input[which(gene_input[[gene.symbol]] == ""), gene.symbol] <- gene_input[which(gene_input[[gene.symbol]] == ""), "gene_name"]
  # CellPhoneDB error in HLA uniprot
  hla.uniprot <- c("P04439", "P01889", "P10321", "P20036", "P04440", "P05538", "P01911")
  names(x = hla.uniprot) <- c("HLAA", "HLAB", "HLAC", "HLADPA1", "HLADPB1", "HLADQB2", "HLADRB1")
  gene_input$uniprot[which(x = gene_input$uniprot %in% names(x = hla.uniprot))] <- hla.uniprot
  hla.match <- interaction_input$partner_a[which(x = interaction_input$partner_a %in% names(x = hla.uniprot))]
  interaction_input$partner_a[which(x = interaction_input$partner_a %in% names(x = hla.uniprot))] <- hla.uniprot[hla.match]
  
  pairs <- interaction_input[,c("partner_a", "partner_b")]
  genename_a <- gene_input$uniprot[match(x = pairs$partner_a, table = gene_input$uniprot)]
  genename_b <- gene_input$uniprot[match(x = pairs$partner_b, table = gene_input$uniprot)]
  complex_a <- complex_input[,grep(pattern = "uniprot", x = colnames(x = complex_input))][match(x = pairs$partner_a[is.na(x = genename_a)], table = complex_input$complex_name),]
  complex_b <- complex_input[,grep(pattern = "uniprot", x = colnames(x = complex_input))][match(x = pairs$partner_b[is.na(x = genename_b)], table = complex_input$complex_name),]
  complex_a <- apply(X = complex_a, MARGIN = 1, FUN = function(x) paste(x[!is.na(x = x) & x!= ""], collapse = ","))
  complex_b <- apply(X = complex_b, MARGIN = 1, FUN = function(x) paste(x[!is.na(x = x) & x!= ""], collapse = ","))
  genename_a[is.na(x = genename_a)] = complex_a
  genename_b[is.na(x = genename_b)] = complex_b
  
  data.comb <- cbind.data.frame(genename_a, genename_b, interaction_input, stringsAsFactors = FALSE)
  split_a <- strsplit(x = data.comb$genename_a, split = ",")
  split_b <- strsplit(x = data.comb$genename_b, split = ",")
  complex_a <- unlist(x = lapply(X = split_a, FUN = length))
  complex_b <- unlist(x = lapply(X = split_b, FUN = length))
  if(max(complex_a) > 1){
    subunit_a <- paste("subunit_a", 1:max(complex_a), sep = "_")
    data.temp <- as.data.frame(x = matrix("", ncol = max(complex_a)), stringsAsFactors = FALSE)
    colnames(x = data.temp) <-  subunit_a
    data.comb <- cbind(data.comb, data.temp)
    data.temp <- lapply(X = which(x = complex_a > 1), FUN = function(x) unlist(x = strsplit(x = data.comb[x, "genename_a"], split = ",")))
    data.comb[which(x = complex_a > 1), subunit_a] <- t(x = data.frame(lapply(X = data.temp, FUN = "length<-", max(lengths(x = data.temp)))))
  }
  if(max(complex_b) > 1){
    subunit_b <- paste("subunit_b", 1:max(complex_b), sep = "_")
    data.temp <- as.data.frame(x = matrix("", ncol = max(complex_b)), stringsAsFactors = FALSE)
    colnames(x = data.temp) <-  subunit_b
    data.comb <- cbind(data.comb, data.temp)
    data.temp <- lapply(X = which(x = complex_b > 1), FUN = function(x) unlist(x = strsplit(x = data.comb[x, "genename_b"], split = ",")))
    data.comb[which(x = complex_b > 1), subunit_b] <- t(x = data.frame(lapply(X = data.temp, FUN = "length<-", max(lengths(x = data.temp)))))
  }
  
  split_b <- strsplit(x = data.comb$genename_b, split = ",")
  data.comb <- data.frame(genename_b = unlist(x = split_b), data.comb[rep(x = 1:nrow(x = data.comb), sapply(X = split_b, FUN = length)), -2], stringsAsFactors = FALSE)
  split_a <- strsplit(x = data.comb$genename_a, split = ",")
  data.comb <- data.frame(genename_a = unlist(x = split_a), data.comb[rep(x = 1:nrow(x = data.comb), sapply(X = split_a, FUN = length)), -2], stringsAsFactors = FALSE)
  data.sort <- t(x = apply(X = data.comb[,1:2], MARGIN = 1, FUN = sort))
  data.comb <- data.comb[!duplicated(x = data.sort), ]
  data.comb$genename_a <- gene_input[[gene.symbol]][match(x = data.comb$genename_a, table = gene_input$uniprot)]
  data.comb$genename_b <- gene_input[[gene.symbol]][match(x = data.comb$genename_b, table = gene_input$uniprot)]
  if(max(complex_a) > 1){
    data.comb[, subunit_a] <- apply(X = data.comb[, subunit_a], MARGIN = 2, FUN = function(x) gene_input[[gene.symbol]][match(x = x, table = gene_input$uniprot)])
  }
  if(max(complex_b) > 1){
    data.comb[, subunit_b] <- apply(X = data.comb[, subunit_b], MARGIN = 2, FUN = function(x) gene_input[[gene.symbol]][match(x = x, table = gene_input$uniprot)])
  }
  data.comb <- data.comb[order(data.comb$genename_a), ]
  if(!is.null(x = annotation.strategy)){
    if(invert){
      data.comb <- data.comb[data.comb[["annotation_strategy"]] != annotation.strategy,]
    } else {
      data.comb <- data.comb[data.comb[["annotation_strategy"]] == annotation.strategy,]
    }
  } 
  rownames(x = data.comb) <- 1:nrow(x = data.comb)
  return(data.comb)
}

#' Data Normalization
#'
#' Normalize count or other data.
#'
#' @param counts Data to normalize. An N x M matrix with N rows of features and M columns of data points.
#' @param normalization Normalization method used. Default is cosine.
#' \itemize{
#'   \item logTPM, logged Transcript Per Million normalization: feature counts for each data point are divided by the total sum of them. Then the data is multiplied by the scale.factor before taking a log-transformed by log(1+x) (Default).
#'   \item cosine, cosine nomalization: feature counts for each data point are divided by their L2.
#'   \item none, additional normalization is not performed.
#'   }
#' @param normalize_factor Normalization factor used with lognorm method. Default is 10000.
#' @param zero_percent Zero-entry percentage threshold. If the number of zero entries in the returned matrices is above this number, a sparse matrix will be returned. Default is 0.7 aka 70\%.
#' @param verbose Whether to display a process bar for cosine normalization. Default is FALSE.
#' @return Returns the normalized data.
#' @export
#' @importFrom Matrix nnzero Matrix
#' @examples \dontrun{
#' normalized.data <- Normalization(data)
#' }
Normalization <- function(counts, normalization = c("logTPM", "cosine", "none"), normalize_factor = 1e4, zero_percent = 0.7, verbose = FALSE){
  normalization <- match.arg(arg = normalization)
  if(any(class(x = counts) %in% c("dgCMatrix", "dsCMatrix"))){
    normalized.data <- switch(normalization,
                              logTPM = logTPM(mat = counts, normalize_factor = normalize_factor),
                              cosine = CosineNormSparse(data = counts, display_progress = verbose),
                              none = counts)
  } else {
    normalized.data <- switch(normalization,
                              logTPM = logTPM(mat = counts, normalize_factor = normalize_factor),
                              cosine = CosineNorm(data = counts, display_progress = verbose),
                              none = counts)
  }
  normalized.data <- asSparse(mat = normalized.data, zero_percent = zero_percent)
  dimnames(x = normalized.data) <- dimnames(x = counts)
  return(normalized.data)
}

#' Data Scaling
#'
#' Scale and center data.
#'
#' @keywords internal
#' @param mat Data matrix.
#' @param scale_max Maximum cutoff for scaled data. Default is 10.
#' @return Returns a scaled and centered data.
#' @examples \dontrun{
#' scaled.data <- Scaling(data)
#' }
Scaling <- function(mat, scale_max = 10){
  if(!any(class(x = mat) == "matrix")){
    mat <- as.matrix(x = mat)
  }
  scale.data <- rowScale(mat = mat, scale_max = scale_max)
  dimnames(x = scale.data) <- dimnames(x = mat)
  return(scale.data)
}

#' Run Geometric Sketching
#'
#' Run Geometric sketching sampling.
#' See \href{https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30152-8}{Hie et al. 2019 Cell Systems} for details.
#'
#' @keywords internal
#' @param data An M x d sparse matrix or data.frame with M rows of data points and d columns of features.
#' @param geom_size Size of geometric sketches to return. Default is 1000.
#' @param n_pca Number of PCA dimensions to use. Default is 10.
#' @param seed Random seed number for PCA. Default is 1.
#' @param which_python Path to python3 used.
#' @param scale Whether to scale the data. Default is TRUE.
#' @return Returns geometric sketch IDs.
#' @importFrom reticulate use_python py_module_available import
#' @importFrom irlba irlba
#' @importFrom Matrix rowMeans
geomSketch <- function(data, geom_size = 1000, n_pca = 10, scale = TRUE, seed = 1, which_python = Sys.which(names = "python3")){
  use_python(python = which_python, required = TRUE)
  if(!py_module_available("geosketch")){
    stop("Cannot find geosketch, please install through pip (e.g. pip install geosketch).")
  }
  message("Preparing")
  feature.means <- rowMeans(x = data)
  feature.sd <- apply(X = data, MARGIN = 1, FUN = sd)
  if(any(feature.means == 0 | feature.sd == 0)){
    message("Removing Missing Features")
    data <- data[-which(feature.means == 0 | feature.sd == 0),]
  } 
  geosketch <- import(module = 'geosketch')
  if(scale){
    message("Scaling")
    data <- Scaling(mat = data, scale_max = Inf)
  }
  message("Run PCA")
  set.seed(seed = seed)
  pca.res <- irlba(A = t(x = data), nv = n_pca)
  message("Run Geometric Sketch")
  sketch.ind <- unlist(x = geosketch$gs(pca.res$u, as.integer(x = geom_size))) + 1
  return(sketch.ind)
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

#' Filter Interaction Data
#'
#' Filter interaction data based on interacting strength and p value cutoffs.
#' 
#' @keywords internal
#' @param interact.strength Interaction strength matrix.
#' @param interact.pval Interaction p-value matrix.
#' @param interaction.meta Interaction metadata with first two columns as interaction partners matching interact.strength and interact.pval data.
#' @param pval.cutoff p value cutoff. Default is 0.05.
#' @param strength.pct Top strength quantile to select. Default is 0.1 aka the top 10\%.
#' @param filter.cells Whether to filter the cells. Default is TRUE.
#' @return Returns a list containing filtered interaction strength, p value matrix and merged interaction data.
#' @importFrom Matrix rowSums colSums
filterInteractions <- function(interact.strength, interact.pval, interaction.meta, pval.cutoff = 0.05, strength.pct = 0.1, filter.cells = TRUE){
  filtered.str <- matrix(data = 0, nrow = nrow(x = interact.strength), ncol = ncol(x = interact.strength), dimnames = dimnames(x = interact.strength))
  filtered.str[as.matrix(x = interact.strength) >= quantile(x = as.matrix(x = interact.strength), probs = 1-strength.pct)] <- 1
  filtered.pval <- matrix(data = 0, nrow = nrow(x = interact.pval), ncol = ncol(x = interact.pval), dimnames = dimnames(x = interact.pval))
  filtered.pval[interact.pval <= pval.cutoff] <- 1
  filtered.id <- filtered.pval * filtered.str
  interact.strength.new <- filtered.id * interact.strength
  interact.pval.new <- filtered.id * interact.pval
  interact.strength.new <- interact.strength.new[, which(x = colSums(x = filtered.id) != 0), drop = FALSE]
  interact.pval.new <- interact.pval.new[, which(x = colSums(x = filtered.id) != 0), drop = FALSE]
  interaction.meta <- interaction.meta[which(x = colSums(x = filtered.id) != 0), , drop = FALSE]
  if(filter.cells){
    interact.strength.new <- interact.strength.new[which(x = rowSums(x = filtered.id) != 0),,drop = FALSE]
    interact.pval.new <- interact.pval.new[which(x = rowSums(x = filtered.id) != 0),,drop = FALSE]
  }
  row.order <- order(rowSums(x = interact.strength.new), decreasing = TRUE)
  col.order <- order(colSums(x = interact.strength.new), decreasing = TRUE)
  interact.strength.new <- interact.strength.new[row.order, col.order, drop = FALSE]
  interact.pval.new <- interact.pval.new[row.order, col.order, drop = FALSE]
  interaction.meta <- interaction.meta[col.order, ,drop = FALSE]
  merge.data <- cbind.data.frame(interaction.meta, t(x = as.matrix(x = interact.strength.new)), row.names = 1:nrow(x = interaction.meta))
  return(list(strength = interact.strength.new, 
              pvalues = interact.pval.new,
              merge.data = merge.data))
}

#' Convert Interaction Matrix
#'
#' Convert interaction matrix into a molten data frame. Interactions with zero strenghth will be removed.
#' 
#' @keywords internal
#' @param interact.strength Interaction strength matrix.
#' @param interact.pval Interaction p-value matrix.
#' @return Returns a molten data frame.
#' @importFrom reshape2 melt
convertData <- function(interact.strength, interact.pval){
  strength.df <- data.frame(melt(data = as.matrix(interact.strength), varnames = c("interactions", "pairs"), value.name = "strength"))
  pval.df <- data.frame(melt(data = as.matrix(interact.pval), varnames = c("interactions", "pairs"), value.name = "pval"))
  merge.df <- cbind.data.frame(strength.df, pval = pval.df$pval)
  merge.df <- merge.df[order(merge.df$strength, merge.df$pval, decreasing = TRUE),]
  merge.df <- merge.df[which(x = merge.df$strength > 0), ]
  rownames(x = merge.df) <- 1:nrow(x = merge.df)
  return(merge.df)
}