#' @include utils.R models.R processing.R 
#'
NULL

#' Screen for Interactions and Initialize Matchmaker Object
#'
#' Screen for interactiins and initialize a Matchmaker object.
#' @param data Expression data. A d x M matrix with d rows of features and M columns of data points (cells).
#' @param annotation Cell type identity in following formats:
#' \itemize{
#'   \item character, cell type identity with either names as cell ids (barcodes) or in the same order matching the column names in data.
#'   \item data.frame, cell ids (barcodes) in row names and cell type identity in first column.
#'   }
#' @param interaction Interaction partner pair data with at least two columns indicating gene-gene interaction partners.
#' @param partner_a Column name or index for interaction partner a. Defaulte is 1, the first column in interaction pair data.
#' @param partner_b Column name or index for interaction partner b. Defaulte is 2, the second column in interaction pairs data.
#' @param project_name Name of the project. 
#' @param path Whether inputs are in directory path format.
#' @param zero_percent Zero-entry percentage threshold. If the number of zero entries in the returned matrices is above this number, a sparse matrix will be returned. Default is 0.7 aka 70\%.
#' @return Returns a Matchmaker object.
#' @importFrom methods new
#' @importFrom utils read.csv packageVersion
#' @export
#' @examples \dontrun{
#' object <- Screening(data = data.use, annotation = idents.use, interaction = interaction.data)
#' }
Screening <- function(data, annotation, interaction, partner_a = 1, partner_b = 2, path = FALSE, project_name = "", zero_percent = 0.7){
  message("Checking for interaction pairs...")
  if(path){
    raw.data <- read.csv(file = data, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    annotation <- read.csv(file = annotation, header = TRUE, stringsAsFactors = TRUE, row.names = 1)[,1]
    interaction <- read.csv(file = interaction, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  }
  if(is.data.frame(x = annotation)){
    cell.id <- rownames(x = annotation)
    idents <- as.character(annotation[,1])
    names(x = idents) <- cell.id
  } else {
    idents <- annotation
    annotation <- data.frame(idents = annotation)
    annotation[,1] <- as.factor(x = annotation[,1])
  }
  raw.data <- asSparse(mat = data, zero_percent = zero_percent)
  idents <- checkIdents(data = raw.data, idents = idents)
  annotation[,1] <- idents
  data.pairs <- matchNames(data = raw.data, pairs = interaction, partner_a = partner_a, partner_b = partner_b)
  object <- new(Class = "Matchmaker", data = data.pairs$data, interaction = data.pairs$pairs, annotation = annotation, project_name = project_name, misc = list(raw_data = raw.data))
  object@command$Screening <- logCommand(pos = 1:3)
  object@command$Version <- packageVersion(pkg = "scMatchmaker")
  return(object)
}

#' Run Matchmaker 
#'
#' Run matchmaker to Identify cell-cell interactions.
#' 
#' @param object Matchmaker object.
#' @param ident1 Cell type identity group 1. Default is NULL and all cell types will be used.
#' @param ident2 Cell type identity group 2. Default is NULL and all cell types will be used.
#' @param stats Statistics used to calculate strength.
#' \itemize{
#'   \item mean, using mean (Default).
#'   \item medain, using median.
#'   }
#' @param pair.fxn Functions used to calculate pair strength.
#' \itemize{
#'   \item +, using sum of the pairs (Default).
#'   \item *, using product of the pairs.
#'   } 
#' @param emd Whethter to run Earth Mover's Distance adjustment. Default is FALSE.
#' @param nbins Number of bins to use for Earth Mover's Distance calculation. Default is 100.
#' @param n_perm Number of random permutations. Default is 100.
#' @param weighted Logic, compute weighted Earth Mover's Distance which will give higher weight/penalty on lowly expressed genes. Default is FALSE.
#' @param p.adjust.method p value adjustment method. Default is "none". See details \code{\link[stats]{p.adjust}}.
#' @param min_emd_pct Earth Mover's Distance below this percentile will be replaced with this percentile value. Default is 0.
#' @param max_emd_pct Earth Mover's Distance above this percentile will be replaced with this percentile value. Default is 1.
#' @param drop_zero Whether to drop all the zeros in the data when calculating EMD. Default is FALSE.
#' @param cutoff.stats Cutoff used to calculate the statistics (mean or median). Default is NA. If set to 0, all 0 values will be removed when calculating mean or median.
#' @param save.perm Whether to save permutation result. Default is FALSE.
#' @param n_cores Number of cores used. Default is to use all cores - 1. See details \code{\link[parallel]{makeCluster}}.
#' @param seed Random seed number for permutation. Default is 1.
#' @return Returns a Matchmaker object.
#' @export
#' @examples \dontrun{
#' object <- Matchmaking(object, stats = "mean", emd = TRUE, n_perm = 100)
#' }
Matchmaking <- function(object, ident1 = NULL, ident2 = NULL, stats = c("mean", "median"), pair.fxn = c("+", "*"), emd = FALSE, nbins = 100, n_perm = 100, weighted = FALSE, 
                        p.adjust.method = "none", min_emd_pct = 0, max_emd_pct = 1, drop_zero = FALSE, cutoff.stats = NA, save.perm = FALSE, n_cores = NULL, seed = 1){
  if(!(is.null(x = ident1) & is.null(x = ident2))){
    if(is.null(x = ident1)) ident1 <- levels(x = object@object@annotation[,1])
    if(is.null(x = ident2)) ident2 <- levels(x = object@object@annotation[,1])
    idx.use <- which(x = object@annotation[,1] %in% unique(x = c(ident1, ident2)))
    data.use <- object@data[, idx.use, drop = FALSE]
    idents.use <- object@object@annotation[idx.use,1]
  } else {
    data.use <- object@data
    idents.use <- object@annotation[,1]
  }
  message("Matchmaking...")
  interactions <- findInteractions(data = data.use, idents = idents.use, ligand = object@interaction[,object@command$Screening[["partner_a"]]], 
                                   receptor = object@interaction[,object@command$Screening[["partner_b"]]], stats = stats, pair.fxn = pair.fxn, n_perm = n_perm, emd = emd, weighted = weighted,
                                   nbins = nbins, p.adjust.method = p.adjust.method, min_pct = min_emd_pct, max_pct = max_emd_pct, drop_zero = drop_zero, cutoff.stats = cutoff.stats, 
                                   return.perm = save.perm, n_cores = n_cores, seed = seed)
  object@strength <- asSparse(mat = interactions$strength, zero_percent = object@command$Screening[["zero_percent"]])
  object@pvalue <- asSparse(mat = interactions$pvalues, zero_percent = object@command$Screening[["zero_percent"]])
  object@misc <- c(object@misc, interactions[-c(1:2)])
  object@command$Matchmaking <- logCommand()
  return(object)
}

#' Filter and Select Interaction Data
#'
#' Filter interaction using interacting strength and p value cutoffs. 
#' Non-significant interactions will be either removed or set to zero.
#' 
#' @param object Matchmaker object.
#' @param pval.cutoff p value cutoff. Default is 0.05.
#' @param strength.pct Top strength quantile to select. Default is 0.1 aka the top 10\%.
#' @param filter.cells Whether to filter the cells. Default is TRUE.
#' @return Returns a Matchmaker object with filtered interactions saved in @selected slot.
#' @export
#' @examples \dontrun{
#' object <- Selecting(object, pval.cutoff = 0.05, strength.pct = 0.1)
#' }
Selecting <- function(object, pval.cutoff = 0.05, strength.pct = 0.1, filter.cells = TRUE){
  filtered.data <- filterInteractions(interact.strength = object@strength, interact.pval = object@pvalue, interaction.meta = object@interaction, 
                                      pval.cutoff = pval.cutoff, strength.pct = strength.pct, filter.cells = filter.cells)
  object@selected$strength <- asSparse(mat = filtered.data$strength, zero_percent = object@command$Screening[["zero_percent"]])
  object@selected$pvalue <- asSparse(mat = filtered.data$pvalues, zero_percent = object@command$Screening[["zero_percent"]])
  object@selected$merged.data <- filtered.data$merge.data
  object@command$Selecting <- logCommand()
  return(object)
}

#' Convert Interaction Data
#'
#' Convert interaction data frame into a molten data frame using \code{\link[reshape2]{melt}}. If there are many interactions the molten data can be very long. 
#' If \code{\link[scMatchmaker]{Selecting}} was run, it will convert the data in @selected slot Interactions with zero strenghth will be removed.
#' 
#' @param object Matchmaker object.
#' @param selected Use selected data if calculated. Default is TRUE.
#' @return Returns a ranked list of interactions.
#' @export
#' @examples \dontrun{
#' convert.data <- Converting(object)
#' }
Converting <- function(object, selected = TRUE){
  if(selected){
    if(length(x = object@selected) > 0){
      message("Using @selected data")
      convert.data <- convertData(interact.strength = object@selected$strength, interact.pval = object@selected$pvalue)
    } else {
      stop("@selected not found, please run Selecting(object, ...) first!", call. = FALSE)
    }
  } else {
    convert.data <- convertData(interact.strength = object@strength, interact.pval = object@pvalue)
  }
  return(convert.data)
}

#' Downsampling and Geometric Sketching
#'
#' Run downsampling or geometric sketching to reduce number of cells for interaction calculation.
#' See \href{https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30152-8}{Hie et al. 2019 Cell Systems} for details.
#' 
#' @param object Matchmaker object.
#' @param cell_id Cell IDs used to subset the data. Default is NULL. If provided the downsampling will be skipped.
#' @param downsampling Perform downsampling based on the idents probability rather than geometric sketching. Default is TRUE.
#' @param size Number of cells to return. Default is 2000. Please note that if downsampling is performed the actual size may vary due to rounding effect.
#' @param geom_pca_dims Number of PCA dimensions to use. Default is 10.
#' @param scale Whether to scale the data. Default is TRUE.
#' @param seed Random seed number for sampling. Default is 1.
#' @param which_python Path to python3 used.
#' @return Returns a Matchmaker object.
#' @export
#' @examples \dontrun{
#' object <- Sketching(object, size = 2000)
#' }
Sketching <- function(object, cell_id = NULL, downsampling = TRUE, size = 2000, geom_pca_dims = 10, scale = TRUE, seed = 1, which_python = Sys.which(names = "python3")){
  if(is.null(x = cell_id)){
    if(downsampling){
      idents.counts <- round(x = (table(object@annotation[,1])/length(x = object@annotation[,1]))*size)
      idents.split <- split(x = 1:ncol(x = object@data), f = object@annotation[,1])
      set.seed(seed = seed)
      sketch.id <- unlist(x = lapply(X = 1:length(x = idents.split), FUN = function(x) sample(x = idents.split[[x]], size = idents.counts[x])))
    } else {
      sketch.id <- geomSketch(data = object@data, geom_size = size, n_pca = geom_pca_dims, scale = scale, seed = seed, which_python = which_python)
    }
  } else {
    sketch.id <- cell_id
  }
  object@data <- object@data[, sketch.id]
  object@annotation <- object@annotation[sketch.id,,drop=FALSE]
  object@misc$sketch_id <- sketch.id
  object@command$Sketching <- logCommand() 
  return(object)
}

#' Subset Interactions
#'
#' Subset the interactions based on cell type identity and interacting partners.
#'  
#' @param object Matchmaker object.
#' @param ident1 Cell type identity group 1 to subset.
#' @param ident2 Cell type identity group 2 to subset. Default is NULL, and all cell types will be used.
#' @param partner_a Interaction partner a (See \code{\link[scMatchmaker]{Screening}}) to subset. Defaulte is NULL, and all interactions will be selected.
#' @param partner_b Interaction partner b (See \code{\link[scMatchmaker]{Screening}}) to subset. Defaulte is NULL, and all interactions will be selected.
#' @param directed Logic, if TRUE subset only ident1 with ligands (partner_a) and ident2 with receptor (partner_b). 
#'                        If FALSE, ident2 with ligands and ident1 with receptor will be included. Default is TRUE.
#' @return Returns a new Matchmaker object.
#' @export
#' @examples \dontrun{
#' # subset all Macrophage related interactions
#' interaction.mac <- Subsetting(object, ident1 = "Macrophage")
#' # subset all CSF1 related interactions
#' interaction.csf <- Subsetting(object, partner_a = "CSF1")
#' }
Subsetting <- function(object, ident1 = NULL, ident2 = NULL, partner_a = NULL, partner_b = NULL, directed = FALSE){
  idents.all <- levels(x = object@annotation[,1])
  if(is.null(x = ident1) & is.null(x = ident2) & is.null(x = partner_a) & is.null(x = partner_b)) 
    stop("At least one of the arguments: ident1, ident2, partner_a or partner_b has to be provided.", call. = FALSE)
  if(is.null(x = ident1))  ident1 <- levels(x = object@annotation[,1])
  if(is.null(x = ident2)) ident2 <- levels(x = object@annotation[,1])
  if(!all(ident1 %in% idents.all)) stop("ident1 not found in @annotation", call. = FALSE)
  if(!all(ident2 %in% idents.all)) stop("ident2 not found in @annotation", call. = FALSE)
  idx.use <- which(x = object@annotation[,1] %in% c(ident1, ident2))
  object@annotation <- object@annotation[idx.use,,drop = FALSE]
  object@data <- object@data[,idx.use]
  object@misc$raw.data <- object@misc$raw.data[,idx.use]
  idents.use <- apply(X = expand.grid(ident1, ident2), MARGIN = 1, FUN = paste, collapse = "|")
  if(!directed) idents.use <- c(idents.use, apply(X = expand.grid(ident2, ident1), MARGIN = 1, FUN = paste, collapse = "|"))
  idx.use1 <- which(x = rownames(x = object@strength) %in% idents.use)
  object@strength <- object@strength[idx.use1, ,drop = FALSE]
  object@pvalue <- object@pvalue[idx.use1, ,drop = FALSE]
  if(!is.null(x = object@selected$strength)){
    idx.use2 <- which(x = rownames(x = object@selected$strength) %in% idents.use)
    object@selected$strength <- object@selected$strength[idx.use1, ,drop = FALSE]
    object@selected$pvalue <- object@selected$pvalue[idx.use1, ,drop = FALSE]
  }
  if(!is.null(x = partner_a) | !is.null(x = partner_b)){
    if(is.null(x = partner_a)) partner_a <- unique(x = object@interaction[,1])
    if(is.null(x = partner_b)) partner_b <- unique(x = object@interaction[,2])
    if(!all(partner_a %in% object@interaction[,1])) stop("partner_a not found in @interaction.", call. = FALSE)
    if(!all(partner_b %in% object@interaction[,2])) stop("partner_b not found in @interaction.", call. = FALSE)
    genes.use <- apply(X = expand.grid(partner_a, partner_b), MARGIN = 1, FUN = paste, collapse = "_")
    genes.use1 <- which(x = colnames(x = object@strength) %in% genes.use)
    object@interaction <- object@interaction[object@interaction[,1] %in% partner_a & object@interaction[,2] %in% partner_b,]
    object@strength <- object@strength[, genes.use1, drop = FALSE]
    object@pvalue <- object@pvalue[, genes.use1, drop = FALSE]
    if(!is.null(x = object@selected$strength)){
      genes.use2 <- which(x = colnames(x = object@selected$strength) %in% genes.use)
      object@selected$strength <- object@selected$strength[, genes.use2, drop = FALSE]
      object@selected$pvalue <- object@selected$pvalue[, genes.use2, drop = FALSE]
    }
  }
  object@command$Subsetting <- logCommand() 
  return(object)
}

#' Save Interaction Data
#'
#' Save the interactions calculated.
#'  
#' @param object Matchmaker object.
#' @param file_name File name prefix to save as.
#' @param dir Directory. Default is current directory (NULL).
#' @param selected Use selected data if calculated. Default is FALSE.
#' @return Return three csv files:
#' \itemize{
#'   \item filename_strength.csv, interaction strength data.
#'   \item filename_pval.csv, interaction p-value data.
#'   \item filename_rank.csv, ranked list of interactions.
#'   }
#' @importFrom utils write.csv
#' @export
#' @examples \dontrun{
#' Saving(object, file_name = "decidua")
#' }
Saving <- function(object, file_name, dir = NULL, selected = FALSE){
  if(selected){
    if(length(x = object@selected) == 0) stop("@selected not found, please run Selecting(object, ...) first!", call. = FALSE)
    strength <- object@selected$strength
    pval <- object@selected$pvalue
    rank.list <- Converting(object = object, selected = selected)
    file_name <- paste(file_name, "selected", sep = "_")
  } else {
    strength <- object@strength
    pval <- object@pvalue
    rank.list <- Converting(object = object, selected = selected)
  }
  gene.use <- paste(object@interaction[,1], object@interaction[,2], sep = "_")
  meta.data <- object@interaction[match(x = colnames(x = strength), table = gene.use),]
  merge.strength <- cbind.data.frame(meta.data, t(x = as.matrix(x = strength)), row.names = 1:nrow(x = meta.data))
  merge.pvalue <- cbind.data.frame(meta.data, t(x = as.matrix(x = pval)), row.names = 1:nrow(x = meta.data))
  if(!is.null(x = dir)){
    write.csv(x = merge.strength, file = paste(dir, "/", paste(file_name, "strength.csv", sep = "_"),sep = ""))
    write.csv(x = merge.pvalue, file = paste(dir, "/", paste(file_name, "pval.csv", sep = "_"),sep = ""))
    write.csv(x = rank.list, file = paste(dir, "/", paste(file_name, "rank.csv", sep = "_"),sep = ""))
  } else {
    write.csv(x = merge.strength, file = paste(file_name, "strength.csv", sep = "_"))
    write.csv(x = merge.pvalue, file = paste(file_name, "pval.csv", sep = "_"))
    write.csv(x = rank.list, file = paste(file_name, "rank.csv", sep = "_"))
  }
}

#' Calculate Complex Interaction 
#'
#' Calculate complex interactions by combining its subunits. 
#' Each complex is named by its subunits saperated by semicolon, such as ITGB1:ITGA2.
#' Therefore COL9A2_ITGB1:ITGA2 represents interactions between COL9A2 and ITGB1:ITGA2 complex.
#'  
#' @param object Matchmaker object.
#' @param subunit_a Column names or indices (of @interactions slot) containing complex subunits information for interaction partner_a. Default is subunit_a_1, subunit_a_2, subunit_a_3.
#' @param subunit_b Column names or indices (of @interactions slot) containing complex subunits information for interaction partner_b. Default is subunit_b_1, subunit_b_2, subunit_b_3.
#' @param strength_comb Method to combine complex interaction strengths from its subunits. Default is min (minimum).
#' \itemize{
#'   \item min, minimum strength of its subunits.
#'   \item average, average strength of its subunits.
#'   \item max, maximum strength of its subunits.
#'   }
#' @param pval_comb Method to combine complex interaction p values from its subunits. Default is max (maximum).
#' \itemize{
#'   \item max, maximum p value of its subunits.
#'   \item average, average p value of its subunits.
#'   \item min, minimum p value of its subunits.
#'   }
#' @return Returns a Matchmaker object.
#' @export
#' @examples \dontrun{
#' # For example, if CellPhoneDB data is used via LoadCellPhoneDB function, 
#' # the subunits are named as subunit_a_1, subunit_a_2 ...
#' object <- Complexing(object, strength_comb = "min", pval_comb = "max")
#' }
Complexing <- function(object, subunit_a = c("subunit_a_1","subunit_a_2","subunit_a_3"), subunit_b = c("subunit_b_1","subunit_b_2","subunit_b_3"), 
                       strength_comb = c("min", "average", "max"), pval_comb = c("max", "average", "min")){
  strength.comb <- match.arg(arg = strength_comb)
  strength.fxn <- switch(strength.comb,
                         min = min,
                         average = mean,
                         max = max)
  pval.comb <- match.arg(arg = pval_comb)
  pval.fxn <- switch(pval.comb, 
                     min = min,
                     average = mean,
                     max = max)
  idx.a <- apply(X = object@interaction[, subunit_a], MARGIN = 2, FUN = function(x) x %in% object@interaction[,object@command$Screening$partner_a])
  idx.b <- apply(X = object@interaction[, subunit_b], MARGIN = 2, FUN = function(x) x %in% object@interaction[,object@command$Screening$partner_b])
  idx.a <- apply(X = idx.a, MARGIN = 1, FUN = any)
  idx.b <- apply(X = idx.b, MARGIN = 1, FUN = any)
  complex.table <- object@interaction[idx.a|idx.b,]
  single.table <- object@interaction[!(idx.a|idx.b),]
  if(nrow(x = complex.table) == 0) stop("No complex detected.", call. = FALSE)
  complex.table$interaction_type <- "complex"
  if(nrow(x = single.table) > 0){
    single.table$interaction_type <- "single"
    single.names <- paste(single.table[,object@command$Screening$partner_a], single.table[,object@command$Screening$partner_b], sep = "_")
    single.strength <- object@strength[,single.names]
    single.pval <- object@pvalue[,single.names]
  } else {
    single.strength <- NULL
    single.pval <- NULL
  }
  complex.names <- apply(X = complex.table, MARGIN = 1, FUN = function(x){
    if(any(!is.na(x = x[subunit_a]) & x[subunit_a] != "")){
      genes.a <- x[subunit_a][!is.na(x = x[subunit_a]) & x[subunit_a] != ""]
      names.a <- paste(genes.a, collapse = ":")
      if(any(!is.na(x = x[subunit_b]) & x[subunit_b] != "")){
        genes.b <- x[subunit_b][!is.na(x = x[subunit_b]) & x[subunit_b] != ""]
        names.b <- paste(genes.b, collapse = ":")
        names.new <- paste(names.a, names.b, sep = "_")
        interaction.names1 <- expand.grid(genes.a, genes.b)
        interaction.names2 <- expand.grid(genes.b, genes.a)
        interaction.names1 <- paste(interaction.names1[,1], interaction.names1[,2], sep = "_")
        interaction.names2 <- paste(interaction.names2[,1], interaction.names2[,2], sep = "_")
        interaction.names <- intersect(x = c(interaction.names1, interaction.names2), y = colnames(object@strength))
        strength.new <- apply(X = object@strength[,interaction.names,drop = FALSE], MARGIN = 1, FUN = strength.fxn)
        pvalues.new <- apply(X = object@pvalue[,interaction.names,drop = FALSE], MARGIN = 1, FUN = pval.fxn)
        data.new <- cbind(strength.new, pvalues.new)
        colnames(x = data.new) <- c(names.new, names.new)
        return(list(data.new))
      } else {
        genes.b <- x[object@command$Screening$partner_b]
        names.new <- paste(names.a, genes.b, sep = "_")
        interaction.names1 <- expand.grid(genes.a, genes.b)
        interaction.names2 <- expand.grid(genes.b, genes.a)
        interaction.names1 <- paste(interaction.names1[,1], interaction.names1[,2], sep = "_")
        interaction.names2 <- paste(interaction.names2[,1], interaction.names2[,2], sep = "_")
        interaction.names <- intersect(x = c(interaction.names1, interaction.names2), y = colnames(x = object@strength))
        strength.new <- apply(X = object@strength[,interaction.names,drop = FALSE], MARGIN = 1, FUN = strength.fxn)
        pvalues.new <- apply(X = object@pvalue[,interaction.names,drop = FALSE], MARGIN = 1, FUN = pval.fxn)
        data.new <- cbind(strength.new, pvalues.new)
        colnames(x = data.new) <- c(names.new, names.new)
        return(list(data.new))
      }
    } else {
      genes.b <- x[subunit_b][!is.na(x = x[subunit_b]) & x[subunit_b] != ""]
      names.b <- paste(genes.b, collapse = ":")
      genes.a <- x[object@command$Screening$partner_a]
      names.new <- paste(genes.a, names.b, sep = "_")
      interaction.names1 <- expand.grid(genes.a, genes.b)
      interaction.names2 <- expand.grid(genes.b, genes.a)
      interaction.names1 <- paste(interaction.names1[,1], interaction.names1[,2], sep = "_")
      interaction.names2 <- paste(interaction.names2[,1], interaction.names2[,2], sep = "_")
      interaction.names <- intersect(x = c(interaction.names1, interaction.names2), y = colnames(x = object@strength))
      strength.new <- apply(X = object@strength[,interaction.names,drop = FALSE], MARGIN = 1, FUN = strength.fxn)
      pvalues.new <- apply(X = object@pvalue[,interaction.names,drop = FALSE], MARGIN = 1, FUN = pval.fxn)
      data.new <- cbind(strength.new, pvalues.new)
      colnames(x = data.new) <- c(names.new, names.new)
      return(list(data.new))
    }
  })
  complex.names <- lapply(X = complex.names, FUN = "[[", ... = 1L)
  complex.strength <- lapply(X = complex.names, FUN = function(x) x[,1,drop=FALSE])
  complex.pval <- lapply(X = complex.names, FUN = function(x) x[,2,drop=FALSE])
  complex.strength <- do.call(cbind, complex.strength)
  complex.pval <- do.call(cbind, complex.pval)
  comb.table <- rbind(single.table, complex.table)
  comb.table$interaction_name <- c(colnames(x = single.strength), colnames(x = complex.strength))
  complex.strength <- complex.strength[, unique(x = colnames(x = complex.strength))]
  complex.pval <- complex.pval[, unique(x = colnames(x = complex.pval))]
  object@misc$single_interaction <- object@interaction
  object@misc$single_strength <- object@strength 
  object@misc$single_pvalue <- object@pvalue
  object@interaction <- comb.table
  object@strength <- cbind(single.strength, complex.strength)
  object@pvalue <- cbind(single.pval, complex.pval)
  object@command$Complexing <- logCommand() 
  message("Complexing finished, re-run Selecting(object, ...) if neccesary.")
  return(object)
}

#' Reset Interactions
#'
#' Reset complex-complex interactions back to single gene-gene interactions.
#' The single gene-gene interactions are stored in the @misc slot with names: single_interaction, single_strength, single_pvalue.
#' 
#' @param object Matchmaker object.
#' @param by Method to reset interaction data. Default is complex.
#' \itemize{
#'   \item complex, reset to single gene-gene interactions before Complexing is called.
#'   \item merge, reset to unmerged (directed) interactions before Merging is called.
#'   }
#' @return Returns a Matchmaker object.
#' @export
#' @examples \dontrun{
#' # reset to single gene-tene interaction
#' object <- Resetting(object, by = "complex")
#' # reset to unmerged interaction
#' object <- Resetting(object, by = "merge")
#' }
Resetting <- function(object, by = c("complex", "merge")){
  reset.by <- match.arg(arg = by)
  if(reset.by == "merge"){
    if(is.null(x = object@misc$unmerged_strength) | is.null(x = object@misc$unmerged_pvalue)) 
      stop("Unmerged interactions not found, please make sure if Merging(object, ...) has been run.", call. = FALSE)
    object@strength <- object@misc$unmerged_strength
    object@pvalue <- object@misc$unmerged_pvalue
    message("Resetting merged interactions finished, @strength, @pvalue have been resetted.")
  }
  if(reset.by == "complex"){
    if(is.null(x = object@misc$single_strength) | is.null(x = object@misc$single_strength) | is.null(x = object@misc$single_pvalue)) 
      stop("Single gene-gene interactions not found, please make sure if Complexing(object, ...) has been run.", call. = FALSE)
    object@interaction <- object@misc$single_interaction
    object@strength <- object@misc$single_strength
    object@pvalue <- object@misc$single_pvalue
    message("Resetting comlex interactions finished, @interaction, @strength, @pvalue have been resetted.")
  }
  object@command$Resetting <- logCommand() 
  message("Re-run Selecting(object, ...) if neccesary.")
  return(object)
}

#' Merge Interactions
#'
#' Merge directed (one-way) interactions into undirected (two-way) interactions.
#' For example, before merging, DC1|dM1 and dM1|DC1 represents two different cell-cell interactions. 
#' After merging, they will be combined as one undirected cell-cell interaction between DC1 and dM1.
#' 
#' @param object Matchmaker object.
#' @param strength_merge Method to combine complex interaction strength from its subunits. Default is max (maximum).
#' \itemize{
#'   \item max, maximum strength of directed interactions.
#'   \item average, average strength of directed interactions.
#'   \item min, minimum strength of directed interaction.
#'   }
#' @param pval_merge Method to combine complex interaction p value from its subunits. Default is min (minimum).
#' \itemize{
#'   \item min, minimum p value of directed interactions.
#'   \item average, average p value of directed interactions.
#'   \item max, maximum p value of directed interactions.
#'   }
#' @return Returns a Matchmaker object.
#' @importFrom matrixStats colMaxs colMins colMeans2
#' @importFrom utils combn
#' @export
#' @examples \dontrun{
#' object <- Merging(object, strength_merge = "average", pval_merge = "max")
#' }
Merging <- function(object, strength_merge = c("max","average", "min"), pval_merge = c("min", "average", "max")){
  strength.merg <- match.arg(arg = strength_merge)
  strength.fxn <- switch(strength.merg,
                         average = colMeans2,
                         min = colMins,
                         max = colMaxs)
  pval.merg <- match.arg(arg = pval_merge)
  pval.fxn <- switch(pval.merg, 
                     max = colMaxs,
                     average = colMeans2,
                     min = colMins)
  all.comb <- combn(x = unique(x = object@annotation[,1]), m = 2)
  self <- paste(unique(x = object@annotation[,1]), unique(x = object@annotation[,1]), sep = "|")
  merged.strength <- c()
  merged.pval <- c()
  names.use <- c()
  for(i in 1:ncol(x = all.comb)){
    idx <- c(paste(all.comb[,i], collapse = "|"), paste(rev(x = all.comb[,i]), collapse = "|"))
    merged.strength <- rbind(merged.strength, strength.fxn(x = object@strength[idx, ]))
    merged.pval <- rbind(merged.pval, pval.fxn(x = object@pvalue[idx, ]))
    names.use <- c(names.use, idx[1])
  }
  rownames(x = merged.strength) <- names.use
  rownames(x = merged.pval) <- names.use
  merged.strength <- rbind(merged.strength, object@strength[self,])
  merged.pval <- rbind(merged.pval, object@pvalue[self,])
  object@misc$unmerged_strength <- object@strength
  object@misc$unmerged_pvalue <- object@pvalue
  object@strength <- merged.strength
  object@pvalue <- merged.pval
  object@command$Merging <- logCommand() 
  message("Merging finished, re-run Selecting(object, ...) if neccesary.")
  return(object)
}

#' Differential Interactions
#'
#' Conduct differential interaction analysis based on permutation similar to \code{\link[scMatchmaker]{Matchmaking}}.
#' To use this feature, save.perm parameter in \code{\link[scMatchmaker]{Matchmaking}} has to be set to TRUE when running Matchmaking.
#'  
#' @param object Matchmaker object.
#' @param interaction1 Cell-Cell interactions of interests.
#' @param interaction2 Cell-Cell interactions to compare against. Default is NULL, all other cell types except the ones in interaction1 will be used.
#' @param diff.thresh Differential interaction strength cutoff, default is 1.
#' @param p.val.cutoff p value cutoff. Default is 0.05
#' @param p.adjust.method p value adjustment method. Default is "none". See details \code{\link[stats]{p.adjust}}.
#' @return A data.frame with interaction differences and p values.
#' @importFrom magrittr %>% divide_by
#' @importFrom Matrix colMeans
#' @export
#' @examples \dontrun{
#' object <- Differing(object, interaction1 = "Macrophage|CAF", interaction2 = "Macrophage|Tumor")
#' }
Differing <- function(object, interaction1, interaction2 = NULL, diff.thresh = 1, p.val.cutoff = 0.05, p.adjust.method = "none"){
  if(!all(interaction1 %in% rownames(x = object@strength))) stop(setdiff(x = interaction1, y = rownames(x = object@strength)), " in interaction1 not found.", call. = FALSE)
  if(is.null(x = interaction2)) interaction2 <- setdiff(x = rownames(x = object@strength), y = interaction1)
  if(!all(interaction2 %in% rownames(x = object@strength))) stop(setdiff(x = interaction1, y = rownames(x = object@strength)), " in interaction2 not found.", call. = FALSE)
  if(is.null(x = object@misc$permute_result)) stop("Permutation results not saved, try setting save.perm=TRUE when running Matchmaking.", call. = FALSE)
  if(is.null(x = object@command$Complexing)){
    strength.null <- object@strength
  } else {
    strength.null <- object@misc$single_strength
  }
  interact.diff.null <- colMeans(x = strength.null[interaction1, ,drop=FALSE]) - colMeans(x = strength.null[interaction2, ,drop=FALSE])
  interact.diff.abs <- abs(x = interact.diff.null)
  interact.diff.pval <-
    lapply(X = object@misc$permute_result, FUN = function(x){
      interact.diff <- abs(x = colMeans(x = x[interaction1, ,drop=FALSE]) - colMeans(x = x[interaction2, ,drop=FALSE]))
      interact.diff >= interact.diff.abs
    }) %>%
    Reduce(f = "+") + 1 %>%
    divide_by(length(x = object@misc$permute_result) + 1) %>%
    p.adjust(method = p.adjust.method)
  diff.data <- cbind.data.frame(Difference = interact.diff.null, p_value = interact.diff.pval, 
                                row.names = make.names(names(x = interact.diff.null), unique = TRUE))
  diff.data <- subset(x = diff.data, abs(Difference) >= diff.thresh & p_value <= p.val.cutoff)
  diff.data <- diff.data[order(diff.data[["Difference"]], decreasing = TRUE),]
  return(diff.data)
}

