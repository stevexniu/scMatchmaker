library(Matrix)
setClassUnion(name = 'matrices', members = c("matrix", "dgCMatrix", "dsCMatrix"))

#' The Matchmaker S4 Class
#'
#' The Matchmaker object with the slot information listed as follow:
#'
#' @slot data Expression data. A d x M matrix with d rows of interaction genes and M columns of data points (cells).
#' @slot annotation Data frame of cell type annotatiuons.
#' \itemize{
#'   \item first column, stores the cell type identities for each cell in the data.
#'   \item ..., additional slots for future use.
#'   }
#' @slot interaction Ligand-receptor interaction pairs. A dataframe with at least two columns for ligand-receptor pair gene names.
#' @slot strength Matrix of relative strengths of interactions.
#' \itemize{
#'   \item row names, interacting cell types seperated by "|". 
#'   \item column names, interacting partner pairs seperated by "_".
#'   }
#' @slot pvalue Matrix of p values of interaction strengths. Row and column naming are the same as the strength matrix.
#' @slot selected A list of selected and filtered interactions:
#' \itemize{
#'   \item stregth, Matrix of relative strengths of interactions.
#'   \item pvalue, Matrix of p values of interaction strengths.
#'   \item merged.data, merged interaction pair data and strength data.
#'   }
#' @slot project_name Name of the project.
#' @slot misc List of miscellaneous:
#' \itemize{
#'   \item raw_data, raw expression data.
#'   \item stats_null, null interaction stength without earth mover's distance adjustment if emd=TRUE option is used.
#'   \item emd, earth mover's distance similarity.
#'   \item sketch_id, geometric sketching ID.
#'   \item graph, graph visualization.
#'   \item single_interaction, single gene-gene interaction table in the format of single interactions if Complexing function is called.
#'   \item single_strength, single gene-gene interaction strengths if Complexing function is called.
#'   \item single_pvalue, single gene-gene interaction p-values of interactions strengths if Complexing function is called.
#'   \item unmerged_strength, unmerged interaction strenghts if Merging function is called.
#'   \item unmerged_pvalue, unmerged interaction p-values if Merging function is called.
#'   \item ..., additional slots for future use.
#'   }
#' @slot command Nested list of commands used.
#' @name Matchmaker-class
#' @rdname Matchmaker-class
#' @exportClass Matchmaker
#'
setClass(Class = "Matchmaker",
         slots = c(
           data = "matrices",
           annotation = "data.frame",
           interaction = "data.frame",
           strength = "matrices",
           pvalue = "matrices",
           selected = "list",
           project_name = "character",
           misc = "list",
           command = "list"
         )
)

setMethod(f = "show",signature = "Matchmaker",
          definition = function(object){
            if(object@project_name == ""){
              cat("A Matchmaker object:", "\n")
            } else {
              cat(object@project_name, "Matchmaker object:", "\n")
            }
            cat(nrow(x = object@interaction), "interaction pairs across", ncol(x = object@data), "cells.", "\n")
            cat(ncol(x = object@strength), "ligand-recptor pairs across", nrow(x = object@strength), "cell-cell interactions calculated.", "\n")
          }
)
