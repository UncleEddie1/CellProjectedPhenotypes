#' Calculate Cell Phenotype Scores
#'
#' This function calculates phenotype scores for each cell based on the association of its neighborhoods with particular phenotypes. 
#' The scores are derived from the log fold change (logFC) values in the `neighborhoodPhenotypeAssociations` data structure and 
#' are adjusted by the distances between cells and neighborhood centers.
#'
#' @param milo_obj A `Milo` class object, containing the neighborhood connectivity matrix (`@nhoods`) and cell metadata (`@colData`).
#' @param neighborhoodPhenotypeAssociations A list where each element corresponds to a phenotype and contains a data frame with neighborhood-level associations. Produced by `calculateNeighborhoodCenterMatrix`.
#' @param cellNeighborhoodDistanceMatrix A matrix representing the Euclidean distances between cells and the neighborhoods in which they are included, output of `cellNeighborhoodDistanceMatrix`.
#' @param col_name The name of the column in `milo_obj@colData` that contains the cell identifiers.
#' @param varsToTest A character vector of phenotypes for which to calculate cell scores.
#'
#' @return A data frame where each column corresponds to a phenotype, and each row corresponds to a cell, containing the calculated phenotype scores for each cell.
#' @export
calculateCellScores <- function(milo_obj, neighborhoodPhenotypeAssociations, cellNeighborhoodDistanceMatrix, col_name, varsToTest) {
  cellNeighborhoodDistanceMatrix <- as(cellNeighborhoodDistanceMatrix, "dgCMatrix")
  cellscores <- lapply(neighborhoodPhenotypeAssociations,function(neighborhoodPhenotypeAssociations){    
    tmp <- as(Matrix::t(milo_obj@nhoods) * neighborhoodPhenotypeAssociations$logFC, "dgCMatrix")
    scores <- col_scoring_ignore0(tmp,cellNeighborhoodDistanceMatrix)
    colInd <- as.character(1:ncol(tmp))
    scores <- scores[colInd]
    return(scores)
  }) |> do.call(what = "cbind") |> as.data.frame()
  colnames(cellscores) <- varsToTest
  rownames(cellscores) <- colnames(milo_obj@assays@data$counts)
  cellscores[[col_name]] <- unlist(milo_obj@colData[[col_name]])
  return(cellscores)
}

#' Calculate Weighted Average of Sparse Matrices
#'
#' This internal function calculates the weighted average of sparse matrices. The weighting is done 
#' using the softmax function applied to the distances between cells and neighborhood centers.
#'
#' @param dgCMat A sparse matrix of phenotype association scores for the neighborhoods.
#' @param dgCMat2 A sparse matrix of distances to the neighborhoods.
#'
#' @return A numeric vector of weighted averages for each column of the input matrices.
#' @keywords internal
col_scoring_ignore0 <- function (dgCMat, dgCMat2) {
  nnz_per_col <- diff(dgCMat@p)
  ColInd <- rep.int(1:ncol(dgCMat), nnz_per_col)              
  matrixSplit1 <- split(dgCMat@x, ColInd)
  matrixSplit2<- split(dgCMat2@x, ColInd)        
  out <- sapply(1:length(matrixSplit1),function(i){            
    associations <- matrixSplit1[[i]]
    distances <- matrixSplit2[[i]]
    weights <- softmax(distances*-1,t = 0.3)           
    sum(associations * weights)            
  })
  names(out) <- names(matrixSplit1)
  return(out)         
}

#' Softmax function
#'
#' This function computes the softmax of a numeric vector. 
#'
#' @param x A numeric vector of values.
#' @param temperature A positive number controlling the sharpness of the output probabilities. 
#'
#' @return A numeric vector the same length as the input vector.
#' @keywords internal
softmax <- function(x, temperature = 1) {
  # Adjust the input by the temperature
  x_adjusted <- x / temperature
  # Subtract the max value for numerical stability
  x_adjusted <- x_adjusted - max(x_adjusted)
  # Compute the exponentials
  exp_x <- exp(x_adjusted)
  # Compute the softmax values
  softmax_values <- exp_x / sum(exp_x)
  return(softmax_values)
}  
