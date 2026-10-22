#' Calculate Cell Neighborhood Distance Matrix
#'
#' This function calculates a sparse matrix representing the Euclidean distance between each cell and the neighborhoods in which it is included.
#'
#' @param NeighborhoodCenters A matrix where each row represent PC coordinates of neighborhood centers. Can be generated with `calculateNeighborhoodCenterMatrix`.
#' @param cellPCs A matrix where each row corresponds to a cell and each column corresponds to a principal component. Can be accessed from `seurat_obj@reductions$pca@cell.embeddings[,1:12]` following standard Seurat PCA. Only the first 12 PCs are used in this function. 
#' @param milo_obj A `Milo` class object containing the neighborhood connectivity matrix (`@nhoods`).
#' @param chunk_size An integer to specify chunk size in Euclidean distance calculations. Memory scales with this parameter.
#'
#' @return A sparse matrix where each element represents the Euclidean distance between a cell and a neighborhood in which it is included.
#'
#' @export
calculateCellNeighborhoodDistanceMatrix <- function(NeighborhoodCenters, cellPCs, milo_obj, chunk_size = 1000) {
  # Calculate the Euclidean distance matrix between cells and neighborhood centers.
  dist_matrix <- euclideanDistanceChunks(NeighborhoodCenters,cellPCs[rownames(milo_obj@nhoods),],chunk_size = 1000)
  dist_matrix <- t(dist_matrix)
  colnames(dist_matrix) <- rownames(NeighborhoodCenters)
  rownames(dist_matrix) <- rownames(milo_obj@nhoods)
  # multiplying the neighborhood connectivity matrix by the distance matrix results in a sparse matrix.
  dist_matrix_Sparse <- Matrix::t(as(dist_matrix * milo_obj@nhoods, "sparseMatrix"))
}

#' Calculate Euclidean distances between rows of two matrices.
#' 
#' euclidean_distances_chunks calculates the Euclidean distances between rows of two matrices.
#' Matrices are processed in smaller chunks for memory efficiency.
#' Note that the number of columns in matrix1 and matrix2 must be the same.
#' 
#' @param matrix1 A float matrix.
#' @param matrix2 A float matrix.
#' @param chunk_size An integer.
#' @importFrom Rfast rowsums
#' @importFrom Rfast mat.mult
#'
#' @return Matrix of pairwise Euclidean distances between the rows of 'matrix1' and 'matrix2'.
euclideanDistanceChunks <- function(matrix1, matrix2, chunk_size = 1000) {
  if (ncol(matrix1) != ncol(matrix2)) {
    stop("Error (in euclidean_distances_chunks): for valid matrix multiplication, matrix1 and matrix2 must have the same number of columns.")
  }
  numRows1 <- nrow(matrix1)
  numRows2 <- nrow(matrix2)
  dists <- matrix(0, numRows1, numRows2) # initialize elements to 0. numRows1 rows, numRows2 columns
  
  chunk_dist <- function(mat1_chunk, mat2_chunk) {
    sum1 <- Rfast::rowsums(mat1_chunk^2)
    sum2 <- Rfast::rowsums(mat2_chunk^2)
    dot_prod <- Rfast::mat.mult(mat1_chunk, t(mat2_chunk)) # dot product between rows of mat1_chunk and mat2_chunk
    sqrt(outer(sum1, sum2, "+") - 2 * dot_prod)
  }
  for (i in seq(1, numRows1, by = chunk_size)) {
    for (j in seq(1, numRows2, by = chunk_size)) {
      mat1_chunk <- matrix1[i:min(i + chunk_size - 1, numRows1), , drop = FALSE]
      mat2_chunk <- matrix2[j:min(j + chunk_size - 1, numRows2), , drop = FALSE]
      dists[i:min(i + chunk_size - 1, numRows1), j:min(j + chunk_size - 1, numRows2)] <- chunk_dist(mat1_chunk, mat2_chunk)
    }
  }
  return(dists)
}
