#' Test neighborhoods for associations with metadata test variables. 
#'
#' @param milo_obj A Milo object for which makeNhoods and countCells has already been run.
#' @param varsToTest A list of variable to test. Ex. c("amyloid","tangles")
#' @param varsToCorrect A list of covariates to account for. Ex. c("fixation_interval","msex","pmi","study")
#' @param col_name The name of the colnames in the nhoodCounts(milo_obj). Ex "projid".
#' @importFrom miloR testNhoods
#' @return A list of data frames, one per test variable. Each data frame will be a table, one row per neighborhood, of enrichment results.
#' @export
runNeighborhoodAssociationTest <- function(milo_obj, varsToTest, varsToCorrect, col_name, exclude = NULL) {
  # set up design matrix 
  designMatrix <- setNames(data.frame(colnames(nhoodCounts(milo_obj))), col_name)
  designMatrix <- cbind(designMatrix,meta[match(designMatrix[[col_name]],meta[[col_name]]),c(varsToTest,varsToCorrect)])
  rownames(designMatrix) <- designMatrix[[col_name]]
  designMatrix <- na.omit(designMatrix) #omit rows with NAs in metadata
  designMatrix <- droplevels(designMatrix[!rownames(designMatrix) %in% as.character(exclude), , drop = FALSE]) # exclude held-out donors
  dfToRun <- data.frame(testVar = varsToTest,
                        covariates = paste0("~", paste(varsToCorrect, collapse = " + "), " + %s"))
  phenoResults <- lapply(1:nrow(dfToRun),function(i){
    message(sprintf(dfToRun[i,"covariates"],dfToRun[i,"testVar"]))
    form <- as.formula(sprintf(dfToRun[i,"covariates"],dfToRun[i,"testVar"]))
    da_results <- testNhoods(milo_obj, design = form, design.df = designMatrix, fdr.weighting = "none")
    return(da_results)
  })
  names(phenoResults) <- varsToTest
  return(phenoResults)
}