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

runFastLeaveOneOutAssociations <- function(milo_obj, varsToTest, varsToCorrect, col_name) {
  designMatrix <- setNames(data.frame(colnames(nhoodCounts(milo_obj))), col_name)
  designMatrix <- cbind(designMatrix,meta[match(designMatrix[[col_name]],meta[[col_name]]),c(varsToTest,varsToCorrect)])
  rownames(designMatrix) <- designMatrix[[col_name]]
  designMatrix <- na.omit(designMatrix) #omit rows with NAs in metadata
  counts <- as.matrix(nhoodCounts(milo_obj)[, rownames(designMatrix)])
  
  # covariates-only model
  X0 <- model.matrix(reformulate(varsToCorrect), designMatrix)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts, lib.size = colSums(counts)), method = "TMM")
  dge <- edgeR::estimateDisp(dge, X0)
  aug <- edgeR::addPriorCount(counts, offset = edgeR::getOffset(dge), prior.count = 0.125)   # edgeR's logFC shrinkage
  mu0 <- edgeR::glmFit(aug$y, X0, offset = aug$offset, dispersion = dge$trended.dispersion, prior.count = 0)$fitted.values
  w <- mu0 / (1 + dge$trended.dispersion * mu0)
  S <- aug$y / mu0 - 1
  S[!is.finite(S)] <- 0 # neighborhoods with no expected cells
  
  # per-trait, donor held-out
  lapply(setNames(varsToTest, varsToTest), function(trait) {
    yt  <- lm.fit(X0, as.numeric(designMatrix[[trait]]))$residuals   # trait with covariates removed
    num <- as.vector((w * S) %*% yt)
    den <- as.vector(w %*% yt^2)
    loo <- (num - sweep(w * S, 2, yt, "*")) / (den - sweep(w, 2, yt^2, "*")) / log(2)
    missing <- setdiff(colnames(nhoodCounts(milo_obj)), colnames(loo)) # for donors missing metadata
    cbind(loo, matrix(num / den / log(2), nrow(loo), length(missing), dimnames = list(NULL, missing)))
  })
}