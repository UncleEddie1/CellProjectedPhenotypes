# CellProjectedPhenotypes

**CellProjectedPhenotypes** (CPP) is a method to calculate phenotype scores for single cells based on their distance to phenotype-associated cell transcriptional neighborhoods. 

## Installation

### Prerequisites

Before installing the package, ensure that the required dependencies are installed:

```
install.packages("remotes")
remotes::install_github("immunogenomics/lisi")
BiocManager::install("sparseMatrixStats")
```

CPP makes extensive use of _Milo_, a method to calculate differential abundance on KNN graphs from single-cell datasets. More information about _Milo_ can be obtained from their [github page](https://github.com/MarioniLab/miloR) and in the associated [manuscript](https://doi.org/10.1038/s41587-021-01033-z).

CPP uses the Local Inverse Simpson's Index (LISI) to remove poorly mixed cells. Details of this method can be gained from [this manuscript](https://www.nature.com/articles/s41592-019-0619-0) and from the corresponding [github page](https://github.com/immunogenomics/LISI).

## Usage
Here we demonstrate basic CPP usage by calculating cell projected phenotype scores
from a Seurat-derived RDS file with corresponding metadata. Both of these files are available in our _data_ directory. 

1. Load essential libraries
```
library(Seurat)
library(SeuratObject)
library(miloR)
library(SingleCellExperiment)
```

2. Set input variables, including the name of the input file, corresponding metadata file, varsToTest (metadata variables to test for neighborhood associations), varsToCorrect(metadata covariates to control for in association testing), and col_name (the column of seurat object corresponding to the donor identifier).
```
file <- 'Vasculature_cells.rds'
meta <- read.csv("Metadata.csv")
varsToTest <- c("amyloid","tangles")
varsToCorrect <- c("fixation_interval","msex","pmi","study")
col_name <- "projid"
```

3. Load single-cell data and filter and subsample individuals with extreme cell counts.
```
raw_cells <- readRDS(file)
seurat_obj <- CreateSeuratObject(counts = raw_cells@assays$RNA@counts)
seurat_obj[[col_name]] <- raw_cells[[col_name]]
seurat_obj$cell_type <- raw_cells$cell_type_high_resolution
seurat_obj <- filterAndSubsampleByCellCount(seurat_obj)
```

4. Make a PCA embedding from the Seurat object and save first 12 principal components for later
```
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
cellPCs <- seurat_obj@reductions$pca@cell.embeddings[,1:12]
seurat_obj <- removePoorlyMixedCells(seurat_obj, cellPCs, col_name, 3)
```

5. Construct a KNN graph, neighborhoods, and neighborhood cell counts with _Milo_
```
sce <- as.SingleCellExperiment(seurat_obj)
milo_obj <- Milo(sce)
milo_obj <- buildGraph(milo_obj, k = 30, d =12) # 12 PCA dimensions, 30 nearest neighbors
milo_obj <- makeNhoods(milo_obj, prop = 0.10, k = 30, d=12, refined = TRUE) # prop is proportion of vertices to sample for neighborhoods.
milo_obj <- countCells(milo_obj, meta.data = data.frame(colData(milo_obj)), samples=col_name)
```

6. calculate neighborhood phenotype associations, the PCA coordinates of neighborhood centers, a matrix of distances between cells and neighborhoods, and a data frame of cell phenotype scores.
```
NeighborhoodCenters <- calculateNeighborhoodCenterMatrix(milo_obj, cellPCs)
cellNeighborhoodDistanceMatrix <- calculateCellNeighborhoodDistanceMatrix(NeighborhoodCenters, cellPCs, milo_obj)

# in-sample
neighborhoodPhenotypeAssociations <- runNeighborhoodAssociationTest(milo_obj, varsToTest, varsToCorrect, col_name)
cellScores <- calculateCellScores(milo_obj, neighborhoodPhenotypeAssociations, cellNeighborhoodDistanceMatrix, col_name, varsToTest)

# leave-one-out (LOO)
donors <- unique(milo_obj[[col_name]])
cellScores <- do.call(rbind, lapply(donors, function(d) {
  assoc_d  <- runNeighborhoodAssociationTest(milo_obj, varsToTest, varsToCorrect, col_name, exclude = d)
  scores_d <- calculateCellScores(milo_obj, assoc_d, cellNeighborhoodDistanceMatrix, col_name, varsToTest)
  scores_d[scores_d[[col_name]] == d, ]
}))
cellScores <- cellScores[colnames(milo_obj), ] #same cell order as the original cellScores

# fast LOO
donors <- unique(milo_obj[[col_name]])
fastLOO <- runFastLeaveOneOutAssociations(milo_obj, varsToTest, varsToCorrect, col_name)
cellScores <- do.call(rbind, lapply(donors, function(d) {
  assoc_d  <- lapply(fastLOO, function(B) data.frame(logFC = B[, as.character(d)]))
  scores_d <- calculateCellScores(milo_obj, assoc_d, cellNeighborhoodDistanceMatrix, col_name, varsToTest)
  scores_d[scores_d[[col_name]] == d, ]
}))
cellScores <- cellScores[colnames(milo_obj), ]
```

## Contributions
We welcome contributions and suggestions from the community. Please report ideas to our issues page.

## License
This package is licensed under the GPL-3 License. 

## Authors
* **Riley J. Mangan** Author and maintainer
* **Gerard Bouland** Author and contributer
* **Manolis Kellis** Author and oversight

## Citation
We ask that the following publication be cited in work using the CPP methodology or code base:
Mangan, Riley J. and Bouland, Gerard A. and Schubert, Tim and Anwer, Danish and Bai, Judy and Ekiz, Gülfem and Pigadas, Nicolas and Dattatri, Amogh and Polster, Annikka and Mahfouz, Ahmed and Reinders, Marcel and Kellis, Manolis, Cell-projected Phenotypes Link Transcriptional and Phenotypic Heterogeneity in Alzheimer’s Disease. Available at SSRN: https://ssrn.com/abstract=6615879 or http://dx.doi.org/10.2139/ssrn.6615879
