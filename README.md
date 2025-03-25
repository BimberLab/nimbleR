![R Build and Checks](https://github.com/BimberLab/nimbleR/workflows/R%20Build%20and%20Checks/badge.svg)

# nimbleR
The nimbleR package is a companion to the [nimble aligner](https://github.com/BimberLab/nimble/). The R package is primarily designed to process and append the [per-cell](https://github.com/BimberLab/nimble/wiki/Quickstart#3-generate-a-feature-by-cell-matrix) alignment data produced by nimble (in TSV format), and append it to a seurat object either as a new assay, or merged with an existing assay. 

It also contains several specialized normalization functions, developed to work with nimble's supplemental alignment data and also MHC/HLA data.  

### <a name="installation">Installation</a>

```{r}
# Install requirements. Other dependencies will be downloaded automatically
install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies = TRUE, ask = FALSE)

# Updating your Rprofile (i.e. ~/.Rprofile), with the following line will ensure install.packages() pulls from Bioconductor repos:
local({options(repos = BiocManager::repositories())})

# Install latest version:
devtools::install_github(repo = 'bimberlab/nimbleR', dependencies = TRUE)
```

### <a name="usage">Basic Usage</a>

Append the output from nimble to a Seurat object:
```{r}
library(nimbleR)

# Path to the TSV file produced by 'nimble report'
nimbleFile <- "nimble_per_cell_data.tsv"

# Append the nimble data as a new assay called 'NimbleCounts'
seuratObj <- AppendNimbleCounts(
  seuratObj = seuratObj,
  nimbleFile = nimbleFile,
  targetAssayName = "NimbleCounts",
  normalizeData = FALSE
)

# Custom library size normalization of the NimbleCounts assay using the RNA assay as the library size reference
seuratObj <- LogNormalizeUsingAlternateAssay(
  seuratObj,
  assayToNormalize = "NimbleCounts",
  assayForLibrarySize = "RNA"
)

# Alternatively, append to an existing assay
seuratObj <- AppendNimbleCounts(
  seuratObj = seuratObj,
  nimbleFile = nimbleFile,
  targetAssayName = "RNA",
  renameConflictingFeatures = TRUE,  # Prevents collisions with existing features
  duplicateFeatureSuffix = ".Nimble",  # Suffix added to any renamed features
)
```
