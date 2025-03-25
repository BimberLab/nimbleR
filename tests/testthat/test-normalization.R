library(SeuratData)

context("Normalization")

getBaseSeuratData <- function(){
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(pbmc3k))
  
  return(seuratObj)
}

test_that("LogNormalizeUsingAlternateAssay works as expected", {
  seuratObj <- getBaseSeuratData()

  assayToAdd <- Seurat::GetAssayData(seuratObj, assay = 'RNA', layer = 'counts')
  assayToAdd <- floor(assayToAdd[c('H6PD', 'H6PD', 'JAK1', 'XAF1', 'TMEM52'),])

  rownames(assayToAdd) <- paste0('Feature', LETTERS[1:nrow(assayToAdd)])

  seuratObj[['Norm']] <- Seurat::CreateAssayObject(assayToAdd)

  seuratObj <- LogNormalizeUsingAlternateAssay(seuratObj, assayToNormalize = 'Norm', assayForLibrarySize = 'RNA', maxLibrarySizeRatio = 0.2)

  nd <- Seurat::GetAssayData(seuratObj, assay = 'Norm', layer = 'data')
  expect_equal(max(nd[,2]), 1.111578, tolerance = 0.000001)
  expect_equal(max(nd[,200]), 1.726143, tolerance = 0.000001)
})