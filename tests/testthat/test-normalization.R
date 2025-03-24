context("Normalization")

test_that("LogNormalizeUsingAlternateAssay works as expected", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))

  assayToAdd <- Seurat::GetAssayData(seuratObj, assay = 'RNA', layer = 'counts')
  assayToAdd <- floor(assayToAdd[1:10,] / 5)

  rownames(assayToAdd) <- paste0('Feature', LETTERS[1:10])

  seuratObj[['Norm']] <- Seurat::CreateAssayObject(assayToAdd)

  seuratObj <- LogNormalizeUsingAlternateAssay(seuratObj, assayToNormalize = 'Norm', assayForLibrarySize = 'RNA')

  nd <- Seurat::GetAssayData(seuratObj, assay = 'Norm', layer = 'data')
  expect_equal(max(nd[,4]), 3.442982, tolerance = 0.000001)
  expect_equal(max(nd[,101]), 2.823479, tolerance = 0.000001)
})