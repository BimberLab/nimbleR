
#' @title LogNormalizeUsingAlternateAssay
#'
#' @param seuratObj The seurat object
#' @param assayToNormalize The name of the assay to normalize
#' @param assayForLibrarySize The name of the assay from which to derive library sizes. This will be added to the library size of assayToNormalize.
#' @param scale.factor A scale factor to be applied in normalization
#' @param maxLibrarySizeRatio This normalization relies on the assumption that the library size of the assay being normalized in negligible relative to the assayForLibrarySize. To verify this holds true, the method will error if librarySize(assayToNormalize)/librarySize(assayForLibrarySize) exceeds this value
#'
#' @import ggplot2
#' @export
LogNormalizeUsingAlternateAssay <- function(seuratObj, assayToNormalize, assayForLibrarySize = 'RNA', scale.factor = 1e4, maxLibrarySizeRatio = 0.01) {
  toNormalize <- Seurat::GetAssayData(seuratObj, assayToNormalize, layer = 'counts')
  assayForLibrarySizeObj <- Seurat::GetAssayData(seuratObj, assay = assayForLibrarySize, layer = 'counts')

  if (any(colnames(toNormalize) != colnames(assayForLibrarySize))) {
    stop(paste0('The assayToNormalize and assayForLibrarySize do not have the same cell names!'))
  }

  if (is.null(maxLibrarySizeRatio) || is.na(maxLibrarySizeRatio)) {
    maxLibrarySizeRatio <- Inf
  }

  margin <- 2
  ncells <- dim(x = toNormalize)[margin]

  start_time <- Sys.time()
  assayForLibrarySizeData <- Matrix::colSums(assayForLibrarySizeObj)
  ratios <- unlist(sapply(seq_len(length.out = ncells), function(i){
    x <- toNormalize[, i]
    sumX <- sum(x)
    librarySize <- sumX + assayForLibrarySizeData[i]

    lsr <- (sumX / librarySize)
    if (lsr > maxLibrarySizeRatio) {
      stop(paste0('The ratio of library sizes was above maxLibrarySizeRatio for cell: ', colnames(assayForLibrarySizeObj)[i], ', on assay: ', assayToNormalize, '. Ratio was: ', lsr, ' (', sumX, ' / ', librarySize, ')'))
    }

    xnorm <- log1p(x = x / librarySize * scale.factor)
    toNormalize[, i] <<- xnorm

    return(lsr)
  }))
  end_time <- Sys.time()

  print('Normalization time:')
  print(end_time - start_time)

  print(ggplot(data.frame(lsr = ratios), aes(x = lsr)) +
          geom_histogram() +
          ggtitle(paste0("Library size ratios: ", assayToNormalize)) +
          egg::theme_article()
  )

  seuratObj <- Seurat::SetAssayData(seuratObj, assay = assayToNormalize, layer = 'data', new.data = toNormalize)

  return(seuratObj)
}