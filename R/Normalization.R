utils::globalVariables(
  names = c('lsr'),
  package = 'nimbleR',
  add = TRUE
)

#' @include Utils.R

#' @title LogNormalizeUsingAlternateAssay
#'
#' @param seuratObj The seurat object
#' @param assayToNormalize The name of the assay to normalize
#' @param assayForLibrarySize The name of the assay from which to derive library sizes. This will be added to the library size of assayToNormalize.
#' @param scale.factor A scale factor to be applied in normalization
#' @param maxLibrarySizeRatio This normalization relies on the assumption that the library size of the assay being normalized in negligible relative to the assayForLibrarySize. To verify this holds true, the method will error if librarySize(assayToNormalize)/librarySize(assayForLibrarySize) exceeds this value
#'
#' @import ggplot2 Matrix
#' @export
LogNormalizeUsingAlternateAssay <- function(seuratObj, assayToNormalize, assayForLibrarySize = 'RNA', scale.factor = 1e4, maxLibrarySizeRatio = 0.01) {
  toNormalize <- Seurat::GetAssayData(seuratObj, assay = assayToNormalize, layer = 'counts')
  toNormalize <- as(toNormalize, 'CsparseMatrix')   # ensure dgCMatrix
  libSizeCounts <- Seurat::GetAssayData(seuratObj, assay = assayForLibrarySize, layer = 'counts')

  if (!setequal(colnames(toNormalize), colnames(libSizeCounts))) {
    stop('The assayToNormalize and assayForLibrarySize do not have the same cell names!')
  }

  if (is.null(maxLibrarySizeRatio) || is.na(maxLibrarySizeRatio)) {
    maxLibrarySizeRatio <- Inf
  }

  start_time <- Sys.time()
  if (.doDebugLogging()) {
    logger::log_info('Calculating cell sizes for nimble nimble normalization')
  }

  sumX <- Matrix::colSums(toNormalize)
  otherSizes <- Matrix::colSums(libSizeCounts)[colnames(toNormalize)]
  librarySize <- sumX + otherSizes

  ratios <- ifelse(librarySize > 0, sumX / librarySize, 0)
  bad <- which(ratios > maxLibrarySizeRatio)
  if (length(bad) > 0) {
    i <- bad[1]
    stop(paste0('The ratio of library sizes was above maxLibrarySizeRatio for ', length(bad),
                ' cell(s), e.g.: ', colnames(toNormalize)[i], ', on assay: ', assayToNormalize,
                '. Ratio was: ', ratios[i], ' (', sumX[i], ' / ', librarySize[i], ')'))
  }

  # Scale each nonzero entry by its column's library size
  normalizedMat <- toNormalize
  colIdx <- rep.int(seq_len(ncol(normalizedMat)), diff(normalizedMat@p))
  normalizedMat@x <- log1p(normalizedMat@x / librarySize[colIdx] * scale.factor)

  end_time <- Sys.time()

  if (all(is.null(rownames(normalizedMat))) || any(rownames(toNormalize) != rownames(normalizedMat))) {
    stop('Row names not equal')
  }

  if (all(is.null(colnames(normalizedMat))) || any(colnames(toNormalize) != colnames(normalizedMat))) {
    stop('Column names not equal')
  }

  if (.doDebugLogging()) {
    logger::log_info('Normalization done')
  }

  print('Normalization time:')
  print(end_time - start_time)

  print(ggplot(data.frame(lsr = ratios), aes(x = lsr)) +
          geom_histogram() +
          ggtitle(paste0("Library size ratios: ", assayToNormalize)) +
          egg::theme_article()
  )

  return(Seurat::SetAssayData(seuratObj, assay = assayToNormalize, layer = 'data', new.data = normalizedMat))
}