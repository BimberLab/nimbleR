utils::globalVariables(
  names = c('V1','V2','V3'),
  package = 'nimbleR',
  add = TRUE
)

# @include Normalization.R

#' @import Seurat
#' @importFrom tidyr pivot_wider

#' @title AppendNimbleCounts
#' @description Reads a given seurat object and a nimble file, and appends the nimble data to the object.
#'
#' @param seuratObj A Seurat object.
#' @param nimbleFile A nimble file, which is a TSV of feature counts created by nimble
#' @param maxAmbiguityAllowed If provided, any features representing more than ths value will be discarded. For example, 'Feat1,Feat2,Feat3' represents 3 features. maxAmbiguityAllowed=1 results in removal of all ambiguous features.
#' @param targetAssayName The target assay. If this assay exists, features will be appended (and an error thrown if there are duplicates). Otherwise a new assay will be created.
#' @param renameConflictingFeatures If true, when appending to an existing assay, any conflicting feature names will be renamed, appending the value of duplicateFeatureSuffix
#' @param duplicateFeatureSuffix If renameConflictingFeatures is true, this string will be appended to duplicated feature names
#' @param normalizeData If true, data will be normalized after appending/creating the assay. This will default to CellMembrane::LogNormalizeUsingAlternateAssay; however, if assayForLibrarySize equals targetAssayName then Seurat::NormalizeData is used.
#' @param performDietSeurat If true, DietSeurat will be run, which removes existing reductions. This may or may not be required based on your usage, but the default is to perform this if the targetAssay exists.
#' @param assayForLibrarySize If normalizeData is true, then this is the assay used for librarySize when normalizing. If assayForLibrarySize equals targetAssayName, Seurat::NormalizeData is used.
#' @param maxLibrarySizeRatio If normalizeData is true, then this is passed to CellMembrane::LogNormalizeUsingAlternateAssay
#' @param doPlot If true, FeaturePlots will be generated for the appended features
#' @param maxFeaturesToPlot If doPlot is true, this is the maximum number of features to plot
#' @param replaceExistingAssayData If true, any existing data in the targetAssay will be deleted
#' @param featureRenameList An optional named list in the format <OLD_NAME> = <NEW_NAME>. If any <OLD_NAME> are present, the will be renamed to <NEW_NAME>. The intention of this is to recover specific ambiguous classes.
#' @return A modified Seurat object.
#'
#' @import dplyr
#' @export
AppendNimbleCounts <- function(seuratObj, nimbleFile, targetAssayName, maxAmbiguityAllowed = 0, renameConflictingFeatures = TRUE, duplicateFeatureSuffix = ".Nimble", normalizeData = TRUE, performDietSeurat = (targetAssayName %in% names(seuratObj@assays)), assayForLibrarySize = 'RNA', maxLibrarySizeRatio = 0.05, doPlot = FALSE, maxFeaturesToPlot = 40, replaceExistingAssayData = TRUE, featureRenameList = NULL) {
  if (!file.exists(nimbleFile)) {
    stop(paste0("Nimble file does not exist: ", nimbleFile))
  }
  
  # Read file and construct df
  df <- NULL
  tryCatch({
    df <- utils::read.table(nimbleFile, sep="\t", header=FALSE)
  }, error = function(e){
    if (conditionMessage(e) != 'no lines available in input') {
      stop(e)
    } else {
      print(paste0('No lines in nimble file: ', nimbleFile))
      df <- data.frame(V1 = character(), V2 = character(), V3 = character())
    }
  })

  # Indicates no data in TSV
  if (all(is.null(df))){
    return(seuratObj)
  }

  if (sum(df$V1 == "") > 0) {
    stop("The nimble data contains blank feature names. This should not occur.")
  }

  if (sum(grepl(df$V1, pattern = "^,")) > 0) {
    stop("The nimble data contains features with leading commas. This should not occur.")
  }

  if (sum(grepl(df$V1, pattern = ",$")) > 0) {
    stop("The nimble data contains features with trailing commas. This should not occur.")
  }

  d <- as.integer(df$V2)
  if (any(is.na(d))){
    stop(paste0('Non-integer count values found, were: ', paste0(utils::head(df$V2[is.na(d)]), collapse = ',')))
  }

  if (is.na(maxAmbiguityAllowed) || is.null(maxAmbiguityAllowed)){
    maxAmbiguityAllowed <- Inf
  } else if (maxAmbiguityAllowed == 0) {
    maxAmbiguityAllowed <- 1
  }

  # Ensure consistent sorting of ambiguous features, and re-group if needed:
  if (any(grepl(df$V1, pattern = ','))) {
    print('Ensuring consistent feature sort within ambiguous features:')
    df$V1 <- unlist(sapply(df$V1, function(y){
      return(paste0(sort(unlist(strsplit(y, split = ','))), collapse = ','))
    }))

    df <- df %>%
      group_by(V1, V3) %>%
      summarize(V2 = sum(V2))

    df <- df[c('V1', 'V2', 'V3')]

    paste0('Distinct features after re-grouping: ', length(unique(df$V1)))
  }

  if (!all(is.null(featureRenameList))) {
    print('Potentially renaming features:')
    df$V1 <- as.character(df$V1)
    totalRenamed  <- 0
    for (featName in names(featureRenameList)) {
      if (featName %in% df$V1) {
        df$V1[df$V1 == featName] <- featureRenameList[[featName]]
        totalRenamed <- totalRenamed + 1
      }
    }

    print(paste0('Total features renamed: ', totalRenamed))
  }

  #Remove ambiguous features
  totalHitsByRow <- sapply(df$V1, function(y){
    return(length(unlist(strsplit(y, split = ','))))
  })

  ambigFeatRows <- totalHitsByRow > maxAmbiguityAllowed
  if (sum(ambigFeatRows) > 0) {
    print(paste0('Dropping ', sum(ambigFeatRows), ' rows with ambiguous features (>', maxAmbiguityAllowed, '), ', sum(ambigFeatRows),' of ', nrow(df)))
    totalUMI <- sum(df$V2)
    x <- df$V1[ambigFeatRows]
    totalHitsByRow <- totalHitsByRow[ambigFeatRows]
    x[totalHitsByRow > 3] <- 'ManyHits'

    x <- sort(table(x), decreasing = T)
    x <- data.frame(Feature = names(x), Total = as.numeric(unname(x)))

    x$Fraction <- x$Total / totalUMI
    x <- x[x$Fraction > 0.005,]

    if (nrow(x) > 0) {
      x$Name <- substr(x$Feature, start = 1, stop = 40)
      x$Name[x$Name != x$Feature] <- paste0(x$Name[x$Name != x$Feature], '..')
      x$Total <- paste0(x$Total, ' (', scales::percent(x$Total / totalUMI, accuracy = 0.001), ')')

      print('Top ambiguous combinations:')
      print(utils::head(x[c('Name', 'Total')], n = 100))
    }

    df <- df[!ambigFeatRows, , drop = F]
    paste0('Distinct features after pruning: ', length(unique(df$V1)))
  }

  if (any(duplicated(df[c('V1','V3')]))) {
    print(paste0('Duplicate cell/features found. Rows at start: ', nrow(df)))
    df <- df %>%
      group_by(V1, V3) %>%
      summarize(V2 = sum(V2))

    df <- df[c('V1', 'V2', 'V3')]

    print(paste0('After re-grouping: ', nrow(df)))
  }

  tryCatch({
    # Group to ensure we have one value per combination:
    d <- as.integer(df$V2)
    if (any(is.na(d))){
        stop(paste0('Non-integer count values found, were: ', paste0(df$V2[is.na(d)], collapse = ',')))
    }
    rm(d)

    # Remove barcodes from nimble that aren't in seurat
    uniqueCells <- unique(df$V3)
    cellsToDrop <- uniqueCells[ ! uniqueCells %in% colnames(seuratObj) ]
    if (length(cellsToDrop) > 0) {
      nKept <- length(uniqueCells[ uniqueCells %in% colnames(seuratObj) ])
      print(paste0('Dropping ', length(cellsToDrop), ' cell barcodes not in the seurat object (out of ', length(uniqueCells), '). Keeping: ', nKept))
      df <- df %>% filter( ! V3 %in% cellsToDrop)
    }

    paste0('Distinct features: ', length(unique(df$V1)))
    paste0('Distinct cells: ', length(unique(df$V3)))

    df <- tidyr::pivot_wider(df, names_from=V3, values_from=V2, values_fill=0)
  }, error = function(e){
    utils::write.table(df, file = 'debug.nimble.txt.gz', sep = '\t', quote = F, row.names = F)

    print(paste0('Error pivoting input data for assay:', targetAssayName, ', results saved to: debug.nimble.txt.gz'))
    print(conditionMessage(e))
    traceback()
    e$message <- paste0('Error pivoting nimble data. target assay: ', targetAssayName)
    stop(e)
  })

  if (replaceExistingAssayData && targetAssayName %in% names(seuratObj@assays)) {
    print('Replacing existing assay')
    seuratObj@assays[[targetAssayName]] <- NULL
  }

  appendToExistingAssay <- targetAssayName %in% names(seuratObj@assays)

  # Fill zeroed barcodes that are in seurat but not in nimble
  zeroedBarcodes <- setdiff(colnames(seuratObj), colnames(df)[-1])
  print(paste0('Total cells lacking nimble data: ', length(zeroedBarcodes), ' of ', ncol(seuratObj), ' cells'))
  for (barcode in zeroedBarcodes) {
    df[barcode] <- 0
  }
  
  # Cast nimble df to matrix
  featureNames <- df$V1
  if (any(duplicated(featureNames))) {
    stop('Error, there were duplicate feature names')
  }

  df <- subset(df, select=-(V1))
  m <- Seurat::as.sparse(df)
  dimnames(m) <- list(featureNames, colnames(df))
  if (is.null(colnames(m))) {
    stop(paste0('Error: no column names in nimble count matrix, size: ', paste0(dim(m), collapse = ' by ')))
  }

  m <- m[,colnames(seuratObj), drop=FALSE] # Ensure column order matches
  if (appendToExistingAssay && ncol(m) != ncol(seuratObj@assays[[targetAssayName]])) {
    stop(paste0('Error parsing nimble data, ncol not equal after subset, was ', ncol(m)))
  }

  if (is.null(colnames(m))) {
    stop(paste0('Error: no column names in matrix after subset, size: ', paste0(dim(m), collapse = ' by ')))
  }

  if (appendToExistingAssay) {
    if (any(rownames(m) %in% rownames(seuratObj@assays[[targetAssayName]]))) {
      conflicting <- rownames(m)[rownames(m) %in% rownames(seuratObj@assays[[targetAssayName]])]

      if (renameConflictingFeatures) {
        print(paste0('The following nimble features have conflicts in the existing assay and will be renamed: ', paste0(conflicting, collapse = ',')))
        newNames <- rownames(m)
        names(newNames) <- newNames
        newNames[conflicting] <- paste0(conflicting, duplicateFeatureSuffix)
        newNames <- unname(newNames)
        rownames(m) <- newNames
      } else {
        stop(paste0('The following nimble features conflict with features in the seuratObj: ', paste0(conflicting, collapse = ',')))
      }
    }

    if (performDietSeurat) {
      print('Running DietSeurat')
      seuratObj <- Seurat::DietSeurat(seuratObj)
    }

    # Append nimble matrix to seurat count matrix
    existingBarcodes <- colnames(Seurat::GetAssayData(seuratObj, assay = targetAssayName, slot = 'counts'))
    if (sum(colnames(m) != existingBarcodes) > 0) {
      stop('cellbarcodes do not match on matrices')
    }

    # If feature source exists, retain it. Otherwise assume these are from cellranger:
    slotName <- .GetAssayMetaSlotName(seuratObj[[targetAssayName]])
    if ('FeatureSource' %in% names(methods::slot(seuratObj@assays[[targetAssayName]], slotName))) {
      fs <- methods::slot(seuratObj@assays[[targetAssayName]], slotName)$FeatureSource
    } else {
      fs <- rep('CellRanger', nrow(seuratObj@assays[[targetAssayName]]))
    }

    fs <- c(fs, rep('Nimble', nrow(m)))

    # perform in two steps to avoid warnings:
    ad <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(rbind(Seurat::GetAssayData(seuratObj, assay = targetAssayName, slot = 'counts'), m)))
    if (targetAssayName != Seurat::DefaultAssay(seuratObj)) {
      seuratObj[[targetAssayName]] <- NULL
    }
    seuratObj[[targetAssayName]] <- ad
    
    names(fs) <- rownames(seuratObj@assays[[targetAssayName]])
    seuratObj@assays[[targetAssayName]] <- Seurat::AddMetaData(seuratObj@assays[[targetAssayName]], metadata = fs, col.name = 'FeatureSource')

    if (sum(colnames(Seurat::GetAssayData(seuratObj, assay = targetAssayName, slot = 'counts')) != existingBarcodes) > 0) {
      stop('cellbarcodes do not match on matrices after assay replacement')
    }
  } else {
    # Add nimble as separate assay
    if (any(duplicated(rownames(m)))) {
      stop('Error: The count matrix had duplicate rownames')
    }

    if (nrow(m) == 0) {
      warning('No features present in incoming data, skipping append')
      return(seuratObj)
    }

    seuratObj[[targetAssayName]] <- Seurat::CreateAssayObject(counts = m, min.features = 0, min.cells = 0)

    fs <- rep('Nimble', nrow(seuratObj@assays[[targetAssayName]]))
    names(fs) <- rownames(seuratObj@assays[[targetAssayName]])
    seuratObj@assays[[targetAssayName]] <- Seurat::AddMetaData(seuratObj@assays[[targetAssayName]], metadata = fs, col.name = 'FeatureSource')
  }

  if (normalizeData) {
    if (targetAssayName == assayForLibrarySize) {
      print('Normalizing using Seurat::NormalizeData')
      seuratObj <- Seurat::NormalizeData(seuratObj, assay = targetAssayName, verbose = FALSE)
    } else {
      print(paste0('Normalizing using LogNormalizeUsingAlternateAssay with ', assayForLibrarySize))
      seuratObj <- LogNormalizeUsingAlternateAssay(seuratObj, assay = targetAssayName, assayForLibrarySize = assayForLibrarySize, maxLibrarySizeRatio = maxLibrarySizeRatio)
    }
  }

  if (doPlot) {
    if (!requireNamespace("RIRA", quietly = TRUE)) {
      warning("The RIRA package must be installed to perform plotting")
    } else {
      print('Plotting features')
      reductions <- intersect(names(seuratObj@reductions), c('umap', 'tsne', 'wnn'))
      if (length(reductions) == 0){
          print('No reductions, cannot plot')
      } else {
        feats <- paste0(seuratObj[[targetAssayName]]@key, rownames(seuratObj[[targetAssayName]]))
        rowSums <- Matrix::rowSums(Seurat::GetAssayData(seuratObj, assay = targetAssayName, layer = 'counts'))
        feats <- feats[rowSums > 0]
        if (length(feats) == 0) {
          print('All features are zero, skipping plot')
        } else {
          if (length(feats) > maxFeaturesToPlot){
              print(paste0('Too many features, will plot the first: ', maxFeaturesToPlot))
              feats <- feats[1:maxFeaturesToPlot]
          }

          RIRA::PlotMarkerSeries(seuratObj, reductions = reductions, features = feats, title = targetAssayName)
        }
      }
    }
  }
  
  return(seuratObj)
}

.GetAssayMetaSlotName <- function(assayObj) {
  slotName <- ifelse('meta.features' %in% methods::slotNames(assayObj), yes = 'meta.features', no = 'meta.data')
  if (! slotName %in% methods::slotNames(assayObj)) {
    stop(paste0('Assay object lacks slot: ', slotName))
  }

  return(slotName)
}