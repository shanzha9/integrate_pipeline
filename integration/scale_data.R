multiline_comments({
  "This module is used for 'normalize the data' based on individual data. Now, the 'NormalizeData' and 'SCTransform' in 'Seurat' are support"
  "1. Normalizethe data"
  "2. Find variable features"
  "3. Choose whether to remove T/BCR genes and ribo genes"
  "4. Normalization, find variable features, and scale data"
})

define_genetoregress <- function(to_regress) {
  if (to_regress == "all") {
    varstoregress <- c("percent.mt", "s.score", "g2m.score", "percent.ribo")
  } else if (to_regress == "difference") {
    object_foo$cc.difference <- object_foo$s.score - object_foo$g2m.score
    object_foo$s.score <- NULL
    object_foo$g2m.score <- NULL
    varstoregress <- c("percent.mt", "cc.difference", "percent.ribo")
  } else if (to_regress == "none") {
    varstoregress <- c("percent.mt", "percent.ribo")
  } else {
    stop("Wrong genes to be regressed out!!! Please check config Yaml file!!!")
  }
  return(varstoregress)
}

calculate_cellcyclesocre <- function(object_foo) {
  # Cell-Cycle Scoring and Regression
  cat_time("Cell-Cycle Scoring...")
  object_foo <- CellCycleScoring(object_foo,
    s.features = species_genes$s.genes,
    g2m.features = species_genes$g2m.genes
  )
  object_foo$s.score <- object_foo$S.Score
  object_foo$S.Score <- NULL
  object_foo$g2m.score <- object_foo$G2M.Score
  object_foo$G2M.Score <- NULL
  object_foo$phase <- object_foo$Phase
  object_foo$Phase <- NULL
  return(object_foo)
}

normalize_data <- function(object_foo, method, toregress1, toregress2, two_step_norm) {
  if (method == "sct") {
    cat_time("Step one sct normalization...")
    # scale data with regress
    varstoregress <- define_genetoregress(toregress1)
    object_foo <- SCTransform(object_foo,
      verbose = FALSE,
      variable.features.n = parameters$integration$normalization$sctransform$nfeatures,
      vars.to.regress = varstoregress
    )

    # scale data with regress
    if (two_step_norm == TRUE) {
      object_foo <- calculate_cellcyclesocre(object_foo)
      cat_time("Step two sct normalization...")
      varstoregress <- define_genetoregress(toregress2)
      object_foo <- SCTransform(object_foo,
        verbose = FALSE,
        variable.features.n = parameters$integration$normalization$sctransform$nfeatures,
        vars.to.regress = varstoregress
      )
    }
  } else {
    cat_time("Step one log normalization...")
    # scale data with regress
    varstoregress <- define_genetoregress(toregress1)
    object_foo <- NormalizeData(object_foo,
      normalization.method = parameters$integration$normalization$lognorlization$method,
      verbose = FALSE
    )
    object_foo <- FindVariableFeatures(object_foo,
      selection.method = "vst",
      nfeatures = parameters$integration$normalization$lognorlization$nfeatures,
      verbose = FALSE
    )

    # scale data with regress
    if (two_step_norm == TRUE) {
      object_foo <- calculate_cellcyclesocre(object_foo)
      cat_time("Step two log normalization...")
      varstoregress <- define_genetoregress(toregress2)
      object_foo <- ScaleData(object_foo,
        vars.to.regress = varstoregress,
        verbose = FALSE
      )
    }
  }

  return(object_foo)
}
