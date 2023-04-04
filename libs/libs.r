suppressPackageStartupMessages(library(yaml))

base_params <- yaml.load_file("./base_configure.yaml")

# long comments
multiline_comments <- function(blah) {
  cat(NULL)
}

# log timestamp output
cat_time <- function(x, level = "info") {
  if (level == "info") {
    cat(format(Sys.time(), "[%H:%M:%S]"), "[INFO]", x, "\n")
  } else if (level == "warning") {
    cat(format(Sys.time(), "[%H:%M:%S]"), "[WARNING]", x, "\n")
  }
}

# cat string generated by sprintf
catf <- function(format, ...) {
  cat_time(sprintf(format, ...))
}

# make directory if it doesn't already exist
make_dir <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

# Function to make PC variation plots based on PC standard deviation
make_pc_plots <- function(object, file, jackstraw = NULL) {
  sdev <- object@reductions$pca@stdev
  var <- sdev^2
  cum_var <- cumsum(var) / sum(var)
  df <- data.frame(PC = 1:length(cum_var), var = var, cum_var = cum_var)

  pdf(file)

  # Make cumulative variation plot
  p1 <- ggplot(df, aes(x = PC, y = cum_var)) +
    geom_point() +
    expand_limits(y = 0)
  p1 <- p1 + labs(x = "PC", y = "Cumulative Proportion of Variation", title = "PC Cumulative Proportion of Variation")
  print(p1)

  # Make variation plot
  p2 <- ggplot(df, aes(x = PC, y = var)) +
    geom_point() +
    expand_limits(y = 0)
  p2 <- p2 + labs(x = "PC", y = "Variance", title = "PC Variance")
  print(p2)

  # Make JackStraw plot
  if (!is.null(jackstraw)) {
    print(JackStrawPlot(object, dims = 1:jackstraw))
  }

  invisible(dev.off())
}

# Write matrix to tsv, works for sparse matrices and adds column name Gene to rownames
matrix_to_tsv <- function(mat, file) {
  write_tsv(rownames_to_column(as.data.frame(as.matrix(mat)), "GENE"), file)
}

generate_file_name <- function(data_type, filt_method=NULL,
                               norm_method=NULL, norm_features=NULL, keeptcr=NULL, 
                               pcs=NULL, 
                               integrate_method=NULL, integrate_pcs=NULL,
                               reduction=NULL, reduction_dims=NULL,
                               resolution=NULL
                               
                               ) {
  data_type <- paste0("{", data_type, "}")
  if (filt_method == "SD") {
    filt <- paste0("{", filt_method,"-", parameters$quality_control_method$standard_deviation$minimum_umi, "-", parameters$quality_control_method$standard_deviation$confidence_interval,  "SD}")
  } else {
    filt <- paste0("{", filt_method, "-gene-", 
                   parameters$quality_control_method$extremum$minimum_gene, "-", 
                   parameters$quality_control_method$extremum$maximum_gene, "-UMI-", 
                   parameters$quality_control_method$extremum$minimum_umi, "-", 
                   parameters$quality_control_method$extremum$maximum_umi,  "}")
  }
  
  if (!is.null(norm_method)) {
    norm <- paste0("{", norm_method, "-", norm_features, "}")
  } else {
    norm <- "{NOT-NORMLIZE}"
  }
  
  if (!is.null(keeptcr)) {
    tcr <- paste0("{", "KEEPTCRGENE", "-", keeptcr, "}")
  } else {
    tcr <- "{NOT-PROCESS-TCR}"
  }
  
  if (!is.null(pcs)) {
    PCA <- paste0("{", "PCA", "-", pcs, "}")
  } else {
    PCA <- "{NOT-PCA}"
  }
  
  if (!is.null(integrate_method)) {
    integrate <- paste0("{", integrate_method, "-", integrate_pcs, "}")
  } else {
    integrate <- "{NOT-INTEGRATED}"
  }
  
  if (!is.null(reduction)) {
    reduc <- paste0("{", reduction, "-", reduction_dims, "}")
  } else {
    reduc <- "{NOT-REDUCTION}"
  }
  
  if (!is.null(resolution)) {
    plot <- paste0("{RESOLUTION-", resolution, "}")
  } else {
    plot <- "{NOT-REDUCTION}"
  }
  
  
  
  file_name <- paste(data_type, filt, norm, tcr, PCA, integrate, reduc, plot, sep = "_")
  return(file_name)
}
