multiline_comments({
  "This module is used for 'qulity control' for single sample according to the <project_configure.yaml>"
  "1. Create 'Seaurat' object"
  "2. Calculate percent of <MT->, <RPL->, <HB->"
  "3. Add meta data"
  "4. Pick out cells according to the 'methods' in <project_configure.yaml>"
  "5. Pick out T/BCR genes and ribo genes in ench dataset and save to "
})

qulity_control <- function(num, addtcr) {
  # make dir for each sample output
  sample_path <- file.path(parameters$project_output_path, "individual", parameters$data[[num]]$orig.ident)
  make_dir(sample_path)
  filted_data_path <- file.path(sample_path, "filted_data")
  plot_path <- file.path(sample_path, "plot")
  log_path <- file.path(sample_path, "log")
  make_dir(filted_data_path)
  make_dir(plot_path)
  make_dir(log_path)

  count_data <- Read10X(parameters$data[[num]]$path)
  object_foo <- CreateSeuratObject(counts = count_data, project = parameters$data[[num]]$orig.ident, min.cells = 3, min.features = 200)

  # Calculate percent mitochondrial, add to seurat, and filter on it
  object_foo <- PercentageFeatureSet(object_foo, pattern = "^[Mm][Tt]-", col.name = "percent.mt")
  object_foo <- PercentageFeatureSet(object_foo, pattern = "^[Rr][Pp][LlSs]", col.name = "percent.ribo") # nolint: line_length_linter.
  hb_match <- match(species_genes$hb.genes, rownames(object_foo@assays$RNA))
  hb_genes <- rownames(object_foo@assays$RNA)[hb_match]
  hb_genes <- na.omit(hb_genes)
  object_foo <- PercentageFeatureSet(object_foo, features = hb_genes, col.name = "percent.hb")

  # doublets
  scr_py <- import("scrublet")
  counts_matrix <- object_foo@assays$RNA@counts
  counts_matrix <- t(counts_matrix)
  scrub <- scr_py$Scrublet(counts_matrix, expected_doublet_rate = parameters$scrub$doublet_rate)
  doublets_result <- scrub$scrub_doublets(verbose = FALSE)
  doublets <- as.data.frame(doublets_result)
  colnames(doublets) <- c("doublet_score", "is_doublet")
  rownames(doublets) <- rownames(object_foo@meta.data)
  write.csv(doublets, file.path(log_path, "scrublet_doublets_score.csv"), row.names = F)

  # add meta data
  object_foo <- AddMetaData(object_foo, doublets)

  # raw data plot
  make_violin_plot(object_foo, file.path(plot_path, "unfilted_violin_plot.png"), features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "doublet_score"))

  for (key in parameters$data$meta_data_key) {
    object_foo[[key]] <- parameters$data[[num]][[key]]
  }

  # pick cells out
  if (parameters$quality_control_method$standard_deviation$use) {
    mean_umi <- mean(object_foo@meta.data$nCount_RNA)
    sd_umi <- sd(object_foo@meta.data$nCount_RNA)
    mean_gene <- mean(object_foo@meta.data$nFeature_RNA)
    sd_gene <- sd(object_foo@meta.data$nFeature_RNA)
    maximum_umi <- mean_umi + parameters$quality_control_method$standard_deviation$confidence_interval * sd_umi
    maximum_gene <- mean_gene + parameters$quality_control_method$standard_deviation$confidence_interval * sd_gene
    minimum_umi <- parameters$quality_control_method$standard_deviation$minimum_umi
    minimum_gene <- parameters$quality_control_method$standard_deviation$minimum_gene
    cat_time(paste0("The mean of gene is: ", mean_gene))
    cat_time(paste0("The mean of umi is: ", mean_umi))
    cat_time(paste0("The SD of gene is: ", sd_gene))
    cat_time(paste0("The SD of umi is: ", sd_umi))
  } else if (parameters$quality_control_method$extremum$use) {
    maximum_umi <- parameters$quality_control_method$extremum$maximum_umi
    maximum_gene <- parameters$quality_control_method$extremum$maximum_gene
    minimum_umi <- parameters$quality_control_method$extremum$minimum_umi
    minimum_gene <- parameters$quality_control_method$extremum$minimum_gene
  } else {
    stop("The quality control method is wrong!!!Plwase check!!!")
  }

  # doublets
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$is_doublet == "TRUE", "Doublet", "Pass")

  # umi
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$nCount_RNA < minimum_umi & object_foo@meta.data$quality == "Pass", "Low_nCount", object_foo@meta.data$quality)
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$nCount_RNA < minimum_umi & object_foo@meta.data$quality != "Pass", paste("Low_nCount", object_foo@meta.data$quality, sep = ","), object_foo@meta.data$quality)
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$nCount_RNA > maximum_umi & object_foo@meta.data$quality == "Pass", "High_nCount", object_foo@meta.data$quality)
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$nCount_RNA > maximum_umi & object_foo@meta.data$quality != "Pass", paste("High_nCount", object_foo@meta.data$quality, sep = ","), object_foo@meta.data$quality)

  # gene
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$nFeature_RNA < minimum_gene & object_foo@meta.data$quality == "Pass", "Low_nFeature", object_foo@meta.data$quality)
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$nFeature_RNA < minimum_gene & object_foo@meta.data$quality != "Pass" & object_foo@meta.data$quality != "Low_nFeature", paste("Low_nFeature", object_foo@meta.data$quality, sep = ","), object_foo@meta.data$quality)
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$nFeature_RNA > maximum_gene & object_foo@meta.data$quality == "Pass", "High_nFeature", object_foo@meta.data$quality)
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$nFeature_RNA > maximum_gene & object_foo@meta.data$quality != "Pass" & object_foo@meta.data$quality != "High_nFeature", paste("High_nFeature", object_foo@meta.data$quality, sep = ","), object_foo@meta.data$quality)

  # MT
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$percent.mt > parameters$quality_control_method$filtering_standard$mt & object_foo@meta.data$quality == "Pass", "High_MT", object_foo@meta.data$quality)
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$percent.mt > parameters$quality_control_method$filtering_standard$mt & object_foo@meta.data$quality != "Pass" & object_foo@meta.data$quality != "High_MT", paste("High_MT", object_foo@meta.data$quality, sep = ","), object_foo@meta.data$quality)

  # HB
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$percent.hb > parameters$quality_control_method$filtering_standard$hb & object_foo@meta.data$quality == "Pass", "High_Erythrocyte", object_foo@meta.data$quality)
  object_foo[["quality"]] <- ifelse(object_foo@meta.data$percent.hb > parameters$quality_control_method$filtering_standard$hb & object_foo@meta.data$quality != "Pass" & object_foo@meta.data$quality != "High_Erythrocyte", paste("High_Erythrocyte", object_foo@meta.data$quality, sep = ","), object_foo@meta.data$quality)

  cell_statistic <- as.data.frame(table(object_foo[["quality"]]))
  write.csv(cell_statistic, file.path(log_path, "cell_statistic.csv"), row.names = F)

  object_foo <- subset(object_foo, subset = quality == "Pass")

  object_foo <- RenameCells(object_foo, new.names = paste(parameters$data[[num]]$orig.ident, colnames(object_foo), sep = "_"))

  # T/BCR add
  if (addtcr) {
      if (file.exists(file.path(log_path, 'BCR.csv'))) {
        file.remove(file.path(log_path, 'BCR.csv'))
      }
      if (file.exists(file.path(log_path, 'TCR.csv'))) {
        file.remove(file.path(log_path, 'TCR.csv'))
      }
      tbcr_filt <- source_python('./libs/tbcr.py')
      count_umi(input_file_path = parameters$data[[num]]$bcrpath, 
                sample_name = parameters$data[[num]]$orig.ident, 
                output_file_name = file.path(log_path, 'BCR.csv'))
      count_umi(input_file_path = parameters$data[[num]]$tcrpath, 
                sample_name = parameters$data[[num]]$orig.ident, 
                output_file_name = file.path(log_path, 'TCR.csv'))
      bcr_csv <- read.csv(file.path(log_path, 'BCR.csv'))
      tcr_csv <- read.csv(file.path(log_path, 'TCR.csv'))
      bcr_csv <- column_to_rownames(bcr_csv, var = "sample_name")
      tcr_csv <- column_to_rownames(tcr_csv, var = "sample_name")
      shared_row_names <- intersect(rownames(bcr_csv), rownames(tcr_csv))
      shared_row_dataframe <- cbind(bcr_csv[shared_row_names,], tcr_csv[shared_row_names,])
      unshared_row_dataframe <-  bind_rows(bcr_csv[setdiff(rownames(bcr_csv), shared_row_names),], 
                                          tcr_csv[setdiff(rownames(tcr_csv), shared_row_names),])  
      tbcr_meta <- rbind(shared_row_dataframe, unshared_row_dataframe)
      object_foo <- AddMetaData(object_foo, tbcr_meta)
    }

  # save filted data
  cat_time("Write individual data...")
  saveRDS(object_foo, file.path(filted_data_path, "filted.rds"))

  # filted data plot
  make_violin_plot(object_foo, file.path(plot_path, "filted_violin_plot.png"), features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "doublet_score"))
  make_scatter_plot(object_foo, file.path(plot_path, "filted_scatter_plot.png"))

  return(object_foo)
}
