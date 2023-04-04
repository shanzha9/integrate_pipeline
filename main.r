source("./libs/libs.r", chdir = FALSE)
source("./preparation/pre-check.R")
source("./preparation/specific_constant.R")
source("./visualization/visualization.R")
source("./qulity_control/QC.R")
source("./integration/scale_data.R")

# load configure file
parameters <- yaml.load_file("./configure.yaml")

# make dir
setwd(parameters$project_root_path)
make_dir(parameters$project_output_path)
integrated_path <- file.path(parameters$project_output_path, "integrated")
individual_path <- file.path(parameters$project_output_path, "individual")
rdsdata_path <- file.path(parameters$project_output_path, "data")
make_dir(individual_path)
make_dir(integrated_path)
make_dir(rdsdata_path)
make_dir(file.path(integrated_path, "pca"))
make_dir(file.path(integrated_path, "harmony"))
make_dir(file.path(integrated_path, "cluster"))
make_dir(file.path(integrated_path, "RPCA"))
make_dir(file.path(integrated_path, "CCA"))
make_dir(file.path(integrated_path, "plot"))

# pre-check
check_packages(base_params)

# load packages
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(tidyr))

# check whether the parameters set is right
pass <- check_parameters(parameters)
if (!pass) {
    stop("Please check config Yaml file!!!See Warning info!!!")
} else {
    cat_time("All samples config Yaml file check pass!!!")
}

# define specific constant
species_genes <- check_species(parameters)


# decide the baseline data of integration
if (parameters$project_baseline$cellranger_data$use) {
    cat_time("Use cellranger data...")
    cat_time("Quality control...")
    if (parameters$quality_control_method$standard_deviation$use) {
      filt_method <- "SD"
    } else {
      filt_method <- "Max-Min"
    }
    filted_data_list <- c()
    for (i in 2:length(parameters$data)) {
        cat_time(paste0(length(parameters$data) - 1, " samples need to be processed, number ", i - 1, " is proccessing...", sep = ""))
        filted_data <- qulity_control(i, addtcr=parameters$addtcttometa)
        filted_data_list <- c(filted_data_list, filted_data)
    }

    # save filted data list
    make_dir(file.path(integrated_path, "QC"))
    cat_time("Write filted_data_list data...")
    file_name <- generate_file_name(data_type = "filted_data_list",
                                    filt_method = filt_method)
    
    saveRDS(filted_data_list, file.path(integrated_path, "QC", paste0(file_name, '.rds')))

    # merge
    seurat_object_merge <- merge(filted_data_list[[1]], filted_data_list[-1], merge.data = TRUE)
    # TCR BCR log(umi + 1)
    seurat_object_merge$bcr_umi <- replace_na(seurat_object_merge$bcr_umi, 0)
    seurat_object_merge$bcr_logumi <- log2(seurat_object_merge$bcr_umi + 1)
    seurat_object_merge$tcr_umi <- replace_na(seurat_object_merge$tcr_umi, 0)
    seurat_object_merge$tcr_logumi <- log2(seurat_object_merge$tcr_umi + 1)
    # keep meta
    old_meta_keep <- seurat_object_merge@meta.data
    write_csv(
      rownames_to_column(as.data.frame(old_meta_keep), "barcodes"),
      file.path(parameters$project_output_path, "old_meta_data.csv")
    )
    
    # write data
    cat_time("Write seurat_object_merge data...")
    file_name <- generate_file_name(data_type = "FILTED_DATA_MERGE", 
                                    filt_method = filt_method)
    saveRDS(seurat_object_merge, file.path(rdsdata_path, paste0(file_name, '.rds')))
    
} else if (parameters$project_baseline$filted_data_list$use) {
    cat_time("Use filted data list...")
    # read data
    cat_time("Load raw data list...")
    filt_method <- "NONE"
    filted_data_list <- readRDS(parameters$project_baseline$filted_data_list$data_path)
    seurat_object_merge <- merge(filted_data_list[[1]], filted_data_list[-1], merge.data = TRUE)
    cat_time("Clean old data...")
    # keep meta
    old_meta_keep <- seurat_object_merge@meta.data
    write_csv(
        rownames_to_column(as.data.frame(old_meta_keep), "barcodes"),
        file.path(parameters$project_output_path, "old_meta_data.csv")
    )
    cat_time("Clean data and calculate percent.X...")
    seurat_object_merge@active.assay <- "RNA"
    # clean old data
    for (i in names(seurat_object_merge@assays)) {
        if (i != "RNA") {
            seurat_object_merge@assays[i] <- NULL
        }
    }
    for (i in names(seurat_object_merge@reductions)) {
        seurat_object_merge@reductions[i] <- NULL
    }
    for (i in names(seurat_object_merge@graphs)) {
        seurat_object_merge@graphs[i] <- NULL
    }
    for (i in names(seurat_object_merge@neighbors)) {
        seurat_object_merge@neighbors[i] <- NULL
    }
    seurat_object_merge <- SetIdent(seurat_object_merge,
        cells = rownames(seurat_object_merge@meta.data),
        value = seurat_object_merge@meta.data$orig.ident
    )
    # clean meta
    seurat_object_merge@meta.data <- seurat_object_merge@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA")]
    # calculate percent.X
    seurat_object_merge <- PercentageFeatureSet(seurat_object_merge, pattern = "^[Mm][Tt]-", col.name = "percent.mt")
    seurat_object_merge <- PercentageFeatureSet(seurat_object_merge, pattern = "^[Rr][Pp][LlSs]", col.name = "percent.ribo") # nolint: line_length_linter.
    hb_match <- match(species_genes$hb.genes, rownames(seurat_object_merge@assays$RNA))
    hb_genes <- rownames(seurat_object_merge@assays$RNA)[hb_match]
    hb_genes <- na.omit(hb_genes)
    seurat_object_merge <- PercentageFeatureSet(seurat_object_merge, features = hb_genes, col.name = "percent.hb")
    # write data
    cat_time("Write clean data...")
    saveRDS(seurat_object_merge, file.path(rdsdata_path, "raw_data_after_clean.rds"))
} else if (parameters$project_baseline$filted_data_merge$use) {
    # read data
    filt_method <- "NONE"
    cat_time("Use filted merge data...")
    cat_time("Load raw data list...")
    cat_time("Clean old data...")
    seurat_object_merge <- readRDS(parameters$project_baseline$filted_data_merge$data_path)
    # keep meta
    old_meta_keep <- seurat_object_merge@meta.data
    write_csv(
        rownames_to_column(as.data.frame(old_meta_keep), "barcodes"),
        file.path(parameters$project_output_path, "old_meta_data.csv")
    )
    cat_time("Clean data and calculate percent.X...")
    seurat_object_merge@active.assay <- "RNA"
    # clean old data
    for (i in names(seurat_object_merge@assays)) {
        if (i != "RNA") {
            seurat_object_merge@assays[i] <- NULL
        }
    }
    for (i in names(seurat_object_merge@reductions)) {
        seurat_object_merge@reductions[i] <- NULL
    }
    for (i in names(seurat_object_merge@graphs)) {
        seurat_object_merge@graphs[i] <- NULL
    }
    for (i in names(seurat_object_merge@neighbors)) {
        seurat_object_merge@neighbors[i] <- NULL
    }
    seurat_object_merge <- SetIdent(seurat_object_merge,
        cells = rownames(seurat_object_merge@meta.data),
        value = seurat_object_merge@meta.data$orig.ident
    )
    # clean meta
    seurat_object_merge@meta.data <- seurat_object_merge@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA")]
    # calculate percent.X
    seurat_object_merge <- PercentageFeatureSet(seurat_object_merge, pattern = "^[Mm][Tt]-", col.name = "percent.mt")
    seurat_object_merge <- PercentageFeatureSet(seurat_object_merge, pattern = "^[Rr][Pp][LlSs]", col.name = "percent.ribo") # nolint: line_length_linter.
    hb_match <- match(species_genes$hb.genes, rownames(seurat_object_merge@assays$RNA))
    hb_genes <- rownames(seurat_object_merge@assays$RNA)[hb_match]
    hb_genes <- na.omit(hb_genes)
    seurat_object_merge <- PercentageFeatureSet(seurat_object_merge, features = hb_genes, col.name = "percent.hb")
    # write data
    cat_time("Write clean data...")
    saveRDS(seurat_object_merge, file.path(rdsdata_path, "raw_data_after_clean.rds"))
} else {
    stop("Wrong baseline data!!!Please check config Yaml file!!!")
}

seurat_object_merge@meta.data <- seurat_object_merge@meta.data[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.hb', 'percent.mt', 'percent.ribo')]
# pick out genes to aviod unique genes in uncertain dataset
cat_time("Match RP- IG- MT- TR- genes...")
species_genes$ribo <- grep(pattern = "^[Rr][Pp][LlSs]", x = rownames(x = seurat_object_merge@assays$RNA), value = TRUE)
species_genes$igh <- grep(pattern = "^[Ii][Gg][HhKkLl][CcVvJjGgDdAa]", x = rownames(x = seurat_object_merge@assays$RNA), value = TRUE)
species_genes$mt <- grep(pattern = "^[Mm][Tt]-", x = rownames(x = seurat_object_merge@assays$RNA), value = TRUE)
species_genes$tr <- grep(pattern = "^[Tt][Rr][AaBbPp][CcVvJjGgDd]", x = rownames(x = seurat_object_merge@assays$RNA), value = TRUE)


# scale data and regress out batch effects
if (parameters$integration$normalization$sctransform$use) {
    cat_time("Use SCTransform to normalize the data!!!")
    norm_method <- "SCT"
    norm_features <- parameters$integration$normalization$sctransform$nfeatures
    seurat_object_merge <- normalize_data(seurat_object_merge,
        method = "sct",
        toregress1 = "none",
        toregress2 = "all",
        two_step_norm = TRUE
    )
} else if (parameters$integration$normalization$lognorlization$use) {
    cat_time("Use Log to normalize the data!!!")
    norm_method <- "LOG"
    norm_features <- parameters$integration$normalization$lognorlization$nfeatures
    seurat_object_merge <- normalize_data(seurat_object_merge,
        method = "log",
        toregress1 = "none",
        toregress2 = "all",
        two_step_norm = TRUE
    )


} else {
    stop("Wrong 'integration norlization' methods!!!Please check config Yaml file!!!")
}

cat_time("Write clean data...")
file_name <- generate_file_name(data_type = "FILTED_DATA_MERGE", 
                                filt_method = filt_method,
                                norm_method = norm_method, norm_features = norm_features)
saveRDS(seurat_object_merge, file.path(rdsdata_path, paste0(file_name, '.rds')))



# SelectIntegrationFeatures
cat_time("Splity object by batch...")
filted_data_list <- SplitObject(seurat_object_merge, split.by = "orig.ident")

cat_time("Sample list proccessing...")
for (i in 1:length(filted_data_list)) {
    cat_time(paste0(length(filted_data_list), " samples need to be processed, number ", i, " is proccessing...", sep = ""))
    if (parameters$integration$normalization$sctransform$use) {
        filted_data_list[[i]] <- normalize_data(filted_data_list[[i]],
            method = "sct",
            toregress1 = "all",
            toregress2 = "NULL",
            two_step_norm = FALSE
        )
    } else if (parameters$integration$normalization$lognorlization$use) {
        filted_data_list[[i]] <- normalize_data(filted_data_list[[i]],
            method = "log",
            toregress1 = "all",
            toregress2 = "NULL",
            two_step_norm = FALSE
        )
    } else {
        stop("Wrong 'integration norlization' methods!!!Please check config Yaml file!!!")
    }
}

# write data.
cat_time("Write clean data...")
file_name <- generate_file_name(data_type = "MERGE_DATA_BATCH_LIST", 
                                  filt_method = filt_method,
                                  norm_method = norm_method, norm_features = norm_features)
saveRDS(filted_data_list, file.path(rdsdata_path, paste0(file_name, ".rds")))


cat_time("Select integration features...")
integrate_features <- SelectIntegrationFeatures(filted_data_list,
    nfeatures = parameters$reduction$nfeaturesforpca
)
markers_todelete <- c(species_genes$g2m.genes, species_genes$s.genes, 
                      species_genes$ribo, species_genes$mt, species_genes$hb.genes)
integrate_features <- setdiff(integrate_features, markers_todelete)
integrate_features_withtcr <- integrate_features
integrate_features_withouttcr <- setdiff(integrate_features, c(species_genes$igh, species_genes$tr))

if (!parameters$keeptcr) {
  cat_time("Remove TCR genes and BCR genes from integration features...")
  integrate_features <- integrate_features_withouttcr
} else {
  integrate_features <- integrate_features_withtcr
}


cat_time("write integration features...")
write_csv(
    as.data.frame(integrate_features),
    file.path(parameters$project_output_path, paste0("integrate_features_", parameters$reduction$nfeaturesforpca, "_keeptcr_", parameters$keeptcr,".csv"))
)


# integrate data
if (parameters$integration$method$harmony$use) {
    cat_time("use Harmony to integrate data!!!")
    # PCA
    cat_time("Merge data...")
    integrate_method <- "Harmony"
    integrate_pcs = parameters$integration$method$harmony$ndims
    seurat_object_merge <- merge(filted_data_list[[1]], filted_data_list[-1], merge.data = TRUE)
    VariableFeatures(seurat_object_merge) <- integrate_features
    cat_time(paste0("RunPCA, use assay ", seurat_object_merge@active.assay, "..."))
    seurat_object_merge <- RunPCA(seurat_object_merge,
        assay.use = seurat_object_merge@active.assay,
        verbose = FALSE,
        npcs = parameters$reduction$pca$npcs,
        features = integrate_features
    )
    cat_time("Write clean data...")
    file_name <- generate_file_name(data_type = "MERGE_DATA_BATCH_LIST", 
                                    filt_method = filt_method,
                                    norm_method = norm_method, norm_features = 1500,
                                    keeptcr = FALSE,
                                    pcs = parameters$reduction$pca$npcs)
    saveRDS(seurat_object_merge, file.path(integrated_path, "pca", paste0(file_name, ".rds")))
    cat_time("Make PCA plot...")
    if (is.null(parameters$reduction$pca$jackstraw)) {
        make_pc_plots(seurat_object_merge, file.path(integrated_path, "pca", paste0(file_name, "_variation.pdf")))
    } else {
        make_pc_plots(seurat_object_merge, file.path(integrated_path, "pca", paste0(file_name, "_variation.pdf"), parameters$reduction$pca$jackstraw))
    }
    write_tsv(
        rownames_to_column(as.data.frame(seurat_object_merge@reductions$pca@feature.loadings), "GENE"),
        file.path(integrated_path, "pca", paste0(file_name, "_loadings.tsv"))
    )
    # define 'group.by.vars' and 'theta' in RunHarmony
    groupby_set <- parameters$integration$method$harmony$groupby
    lamda_set <- parameters$integration$method$harmony$lamda

    # define ndims to run
    cat_time("RunHarmony...")
    for (i in 1:length(parameters$integration$method$harmony$ndims)) {
        ndims <- parameters$integration$method$harmony$ndims[[i]]
        cat_time(paste0(length(parameters$integration$method$harmony$ndims), " harmony dims need to be processed, ", ndims, " is proccessing...", sep = ""))
        reduction_key <- paste0("harmonyPC", ndims, spe = "")
        seurat_object_merge <- RunHarmony(seurat_object_merge,
            assay.use = seurat_object_merge@active.assay,
            reduction = "pca",
            dims.use = 1:ndims,
            group.by.vars = 'orig.ident',
            plot_convergence = FALSE,
            reduction.save = reduction_key,
            verbose = FALSE
        )
    }
    cat_time("Write clean data...")
    file_name <- generate_file_name(data_type = "MERGE_DATA_BATCH_LIST", 
                                           filt_method = filt_method,
                                           norm_method = norm_method, norm_features = norm_features,
                                           keeptcr = parameters$keeptcr,
                                           pcs = parameters$reduction$pca$npcs,
                                           integrate_method = integrate_method, integrate_pcs = parameters$integration$method$harmony$ndims
                                           )
    saveRDS(seurat_object_merge, file.path(integrated_path, "harmony", paste0(file_name, ".rds")))
} else if (parameters$integration$method$RPCA$use) {
    cat_time("Run RPCA...")
    integrate_method <- "RPCA"
    if (parameters$integration$normalization$sctransform$use) {
        cat_time("Run PrepSCTIntegration...")
        filted_data_list <- PrepSCTIntegration(
            object.list = filted_data_list,
            anchor.features = integrate_features
        )
        normalization_method <- "SCT"
    } else {
        normalization_method <- "LogNormalize"
    }
    filted_data_list <- lapply(filted_data_list,
                               Seurat::RunPCA,
                               features = integrate_features,
                               verbose = FALSE
    )
    cat_time("Run FindIntegrationAnchors...")
    seurat_anchors <- FindIntegrationAnchors(
        object.list = filted_data_list,
        reduction = "rpca",
        normalization.method = normalization_method,
        anchor.features = integrate_features,
        k.anchor = 20,
        verbose = F
        
    )
    cat_time("Run IntegrateData...")
    seurat_object_merge <- IntegrateData(
        anchorset = seurat_anchors,
        normalization.method = normalization_method,
        dims = 1:parameters$integration$method$RPCA$ndims,
        verbose = F
    )
    integrate_pcs <- parameters$integration$method$RPCA$ndims
    cat_time("Write clean data...")
    file_name <- generate_file_name(data_type = "MERGE_DATA_BATCH_LIST", 
                             filt_method = filt_method,
                             norm_method = norm_method, norm_features = norm_features,
                             keeptcr = parameters$keeptcr,
                             pcs = parameters$reduction$pca$npcs,
                             integrate_method = integrate_method, 
                             integrate_pcs = parameters$integration$method$RPCA$ndims)
    saveRDS(seurat_object_merge, file.path(integrated_path, "RPCA", paste0(file_name, ".rds")))
} else if (parameters$integration$method$CCA$use) {
    if (parameters$integration$normalization$sctransform$use) {
        cat_time("Run PrepSCTIntegration...")
        filted_data_list <- PrepSCTIntegration(
            object.list = filted_data_list,
            anchor.features = integrate_features
        )
        normalization_method <- "SCT"
    } else {
        normalization_method <- "LogNormalize"
    }
    cat_time("Run FindIntegrationAnchors...")
    seurat_anchors <- FindIntegrationAnchors(
        object.list = filted_data_list,
        normalization.method = normalization_method,
        anchor.features = integrate_features
    )
    cat_time("Run IntegrateData...")
    seurat_object_merge <- IntegrateData(
        anchorset = seurat_anchors,
        normalization.method = normalization_method,
        dims = 1:parameters$integration$method$CCA$ndims
    )
} else {
    NULL
}

# Run UMAp, tSNE and build SNN
if (parameters$integration$method$harmony$use) {
    for (i in names(seurat_object_merge@reductions)) {
        if (grepl("harmony", i)) {
            cat_time(paste(i, "is run FindNeighbors..."))
            seurat_object_merge <- FindNeighbors(seurat_object_merge,
                reduction = i,
                assay = seurat_object_merge@active.assay,
                dims = 1:parameters$cluster$ndims, verbose = FALSE
            )
            cat_time(paste(i, "is run FindClusters..."))
            seurat_object_merge <- FindClusters(seurat_object_merge,
                resolution = c(0.2, 0.5),
                verbose = FALSE
            )
            cat_time(paste(i, "is RunUMAP..."))
            seurat_object_merge <- RunUMAP(seurat_object_merge,
                reduction = i,
                dims = 1:parameters$cluster$ndims,
                verbose = FALSE,
                assay = seurat_object_merge@active.assay,
                umap.method = "umap-learn"
            )
            # cat_time(paste(i, "is RunTSNE..."))
            # seurat_object_merge <- RunTSNE(seurat_object_merge,
            #     reduction = i,
            #     dims = 1:parameters$cluster$ndims,
            #     assay = seurat_object_merge@active.assay,
            #     verbose = FALSE
            # )
            file_name <- generate_file_name(data_type = "MERGE_DATA_BATCH_LIST", 
                                            filt_method = filt_method,
                                            norm_method = norm_method, norm_features = norm_features,
                                            keeptcr = FALSE,
                                            pcs = parameters$reduction$pca$npcs,
                                            integrate_method = integrate_method, integrate_pcs = integrate_pcs,
                                            reduction = "UMAP-TSNE", reduction_dims = parameters$cluster$ndims
            )
            cat_time("Write clean data...")
            saveRDS(seurat_object_merge, file.path(integrated_path, "cluster", paste0(file_name, ".rds")))
        }
    }
} else {
    # original unmodified data still resides in the 'RNA' assay
    seurat_object_merge@active.assay <- "integrated"
    # Run the standard workflow for visualization and clustering
    cat_time("Run PCA...")
    # seurat_object_merge  <- ScaleData(object = seurat_object_merge)
    # seurat_object_merge <- FindVariableFeatures(seurat_object_merge,
    #                                             assay="integrated",
    #                                             nFeatures=1500)
    seurat_object_merge <- RunPCA(seurat_object_merge,
        npcs = parameters$reduction$pca$npcs,
        verbose = FALSE
        # features = integrate_features
    )
    cat_time("RunUMAP...")
    seurat_object_merge <- RunUMAP(seurat_object_merge,
        reduction = "pca",
        dims = 1:parameters$cluster$ndims,
        verbose = FALSE,
        umap.method = 'umap-learn',
        metric = 'correlation',
        min.dist = 0.2
    )
    cat_time("RunTSNE...")
    # seurat_object_merge <- RunTSNE(seurat_object_merge,
    #     reduction = "pca",
    #     dims = 1:parameters$cluster$ndims,
    #     verbose = FALSE
    # )
    cat_time("FindNeighbors...")
    seurat_object_merge <- FindNeighbors(seurat_object_merge,
        reduction = "pca",
        dims = 1:parameters$cluster$ndims,
        verbose = FALSE
    )
    cat_time("FindClusters...")
    seurat_object_merge <- FindClusters(seurat_object_merge,
        resolution = c(0.2, 0.7),
        verbose = FALSE
    )
    file_name <- generate_file_name(data_type = "MERGE_DATA_BATCH_LIST", 
                                    filt_method = filt_method,
                                    norm_method = norm_method, norm_features = norm_features,
                                    keeptcr = FALSE,
                                    pcs = parameters$reduction$pca$npcs,
                                    integrate_method = integrate_method, integrate_pcs = integrate_pcs,
                                    reduction = "UMAP-TSNE", reduction_dims = parameters$cluster$ndims,
                                    resolution = "0.2-0.7"
    )
    cat_time("Write clean data...")
    seurat_object_merge <- AddMetaData(seurat_object_merge, old_meta_keep)
    saveRDS(seurat_object_merge, file.path(integrated_path, "cluster", paste0(file_name, ".rds")))
 }



# plot and calculate marker genes
multiline_comments({
    "Visualization of the results:"
    "1. Dimplot of umap and tsne"
    "2. VlnPlot of all features"
    "3. FeaturePlot of selected features"
    "4. Calculate marker genes"
    "Step:"
    "1. merge old data: keep the key in new data"
    "2. merge the define marker genes"
    "3. loop cluster and plot"
    "4. Calculate marker genes"
})

# step1
cat_time("Plot given marker genes...")
if (!is.null(old_meta_keep)) {
  key_to_remove_in_old_meta <- intersect(colnames(old_meta_keep), colnames(seurat_object_merge@meta.data))
  old_meta_keep <- select(old_meta_keep, select = -key_to_remove_in_old_meta)
  seurat_object_merge <- AddMetaData(old_meta_keep)
}


# step2
vlnplot_markers <- c()
featureplt_markers <- c()
for (i in names(parameters$markergenes)) {
    if (i != "featureplot") {
        vlnplot_markers <- unique(c(vlnplot_markers, parameters$markergenes[[i]]))
    } else {
        featureplt_markers <- unique(parameters$markergenes[["featureplot"]])
    }
}

# step3
if (integrate_method == "SCT") {
  ident_cluster <- c("SCT_snn_res.0.1", "SCT_snn_res.0.2", "SCT_snn_res.0.3",
                     "SCT_snn_res.0.4", "SCT_snn_res.0.5", "SCT_snn_res.0.6",
                     "SCT_snn_res.0.7", "SCT_snn_res.0.8", "SCT_snn_res.0.9"
                     )
}

for (i in ident_cluster) {
  if (i %in% colnames(seurat_object_merge@meta.data)) {
    cat_time(paste0(i, " is processing..."))
    seurat_object_merge <- SetIdent(seurat_object_merge, cells = rownames(seurat_object_merge@meta.data), value = seurat_object_merge@meta.data[[i]])
    file_name <- generate_file_name(data_type = "merge_data_batch_list", 
                                    filt_method = filt_method,
                                    norm_method = norm_method, norm_features = norm_features,
                                    keeptcr = FALSE,
                                    pcs = parameters$reduction$pca$npcs,
                                    integrate_method = integrate_method, integrate_pcs = integrate_pcs,
                                    reduction = "UMAP-TSNE", reduction_dims = parameters$cluster$ndims,
                                    resolution = i)
                                    
    marker_gene_dim_plot(seurat_object_merge, 
                     file.path(integrated_path, "plot", paste0(file_name, ".pdf")),
                     vlnplot_markers)
    
    if (!is.null(featureplt_markers)) {
      marker_gene_plot(seurat_object_merge, 
                       file.path(integrated_path, "plot", paste0(file_name, ".pdf")),
                       vlnplot_markers)
    }
    
    # marker gene
    if (!is.null()) {
      cat_time("Calculate marker genes...")
      seurat_object_merge <- PrepSCTFindMarkers(seurat_object_merge, verbose = T)
      markers_all <- FindAllMarkers(seurat_object_merge, assay = "SCT", only.pos = T)
      
      markers_top <- markers_all %>%
        group_by(cluster) %>%
        top_n(n = 15, wt = avg_log2FC)
      
      write_csv(
        markers_all,
        file.path(parameters$project_output_path, paste0(file_name, "_Allmarkers.csv"))
      )
      
      write_csv(
        markers_top,
        file.path(parameters$project_output_path, paste0(file_name, "_AllmarkersTop.csv"))
      )
    }
  }
}


# label major cluster 
cat_time("Label cell clusters...")
label_cluster_json <- fromJSON(file = file.path(parameters$project_root_path, 'label.json'))
seurat_object_merge@meta.data[["majorType"]] <- "Undefined"
for (i in names(label_cluster_json)) {
  cell_cluster <- str_split(i, pattern = "_")[[1]][2]
  seurat_object_merge@meta.data[["majorType"]] <- ifelse(seurat_object_merge@meta.data[["SCT_snn_res.0.1"]] == cell_cluster, 
                                                       label_cluster_json[[i]], 
                                                       seurat_object_merge$majorType)
}
seurat_object_merge <- SetIdent(seurat_object_merge,
                    cells = rownames(seurat_object_merge@meta.data), 
                    value = seurat_object_merge@meta.data["majorType"])
