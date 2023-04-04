# violin plot for QC
make_violin_plot <- function(object, file, features, pt=0.1) {
  
  # ggsave(file, VlnPlot(object, features = features, ncol = floor(length(features) / 2)))
  
  png(file)
  p1 <- VlnPlot(object, features = features, ncol = floor(length(features) / 2))
  print(p1)
  invisible(dev.off())
}

# scatter plot for QC
make_scatter_plot <- function(object, file) {
  
  # ggsave(file, VlnPlot(object, features = features, ncol = floor(length(features) / 2)))
  
  png(file)
  p1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p1 + p2)
  invisible(dev.off())
}

# Function to make PC variation plots based on PC standard deviation
make_pc_plots <- function(object, file, jackstraw = NULL) {
  sdev <- object@reductions$pca@stdev
  var <- sdev^2
  cum_var <- cumsum(var)/sum(var)
  df <- data.frame(PC = 1:length(cum_var), var = var, cum_var = cum_var)
  
  pdf(file)
  
  # Make cumulative variation plot
  p1 <- ggplot(df, aes(x = PC, y = cum_var)) + geom_point() + expand_limits(y = 0)
  p1 <- p1 + labs(x = "PC", y = "Cumulative Proportion of Variation", title = "PC Cumulative Proportion of Variation")
  print(p1)
  
  # Make variation plot
  p2 <- ggplot(df, aes(x = PC, y = var)) + geom_point() + expand_limits(y = 0)
  p2 <- p2 + labs(x = "PC", y = "Variance", title = "PC Variance")
  print(p2)
  
  # Make JackStraw plot
  if (!is.null(jackstraw))
  {
    print(JackStrawPlot(object,dims=1:jackstraw))
  }
  
  invisible(dev.off())
}

marker_gene_dim_plot <- function(object, file, features) {
  pdf(file)
  p1 <- DimPlot(object, reduction = "umap", label = T, pt.size =0.2)
  print(p1)
  
  p2 <- VlnPlot(object, 
                features = features,
                pt.size = 0, stack = T, sort = F, flip = T) + NoLegend() +
    scale_x_discrete(limits = sort(unique(object@active.ident)))
  print(p2)
  
  p3 <- DotPlot(object, features = features) + coord_flip()
  print(p3)
  
  invisible(dev.off())

}

marker_gene_feature_plot <- function(object, file, features) {
  pdf(file)
  p1 <- FeaturePlot(object, features)
  print(p1)
  
  
  invisible(dev.off())
  
}