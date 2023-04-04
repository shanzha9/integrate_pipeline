multiline_comments({
  "1. Check R packages"
  "2. Check species"
  "3. Check whether the 'keys in 'meta_data_key' are defined in each sample"
  "4. Check whether the 'sample path' is Ture"
})

# check packages
check_packages <- function(params) {
  package_list <- params$reqirement
  to_be_install_list <- c()
  for (i in package_list) {
    if (!requireNamespace(i, quietly = TRUE)) {
      to_be_install_list <- c(to_be_install_list, i)
    }
  }
  if (length(to_be_install_list) > 0) {
    cat_time(paste0(length(to_be_install_list), " R packages are needed, please install. They are:\n", sep = ""), level = "info")
    print(to_be_install_list)
    stop("R packages are needed!!!")
  } else {
    cat_time("R packages check pass!!!")
  }
}

# check parameters
check_parameters <- function(params) {
  # filting method check
  assert_that(parameters$quality_control_method$standard_deviation$use != parameters$quality_control_method$extremum$use, msg = "Only one method of <quality control method> can be used!!! Please check!!!")

  # data path check
  assert_that((parameters$project_baseline$cellranger_data$use != parameters$project_baseline$filted_data_list$use) ||
    (parameters$project_baseline$cellranger_data$use != parameters$project_baseline$filted_data_merge$use) ||
    (parameters$project_baseline$filted_data_merge$use != parameters$project_baseline$filted_data_merge$use), msg = "Only one 'project_baseline' can be used!!! Please check!!!\n")

  # normalization methods check
  assert_that(parameters$integration$normalization$sctransform$use != parameters$integration$normalization$lognorlization$use, msg = "Only one normalization method of <integration> can be used!!! Please check!!!")

  # Harmony "groupby" and "theta"
  assert_that(length(parameters$integration$method$harmony$groupby) == length(parameters$integration$method$harmony$lamda), msg = "The length of Harmony <groupby> and <lamda> should be same!!! Please check!!!")

  if (parameters$project_baseline$cellranger_data$use) {
    meta_keys <- parameters$data$meta_data_key
    for (i in 1:length(parameters$data)) {
      if (names(parameters$data[i]) != "meta_data_key") {
        # check keys
        if (any(!meta_keys %in% names(parameters$data[[i]]))) {
          missing <- which(!meta_keys %in% names(parameters$data[[i]]))
          missing <- meta_keys[missing]
          warning(sprintf("Could not find %s in <%s>\n", missing, names(parameters$data[i])))
          return(FALSE)
        }

        # check sample path
        path <- parameters$data[[i]]$path
        if (!file.exists(parameters$data[[i]]$path)) {
          warning(sprintf("The count file path of <%s> is wrong!!! Please check!!!\n", names(parameters$data[i])))
        }
      }
    }
  } else if (parameters$project_baseline$filted_data_list$use) {
    if (!file.exists(parameters$project_baseline$filted_data_list$data_path)) {
      warning(sprintf("The 'filted_data_list' path of is wrong!!! Please check!!!\n"))
      return(FALSE)
    }
  } else if (parameters$project_baseline$filted_data_merge$use) {
    if (!file.exists(parameters$project_baseline$filted_data_merge$data_path)) {
      warning(sprintf("The 'filted_data_merge' path of is wrong!!! Please check!!!\n"))
      return(FALSE)
    }
  } else {
    warning(sprintf("The 'project_baseline' data path of is wrong!!! Please check!!!\n"))
  }

  return(TRUE)
}
