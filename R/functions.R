#' Infer ADAR activity with the dataset-specific regulon
#'
#' This function is used to get ADAR activity scores starting from an expression
#' matrix with genes as rows and samples as columns, and the gene regulation
#' network obtained by the expression matrix.
#'
#' @usage yourSignature(mat, net, organism = "")
#'
#' @param mat Normalized expression matrix used by ARACNe to compute the
#' network
#' @param net Network obtained from ARACNe
#' @param organism "hsapiens" or "mmusculus"
#'
#' @return Matrix with the activity scores of the genes in the reg_list
#'
#' @examples
#' data(example_network)
#' data(matrix)
#' res <- yourSignature(matrix, example_network, organism = "hsapiens")
#'
#' @import utils
#' @importFrom igraph graph_from_adjacency_matrix as_data_frame
#' @importFrom dplyr filter
#' @importFrom viper aracne2regulon viper
#'
#' @export

yourSignature <- function(mat, net, organism = "") {

  g  <- graph_from_adjacency_matrix(net, weighted = TRUE, mode = "undirected")
  dataf <- as_data_frame(g)

  if (organism == "hsapiens") {
    genes_selected <- c(regulators_list, "ADAR", "ADARB1")
  }
  else {genes_selected <- c(regulators_list_Mm, "Adar", "Adarb1")}

  df_new <- data.frame()
  for (g in genes_selected) {
    df1 <- filter(dataf, dataf[, 1] == g)
    df2 <- filter(dataf, dataf[, 2] == g)
    df2 <- df2[, c(2, 1, 3)]
    colnames(df2) <- c("from", "to", "weight")
    df3 <- rbind(df1, df2)
    df_new <- rbind(df_new, df3)
  }

  df_list <- split(df_new, df_new$from)
  df_list <- lapply(df_list, function(dataf) {dataf[, -1]})

  file_adj_path <- file.path(getwd(), "interactors_scores.adj")

  if (file.exists(file_adj_path)) {
    file.remove(file_adj_path)
  }

  i <- 1
  while (i < length(df_list)) {
    a <- c(names(df_list[i]), paste(df_list[[i]]$to, df_list[[i]]$weight))
    write.table(x = rbind(a), file = file_adj_path, row.names = FALSE,
                col.names = FALSE, quote = FALSE, append = TRUE, sep = " ")
    i <- i + 1
  }

  adj_file <- readLines(file_adj_path)
  modified_file <- gsub(" ", "\t", adj_file)
  writeLines(modified_file, file_adj_path)

  regulon <- aracne2regulon(file_adj_path, mat, verbose = FALSE)
  file.remove(file_adj_path)
  res <- viper(mat, regulon[grep("ADAR", ignore.case = TRUE, names(regulon))],
               verbose = FALSE, minsize = 0)

  res[which(tolower(rownames(res)) == "adar"), ] <-
    1 / (1 + exp(- res[which(tolower(rownames(res)) == "adar"), ]))

  res[which(tolower(rownames(res)) == "adarb1"), ] <-
    1 / (1 + exp(- res[which(tolower(rownames(res)) == "adarb1"), ]))

  res

}

#' Compute ADAR activity with the cancer-related signature
#'
#' The function returns ADAR activity scores for each sample/spot/cell depending
#' on the starting expression matrix given as input. It employs a signature
#' derived from the union of two datasets of A549 cell line samples, which is
#' ideal to measure ADAR activity in a cancer-related context.
#'
#' @usage tumoralSignature(data, object)
#'
#' @param data Either a numeric expression matrix with genes as rows and samples
#' as columns or an ExpressionSet from a \code{Seurat} object,
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} or
#' \code{\link[SpatialExperiment]{SpatialExperiment}}
#' @param object Character string indicating the type of object, either bulk or
#' eset
#'
#' @return Matrix with ADAR activity scores or eset object updated
#'
#' @examples
#' data(tumoral_regulon)
#' data(example_matrix)
#' res <- tumoralSignature(example_matrix, object = "bulk")
#'
#' @importFrom viper viper
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom GSEABase GeneSetCollection GeneSet
#' @importFrom SummarizedExperiment assay
#'
#' @export

tumoralSignature <- function(data, object = c("bulk", "eset")) {

  if (object == "bulk") {
    regulon <- tumoral_regulon
    res <- viper(data, regulon, verbose = FALSE, minsize = 0)
    res["ADAR", ] <- 1 / (1 + exp(- res["ADAR", ]))
    res["ADARB1", ] <- 1 / (1 + exp(- res["ADARB1", ]))

    return(res)
  }

  else {

    if (class(data) == "Seurat") {
      counts <- data@assays$RNA$data
    }
    else {counts <- assay(data, "logcounts")}

    tumoral_adar1_sign_genes <- names(tumoral_regulon$ADAR$tfmode)
    tumoral_adar1_gs <- GeneSet(tumoral_adar1_sign_genes,
                                setName = "tumoral_adar1")
    tumoral_adar2_sign_genes <- names(tumoral_regulon$ADARB1$tfmode)
    tumoral_adar2_gs <- GeneSet(tumoral_adar2_sign_genes,
                                setName = "tumoral_adar2")

    all_sets <- GeneSetCollection(list(tumoral_adar1_gs, tumoral_adar2_gs))
    ssgsea_obj <- ssgseaParam(counts, geneSets = all_sets)
    ssgsea <- gsva(ssgsea_obj)

    data[["adar1_activity"]] <- ssgsea[1, ]
    data[["adar2_activity"]] <- ssgsea[2, ]

    return(data)
  }

}

#' Compute ADAR activity with the neuronal signature
#'
#' The function returns ADAR activity scores for each sample/spot/cell depending
#' on the starting expression matrix given as input. It employs a signature
#' derived from a dataset of hESC h9 cell line samples, which is ideal to
#' measure ADAR activity in a human neuronal context.
#'
#' @usage neuronalSignature(data, object)
#'
#' @param data Either a numeric expression matrix with genes as rows and samples
#' as columns or an ExpressionSet from a \code{Seurat} object,
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} or
#' \code{\link[SpatialExperiment]{SpatialExperiment}}
#' @param object Character string indicating the type of object, either bulk or
#' eset
#'
#' @return Matrix with ADAR activity scores or eset object updated
#'
#' @examples
#' data(neuronal_regulon)
#' data(example_matrix)
#' res <- neuronalSignature(example_matrix, object = "bulk")
#'
#' @importFrom viper viper
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom GSEABase GeneSetCollection GeneSet
#' @importFrom SummarizedExperiment assay
#'
#' @export

neuronalSignature <- function(data, object = c("bulk", "eset")) {

  if (object == "bulk") {
    regulon <- neuronal_regulon
    res <- viper(data, regulon, verbose = FALSE, minsize = 0)
    res["ADAR", ] <- 1 / (1 + exp(- res["ADAR", ]))
    res["ADARB1", ] <- 1 / (1 + exp(- res["ADARB1", ]))
    res
  }

  else {

    if (class(data) == "Seurat") {
      counts <- data@assays$RNA$data
    }
    else {counts <- assay(data, "logcounts")}

    neuronal_adar1_sign_genes <- names(neuronal_regulon$ADAR$tfmode)
    neuronal_adar1_gs <- GeneSet(neuronal_adar1_sign_genes,
                                 setName = "neuronal_adar1")

    ssgsea_obj <- ssgseaParam(counts, geneSets = neuronal_adar1_gs)
    ssgsea <- gsva(ssgsea_obj)

    data[["adar1_activity"]] <- ssgsea[1, ]

    return(data)

  }

}

#' Compute ADAR activity with the mouse-specific neuronal signature
#'
#' The function returns ADAR activity scores for each sample/spot/cell depending
#' on the starting expression matrix given as input. It employs a signature
#' derived from a dataset of neurons and whole brain of mice, which is
#' ideal to measure ADAR activity in a mouse neuronal context.
#'
#' @usage neuronalMmSignature(data, object)
#'
#' @param data Either a numeric expression matrix with genes as rows and samples
#' as columns or an ExpressionSet from a \code{Seurat} object,
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} or
#' \code{\link[SpatialExperiment]{SpatialExperiment}}
#' @param object Character string indicating the type of object, either bulk or
#' eset
#'
#' @return Matrix with ADAR activity scores or eset object updated
#'
#' @examples
#' data(mouse_neuronal_regulon)
#' data(example_matrix_Mm)
#' res <- neuronalMmSignature(example_matrix_Mm, object = "bulk")
#'
#' @importFrom viper viper
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom GSEABase GeneSetCollection GeneSet
#' @importFrom SummarizedExperiment assay
#'
#' @export

neuronalMmSignature <- function(data, object = c("bulk", "eset")) {

  if (object == "bulk") {
    regulon <- mouse_neuronal_regulon
    res <- viper(data, regulon, verbose = FALSE, minsize = 0)
    res["Adar", ] <- 1 / (1 + exp(- res["Adar", ]))
    res["Adarb1", ] <- 1 / (1 + exp(- res["Adarb1", ]))
    res["Adarb2", ] <- 1 / (1 + exp(- res["Adarb2", ]))
    res
  }

  else {

    if (class(data) == "Seurat") {
      counts <- data@assays$RNA$data
    }
    else {counts <- assay(data, "logcounts")}

    Mm_adar1_sign_genes <- names(mouse_neuronal_regulon$Adar$tfmode)
    Mm_adar1_gs <- GeneSet(Mm_adar1_sign_genes,
                                setName = "neuronalMm_adar1")
    Mm_adar2_sign_genes <- names(mouse_neuronal_regulon$Adarb1$tfmode)
    Mm_adar2_gs <- GeneSet(Mm_adar2_sign_genes,
                                setName = "neuronalMm_adar2")

    all_sets <- GeneSetCollection(list(Mm_adar1_gs, Mm_adar2_gs))
    ssgsea_obj <- ssgseaParam(counts, geneSets = all_sets)
    ssgsea <- gsva(ssgsea_obj)

    data[["adar1_activity"]] <- ssgsea[1, ]
    data[["adar2_activity"]] <- ssgsea[2, ]

    return(data)

  }

}

#' Visualize ADAR CAS with boxplots
#'
#' The function generates a boxplot showing the distribution of ADAR CAS in the
#' cells or spots of a \code{Seurat} or \linkS4class{SummarizedExperiment}
#' object across user-specified categories.
#'
#' @usage CASboxplots(data, by = "", color = "", adar =  c("ADAR1", "ADAR2"))
#'
#' @param data Either a \code{Seurat} object,
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} or
#' \code{\link[SpatialExperiment]{SpatialExperiment}}
#' @param by Character string indicating the name of the categorical variable
#' present in the metadata to be used to divide the samples into the boxplots
#' (e.g. cluster, condition, group)
#' @param color Character string indicating the name of the variable in the
#' metadata to use for coloring boxplots
#' @param adar Character string specifying the CAS to visualize, between ADAR1
#' and ADAR2
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'
#' @examples
#' data(example_sce)
#' CASboxplots(example_sce, by = "Group", color = "Group", adar = "ADAR1")
#'
#' @import ggplot2
#'
#' @export

CASboxplots <- function (data, by = "", color = "",
                         adar =  c("ADAR1", "ADAR2")) {

  if (class(data) == "Seurat") {
    meta <- data@meta.data
  } else {
    meta <- data@colData
  }

  if (adar == "ADAR1") {
    cas <- meta$adar1_activity
  } else {
    cas <- meta$adar2_activity
  }

  gg_df <- data.frame(adar_cas = cas,
                      groups = meta[[by]],
                      color = meta[[color]])

  ggplot(gg_df, aes(x = "groups", y = "adar_cas", fill = "color")) +
    geom_boxplot() +
    labs(title = NULL, x = NULL, y = "CAS") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90),
          legend.position = "none")

}

#' Dimensional visualization of ADAR CAS
#'
#' The function shows the distribution of ADAR CAS on a dimensional reduction
#' plot using \code{dittoDimPlot}.
#'
#' @usage CASdittoplot(data, reduction = "", adar =  c("ADAR1", "ADAR2"))
#'
#' @param data Either a \code{Seurat} object,
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} or
#' \code{\link[SpatialExperiment]{SpatialExperiment}}
#' @param reduction Character string indicating the name of a dimensionality
#' reduction (e.g. umap, pca) slot within the object
#' @param adar Character string specifying the CAS to visualize, between ADAR1
#' and ADAR2
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'
#' @examples
#' data(example_sce)
#' CASdittoplot(example_sce, reduction = "PCA", adar = "ADAR1")
#'
#' @import ggplot2
#' @importFrom dittoSeq dittoDimPlot
#' @importFrom scales rescale
#'
#' @export

CASdittoplot <- function(data, reduction = "", adar = c("ADAR1", "ADAR2")) {

  if (adar == "ADAR1") {
    var <- "adar1_activity"
  } else {
    var <- "adar2_activity"
  }

  dittoDimPlot(data, var = var, reduction.use = reduction) +
    theme(legend.position = "bottom") +
    scale_color_gradientn(colors = c("#2D1160FF", "#721F81FF", "white",
                                     "#ED7953FF", "#FDB32FFF"),
                          values = rescale(seq(from = min(data[[var]]),
                                               to = max(data[[var]]),
                                               length.out = 5)),
                          limits = c(min(data[[var]]), max(data[[var]])),
                          breaks = c(min(data[[var]]), max(data[[var]])),
                          labels = c("low activity", "high activity")) +
    labs(color = "CAS", title = NULL)

}
