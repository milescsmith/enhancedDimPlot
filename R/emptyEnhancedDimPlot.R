#' @title emptyEnhancedDimPlot
#'
#' @description Create a scatter plot where the points are not displayed but labels
#' for each group are.
#'
#' @param object scRNAseq data object
#' @param reduction dimensional reduction to use. Default: tsne
#' @param grouping_var Metadata variable to use in grouping data. If none is provided,
#' the current ident will be used.  Default: NULL
#' @param group_plot If provided, only this identity group will be displayed. Default: NULL
#' @param dim_1 Dimension to display along the x-axis. Default: 1
#' @param dim_2 Dimension to display along the y-axis. Default: 2
#' @param pt_size Point size. Default: 1
#' @param alpha Point alpha. Default: 1 (i.e. fully opaque)
#' @param force Force parameter to pass to ggrepel.  See \code{\link{geom_label_repel}}
#' @param label Should labels be shown? Default: TRUE
#' @param label_size Label font size. Default: 3
#' @param label_text_color Label font color. Default: black
#' @param ... Additional parameters
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr inner_join select filter group_by summarise
#' @importFrom stats median
#' @importFrom ggplot2 ggplot aes geom_point theme
#' @importFrom ggrepel geom_label_repel
#' @importFrom glue glue
#'
#' @return
#' @export
#'
#' @examples
#'
emptyEnhancedDimPlot <- function(object, ...){
  UseMethod("emptyEnhancedDimPlot")
}

#' @rdname emptyEnhancedDimPlot
#' @method emptyEnhancedDimPlot Seurat
#' @importFrom Seurat Embeddings
#' @export
#' @return
emptyEnhancedDimPlot.Seurat <- function(object,
                                        reduction = "tsne",
                                        grouping_var = NULL,
                                        group_plot = NULL,
                                        dim_1 = 1,
                                        dim_2 = 2,
                                        pt_size = 1,
                                        alpha = 1,
                                        force = 1,
                                        label = TRUE,
                                        label_size = 3,
                                        label_text_color = 'black',
                                        ...){
  try(
    if (is.null(grouping_var)){
      grouping_var <- "ident"
    }, silent = TRUE)
  try(
    if (is.character(grouping_var)) {
      grouping_var <- as.name(substitute(grouping_var))
    }, silent = TRUE
  )
  grouping_var <- enquo(grouping_var)
  varlist <- c(quo_name(grouping_var))

  if (!reduction %in% names(object)) {
    stop(glue("{reduction} coordinates were not found in object"))
  }
  dimData <- Embeddings(object = object,
                        reduction = reduction)

  metaData <- FetchData(object = object,
                        vars = varlist) %>%
    rownames_to_column('cell')

  dimNames <- colnames(dimData)

  dim_1 <- dimNames[[1]]
  dim_2 <- dimNames[[2]]

  plot.data <- dimData %>%
    as.data.frame() %>%
    rownames_to_column('cell') %>%
    inner_join(metaData,
               by = 'cell')

  if (!is.null(group_plot)){
    plot.data %<>% filter(!!grouping_var == group_plot)
  }

  plot.data[["x"]] <- plot.data[, dim_1]
  plot.data[["y"]] <- plot.data[, dim_2]

  centers <- plot.data %>%
    group_by(!!grouping_var) %>%
    summarise(x = median(x = x),
              y = median(x = y))

    p2 <- centers %>%
      ggplot(aes(x = x,
                 y = y),
             size = 0,
             alpha = 0) +
      geom_point(data = centers,
                 mapping = aes(x = x,
                               y = y),
                 size = 0,
                 alpha = 0) +
      geom_label_repel(data = centers,
                       mapping = aes(label = !!grouping_var,
                                     fill = !!grouping_var),
                       size = label_size,
                       color = label_text_color,
                       box.padding = 1,
                       force = force)

  p2 + theme(legend.position = "none")
}
