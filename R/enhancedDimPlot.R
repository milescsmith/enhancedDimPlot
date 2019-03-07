#' @title enhancedDimPlot
#'
#' @description A customized version of Seurat's DimPlot. Used for
#' plotting dimensional reductions.
#'
#' @param object scRNAseq data object
#' @param reduction dimensional reduction to use. Default: tsne
#' @param group_by Metadata variable to use in grouping data. If none is provided,
#' the current ident will be used.  Default: NULL
#' @param split_by Metadata variable to use in faceting the plots. Default: NULL
#' @param group_plot If provided, only this identity group will be displayed. Default: NULL
#' @param dim_1 Dimension to display along the x-axis. Default: 1
#' @param dim_2 Dimension to display along the y-axis. Default: 2
#' @param pt_size Point size. Default: 1
#' @param alpha Point alpha. Default: 1 (i.e. fully opaque)
#' @param force Force parameter to pass to ggrepel.  See \code{ggrepel::\link[geom_label_repel]{geom_label_repel}}
#' @param label Should labels be shown? Default: TRUE
#' @param label_size Label font size. Default: 3
#' @param label_text_color Label font color. Default: black
#'
#' @importFrom dplyr enquo quo_name inner_join filter group_by summarise
#' @importFrom tibble rownames_to_column
#' @importFrom stats median
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap theme
#' @importFrom ggrepel geom_label_repel
#'
#' @return
#' @export
#'
#' @examples
#'
enhancedDimPlot <- function(object, ...){
  UseMethod('enhancedDimPlot')
}

#' @rdname enhancedDimPlot
#' @method enhancedDimPlot Seurat
#' @importFrom Seurat Embeddings FetchData
#' @export
#' @return
enhancedDimPlot.Seurat <- function(object,
                                   reduction = "tsne",
                                   group_by = NULL,
                                   group_plot = NULL,
                                   split_by = NULL,
                                   dim_1 = 1,
                                   dim_2 = 2,
                                   pt.size = 1,
                                   alpha = 1,
                                   force = 1,
                                   label = TRUE,
                                   label_size = 3,
                                   label_text_color = 'black'){
  try(
    if (is.null(group_by)){
      group_by <- "ident"
    }, silent = TRUE)
  try(
    if (is.character(group_by)) {
      group_by <- as.name(substitute(group_by))
    }, silent = TRUE
  )
  group_by <- enquo(group_by)

  varlist <- c(quo_name(group_by))

  try(
    if (is.character(split_by)) {
      split_by <- as.name(substitute(split_by))
    }, silent = TRUE
  )
  if(is.null(split_by)){
    faceting <- FALSE
  } else {
    faceting <- TRUE
    varlist <- c(varlist, quo_name(split_by))
  }
  split_by <- enquo(split_by)

  if (!reduction %in% names(object)) {
    stop(glue("{} coordinates were not found in {object}"))
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
    plot.data %<>% filter(!!group_by == group_plot)
  }

  plot.data[["x"]] <- plot.data[, dim_1]
  plot.data[["y"]] <- plot.data[, dim_2]

  centers <- plot.data %>%
    group_by(!!group_by) %>%
    summarise(x = median(x = x),
              y = median(x = y))

  p1 <- plot.data %>%
    ggplot(aes(x = x,
               y = y,
               color = !!group_by,
               fill = !!group_by,
               label = !!group_by)) +
    geom_point(size = pt.size,
               alpha = alpha,
               shape = 21)

  if (label){
    p2 <- p1 +
      geom_point(data = centers,
                 mapping = aes(x = x,
                               y = y),
                 size = 0,
                 alpha = 0) +
      geom_label_repel(data = centers,
                       mapping = aes(label = !!group_by,
                                     fill = !!group_by),
                       size = label_size,
                       color = label_text_color,
                       box.padding = 1,
                       force = force)
  } else {
    p2 <- p1
  }

  if (isTRUE(faceting)){
    p2 <- p2 + facet_wrap(quo_name(split_by))
  }
  p2 + theme(legend.position = "none")
}
