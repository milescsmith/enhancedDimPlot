#' @title HighlightGroupDimPlot
#'
#' @description
#'
#' @param object scRNAseq data object
#' @param reduction dimensional reduction to use. Default: tsne
#' @param grouping_var Metadata variable to use in grouping data. If none is provided,
#' the current ident will be used.  Default: NULL
#' @param highlight_group Particular group to highlight. Default: NULL
#' @param highlight_color Color to highlight chosen group with. Default: red
#' @param contrast_color Color to use for all other groups. Default: grey
#' @param dim_1 Dimension to display along the x-axis. Default: 1
#' @param dim_2 Dimension to display along the y-axis. Default: 2
#' @param pt_size Point size. Default: 1
#' @param alpha Point alpha. Default: 1 (i.e. fully opaque)
#' @param force Force parameter to pass to ggrepel.  See \code{ggrepel::\link[geom_label_repel]{geom_label_repel}}
#' @param label Should labels be shown? Default: TRUE
#' @param label_size Label font size. Default: 3
#' @param label_text_color Label font color. Default: black
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter select inner_join group_by summarise
#' @importFrom stats median
#' @importFrom ggplot2 ggplot theme aes geom_point
#' @importFrom ggrepel geom_label_repel
#'
#' @return
#' @export
#'
#' @examples
HighlightGroupDimPlot <- function(object, ...){
  UseMethod("HighlightGroupDimPlot")
}

#' @rdname HighlightGroupDimPlot
#' @method HighlightGroupDimPlot Seurat
#' @importFrom Seurat Embeddings
#' @export
#' @return
HighlightGroupDimPlot.Seurat <- function(object,
                                         reduction = "tsne",
                                         grouping_var = NULL,
                                         highlight_group = NULL,
                                         highlight_color = "#E41A1C",
                                         contrast_color = "#AAAAAA",
                                         dim_1 = 1,
                                         dim_2 = 2,
                                         pt_size = 1,
                                         highlight_alpha = 1,
                                         contrast_alpha = 0.75,
                                         force = 1,
                                         label = TRUE,
                                         label_size = 3,
                                         label_text_color = 'black'){
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

  try(
    if (is.character(highlight_group)) {
      highlight_group <- as.name(substitute(highlight_group))
    }, silent = TRUE
  )
  if(is.null(highlight_group)){
    highlight <- FALSE
  } else {
    highlight <- TRUE
  }

  if (!reduction %in% names(object)) {
    stop(glue("{} coordinates were not found in {object}"))
  }
  dimData <- Embeddings(object = object,
                        reduction = reduction)

  dimNames <- colnames(dimData)

  dim_1 <- dimNames[[1]]
  dim_2 <- dimNames[[2]]

  plot.data <- dimData %>%
    as.data.frame() %>%
    rownames_to_column('cell') %>%
    inner_join(object@meta.data %>%
                 rownames_to_column('cell') %>%
                 select(cell,
                        !!grouping_var),
               by = 'cell')

  plot.data[["x"]] <- plot.data[, dim_1]
  plot.data[["y"]] <- plot.data[, dim_2]

  if (is.true(highlight)) {
    filtered.plot.data <- plot.data %>% filter(!!grouping_var == highlight_group)
    remaining.plot.data <- plot.data %>% filter(!(!!grouping_var) == highlight_group)
  }

  centers <- filtered.plot.data %>%
    group_by(!!grouping_var) %>%
    summarise(x = median(x = x),
              y = median(x = y))

  p1 <- remaining.plot.data %>%
    ggplot(aes(x = x,
               y = y)) +
    geom_point(size = pt_size,
               alpha = contrast_alpha,
               color = contrast_color) +
    geom_point(data = filtered.plot.data,
               mapping = aes(x = x,
                             y = y),
               color = highlight_color,
               size = pt_size,
               alpha = highlight_alpha)

  if (label){
    p2 <- p1 +
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
  } else {
    p2 <- p1
  }

  p2 + theme(legend.position = "none")
}
