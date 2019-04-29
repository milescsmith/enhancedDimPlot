#' @title enhancedFeaturePlot
#'
#' @description A customized version of Seurat's DimPlot. Used for
#' plotting dimensional reductions.
#'
#' @param object scRNAseq data object
#' @param reduction dimensional reduction to use. Default: tsne
#' @param feature Feature to use for coloring.  Anything FetchData can extract.
#' @param grouping_var Metadata variable to use in grouping data. If none is provided,
#' the current ident will be used.  Default: NULL
#' @param subset_labels A list of the groups in \code{grouping_var} to display.  By default,
#' all are displayed.  Default: NULL
#' @param group_plot If provided, only this identity group will be displayed. Default: NULL
#' @param split_by Metadata variable to use in faceting the plots. Default: NULL
#' @param lower.cutoff Remove expression values below this quantile. Expressed as decimal. Default: 0.05
#' @param upper.cutoff Remove expression values above this quantile. Expressed as decimal. Default: 0.95
#' @param dim_1 Dimension to display along the x-axis. Default: 1
#' @param dim_2 Dimension to display along the y-axis. Default: 2
#' @param pt_size Point size. Default: 1
#' @param alpha Point alpha. Default: 1 (i.e. fully opaque)
#' @param force Force parameter to pass to ggrepel.  See \code{\link{geom_label_repel}}
#' @param label Should labels be shown? Default: TRUE
#' @param label_size Label font size. Default: 3
#' @param label_text_color Label font color. Default: black
#' @param scale_colors  Colors to use for expression gradient.
#' Takes a list of \code{c(low, high)}.  Default: c("grey90","red")
#' @param fill_palette Palette to use for point fill.  By default, this matches \code{color_palette}
#' @param ... Additional parameters
#'
#' @importFrom dplyr enquo quo_name inner_join filter group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom stats median
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap theme
#' @importFrom ggrepel geom_label_repel
#' @importFrom glue glue
#'
#' @return
#' @export
#'
#' @examples
#'
enhancedFeaturePlot <- function(object, ...){
  UseMethod('enhancedFeaturePlot')
}

#' @rdname enhancedFeaturePlot
#' @method enhancedFeaturePlot Seurat
#' @importFrom Seurat Embeddings FetchData
#' @export
#' @return
enhancedFeaturePlot.Seurat <- function(object,
                                       feature,
                                       reduction = "tsne",
                                       grouping_var = NULL,
                                       subset_labels = NULL,
                                       group_plot = NULL,
                                       split_by = NULL,
                                       lower.cutoff = 0.05,
                                       upper.cutoff = 0.95,
                                       dim_1 = 1,
                                       dim_2 = 2,
                                       pt_size = 1,
                                       alpha = 0.5,
                                       force = 1,
                                       label = TRUE,
                                       label_size = 3,
                                       label_text_color = 'black',
                                       scale_colors = c("grey90","red"),
                                       fill_palette = NULL,
                                       ...){
  if(missing(feature)){
    stop("No feature provided.")
  }
  try(
    if (is.null(grouping_var)){
      grouping_var <- "ident"
    }, silent = TRUE)
  try(
    if (is.character(grouping_var)) {
      grouping_var <- as.name(substitute(grouping_var))
    }, silent = TRUE)
  try(
    if (is.character(feature)) {
      feature <- as.name(substitute(feature))
    }, silent = TRUE
  )
  grouping_var <- enquo(grouping_var)
  varlist <- c(quo_name(grouping_var))
  feature <- enquo(feature)

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
    stop(glue("{reduction} coordinates were not found in object"))
  }
  dimData <- Embeddings(object = object,
                        reduction = reduction)
  colnames(dimData)[[dim_1]] <- "x"
  colnames(dimData)[[dim_2]] <- "y"

  dimData %<>% as_tibble(rownames = "cell")

  metaData <- FetchData(object = object,
                        vars = c(quo_name(feature), unique(varlist))) %>%
    as_tibble(rownames = "cell")

  feature_levels <- metaData[[quo_name(feature)]]
  lowend <- quantile(feature_levels[feature_levels > 0], lower.cutoff)
  highend <- quantile(feature_levels[feature_levels > 0], upper.cutoff)
  feature_levels[feature_levels < lowend] <- lowend
  feature_levels[feature_levels > highend] <- highend
  metaData[[quo_name(feature)]] <- feature_levels

  plotData <- dimData %>%
    inner_join(metaData, by = 'cell')

  if (!is.null(group_plot)){
    plotData %<>% filter(!!grouping_var == group_plot)
  }

  centers <- plotData %>%
    group_by(!!grouping_var) %>%
    summarise(x = median(x = x),
              y = median(x = y))

  if(!is.null(subset_labels)){
    centers %<>% filter(!!grouping_var %in% subset_labels)
  }

  pl <- plotData %>%
    ggplot(aes(x = x,
               y = y,
               color = !!feature)) +
    geom_point(size = pt_size,
               alpha = alpha,
               shape = 16,
               stroke = 0)

  if (isTRUE(faceting)){
    pl <- pl + facet_wrap(quo_name(split_by))
  }

  if (label){
    if (isTRUE(faceting)){
      pl <- pl + geom_label_repel(data = centers,
                                  mapping = aes( x = x,
                                                 y = y,
                                                 label = !!grouping_var),
                                  fill = "grey75",
                                  size = label_size,
                                  color = label_text_color,
                                  box.padding = 1,
                                  force = force,
                                  inherit.aes = FALSE)
    } else {
      pl <- pl + geom_label_repel(data = centers,
                                  mapping = aes( x = x,
                                                 y = y,
                                                 label = !!grouping_var),
                                  fill = "grey75",
                                  size = label_size,
                                  color = label_text_color,
                                  box.padding = 1,
                                  force = force,
                                  inherit.aes = FALSE)
    }
  }

  pl + theme(legend.position = "none") +
    scale_color_gradient(low = scale_colors[[1]],
                         high = scale_colors[[2]]) +
    labs(title = quo_name(feature),
         x = glue("{reduction}_{dim_1}"),
         y = glue("{reduction}_{dim_2}"))
}
