#' @title enhancedDimPlot
#'
#' @description A customized version of Seurat's DimPlot. Used for
#' plotting dimensional reductions.
#'
#' @param object scRNAseq data object
#' @param reduction dimensional reduction to use. Default: tsne
#' @param grouping_var Metadata variable to use in grouping data. If none is provided,
#' the current ident will be used.  Default: NULL
#' @param split_by Metadata variable to use in faceting the plots. Default: NULL
#' @param group_plot If provided, only this identity group will be displayed. Default: NULL
#' @param dim_1 Dimension to display along the x-axis. Default: 1
#' @param dim_2 Dimension to display along the y-axis. Default: 2
#' @param pt_size Point size. Default: 1
#' @param alpha Point alpha. Default: 1 (i.e. fully opaque)
#' @param force Force parameter to pass to ggrepel.  See \code{\link{geom_label_repel}}
#' @param label Should labels be shown? Default: TRUE
#' @param label_size Label font size. Default: 3
#' @param label_text_color Label font color. Default: black
#' @param color_palette Palette to use.  If \code{fill_palette} is provided, this controls
#' the color surrounding objects and the background of labels. Default: palette "bold_color" from package "rcartocolor"
#' @param fill_palette Palette to use for point fill.  By default, this matches \code{color_palette}
#' @param ... Additional parameters
#'
#' @importFrom dplyr enquo quo_name inner_join filter group_by summarise
#' @importFrom tibble rownames_to_column
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
                                   grouping_var = NULL,
                                   group_plot = NULL,
                                   split_by = NULL,
                                   dim_1 = 1,
                                   dim_2 = 2,
                                   pt_size = 1,
                                   alpha = 1,
                                   force = 1,
                                   label = TRUE,
                                   label_size = 3,
                                   label_text_color = 'black',
                                   color_palette = NULL,
                                   fill_palette = NULL,
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

  metaData <- FetchData(object = object,
                        vars = unique(varlist)) %>%
    rownames_to_column('cell')

  dimNames <- colnames(dimData)

  dim_1 <- dimNames[[dim_1]]
  dim_2 <- dimNames[[dim_2]]

  plotData <- dimData %>%
    as.data.frame() %>%
    rownames_to_column('cell') %>%
    inner_join(metaData,
               by = 'cell')

  if (!is.null(group_plot)){
    plotData %<>% filter(!!grouping_var == group_plot)
  }

  plotData[["x"]] <- plotData[, dim_1]
  plotData[["y"]] <- plotData[, dim_2]

  centers <- plotData %>%
    group_by(!!grouping_var) %>%
    summarise(x = median(x = x),
              y = median(x = y))

  pl <- plotData %>%
    ggplot(aes(x = x,
               y = y,
               color = !!grouping_var,
               fill = !!grouping_var)) +
    geom_point(size = pt_size,
               alpha = alpha,
               shape = 21)

  if (isTRUE(faceting)){
    pl <- pl + facet_wrap(quo_name(split_by))
  }

  if (label){
    pl <- pl +
      geom_point(data = centers,
                 mapping = aes(x = x,
                               y = y),
                 size = 0,
                 alpha = 0)
    if (isTRUE(faceting)){
      pl <- pl + geom_label_repel(data = centers,
                                  mapping = aes( x = x,
                                                 y = y,
                                                 label = !!grouping_var,
                                                 fill = !!grouping_var),

                                  size = label_size,
                                  color = label_text_color,
                                  box.padding = 1,
                                  force = force,
                                  inherit.aes = FALSE)
    } else {
      pl <- pl + geom_label_repel(data = centers,
                                  mapping = aes( x = x,
                                                 y = y,
                                                 label = !!grouping_var,
                                                 fill = !!grouping_var),
                                  size = label_size,
                                  color = label_text_color,
                                  box.padding = 1,
                                  force = force,
                                  inherit.aes = FALSE)
    }
  }

  if (!is.null(color_palette)){
    color_palette <- PrepPalette(palette_use = color_palette,
                                 bins = length(unique(plotData[[quo_name(grouping_var)]])))
    if (is.null(fill_palette)){
      fill_palette <- color_palette
    } else {
      color_palette <- PrepPalette(palette_use = fill_palette,
                                   bins = length(unique(plotData[[quo_name(grouping_var)]])))
    }
    pl <- pl + scale_color_manual(values = color_palette) + scale_fill_manual(values = color_palette)
  }

  pl + theme(legend.position = "none")
}
