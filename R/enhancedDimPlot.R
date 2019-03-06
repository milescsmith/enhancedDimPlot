enhancedDimPlot <- function(object, 
                            dim, 
                            group_by,
                            group.plot = NULL,
                            split_by = NULL,
                            dim.1 = 1, 
                            dim.2 = 2, 
                            pt.size = 1,
                            alpha = 1,
                            force = 1,
                            do.label = TRUE,
                            label.size = 3,
                            label.text.color = 'black')
{
  try(
    if (is.character(group_by)) {
      group_by <- as.name(substitute(group_by))
    }, silent = TRUE
  )
  group_by <- enquo(group_by)

  try(
    if (is.character(split_by)) {
      split_by <- as.name(substitute(split_by))
    }, silent = TRUE
  )
  if(is.null(split_by)){
    faceting <- FALSE
  } else {
    faceting <- TRUE
  }
  split_by <- enquo(split_by)
  
  dimData <- Embeddings(object = object,
                        reduction = dim
  )

  metaData <- FetchData(object = object, vars = c(quo_name(group_by), quo_name(split_by))) %>% rownames_to_column('cell')
  
  dimNames <- colnames(dimData)
  
  dim.1 <- dimNames[[1]]
  dim.2 <- dimNames[[2]]
  
  plot.data <- dimData %>%
    as.data.frame() %>%
    rownames_to_column('cell') %>%
    inner_join(metaData,
               by = 'cell')
  
  if (!is.null(group.plot)){
    plot.data %<>% filter(!!group_by == group.plot)
  }
  
  plot.data$x <- plot.data[, dim.1]
  plot.data$y <- plot.data[, dim.2]
  
  centers <- plot.data %>%
    dplyr::group_by(!!group_by) %>%
    dplyr::summarise(x = median(x = x), 
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
  
  if (do.label){
    p2 <- p1 +
      geom_point(data = centers, 
                 mapping = aes(x = x, 
                               y = y),
                 size = 0, 
                 alpha = 0) +
      geom_label_repel(data = centers, 
                       mapping = aes(label = !!group_by,
                                     fill = !!group_by), 
                       size = label.size, 
                       color = label.text.color,
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

DimPlotHighlightGroup <- function(object, 
                            dim, 
                            group_by,
                            highlight.group = NULL,
                            highlight.color = "#E41A1C",
                            contrast.color = "#AAAAAA",
                            dim.1 = 1, 
                            dim.2 = 2, 
                            pt.size = 1,
                            highlight.alpha = 1,
                            contrast.alpha = 0.75,
                            force = 1,
                            do.label = TRUE,
                            label.size = 3,
                            label.text.color = 'black')
{
  group_by <- enquo(group_by)
  
  dimData <- Embeddings(object = object,
                        reduction = dim
  )
  
  dimNames <- colnames(dimData)
  
  dim.1 <- dimNames[[1]]
  dim.2 <- dimNames[[2]]
  
  plot.data <- dimData %>%
    as.data.frame() %>%
    rownames_to_column('cell') %>%
    inner_join(object@meta.data %>% 
                 rownames_to_column('cell') %>%
                 select(cell, 
                        !!group_by), 
               by = 'cell')
  
  plot.data$x <- plot.data[, dim.1]
  plot.data$y <- plot.data[, dim.2]
  
  if (!is.null(highlight.group)) {
    filtered.plot.data <- plot.data %>% filter(!!group_by == highlight.group)
    remaining.plot.data <- plot.data %>% filter(!(!!group_by) == highlight.group)
  }
  
  centers <- filtered.plot.data %>%
    dplyr::group_by(!!group_by) %>%
    dplyr::summarise(x = median(x = x), 
              y = median(x = y))
  
  p1 <- remaining.plot.data %>%
    ggplot(aes(x = x,
               y = y)) +
    geom_point(size = pt.size,
               alpha = contrast.alpha,
               color = contrast.color) +
    geom_point(data = filtered.plot.data,
               mapping = aes(x = x,
                             y = y),
               color = highlight.color,
               size = pt.size,
               alpha = highlight.alpha)
  
  if (do.label){
    p2 <- p1 +
      geom_point(data = centers, 
                 mapping = aes(x = x, 
                               y = y),
                 size = 0, 
                 alpha = 0) +
      geom_label_repel(data = centers, 
                       mapping = aes(label = !!group_by,
                                     fill = !!group_by), 
                       size = label.size, 
                       color = label.text.color,
                       box.padding = 1,
                       force = force)
  } else {
    p2 <- p1
  }

  p2 + theme(legend.position = "none")
}

emptyEnhancedDimPlot <- function(object, 
                            dim, 
                            group_by,
                            group.plot = NULL,
                            dim.1 = 1, 
                            dim.2 = 2, 
                            pt.size = 1,
                            alpha = 1,
                            force = 1,
                            do.label = TRUE,
                            label.size = 3,
                            label.text.color = 'black')
{
  group_by <- enquo(group_by)
  
  dimData <- Embeddings(object = object,
                        reduction = dim
  )
  
  dimNames <- colnames(dimData)
  
  dim.1 <- dimNames[[1]]
  dim.2 <- dimNames[[2]]
  
  plot.data <- dimData %>%
    as.data.frame() %>%
    rownames_to_column('cell') %>%
    inner_join(object@meta.data %>% 
                 rownames_to_column('cell') %>%
                 select(cell, 
                        !!group_by), 
               by = 'cell')
  
  if (!is.null(group.plot)){
    plot.data %<>% filter(!!group_by == group.plot)
  }
  
  plot.data$x <- plot.data[, dim.1]
  plot.data$y <- plot.data[, dim.2]
  
  centers <- plot.data %>%
    dplyr::group_by(!!group_by) %>%
    dplyr::summarise(x = median(x = x), 
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
                       mapping = aes(label = !!group_by,
                                     fill = !!group_by), 
                       size = label.size, 
                       color = label.text.color,
                       box.padding = 1,
                       force = force)
  
  p2 + theme(legend.position = "none")
}

