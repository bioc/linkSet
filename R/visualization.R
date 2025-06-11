#' @rdname geom_linkset
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame start end
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom utils read.table
#' @importFrom scales rescale
#' @importFrom ggplot2 ggplot_add ggplot aes_string scale_color_gradientn
#'   labs theme_classic theme element_blank element_rect
#'   element_text margin scale_y_continuous scale_x_continuous expansion
#'   coord_cartesian geom_curve
#' @importFrom patchwork wrap_plots
#' @export
#'
#'
setMethod("geom_linkset", "linkSet", function(linkSet,
                      score.col = "count",
                      score.threshold = NULL,
                      score.color = c("grey70", "#56B1F7" , "#132B43"),
                      scale.range = 10,
                      plot.space = 0.1,
                      plot.height = 0.2,
                      arrow.size = 0.05,
                      remove_x_axis = FALSE,
                      link_plot_on_top = FALSE,
                      extend.base = 1000000,
                      show.rect = FALSE,
                      x.range = NULL,
                      log.scale = TRUE) {
  structure(
    list(
      linkSet = linkSet,
      score.col = score.col,
      score.threshold = score.threshold,
      score.color = score.color,
      scale.range = scale.range,
      plot.space = plot.space,
      plot.height = plot.height,
      show.rect = show.rect,
      arrow.size = arrow.size,
      remove_x_axis = remove_x_axis,
      link_plot_on_top = link_plot_on_top,
      extend.base = extend.base,
      x.range = x.range,
      log.scale = log.scale
    ),
    class = "interSet"
  )
})



#' @export
ggplot_add.interSet <- function(object, plot, object_name) {
  # get plot data (handle cases where plot may not have layers)
  track.data <- NULL
  if ("patchwork" %in% class(plot)) {
    if (length(plot) > 0 && length(plot[[1]]$layers) > 0) {
      track.data <- plot[[1]]$layers[[1]]$data
    }
  } else {
    if (length(plot$layers) > 0) {
      track.data <- plot$layers[[1]]$data
    }
  }

  # get parameters
  linkSet <- object$linkSet
  score.col <- object$score.col
  score.threshold <- object$score.threshold
  score.color <- object$score.color
  scale.range <- object$scale.range
  plot.curve <- object$plot.curve
  plot.space <- object$plot.space
  plot.height <- object$plot.height
  show.rect <- object$show.rect
  arrow.size <- object$arrow.size
  remove_x_axis <- object$remove_x_axis
  link_plot_on_top <- object$link_plot_on_top
  flip_arrow <- link_plot_on_top
  top_margin <- bottom_margin <- plot.space
  extend.base <- object$extend.base
  x.range <- object$x.range
  log.scale <- object$log.scale


  # prepare plot range
  plot.range.chr <- as.character(seqnames(regionsBait(object$linkSet))[1])
  plot.range.start <- min(start(regions(object$linkSet))) - extend.base
  plot.range.end <- max(end(regions(object$linkSet))) + extend.base
  if (is.null(x.range)) {
    x.range <- c(plot.range.start, plot.range.end)
  }
  # prepare dataframe
  link.point.df <- data.frame(
    chr = as.character(seqnames(regionsBait(linkSet))),
    start = start(oe(linkSet)),
    end = start(regionsBait(linkSet))
  )

  # add score
  if (score.col %in% colnames(mcols(linkSet))) {
    link.point.df$score <- mcols(linkSet)[[score.col]]
    if (!is.null(score.threshold)) {
      link.point.df <- link.point.df[link.point.df$score > score.threshold, ]
    }
    if (log.scale) {
      link.point.df$score <- log1p(link.point.df$score)
    }
  }
  # filter link gr
  link.point.df <- link.point.df[link.point.df$start >= x.range[1] &
                                link.point.df$end >= x.range[1] &
                                link.point.df$end <= x.range[2]&
                                link.point.df$start <= x.range[2], ]
  if(nrow(link.point.df) == 0){
    warning("There are no valid links in the given region!")
    return(NULL)
  }
  rownames(link.point.df) <- 1:nrow(link.point.df)
  # check dataframe
  if (nrow(link.point.df) < 1) {
    warning("There are no valid links in the given region!")
    # create empty plot
    link.basic.plot <- ggplot2::ggplot(data = link.point.df)
  } else {
    # prepare plot dataframe
    link.point.df$group <- seq_len(length.out = nrow(link.point.df))
    link.point.plot <- link.point.df
    link.point.plot$width <- link.point.df$end - link.point.df$start
    #browser()
    # scale width to range
    link.point.plot$rw <- scales::rescale(link.point.plot$width, to = c(1, scale.range))

    if ("score" %in% colnames(link.point.plot)) {
      group_color <- "score"
      scale_color <- ggplot2::scale_color_gradientn(
        colors = score.color,
        limits = c(0, max(link.point.plot$score))
      )
    } else {
      group_color <- NULL
      scale_color <- ggplot2::scale_color_manual()
    }

    y_limit <- ifelse(flip_arrow, 0, 1)
    link.point.plot.pos <- link.point.plot[link.point.plot$width > 0,]
    link.point.plot.neg <- link.point.plot[link.point.plot$width < 0,]
    link.basic.plot <-
      ggplot2::ggplot(data = link.point.plot) +
      ggplot2::geom_curve(
        data = link.point.plot.pos,
        ggplot2::aes_string(
          x = "start",
          xend = "end",
          y = y_limit,
          yend = y_limit,
          color = group_color,
          size = "score"
        ),
        curvature = ifelse(flip_arrow, -0.2, 0.2),
        angle = 90,
        ncp = 15,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow.size, "npc"))
      ) +
      ggplot2::geom_curve(
        data = link.point.plot.neg,
        ggplot2::aes_string(
          x = "start",
          xend = "end",
          y = y_limit,
          yend = y_limit,
          color = group_color,
          size = "score"
        ),
        curvature = ifelse(flip_arrow, 0.2, -0.2),
        angle = 90,
        ncp = 15,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow.size, "npc"))
      ) +
      scale_color +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::scale_size_continuous(range = c(0.5, 2))
  }

  # create plot
  link.plot <-
    link.basic.plot +
    ggplot2::labs(y = "Links") +
    themeLinkset(
      x.range = x.range,
      margin.len = plot.space,
      show.rect = show.rect
    ) +
    ggplot2::guides(size = "none")  # Remove legend for arrow size

  # Add chromosome name to the side
  link.plot <- link.plot +
    ggplot2::annotate("text", x = min(x.range), y = 1,
                      label = plot.range.chr, hjust = 0, vjust = 1)

  # assemble plot
  patchwork::wrap_plots(
    plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
    link.plot,
    ncol = 1,
    heights = c(1, plot.height)
  )
  # Create a function to adjust plot margins and remove x-axis elements
  adjust_plot <- function(p, top_margin, bottom_margin, remove_x_axis = FALSE) {
    p <- p + ggplot2::theme(plot.margin = ggplot2::margin(t = top_margin, b = bottom_margin))
    if (remove_x_axis) {
      p <- p + ggplot2::theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    return(p)
  }

  # Adjust plots based on user preferences
  adjusted_link_plot <- adjust_plot(link.plot, top_margin = top_margin, bottom_margin = bottom_margin)
  adjusted_plot <- adjust_plot(plot, top_margin = top_margin, bottom_margin = bottom_margin, remove_x_axis = object$remove_x_axis)

  # Determine plot order and create list of plots
  plot_list <- if (object$link_plot_on_top) {
    list(adjusted_link_plot, adjusted_plot)
  } else {
    list(adjusted_plot, adjusted_link_plot)
  }

  # Calculate heights based on plot order and overlap
  total_height <- 1 + object$plot.height
  heights <- if (object$link_plot_on_top) {
    c(object$plot.height, total_height - object$plot.height)
  } else {
    c(total_height - object$plot.height, object$plot.height)
  }

  # Combine plots
  combined_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = 1,
    heights = heights
  ) +
    patchwork::plot_layout(guides = "collect")
  return(combined_plot)
}

#' Plot genomic ranges
#'
#' `geomRange()` and `geom_half_range()` draw tiles that are designed to
#' represent range-based genomic features, such as exons. In combination with
#' `geom_intron()`, these geoms form the core components for visualizing
#' transcript structures.
#'
#' @param mapping Set of aesthetic mappings created by [aes()]. If specified and
#'   `inherit.aes = TRUE` (the default), it is combined with the default mapping
#'   at the top level of the plot. You must supply `mapping` if there is no plot
#'   mapping.
#' @param data The data to be displayed in this layer. There are three
#'    options:
#'
#'    If `NULL`, the default, the data is inherited from the plot
#'    data as specified in the call to [ggplot()].
#'
#'    A `data.frame`, or other object, will override the plot
#'    data. All objects will be fortified to produce a data frame. See
#'    [fortify()] for which variables will be created.
#'
#'    A `function` will be called with a single argument,
#'    the plot data. The return value must be a `data.frame`, and
#'    will be used as the layer data. A `function` can be created
#'    from a `formula` (e.g. `~ head(.x, 10)`).
#' @param stat The statistical transformation to use on the data for this
#'    layer, as a string.
#' @param position Position adjustment, either as a string, or the result of
#'  a call to a position adjustment function.
#' @param ... Other arguments passed on to [layer()]. These are
#'   often aesthetics, used to set an aesthetic to a fixed value, like
#'   `colour = "red"` or `size = 3`. They may also be parameters
#'   to the paired geom/stat.
#' @param na.rm If `FALSE`, the default, missing values are removed with
#'   a warning. If `TRUE`, missing values are silently removed.
#' @param show.legend logical. Should this layer be included in the legends?
#'   `NA`, the default, includes if any aesthetics are mapped.
#'   `FALSE` never includes, and `TRUE` always includes.
#'   It can also be a named logical vector to finely select the aesthetics to
#'   display.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics,
#'   rather than combining with them. This is most useful for helper functions
#'   that define both data and aesthetics and shouldn't inherit behaviour from
#'   the default plot specification, e.g. [borders()].
#'
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @examples
#' library(ggplot2)
#'
#' # Create some example data
#' df <- data.frame(
#'   start = c(100, 200, 300, 400),
#'   end = c(150, 250, 350, 450),
#'   strand = c("+", "+", "-", "-"),
#'   gene = c("Gene1", "Gene1", "Gene2", "Gene2")
#' )
#'
#' # Basic usage
#' ggplot(df) +
#'   geomRange(aes(xmin = start, xmax = end, y = gene))
#'
geomRange <- function(mapping = NULL, data = NULL,
                    stat = "identity", position = "identity",
                    ...,
                    na.rm = FALSE,
                    show.legend = NA,
                    inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = ggplot2::GeomRect,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

#' Plot genomic ranges for linkSet objects
#'
#' This function visualizes the genomic interactions in a linkSet object,
#' showing the bait and other end regions as well as the links between them.
#'
#' @param linkset A linkSet object to plot
#' @param showBait Vector of bait regions to display (NULL for all)
#' @param showOE Vector of other end regions to display (NULL for all)
#' @param x.range Range of x-axis to display
#' @param score.col Column name for coloring links
#' @param show.rect Whether to show rectangles for regions
#' @param extend.base Base pairs to extend the plot
#' @param ... Additional arguments
#' @param bait_col Color for bait regions
#' @param oe_col Color for other end regions
#' @param default_col Default color
#' @param vjust Vertical justification
#' @param linejoin Line join style
#' @param na.rm Whether to remove NA values
#' @param minimal_width Minimal width for regions
#' @param show.legend Whether to show legend
#' @param inherit.aes Whether to inherit aesthetics
#' @param link_plot_on_top Whether to draw links on top
#' @param arrow.size Size of arrows
#' @param remove_x_axis Whether to remove x axis
#' @param plot.height Height of the plot
#' @param plot.space Space between plots
#' @param log.scale Whether to use log scale for colors
#'
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes_string geom_curve labs scale_color_gradientn scale_x_continuous scale_y_continuous theme element_text element_blank element_rect expansion margin
#' @importFrom scales rescale
#' @importFrom patchwork wrap_plots
#' @export
#' @aliases plotGenomicRanges,linkSet-method
#' @aliases plotGenomicRanges
#'
#' @examples
#' data(linkExample)
#' plotGenomicRanges(linkExample, extend.base = 10)
setMethod("plotGenomicRanges", "linkSet", function(linkset, showBait = NULL,
                                                  showOE = NULL,
                                                  x.range = NULL,
                                                  score.col = "count",
                                                  show.rect = TRUE,
                                                  extend.base = 10000,
                                                  ...,
                                                  bait_col = "red",
                                                  oe_col = "DeepSkyBlue3",
                                                  default_col = "grey",
                                                  vjust = NULL,
                                                  linejoin = "mitre",
                                                  na.rm = FALSE,
                                                  minimal_width = 0.01,
                                                  show.legend = NA,
                                                  inherit.aes = TRUE,
                                                  link_plot_on_top = FALSE,
                                                  arrow.size = 0.05, remove_x_axis = FALSE,
                                                  plot.height = 0.4, plot.space = 0.1,
                                                  log.scale = TRUE) {
  # Get the bait and oe plots
  plot_oe <- plot_bait <- NULL
  
  # Extract bait and OE regions
  bait_regions <- regionsBait(linkset)
  if (is.null(bait_regions)) {
    warning("No bait regions found. Creating mock bait regions from unique bait names.")
    unique_baits <- unique(bait(linkset))
    bait_regions <- GRanges(
      seqnames = "chr1",
      ranges = IRanges(start = seq_along(unique_baits) * 1000, width = 1000),
      names = unique_baits
    )
    names(bait_regions) <- unique_baits
  }
  bait_ids <- seq_along(bait_regions)
  
  # If showBait or showOE are provided, filter the linkSet
  if (!is.null(showBait)) {
    if (is.character(showBait)) {
      if (!all(showBait %in% names(bait_regions))) {
        stop("Specified showBait names not found in linkSet")
      }
      bait_ids <- which(names(bait_regions) %in% showBait)
    } else if (is.numeric(showBait)) {
      if (!all(showBait %in% seq_along(bait_regions))) {
        stop("Specified showBait indices out of range")
      }
      bait_ids <- showBait
    } else {
      stop("showBait must be a character or numeric vector")
    }
  }
  
  # Get the regions in these bait IDs
  bait_regions_subset <- bait_regions[bait_ids]
  
  # Function to create range plot
  create_range_plot <- function(gr, x.range = NULL, region_color, title = NULL) {
    # Calculate mid point of each range
    start_pos <- GenomicRanges::start(gr)
    end_pos <- GenomicRanges::end(gr)
    mid_pos <- (start_pos + end_pos) / 2
    
    # Determine x-axis range
    if (is.null(x.range)) {
      min_pos <- min(start_pos) - extend.base
      max_pos <- max(end_pos) + extend.base
    } else {
      min_pos <- x.range[1]
      max_pos <- x.range[2]
    }
    
    # Create data frame for plotting
    plot_data <- data.frame(
      chromosome = as.character(GenomicRanges::seqnames(gr)),
      start = start_pos,
      end = end_pos,
      name = if (is.null(names(gr))) paste0("region_", seq_along(gr)) else names(gr)
    )
    
    # Check if there's only one chromosome
    unique_chromosomes <- unique(plot_data$chromosome)
    if (length(unique_chromosomes) > 1) {
      warning("Multiple chromosomes detected. Using the first chromosome for plotting.")
      plot_data <- plot_data[plot_data$chromosome == unique_chromosomes[1], ]
    }
    
    # Create ggplot  
    p <- ggplot(plot_data, aes_string(xmin = "start", xmax = "end", ymin = -0.1, ymax = 0.1)) +
      geomRange(color = region_color, size = 3) +
      scale_x_continuous(limits = c(min_pos, max_pos), 
                         labels = function(x) paste0(x / 1000, "kb"),
                         expand = ggplot2::expansion(mult = 0.01)) +
      themeRange(x.range = c(min_pos, max_pos), show.rect = TRUE) +
      labs(title = title, x = "Position", y = "")
    
    if (remove_x_axis) {
      p <- p + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
    }
    
    return(p)
  }
  
  # Create bait plot
  plot_bait <- create_range_plot(bait_regions_subset, x.range, bait_col, "Bait Regions")
  
  # Process links
  if (!is.null(showBait)) {
    linkset_subset <- linkset[anchor1(linkset) %in% bait_ids]
  } else {
    linkset_subset <- linkset
  }
  
  # If showOE is provided, further filter
  if (!is.null(showOE)) {
    # Implement logic to filter by OE regions
    # This would depend on how OE regions are identified in your linkSet
  }
  
  # Get interaction data
  if (length(linkset_subset) == 0) {
    warning("No interactions to plot")
    return(NULL)
  }
  
  # Check if anchor indices are valid
  anchor1_indices <- anchor1(linkset_subset)
  anchor2_indices <- anchor2(linkset_subset)
  
  # For mock bait regions, map bait names to indices
  if (is.null(regionsBait(linkset))) {
    bait_names <- bait(linkset_subset)
    anchor1_indices <- match(bait_names, names(bait_regions))
    if (any(is.na(anchor1_indices))) {
      warning("Some bait names not found in bait regions")
      anchor1_indices[is.na(anchor1_indices)] <- 1
    }
  }
  
  # Validate anchor1_indices are within bounds
  valid_anchor1 <- anchor1_indices >= 1 & anchor1_indices <= length(bait_regions)
  if (!all(valid_anchor1)) {
    warning("Some anchor1 indices are out of bounds, using first bait region as fallback")
    anchor1_indices[!valid_anchor1] <- 1
  }
  
  # Validate anchor2_indices are within bounds
  valid_anchor2 <- anchor2_indices >= 1 & anchor2_indices <= length(regions(linkset_subset))
  if (!all(valid_anchor2)) {
    warning("Some anchor2 indices are out of bounds, filtering them out")
    # Keep only valid interactions
    valid_interactions <- valid_anchor1 & valid_anchor2
    if (sum(valid_interactions) == 0) {
      warning("No valid interactions found")
      return(NULL)
    }
    anchor1_indices <- anchor1_indices[valid_interactions]
    anchor2_indices <- anchor2_indices[valid_interactions]
    linkset_subset <- linkset_subset[valid_interactions]
  }
  
  anchor1_pos <- GenomicRanges::start(bait_regions[anchor1_indices]) + 
                  (GenomicRanges::end(bait_regions[anchor1_indices]) - 
                     GenomicRanges::start(bait_regions[anchor1_indices])) / 2
  
  anchor2_pos <- GenomicRanges::start(regions(linkset_subset)[anchor2_indices]) + 
                  (GenomicRanges::end(regions(linkset_subset)[anchor2_indices]) - 
                     GenomicRanges::start(regions(linkset_subset)[anchor2_indices])) / 2
  
  # Create link data frame
  link_data <- data.frame(
    x = anchor1_pos,
    y = rep(0, length(anchor1_pos)),
    xend = anchor2_pos,
    yend = rep(1, length(anchor2_pos))
  )
  
  # Add score column if available
  if (score.col %in% colnames(mcols(linkset_subset))) {
    link_data$score <- mcols(linkset_subset)[[score.col]]
    
    # Log transform if requested
    if (log.scale && any(link_data$score > 0)) {
      link_data$score[link_data$score > 0] <- log10(link_data$score[link_data$score > 0])
    }
  } else {
    link_data$score <- 1  # Default score
  }
  
  # Create link plot
  if (link_plot_on_top) {
    link_plot <- ggplot(link_data, aes_string(x = "x", y = "y", xend = "xend", yend = "yend", color = "score")) +
      geom_curve(curvature = 0.2, arrow = arrow(length = unit(arrow.size, "inches"), type = "closed")) +
      scale_color_gradientn(colors = c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")) +
      theme_void() +
      coord_cartesian(xlim = c(min(link_data$x, link_data$xend), max(link_data$x, link_data$xend)))
    
    # Combine plots with patchwork
    final_plot <- plot_bait + link_plot + 
      plot_layout(heights = c(plot.height, 1 - plot.height - plot.space))
  } else {
    # Logic for when links are not on top
    # You can implement this based on your requirements
    final_plot <- plot_bait
  }
  
  return(final_plot)
})

#' @rdname plotGenomicRanges
#' @export
setMethod("plot_genomic_ranges", "linkSet", function(linkset, showBait = NULL,
                                                    showOE = NULL,
                                                    x.range = NULL,
                                                    score.col = "count",
                                                    show.rect = TRUE,
                                                    extend.base = 10000,
                                                    ...,
                                                    bait_col = "red",
                                                    oe_col = "DeepSkyBlue3",
                                                    default_col = "grey",
                                                    vjust = NULL,
                                                    linejoin = "mitre",
                                                    na.rm = FALSE,
                                                    minimal_width = 0.01,
                                                    show.legend = NA,
                                                    inherit.aes = TRUE,
                                                    link_plot_on_top = FALSE,
                                                    arrow.size = 0.05, remove_x_axis = FALSE,
                                                    plot.height = 0.4, plot.space = 0.1,
                                                    log.scale = TRUE) {
  # Call the new function for backward compatibility
  plotGenomicRanges(linkset = linkset, 
                   showBait = showBait, 
                   showOE = showOE,
                   x.range = x.range,
                   score.col = score.col,
                   show.rect = show.rect,
                   extend.base = extend.base,
                   ...,
                   bait_col = bait_col,
                   oe_col = oe_col,
                   default_col = default_col,
                   vjust = vjust,
                   linejoin = linejoin,
                   na.rm = na.rm,
                   minimal_width = minimal_width,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   link_plot_on_top = link_plot_on_top,
                   arrow.size = arrow.size,
                   remove_x_axis = remove_x_axis,
                   plot.height = plot.height,
                   plot.space = plot.space,
                   log.scale = log.scale)
})

#' Extract data from linkSet for plotting
#'
#' @param linkset A linkSet object
#' @return A data.frame with extracted data
#' @keywords internal
#' @noRd
extractDataFromLinkset <- function(linkset) {
  # Extract regions
  regions_data <- as.data.frame(regions(linkset))
  regions_data$type <- "region"
  
  # Extract bait regions
  bait_regions_data <- as.data.frame(regionsBait(linkset))
  bait_regions_data$type <- "bait"
  
  # Combine all data
  all_data <- rbind(regions_data, bait_regions_data)
  
  # Add unique identifiers
  all_data$id <- 1:nrow(all_data)
  
  return(all_data)
}

#' Theme for linkSet plots
#'
#' @param x.range The x-axis range
#' @param margin.len Margin length
#' @param show.rect Whether to show rectangle
#' @return A ggplot2 theme
#' @export
themeLinkset <- function(x.range, margin.len, show.rect) {
  theme <- ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA),
      plot.background = ggplot2::element_rect(fill = NA, color = NA),
      panel.border = if (show.rect) {
        ggplot2::element_rect(fill = NA, color = "black")
      } else {
        ggplot2::element_blank()
      },
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = margin.len, r = margin.len, b = margin.len, l = margin.len)
    )
  
  return(theme)
}

#' Theme for genomic range plots
#'
#' @param x.range The x-axis range
#' @param show.rect Whether to show rectangle
#' @return A ggplot2 theme
#' @export
themeRange <- function(x.range, show.rect) {
  theme <- ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA),
      plot.background = ggplot2::element_rect(fill = NA, color = NA),
      panel.border = if (show.rect) {
        ggplot2::element_rect(fill = NA, color = "black")
      } else {
        ggplot2::element_blank()
      },
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
      axis.ticks.y = ggplot2::element_blank()
    )
  
  return(theme)
}

#' Plot baits in a linkSet object
#' @title Plot Baits
#' @importFrom rlang .data
#' @importFrom GenomicRanges strand
#' @name plotBaits
#' @param linkset A linkSet object
#' @param scoreCol Column name containing scores for coloring points
#' @param countCol Column name containing counts for y-axis values
#' @param n Number of random baits to plot if baits parameter is NULL
#' @param baits Vector of specific baits to plot. If NULL, n random baits are selected
#' @param plotBaitNames Logical indicating whether to show bait names in plot titles
#' @param plevel1 Upper threshold for score coloring (red)
#' @param plevel2 Lower threshold for score coloring (blue)
#' @param outfile Output file path. If NULL, plot is displayed rather than saved
#' @param width Width of output plot in inches
#' @param height Height of output plot in inches
#' @param extend.base Base pairs to extend view range on either side of bait
#' @param bgCol Color for points below plevel2 threshold
#' @param lev2Col Color for points between plevel2 and plevel1 thresholds
#' @param lev1Col Color for points above plevel1 threshold
#' @param ... Additional plotting parameters
#' @return A ggplot object
#' @export
plotBaits <- function(linkset, scoreCol = "score", countCol = "count", n = 4, baits = NULL, plotBaitNames = TRUE, 
                      plevel1 = 5, plevel2 = 3, outfile = NULL,
                      width = 20, height = 20, extend.base = 1e6, bgCol = "black", lev2Col = "blue", 
                      lev1Col = "red", ...) {

  if (!is(linkset, "linkSet")) {
    stop("Input must be a linkSet object")
  }
  if (is.null(baits)) {
    baits <- sample(bait(linkset), n)
  } else {
    n <- length(baits)
  }
  
  # Pre-compute color vector
  color_vector <- c(bgCol, lev2Col, lev1Col)
  
  plot_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    bait <- baits[i]
    this <- subsetBait(linkset, bait)
    
    if (!is.null(extend.base)) {
      baitGr <- regionsBait(this)
      baitGr<- unique(baitGr)
      new_start <- max(0, start(baitGr) - extend.base)
      new_end <- end(baitGr) + extend.base
      
      expandGr <- GRanges(
        seqnames = seqnames(baitGr),
        ranges = IRanges(start = new_start, end = new_end),
        strand = GenomicRanges::strand(baitGr)
      )
      this <- subsetOE(this, expandGr)
    }
    if (length(this) == 0) {
      warning(paste0("No interactions found for bait ", bait))
      next
    }

    plotDf <- as.data.frame(this)
    plotDf$oe_middle <- start(oe(this)) + (end(oe(this)) - start(oe(this))) / 2
    bait_middle <- start(baitGr) + (end(baitGr) - start(baitGr)) / 2
    
    # Compute color factor
    plotDf$color_factor <- cut(plotDf[[scoreCol]], 
                               breaks = c(-Inf, plevel2, plevel1, Inf), 
                               labels = c(1, 2, 3))
    
    title <- if (plotBaitNames) {
      baitName <- bait
      if (grepl(",", baitName)) {
        baitName <- sub("(\\S+,).+", "\\1...", baitName)
      }
      paste0(baitName, " (", as.character(bait), ")")
    } else {
      as.character(bait)
    }
    p <- ggplot2::ggplot(plotDf, 
                         ggplot2::aes(x = .data$oe_middle, y = .data[[countCol]], color = .data$color_factor)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = bait_middle, color = "grey", linetype = "dashed") +
      ggplot2::labs(title = title, x = "Distance from viewpoint", y = countCol) +
      ggplot2::theme_minimal() +
      ggplot2::scale_y_continuous(limits = c(0, NA)) +
      ggplot2::scale_color_manual(values = color_vector, guide = "none")
    
    plot_list[[i]] <- p
  }

  combined_plot <- patchwork::wrap_plots(plotlist = plot_list, ncol = n)

  if (!is.null(outfile)) {
    ggplot2::ggsave(outfile, combined_plot, width = width, height = height)
  } else {
    print(combined_plot)
  }
  invisible(baits)
}