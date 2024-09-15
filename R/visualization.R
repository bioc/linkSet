
#' Add Genome Links to Coverage Plot.
#'
#' @param link.file File contains region link information.
#' @param file.type The type of \code{link.file}, choose from bedpe, pairs. Default: bedpe.
#' @param score.col Column index that contains score information, used when \code{file.type} is bedpe. Default: NULL.
#' @param score.threshold The score threshold, used when \code{score.col} is not NULL. Default: NULL.
#' @param score.color The score color vector. Default: c("grey70", "#56B1F7", "#132B43").
#' @param scale.range Scale the height of links according to width, should be greater than or equal to 1 (not scale). Default: 10.
#' @param plot.curve One of 'curve' or 'bezier', for the latter it is required to install package \code{ggforce}. Default: 'curve'.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of link to coverage plot. Default: 0.2.
#' @param show.rect Logical value, whether to add rect border to the plot. Default: FALSE.
#'
#' @return Plot.
#' @importFrom dplyr %>%
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame start end
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom utils read.table
#' @importFrom scales rescale
#' @importFrom ggplot2 ggplot_add ggplot aes_string scale_color_gradientn
#'   labs theme_classic theme element_blank element_rect
#'   element_text margin scale_y_continuous scale_x_continuous expansion
#'   coord_cartesian geom_curve
#' @importFrom patchwork wrap_plots
#' @references \url{https://stuartlab.org/signac/articles/cicero.html}
#' @export
#'
#' @examples
#' library(ggcoverage)
#' # create random test data
#' # use seed to obtain same result every time
#' set.seed(123)
#'
#' df <- data.frame(
#'   seqnames = "chr2L",
#'   start = seq(from = 8000000, to = 8300000, by = 1000),
#'   end = seq(from = 8001000, to = 8301000, by = 1000),
#'   score = sample(1:100, 301, replace = TRUE),
#'   Type = "Example", Group = "Example"
#' )
#' # get links
#' link.file <- system.file(
#'   "extdata", "HiC", "HiC_link.bedpe",
#'   package = "ggcoverage"
#' )
#'
#' # create plot
#' ggcoverage(
#'   data = df, color = "grey",
#'   mark.region = NULL, range.position = "out"
#' ) +
#'   geom_link(link.file = link.file, file.type = "bedpe", show.rect = TRUE)
#'
#' 
setMethod("geom_linkset", "linkSet", function(linkSet,
                      score.col = "count",
                      score.threshold = NULL,
                      score.color = c("grey70", "#56B1F7", "#132B43"),
                      scale.range = 10,
                      plot.space = 0.1,
                      plot.height = 0.2,
                      arrow.size = 0.2,
                      remove_x_axis = FALSE,
                      link_plot_on_top = FALSE,
                      extend.base = 10000,
                      show.rect = FALSE) {
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
      extend.base = extend.base
    ),
    class = "interSet"
  )
})



#' @export
ggplot_add.interSet <- function(object, plot, object_name) {
  # get plot data
  if ("patchwork" %in% class(plot)) {
    track.data <- plot[[1]]$layers[[1]]$data
  } else {
    track.data <- plot$layers[[1]]$data
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


  # prepare plot range
  plot.range.chr <- as.character(seqnames(regionsBait(object$linkSet))[1])
  plot.range.start <- min(start(regions(object$linkSet))) - extend.base
  plot.range.end <- max(end(regions(object$linkSet))) + extend.base

  # prepare dataframe
  link.point.df <- data.frame(
    chr = as.character(seqnames(regionsBait(linkSet))),
    start = start(regionsBait(linkSet)),
    end = start(oe(linkSet))
  )

  # add score
  if (score.col %in% colnames(mcols(linkSet))) {
    link.point.df$score <- mcols(linkSet)[[score.col]]
    if (!is.null(score.threshold)) {
      link.point.df <- link.point.df[link.point.df$score > score.threshold, ]
    }
  }

  # filter link gr
  link.point.df <- link.point.df[link.point.df$start >= plot.range.start &
                                   link.point.df$end <= plot.range.end, ]
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
    # scale width to range
    link.point.plot$rw <- scales::rescale(link.point.plot$width, to = c(1, scale.range))

    if ("score" %in% colnames(link.point.plot)) {
      group_color <- "score"
      scale_color <- ggplot2::scale_color_gradientn(
        colors = score.color,
        limits = range(link.point.plot$score)
      )
    } else {
      group_color <- NULL
      scale_color <- ggplot2::scale_color_manual()
    }

    #scale_y_limit <- ifelse(flip_arrow, c(0,1), c(0,0))
    y_limit <- ifelse(flip_arrow, 0,1)
    link.basic.plot <-
      ggplot2::ggplot(data = link.point.plot) +
      ggplot2::geom_curve(
        ggplot2::aes_string(
          x = "start",
          xend = "end",
          y = y_limit,
          yend = y_limit,
          color = group_color
        ),
        curvature = ifelse(flip_arrow, -0.2, 0.2),
        angle = 90,
        ncp = 15,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow.size, "npc"))
      ) +
      scale_color+
      ggplot2::scale_y_continuous(limits = c(0,1)) 
  }
  # create plot
  link.plot <-
    link.basic.plot +
    ggplot2::labs(y = "Links") +
    theme_linkset(
      x.range = c(plot.range.start, plot.range.end),
      margin.len = plot.space,
      show.rect = show.rect
    )
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
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")

  return(combined_plot)
}




#' Plot genomic ranges
#'
#' `geom_range()` and `geom_half_range()` draw tiles that are designed to
#' represent range-based genomic features, such as exons. In combination with
#' `geom_intron()`, these geoms form the core components for visualizing
#' transcript structures.
#'
#' `geom_range()` and `geom_half_range()` require the following `aes()`;
#' `xstart`, `xend` and `y` (e.g. transcript name). `geom_half_range()` takes
#' advantage of the vertical symmetry of transcript annotation by plotting only
#' half of a range on the top or bottom of a transcript structure. This can be
#' useful for comparing between two transcripts or free up plotting space for
#' other transcript annotations (e.g. `geom_junction()`).
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @inheritParams ggplot2::geom_tile
#' @inheritParams ggplot2::geom_segment
#' @inheritParams grid::rectGrob
#'
#' @return the return value of a `geom_*` function is not intended to be
#'   directly handled by users. Therefore, `geom_*` functions should never be
#'   executed in isolation, rather used in combination with a
#'   `ggplot2::ggplot()` call.
#'
#' @export
geom_range <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       bait_col = "red",
                       oe_col = "DeepSkyBlue3",
                       default_col = "grey",
                       ...,
                       vjust = NULL,
                       linejoin = "mitre",
                       na.rm = FALSE,
                       minimal_width = 0.02,
                       show.legend = NA,
                       inherit.aes = TRUE) {
    ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomRange,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            minimal_width = minimal_width,
            vjust = vjust,
            linejoin = linejoin,
            na.rm = na.rm,
            bait_col = bait_col,
            oe_col = oe_col,
            default_col = default_col,
            ...
        )
    )
}

#' `GeomRange` is `ggplot2::GeomTile` with modified `aes` to match genetic
#' nomenclature (`xstart`/`xend`)
#' @keywords internal
#' @noRd
GeomRange <- ggplot2::ggproto("GeomRange", ggplot2::GeomTile,
    required_aes = c("xstart", "xend","region"),
    default_aes = ggplot2::aes(
        fill = "grey",
        colour = "black",
        linewidth = 0.25,
        linetype = 1,
        alpha = NA,
        height = NA
    ),
    setup_data = function(data, params) {
        # modified from ggplot2::GeomTile

        data$height <- data$height %||% params$height %||% 0.5
        data$fill <- with(data, case_when(
            region == "bait" ~ params$bait_col,
            region == "oe" ~ params$oe_col,
            TRUE ~ params$default_col  # Default color if not bait or oe
        ))

        transform(
            data,
            xmin = xstart,
            xmax = xend,
            ymin =  - data$height / 200,
            ymax =  data$height / 200,
            height = NULL
        )
    },
    draw_panel = function(self,
                          data,
                          panel_params,
                          coord,
                          vjust = NULL,
                          minimal_width = 0.02,
                          bait_col = "red",
                          oe_col = "DeepSkyBlue3",
                          default_col = "grey",
                          lineend = "butt",
                          linejoin = "mitre") {
        if (!coord$is_linear()) {
            # prefer to match geom_curve and warn
            # rather than copy the implementation from GeomRect for simplicity
            # also don'think geom_range would be used for non-linear coords
            warn("geom_ is not implemented for non-linear coordinates")
        }

        coords <- coord$transform(data, panel_params)
        grid::rectGrob(
            coords$xmin, coords$ymax,
            width = max(coords$xmax - coords$xmin, minimal_width),
            height = coords$ymax - coords$ymin,
            default.units = "native",
            just = c("left", "top"),
            vjust = vjust,
            gp = grid::gpar(
                col = coords$colour,
                fill = ggplot2::alpha(coords$fill, coords$alpha),
                lwd = coords$linewidth * ggplot2::.pt,
                lty = coords$linetype,
                linejoin = linejoin,
                lineend = lineend
            )
        )
    }
)

# plot genomic ranges
#' @export
plot_genomic_ranges <- function(linkset, x.range = NULL, show.rect = TRUE,extend.base = 10000,
                                ...,
                                bait_col = "red",
                                oe_col = "DeepSkyBlue3",
                                default_col = "grey",
                                vjust = NULL,
                                linejoin = "mitre",
                                na.rm = FALSE,
                                minimal_width = 0.02,
                                show.legend = NA,
                                inherit.aes = TRUE) {
    # Extract data from linkset object
    if (is.null(x.range)) {
        plot.range.start <- min(start(regions(linkset))) - extend.base
        plot.range.end <- max(end(regions(linkset))) + extend.base
        x.range <- c(plot.range.start, plot.range.end)
    }
    data <- extract_data_from_linkset(linkset)
    
    
    # Create the base plot
    p <- ggplot2::ggplot(data, aes(xstart =  xstart, xend = xend, region = region)) +
        geom_range(
            minimal_width = minimal_width,
            bait_col = bait_col,
            oe_col = oe_col,
            default_col = default_col,
            ...
        )

    # Apply the theme_linkset

    p <- p + ggplot2::labs(y = "Ranges") +theme_range(x.range, show.rect)

    return(p)
}

#' Extract data from linkSet object
#' @keywords internal
#' @noRd
extract_data_from_linkset <- function(linkset) {
    # Extract data from linkset object

    region_bait <- regionsBait(linkset)  %>% as.data.frame()
    region_oe <- oe(linkset) %>% as.data.frame()

    region_bait$region <- "bait"
    region_oe$region <- "oe"
    regionDf <- rbind(region_bait,region_oe)
    regionDf <- unique(regionDf)
    regionDf <- regionDf[,c("start","end","region")]
    colnames(regionDf) <- c("xstart","xend","region")

    return(regionDf)
} 



#' @export
theme_linkset <- function(x.range, margin.len, show.rect) {
  if (show.rect) {
    list(
      ggplot2::theme_classic(),
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.margin = margin(t = margin.len, b = margin.len)
      ),
      #ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0.5, 0)), position = "right"),
      ggplot2::scale_x_continuous(expand = c(0, 0)),
      ggplot2::coord_cartesian(xlim = x.range)
    )
  } else {
    list(
      ggplot2::theme_classic(),
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = ggplot2::element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = margin.len, b = margin.len)
      ),
      ggplot2::scale_x_continuous(expand = c(0, 0)),
      ggplot2::coord_cartesian(xlim = x.range)
    )
  }
}

#' @export
theme_range <- function(x.range, show.rect) {
  if (show.rect) {
    list(
      ggplot2::theme_classic(),
      ggplot2::theme(
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        #plot.margin = margin(t = margin.len, b = margin.len)
      ),
      #ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0.5, 0)), position = "right"),
      ggplot2::scale_x_continuous(expand = c(0, 0)),
      ggplot2::coord_cartesian(xlim = x.range)
    )
  } else {
    list(
      ggplot2::theme_classic(),
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = ggplot2::element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        #plot.margin = ggplot2::margin(t = margin.len, b = margin.len)
      ),
      ggplot2::scale_x_continuous(expand = c(0, 0)),
      ggplot2::coord_cartesian(xlim = x.range)
    )
  }
}