#' Make Repeated Interactions Unique and Count Repetitions
#'
#' This function takes a linkSet object, identifies repeated interactions
#' (where both bait and other end are repeated), makes them unique, and
#' counts the number of repetitions for each unique interaction.
#'
#' @param x A linkSet object
#'
#' @return A list containing:
#'   \item{unique_linkSet}{A new linkSet object with unique interactions}
#'   \item{interaction_counts}{A data frame with counts for each unique interaction}
#'
#' @export
#'
#' @examples
#' ls = data(linkExample)
#' ls = c(ls,ls)
#' result <- countInteractions(ls)
#' result
#'
setMethod("countInteractions", "linkSet", function(x, baitRegions = TRUE) {
  # Get the regions and bait names
  reg <- regions(x)
  bait_names <- bait(x)
  bait_regions <- paste(regionsBait(x))
  oe_regions <- paste(oe(x))

  if (baitRegions) {
    interaction_ids <- paste(bait_names, bait_regions, oe_regions, sep = "___")
  } else {
    interaction_ids <- paste(bait_names, oe_regions, sep = "___")
  }

  # Count occurrences of each unique interaction
  interaction_counts <- table(interaction_ids)

  # Identify unique interactions
  unique_indices <- !duplicated(interaction_ids)

  # Create new linkSet with unique interactions
  unique_linkSet <- x[unique_indices]
  unique_linkSet <- clean_unused_regions(unique_linkSet)
  # Add interaction counts as metadata
  indiceName <- interaction_ids[unique_indices]
  mcols(unique_linkSet)$count <- as.vector(interaction_counts[indiceName])

  # Return the results
  return(unique_linkSet)
})


#' Filter links for further analysis
#' @export
setMethod("filterLinks", "linkSet", function(x, filter_intra = TRUE,
                                            filter_unannotate = TRUE,
                                            distance = NULL) {
  # Ensure inter_type is annotated
  if (!.exist_inter(x)) {
    x <- annotateInter(x)
  }

  # Ensure distance is calculated
  if (is.null(x$distance)) {
    x <- pairdist(x)
  }

  # Initialize a logical vector for filtering
  keep <- rep(TRUE, length(x))

  if (filter_intra) {
    keep <- keep & (x$inter_type == "inter")
  }

  if (filter_unannotate) {
    bait_regions <- regionsBait(x)
    keep <- keep & (as.character(seqnames(bait_regions)) != "chrNULL")
  }

  if (!is.null(distance)) {
    keep <- keep & (!is.na(x$distance) & x$distance <= distance)
  }

  # Apply the filter
  filtered_x <- x[keep]

  # Clean unused regions
  filtered_x <- clean_unused_regions(filtered_x)

  return(filtered_x)
})

#' Cross gene enhancer
#' @export
setMethod("crossGeneEnhancer", "linkSet", function(x, score_threshold = NULL) {
  if (!is.null(score_threshold) & !"score" %in% colnames(mcols(x))) {
    warning("score column not found.")
    score_threshold = NULL
  }
  else if (!is.null(score_threshold)) {
    x <- x[mcols(x)$score >= score_threshold]
  }

  # Get bait and oe regions
  baitName <- bait(x)
  oeName <- paste(oe(x))

  baitDf <- as.data.frame(table(oeName, baitName))
  baitDf <- baitDf[baitDf$Freq > 0, ]
  enhancerDf <- as.data.frame(table(baitDf$oeName))
  
  # Assign the Freq column to the mcol of x
  mcols(x)$crossFreq <- enhancerDf$Freq[match(oeName, enhancerDf[,1])]

  return(x)
})


#' Order linkSet by mcols
#' @export
setMethod("orderLinks", "linkSet", function(x, by = "count", decreasing = TRUE) {
  if (!by %in% colnames(mcols(x))) {
    stop("Column '", by, "' not found in mcols of linkSet object.")
  }
  order_indices <- order(mcols(x)[[by]], decreasing = decreasing)
  return(x[order_indices])
})
