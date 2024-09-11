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
#' result <- make_interactions_unique(ls)
#' result
#' 
setMethod("countInteractions", "linkSet", function(x, baitRegions = TRUE) {
  # Get the regions and bait names
  reg <- regions(x)
  bait_names <- bait(x)
  bait_regions <- regionsBait(x) %>% paste
  oe_regions <- oe(x) %>% paste
  
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