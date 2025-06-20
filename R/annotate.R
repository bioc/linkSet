# Annotation--------------
# Create a new environment to store database connections
dbCache <- new.env(parent = emptyenv())

# Helper function to get or create database connection
getDbConnection <- function(genome) {
  # Create a key for the connection
  conn_key <- paste0("conn_", genome)

  # Check if connection exists and is valid
  if (exists(conn_key, envir = dbCache)) {
    src <- get(conn_key, envir = dbCache)
    if (DBI::dbIsValid(src$con)) {
      return(src)
    } else {
      # Remove invalid connection
      rm(list = conn_key, envir = dbCache)
    }
  }

  # Create new connection
  src <- if (genome == "hg38") {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      stop("Package TxDb.Hsapiens.UCSC.hg38.knownGene needed for this function")
    }
    Organism.dplyr::src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
  } else if (genome == "hg19") {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
      stop("Package TxDb.Hsapiens.UCSC.hg19.knownGene needed for this function")
    }
    Organism.dplyr::src_organism("TxDb.Hsapiens.UCSC.hg19.knownGene")
  } else if (genome == "mm10") {
    if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
      stop("Package TxDb.Mmusculus.UCSC.mm10.knownGene needed for this function")
    }
    Organism.dplyr::src_organism("TxDb.Mmusculus.UCSC.mm10.knownGene")
  } else {
    stop("Unsupported genome. Please use 'hg38', 'hg19' or 'mm10'.")
  }

  # Store the connection
  assign(conn_key, src, envir = dbCache)

  return(src)
}

# Function to cleanup connections
cleanupConnections <- function() {
  for (conn_key in ls(dbCache)) {
    src <- get(conn_key, envir = dbCache)
    if (inherits(src, "src_dbi") && DBI::dbIsValid(src$con)) {
      tryCatch({
        DBI::dbDisconnect(src$con)
      }, error = function(e) {
        warning(paste("Failed to disconnect from database:", e$message))
      })
    }
  }
  rm(list = ls(dbCache), envir = dbCache)
}


#' Annotate the link set with txDb. Give a gene list, and return a
#'
#' @aliases annotatePromoter
#' @param x linkSet
#' @param genome the genome you want to annotate
#' @param keyType the key type. You can check with AnnotationDbi::keytypes
#' @param upstream The upstream base from the gene
#' @param overwrite Whether to overwrite the regionsBait if it already exists
#'
#' @return linkSet object
#' @export
#'
#' @examples
#'   gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                 ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
#'                 strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
#'   gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                 ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
#'                 strand = "+")
#'   linkset_obj <- linkSet(gr1, gr2, specificCol = "symbol")
#'
#'   # Test annotatePromoter
#'   annotated_linkset <- suppressWarnings(annotatePromoter(linkset_obj, genome = "hg38", upstream = 500,overwrite = TRUE))
#'
#'
# Modified annotatePromoter method
setMethod("annotatePromoter", "linkSet", function(x, genome = "hg38",
                                                  keyType = "symbol", upstream = 5000,
                                                  overwrite = FALSE) {
  if (!is.null(regionsBait(x)) && !overwrite) {
    warning("regionsBait already exists, set overwrite = TRUE to overwrite it")
    return(x)
  }
  if (!is.null(regionsBait(x)) && overwrite) {
    warning("regionsBait is being overwritten")
  }

  tryCatch({
    # Get cached or new connection
    src <- getDbConnection(genome)

    genes <- bait(x)
    geneGr <- Organism.dplyr::genes(src, filter = ~(symbol %in% genes))
    promoterGr <- IRanges::promoters(geneGr, upstream = upstream)

    index <- match(genes, geneGr$symbol)
    if (any(is.na(index))) {
      warning("Some genes are not found in the txDb, they will be set to chrNULL")
      newIndex <- which(!is.na(index))
      grMatch <- promoterGr[index[newIndex]]
      gr <- GRanges(seqnames = rep("chrNULL", length(genes)),
                    ranges = IRanges(rep(0, length(genes)), rep(0, length(genes))))
      mcols(grMatch) <- NULL
      gr[newIndex] <- grMatch
    } else {
      gr <- promoterGr[index]
    }

    regionsBait(x) <- gr
    return(x)
  }, error = function(e) {
    warning(paste("An error occurred:", e$message))
    return(x)
  })
})

# Register cleanup function to be called when R exits
reg.finalizer(dbCache, function(e) {
  cleanupConnections()
}, onexit = TRUE)

# Create a namespace environment to avoid conflicts
.linkset_env <- new.env(parent = emptyenv())

# Custom genes function that handles both Organism.dplyr and fallback cases
.linkset_env$genes <- function(x, columns = NULL, filter = NULL, ...) {
  # If this is called from withTxDb context and Organism.dplyr failed
  if (inherits(x, "mock_organism_src")) {
    # This is our fallback case
    if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE) ||
        !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Required packages not available for mm10 annotation")
    }
    
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    gene_ranges <- GenomicFeatures::genes(txdb)
    
    if ("symbol" %in% columns) {
      gene_ids <- names(gene_ranges)
      symbols <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, 
                                     keys = gene_ids, 
                                     columns = "SYMBOL", 
                                     keytype = "ENTREZID")
      symbol_map <- setNames(symbols$SYMBOL, symbols$ENTREZID)
      mcols(gene_ranges)$symbol <- symbol_map[names(gene_ranges)]
      gene_ranges <- gene_ranges[!is.na(mcols(gene_ranges)$symbol)]
    }
    
    return(gene_ranges)
  } else {
    # Try to use Organism.dplyr if available
    if (requireNamespace("Organism.dplyr", quietly = TRUE)) {
      return(Organism.dplyr::genes(x, columns = columns, filter = filter, ...))
    } else {
      stop("Organism.dplyr package not available")
    }
  }
}

#' @rdname withTxDb
#' @param x Character string specifying the genome ("hg38", "hg19", or "mm10")
#' @param expr Function to execute with database connection
#' @param ... Additional arguments passed to expr
#' @importFrom methods setMethod
#' @export
setMethod("withTxDb", signature(x = "character", expr = "function"),
  function(x, expr, ...) {
    if (!x %in% c("mm10", "hg38", "hg19")) {
      stop("Unsupported genome. Please use 'hg38', 'hg19' or 'mm10'.")
    }
    
    # Use the original approach with proper error handling
    tryCatch({
      src <- getDbConnection(x)
      result <- expr(src, ...)
      return(result)
    }, error = function(e) {
      # If Organism.dplyr fails, try fallback approach
      warning(paste("Organism.dplyr approach failed:", e$message, ". Using fallback."))
      
      if (x == "mm10") {
        # Create a fallback implementation that works directly with the expected call
        if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE) ||
            !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
          stop("Required packages not available for mm10 annotation")
        }
        
        # Create mock source object with its own genes method
        mock_src <- list()
        class(mock_src) <- "mock_organism_src"
        
        # Create custom environment that intercepts Organism.dplyr::genes calls
        custom_env <- new.env(parent = globalenv())
        
        # Create mock Organism.dplyr namespace
        mock_organism_dplyr <- new.env()
        mock_organism_dplyr$genes <- function(src, columns = NULL, filter = NULL, ...) {
          txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
          gene_ranges <- GenomicFeatures::genes(txdb)
          
          if ("symbol" %in% columns) {
            gene_ids <- names(gene_ranges)
            symbols <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, 
                                           keys = gene_ids, 
                                           columns = "SYMBOL", 
                                           keytype = "ENTREZID")
            symbol_map <- setNames(symbols$SYMBOL, symbols$ENTREZID)
            mcols(gene_ranges)$symbol <- symbol_map[names(gene_ranges)]
            gene_ranges <- gene_ranges[!is.na(mcols(gene_ranges)$symbol)]
          }
          
          return(gene_ranges)
        }
        
        # Assign the mock to the custom environment
        assign("Organism.dplyr", mock_organism_dplyr, envir = custom_env)
        
        # Set the expression's environment to our custom one
        environment(expr) <- custom_env
        
        result <- expr(mock_src, ...)
        return(result)
      } else {
        stop("Fallback not implemented for genome: ", x)
      }
    })
  }
)
