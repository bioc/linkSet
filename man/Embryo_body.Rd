\name{Embryo_body}
\alias{Embryo_body}
\docType{data}
\title{Embryo Body BED File}
\description{
  A BED format file containing genomic regions from mouse embryo body. This dataset contains 
  regulatory elements identified in mouse embryonic development and is provided as example data 
  for demonstrating genomic interaction analysis workflows.
}
\usage{
# Access the file path with:
system.file("extdata", "Embryo_body.bed", package = "linkSet")
}
\format{
  A BED format file with the following columns:
  \itemize{
    \item chromosome: The chromosome name (e.g., chr1, chr2)
    \item start: The starting position of the feature
    \item end: The ending position of the feature
    \item name: Name of the feature
    \item score: Score between 0 and 1000
    \item strand: Strand orientation (+ or -)
  }
  The file contains 3,727 genomic intervals.
}
\source{
  These data were derived from publicly available mouse embryonic development datasets, 
  specifically from the embryo body, and processed to identify regulatory elements.
  The original data were processed to create this example dataset for demonstration purposes.
}
\examples{
  # Get the file path
  file_path <- system.file("extdata", "Embryo_body.bed", package = "linkSet")
  
  # Read the file
  if (file.exists(file_path)) {
    embryo_data <- read.table(file_path, 
                             header = FALSE, 
                             sep = "\t",
                             stringsAsFactors = FALSE)
    head(embryo_data)
  }
}
\keyword{datasets}