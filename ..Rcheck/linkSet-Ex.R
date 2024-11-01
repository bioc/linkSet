pkgname <- "linkSet"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('linkSet')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("annotateInter-linkSet-method")
### * annotateInter-linkSet-method

flush(stderr()); flush(stdout())

### Name: annotateInter,linkSet-method
### Title: Annotate linkSet with inter/intra chromosome interactions
### Aliases: annotateInter,linkSet-method annotateInter

### ** Examples

data(linkExample)
linkExample <- annotateInter(linkExample)



cleanEx()
nameEx("countInteractibility-linkSet-method")
### * countInteractibility-linkSet-method

flush(stderr()); flush(stdout())

### Name: countInteractibility,linkSet-method
### Title: Count bait and oe
### Aliases: countInteractibility,linkSet-method countBaitOe

### ** Examples

data(linkExample)
linkSet = c(linkExample,linkExample)
linkSet = countInteractions(linkSet)
linkSet = countInteractibility(linkSet)



cleanEx()
nameEx("countInteractions-linkSet-method")
### * countInteractions-linkSet-method

flush(stderr()); flush(stdout())

### Name: countInteractions,linkSet-method
### Title: Make Repeated Interactions Unique and Count Repetitions
### Aliases: countInteractions,linkSet-method

### ** Examples

data(linkExample)
linkSet = c(linkExample,linkExample)
result <- countInteractions(linkSet)
result




cleanEx()
nameEx("crossGeneEnhancer-linkSet-method")
### * crossGeneEnhancer-linkSet-method

flush(stderr()); flush(stdout())

### Name: crossGeneEnhancer,linkSet-method
### Title: Cross gene enhancer
### Aliases: crossGeneEnhancer,linkSet-method crossGeneEnhancer

### ** Examples

data(linkExample)
linkSet = c(linkExample,linkExample)
linkSet = countInteractions(linkSet)
linkSet = filterLinks(linkSet, filter_intra = FALSE, filter_unannotate = FALSE, distance = 100000)
linkSet = crossGeneEnhancer(linkSet, score_threshold = 10)



cleanEx()
nameEx("dot-dbCache")
### * dot-dbCache

flush(stderr()); flush(stdout())

### Name: .dbCache
### Title: Annotate the link set with txDb. Give a gene list, and return a
### Aliases: .dbCache annotatePromoter
### Keywords: datasets

### ** Examples

  gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
                ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
                strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
  gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
                ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
                strand = "+")
  ls <- linkSet(gr1, gr2, specificCol = "symbol")

  # Test annotatePromoter
  annotated_ls <- suppressWarnings(annotatePromoter(ls, genome = "hg38", upstream = 500,overwrite = TRUE))





cleanEx()
nameEx("filterLinks-linkSet-method")
### * filterLinks-linkSet-method

flush(stderr()); flush(stdout())

### Name: filterLinks,linkSet-method
### Title: Filter links for further analysis
### Aliases: filterLinks,linkSet-method filterLinks

### ** Examples

data(linkExample)
linkSet = c(linkExample,linkExample)
linkSet = countInteractions(linkSet)
linkSet = filterLinks(linkSet, filter_intra = FALSE, filter_unannotate = FALSE, distance = 100000)



cleanEx()
nameEx("linkSet-GRange-Methods")
### * linkSet-GRange-Methods

flush(stderr()); flush(stdout())

### Name: trim,linkSet-method
### Title: linkSet-GRange-Methods
### Aliases: trim,linkSet-method trim resize,linkSet-method resize
###   resizeRegions,linkSet-method resizeRegions narrow,linkSet-method
###   narrow narrowRegions,linkSet-method narrowRegions
###   shift,linkSet-method shift shiftRegions,linkSet-method shiftRegions
###   flank,linkSet-method flank flankRegions,linkSet-method flankRegions
###   promoters,linkSet-method promoters promoterRegions,linkSet-method
###   promoterRegions width,linkSet-method width reduce,linkSet-method
###   reduce

### ** Examples

data(linkExample)
resize_bait <- resizeRegions(linkExample, width = 75, fix = "start", region = "bait")
resize_bait

narrow_bait <- narrowRegions(linkExample, start = 1, width = 5, region = "bait")
narrow_bait

shift_oe <- shiftRegions(linkExample, shift = 10, region = "oe")
shift_oe

flank_bait <- flankRegions(linkExample, width = 100, start = TRUE, both = FALSE, use.names = TRUE, ignore.strand = FALSE, region = "bait")
flank_bait

width(linkExample)




cleanEx()
nameEx("linkSet-class")
### * linkSet-class

flush(stderr()); flush(stdout())

### Name: linkSet-class
### Title: LinkSet object
### Aliases: linkSet-class linkSet LinkSet

### ** Examples

showClass("linkSet")  # shows the known subclasses

set.seed(7000)
N <- 40
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))
genes = c(rep("SP7",4),rep("ASPN",10),rep("XBP1",6))
Np <- 20
all.anchor1 <- sample(N, Np)
gr1 <- all.regions[all.anchor1]
gr1$symbol <- genes
all.anchor2 <- setdiff(1:40,all.anchor1)
gr2 <- all.regions[all.anchor2]
x <- linkSet(gr1, gr2,specificCol = "symbol")
x
x2 <- linkSet(genes, gr2)
x2




cleanEx()
nameEx("linkSet")
### * linkSet

flush(stderr()); flush(stdout())

### Name: linkSet
### Title: linkSet: Base Classes for Storing Genomic Link Data
### Aliases: linkSet

### ** Examples

# Basic usage example
library(linkSet)
# Add a simple example here




cleanEx()
nameEx("orderLinks-linkSet-method")
### * orderLinks-linkSet-method

flush(stderr()); flush(stdout())

### Name: orderLinks,linkSet-method
### Title: Order linkSet by mcols
### Aliases: orderLinks,linkSet-method orderLinks

### ** Examples

data(linkExample)
linkSet = c(linkExample,linkExample)
linkSet = countInteractions(linkSet)
linkSet = filterLinks(linkSet, filter_intra = FALSE, filter_unannotate = FALSE, distance = 100000)
linkSet = orderLinks(linkSet, by = "count", decreasing = TRUE)



cleanEx()
nameEx("pairdist-linkSet-method")
### * pairdist-linkSet-method

flush(stderr()); flush(stdout())

### Name: pairdist,linkSet-method
### Title: Calculate the distance between bait and the other end
### Aliases: pairdist,linkSet-method pairdist

### ** Examples

data(linkExample)
linkExample <- pairdist(linkExample, type="mid")




cleanEx()
nameEx("plot_genomic_ranges-linkSet-method")
### * plot_genomic_ranges-linkSet-method

flush(stderr()); flush(stdout())

### Name: plot_genomic_ranges,linkSet-method
### Title: plot genomic ranges and links
### Aliases: plot_genomic_ranges,linkSet-method plot_genomic_ranges

### ** Examples

data(linkExample)
plot_genomic_ranges(linkExample, extend.base = 10)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
