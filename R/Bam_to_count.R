# library("Rsamtools")
# library("GenomicAlignments")
# library("GenomicRanges")
# library("GenomeInfoDb")
 
#' convert bam files to read count matrix with user defined bin size.
#' this function perform converting bam files to read count matrix
#' @param bam_files the full names of the bam files, inlcude the path and file names
#' @param bin_size the number of bases within a bin 
#' @param chromosomes the chromosomes kept in the output matrix
#' @return bins the bin information
#' @return counts the read count matrix, with bins being the rows and cells being the columns
#' @export

bam_to_counts <- function(bam_files, bin_size = 100000,
                          chromosomes = paste0("chr", c(1:22, "X", "Y"))) {
  
  # read BAM header
  header <- scanBamHeader(bam_files[1])[[1]]$targets
  
  # remove empty names
  header <- header[!is.na(names(header)) & names(header) != "" & 
                     ( names(header) %in% chromosomes | paste0("chr",names(header)) %in% chromosomes)]
  
  # create bins
  bins <- tileGenome(
    seqlengths = header,
    tilewidth = bin_size,
    cut.last.tile.in.chrom = TRUE
  )
  
  # initialize matrix
  count_matrix <- matrix(
    0,
    nrow = length(bins),
    ncol = length(bam_files)
  )
  
  colnames(count_matrix) <- basename(bam_files)
  
  rownames(count_matrix) <- paste0(
    seqnames(bins), ":", start(bins), "-", end(bins)
  )
  
  # count reads
  for (i in seq_along(bam_files)) {
    
    aln <- readGAlignments(bam_files[i])
    
    count_matrix[, i] <- countOverlaps(bins, aln)
  }
  
  return(list(
    bins = bins,
    counts = count_matrix
  ))
}



# bam_files <- list.files(
#   "/orange/feifeixiao/xuanxuan.yu/Bam_to_count/bam_files/",
#   pattern = ".bam$",
#   full.names = TRUE
# )[1]
# result <- bam_to_counts( bam_files, bin_size = 100000)
# count_matrix <- result$counts

