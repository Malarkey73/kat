#' Reading data from a MAF or VCF type file.
#' This is a wrapper for read.delim that specifies the columns that will go to make a VCF data.frame/tbl
#' @param file a delimited file with columns containing at least "chr", "start.position", and "end.position" (a data.frame)
#' @param chr name of chromosome column (character, default "chr")
#' @param start.position name of start.position column (character, default "start.position")
#' @param end.position name of end position column (character, default "end.position")
#' @param strand name of starnd column(character, default "strand")
#' @param WT name of WT sequence column (character, default "WT")
#' @param MUT name of Mutant sequence column (character, default "MUT")
#' @param sampleID a column of sample IDs or names (character, default "sample")
#' @keywords rainfall, mutants, kataegesis
#' @examples
#' # a minimal example
#' read.mutations(file="1317049.tsv", chr="Chromosome", start.position="Genome.start", end.position="Genome.stop")
#'# a more fully specified VCF
#' VCF <- read.mutations(file="1317049.tsv", chr="Chromosome", start.position="Genome.start", end.position="Genome.stop", strand= "Strand", WT="WT.seq", MUT="Mut.seq", sampleID="Sample.Name")


read.mutations <- function(file, chr="chr", start.position="start.position", end.position="end.position", strand="strand", WT="WT", MUT="MUT", sampleID= "sampleID", skip=0, other=F, ...){

  temp = read.delim(file, as.is=T, ...)

  # check that at least these 3 columns exist - fail if not.
  if(all(c(chr, start.position, end.position) %in% colnames(temp)))
  {
    VCF <- data.frame(chr= temp[,chr], start.position= temp[,start.position], end.position= temp[,end.position])
  }
  else
  {
    stop(" A valid VCF needs at least chr, start.position, and end.position columns - check your file has these correctly named columns")
  }

  # using the dplyr tbl throughout this package
  library(dplyr)
  VCF=as.tbl(VCF)

  # add these columns if they exist
  if(strand %in% colnames(temp))
    VCF <- mutate(VCF, strand = temp[,strand])
  if(WT %in% colnames(temp))
    VCF <- mutate(VCF, WT = temp[,WT])
  if(MUT %in% colnames(temp))
    VCF <- mutate(VCF, MUT = temp[,MUT])
  if(sampleID %in% colnames(temp))
    VCF <- mutate(VCF, sampleID = temp[,sampleID])

  # if you want the other columns in the VCF then set "other = TRUE"
  if(other==TRUE)
  {
    other.cols<- setdiff(colnames(temp), c(chr,start.position,end.position, strand, WT, MUT, sampleID))
    VCF <- as.tbl(cbind(VCF, temp[,other]))
  }


  #short summary message and return a chr/pos sorted dplyr style tbl/data.frame
  cat(c("VCF contains: ", paste0(colnames(VCF), ", "), "\n\n"))
  VCF <- VCF %>%
    arrange(chr, start.position)

  return(VCF)

}
