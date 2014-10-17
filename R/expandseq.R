#' Get surrounding sequence context
#' uses the Bioconductor BSGenome framework to find
#' @param VCF a VCF tbl with columns containing at least "chr", "start.position", and "end.position" 
#' @param expand The numberof upstream and downstream bases to get
#' @param genomeseq The Bioconductor genome annotation package to retireve sequence information
#' @keywords Poisson, mutants, kataegesis
#' @export
#' @examples
#' # This will add a column of the triplet context surrounding point mutations
#' expandseq(VCF, expand=1)

expandseq<- function(VCF, expand=1, genomeseq = BSgenome.Hsapiens.UCSC.hg19)
{
  if (!("VCF" %in% class(VCF)))
    stop(" A valid VCF needs at least chr, start.position, and end.position columns - you can read data files using read.mutations() or convert your data into this format using as.VCF()")
  
  require(BSgenome)
  #NSE
    genomeseq=deparse(substitute(genomeseq))
  
  
  # This is the format of the UCSC.hg19 data "chr19" NOT "19", "chrX" not "chr23"
  VCF <- VCF %>%
    mutate(adj.chr = ifelse(chr %in% 1:24, yes=paste0("chr", chr), no=chr)) %>%
    mutate(adj.chr = ifelse(adj.chr=="chr23", yes="chrX", no=adj.chr))

  # load the right library
  library(genomeseq, character.only=T)
  message("This may take a few seconds.")
  # if that library is human then use this to access the sequence data.
  if(grep("Hsapiens", genomeseq))
    gs <- getSeq(Hsapiens, VCF$adj.chr, start=(VCF$start.position)-1, end=(VCF$start.position)+1)

  VCF<- VCF %>% 
    mutate(triplet = as.character(gs)) %>%
    select(-adj.chr)
  
  class(VCF)<- c(class(VCF), "VCF")
  return(VCF)
}
