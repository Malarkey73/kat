#' Get surrounding sequence context
#' uses the Bioconductor BSGenome framework to find
#' @param VCF a VCF tbl with columns containing at least "chr", "start.position", and "end.position" 
#' @param start.position the start.positions column
#' @param expand The numberof upstream and downstream bases to get
#' @param genomeseq The Bioconductor genome annotation package to retireve sequence information
#' @keywords Poisson, mutants, kataegesis
#' @export
#' @examples
#' qqmutation(VCF, start.position= start.position, expand=1)

expandseq<- function(VCF, start.position= start.position, expand=1, genomeseq = BSgenome.Hsapiens.UCSC.hg19)
{
  require(BSgenome)
  #NSE
    start.position=deparse(substitute(start.position))
    genomeseq=deparse(substitute(genomeseq))
  
  
  # This is the format of the UCSC.hg19 data "chr19" NOT "19", "chrX" not "chr23"
  VCF <- VCF %>%
    mutate(chr = ifelse(chr %in% 1:24, yes=paste0("chr", chr), no=chr)) %>%
    mutate(chr = ifelse(chr=="chr23", yes="chrX", no=chr))

  # load the right library
  library(genomeseq, character.only=T)
  message("This may take a few seconds.")
  # if that library is human then use this to access the sequence data.
  if(grep("Hsapiens", genomeseq))
    gs <- getSeq(Hsapiens, VCF$chr, start=(VCF$start.position)-1, end=(VCF$start.position)+1)

  return(as.character(gs))
}
