#' Get surrounding sequence context
#' uses the Bioconductor BSGenome framework to find
#' @param VCF a VCF tbl with columns containing at least "chr", "start.position", and "end.position" (a data.frame)
#' @param expand The numberof upstream and downstream bases to get (numeric, default = 1)
#' @keywords Poisson, mutants, kataegesis
#' @export
#' @examples
#' qqmutation(VCF, start= "start.position", around=1)

expandseq<- function(VCF, start= start.position, expand=1)
{
  #NSE
  start.position=deparse(substitute(start.position))
  
  # This is the format of the UCSC.hg19 data "chr19" NOT "19", "chrX" not "chr23"
  VCF <- VCF %>%
    mutate(VCF, chr = ifelse(chr %in% 1:24, yes=paste0("chr", chr), no=chr)) %>%
    mutate(VCF, chr = ifelse(chr=="chr23", yes="chrX", no=chr))

library(BSgenome.Hsapiens.UCSC.hg19)
  gs <- getSeq(Hsapiens, VCF$chr, start=(VCF$start.position)-1, end=(VCF$start.position)+1)

  return(as.character(gs))
}
