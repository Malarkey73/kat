#' Count of dinucleotide mutations (e.g. CC>AA, TT>CC)#'
#' This function counts the different types of dinucleotide substitutions (e.g. CC to AA). As with the simple pointmut.count function it also counts the reverse complement... so CCtoAA is really the count of both CCtoAA and GGtoTT. Most NGS technologies (not all) are not stranded so you cannot tell if a CCtoAA on positive strand was caused by a GGtoTT on the negative starnd (or vice versa).
#' @param VCF a delimited file with columns "WT" and "MUT", containing values including "A", "T", "G" or "C" (a data.frame)
#' @keywords monte carlo, mutants, substitutions, point mutations
#' @examples
#' res <- dinucmut.count(VCF)

dinucmut.count<- function(VCF)
{
  require(dplyr)
  
  if (!("VCF" %in% class(VCF)))
    stop(" A valid VCF needs at least chr, start.position, end.position, WT and MUT columns - you can read data files using read.mutations() or convert your data into this format using as.VCF()")
  
  # filter down to second substitutions (dinucleotides, dn) same as the one before
  dn <- VCF %>%
    arrange(chr,start.position) %>%
    mutate(im.lag = start.position-lag(start.position)) %>%
    filter(im.lag ==1 & WT == lag(WT) & MUT ==lag(MUT)) 
    
  # this function can be used alone but here it counts the dinucs - crucially ignoring any indels         
  if(nrow(dn)>0)
    {
    res<- pointmut.count(dn)
    }else
    {
      res= rep(0,6)
    }

  names(res)=c("CCtoAA", "CCtoGG", "CCtoTT", "TTtoAA", "TTtoCC", "TTtoGG")
  return(res)

}
