#' Count of dinucleotide mutations (e.g. CC>AA, TT>CC)#'
#' This function counts the different types dinucleotide substitutions (e.g. CC to AA).
#' @param VCF a delimited file with columns "WT" and "MUT", containing values including "A", "T", "G" or "C" (a data.frame)
#' @keywords monte carlo, mutants, substitutions, point mutations
#' @examples
#' res <- dinucmut.count(VCF)

dinucmut.count<- function(VCF)
{
  
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
