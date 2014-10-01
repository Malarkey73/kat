#' Frequency of dinucleotide mutations (e.g. CC>AA, TT>CC)
#'
#' In Nik-Zainal(2012) they describe a method to asees the numberof dinucleotide mutaions by first measuring the frequency of single point mutaion types (e.g A>C, T>C)
#' for each chromosome. These frequencies are then used to simulate many 1000s of  times the expected dinucleotide frequencies given
#' @param VCF a delimited file with columns containing at least "chr", "start.position", and "end.position" (a data.frame)
#' @keywords monte carlo, mutants, dinucleotide
#' @examples
#' res <- dinuc.freq.table(VCF)

dinucmut.count<- function(VCF)
{

  dn <- VCF %>%
    group_by(chr) %>%
    arrange(chr,start.position) %>%
    mutate(intermutation.lag = start.position-lag(start.position), intermutation.lead = lead(start.position)-start.position) %>%
    filter((intermutation.lead ==1 & nchar(WT)==1 & nchar(MUT)==1 & lead(WT)==WT & lead(MUT)==MUT)|
             (intermutation.lag ==1 & nchar(WT)==1 & nchar(MUT)==1 & lag(WT)==WT & lag(MUT)==MUT)) %>%
    mutate(dinucWT = ifelse(intermutation.lead==1, paste0(WT,lead(WT)), NA), dinucMUT = ifelse(intermutation.lead==1, paste0(MUT,lead(MUT)),NA)) %>%
    filter(!is.na(dinucWT)) %>%
    ungroup() %>%
    pointmut.count()


  names(dn)=c("CCtoAA", "CCtoGG", "CCtoTT", "TTtoAA", "TTtoCC", "TTtoGG")
  return(dn)

}
