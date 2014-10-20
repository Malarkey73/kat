#' Dinucleotide Monte Carlo Simulation'
#' In Nik-Zainal et al. 2012 paper "Mutational Processes Molding the Genomes of 21 Breast Cancers", they note a higher rate of dinucleotide mutations than expected by chance. They describe this as estimated by Monte Carlo simulations. This function is my version of such a simulation. I scramble the existing data n=100 times, assuming the same number of mutations on each chromosoem Poisson distributed across the observed range e.g. between the first to last actual mutation. The function then calculates permuted counts (Expected Counts :E) for each dinucleotide type each time. The fraction of E > O (Observed Counts: O) is a Monte Carlo estimate of P (the null hypothesis).   
#' @param VCF a delimited file with columns "WT" and "MUT", containing values including "A", "T", "G" or "C" (a data.frame)
#' @param n.mc the number of Monte Carlo simulations (permutations) of the expected dincleotide rate E. 
#' @keywords monte carlo, mutants, substitutions, point mutations
#' @examples
#' res <- dinucmut.mc(VCF, n=100)

dinucmut.mc <- function(VCF, n.mc=100)
  {
  if (!("VCF" %in% class(VCF)))
    stop(" A valid VCF needs at least chr, start.position, end.position, WT and MUT columns. 
         You can read data files using read.mutations() or convert your data into this format using as.VCF()")
  
  # the observed dinucleotide counts
  Observed <- dinucmut.count(VCF)
  
  # n=100 rows by default, a column per dinuc type, same as Observed 
  Expected=as.data.frame(matrix(NA, n.mc, 6))
  names(Expected)=c("CCtoAA", "CCtoGG", "CCtoTT", "TTtoAA", "TTtoCC", "TTtoGG")
  
  VCF<- group_by(VCF, chr)
  
  VCF.summ<- VCF %>%
    summarise(first.position=first(start.position), last.position=last(start.position), n=n())
  
  
  
  for(i in 1:n.mc)
  {
    
    # This essentially keeps the same mutations of same types on same chromosomes but just moves them to random positions 
    VCF.scramble<- VCF %>%
      mutate(start.position = unlist(apply(VCF.summ, 1, function(x)floor(runif(x[4], x[2], x[3]))))) %>%
      arrange(start.position)
    
    class(VCF.scramble)<- c(class(VCF.scramble), "VCF")
    Expected[i,] <- dinucmut.count(VCF.scramble)
    message(i, " scrambles")
  }
  
  P.mc = rowMeans(apply(Expected,1, function(x) x > Observed))
  names(P.mc) = names(Expected)
  return(P.mc)
}
