#' Dinucleotide Monte Carlo Simulation'
#' In Nik-Zainal et al. 2012 paper "Mutational Processes Molding the Genomes of 21 Breast Cancers", they note a higher rate of dinucleotide mutations than expected by chance. They describe this as estimated by Monte carlo Simulations. This function is my version of such a simulation. I permute the existing data n=100 times (within the actual chromosomes) calculating permuted counts (Expected Counts :E) for each dinucleotide type each time. The fraction of E > O (Observed Counts: O) is a Monte Carlo estimate of P (the null hypothesis).   
#' @param VCF a delimited file with columns "WT" and "MUT", containing values including "A", "T", "G" or "C" (a data.frame)
#' @keywords monte carlo, mutants, substitutions, point mutations
#' @examples
#' res <- dinucmut.mc(VCF)

dinucmut.mc <- function(VCF)
  {
  
  }
