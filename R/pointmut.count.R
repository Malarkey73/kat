#' Count of different point mutations(e.g. C>A, T>C)#'
#' Ignoring Indels (e.g. ATGCGC > A) this function simply counts the different types of point mutation. Note that the CtoA category  contains both CA and the reverse complement - as in most sequencing experiments we cannot distinguish which strand was actually mutated.
#' @param VCF a delimited file with columns "WT" and "MUT", containing values including "A", "T", "G" or "C" (a data.frame)
#' @keywords monte carlo, mutants, substitutions, point mutations
#' @examples
#' res <- point_mutant_count(VCF)
#' 
pointmut.count= function(VCF)
{
  
  #C>A is the same as G>T on the -strand 
  CtoA<- VCF %>% summarise(CA=sum((WT=="C" & MUT=="A")|(WT=="G" & MUT=="T")))
  #C>G is the same as G>C on the -strand 
  CtoG<- VCF %>% summarise(CG=sum((WT=="C" & MUT=="G")|(WT=="G" & MUT=="C")))
  #C>T is the same as G>A on the -strand 
  CtoT<- VCF %>% summarise(CT=sum((WT=="C" & MUT=="T")|(WT=="G" & MUT=="A")))
  
  #T>A is the same as A>T on the -strand 
  TtoA<- VCF %>% summarise(TA=sum((WT=="T" & MUT=="A")|(WT=="A" & MUT=="T")))
  #T>C is the same as A>G on the -strand 
  TtoC<- VCF %>% summarise(TC=sum((WT=="T" & MUT=="C")|(WT=="A" & MUT=="G")))
  #T>G is the same as A>C on the -strand 
  TtoG<- VCF %>% summarise(TG=sum((WT=="T" & MUT=="G")|(WT=="A" & MUT=="C")))
  
  message("NB A count of  e.g. CtoA mutations also includes reverse complement GtoT mutations.")
  return(cbind(CtoA,CtoG,CtoT,TtoA,TtoC,TtoG))
}