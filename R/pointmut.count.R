
pointmut.count= function(VCF)
{
  
  #C>A is the same as G>T on the -strand 
  CA<- VCF %>% summarise(CA=sum((WT=="C" & MUT=="A")|(WT=="G" & MUT=="T")))
  #C>G is the same as G>C on the -strand 
  CG<- VCF %>% summarise(CG=sum((WT=="C" & MUT=="G")|(WT=="G" & MUT=="C")))
  #C>T is the same as G>A on the -strand 
  CT<- VCF %>% summarise(CT=sum((WT=="C" & MUT=="T")|(WT=="G" & MUT=="A")))
  
  #T>A is the same as A>T on the -strand 
  TA<- VCF %>% summarise(TA=sum((WT=="T" & MUT=="A")|(WT=="A" & MUT=="T")))
  #T>C is the same as A>G on the -strand 
  TC<- VCF %>% summarise(TC=sum((WT=="T" & MUT=="C")|(WT=="A" & MUT=="G")))
  #T>G is the same as A>C on the -strand 
  TG<- VCF %>% summarise(TG=sum((WT=="T" & MUT=="G")|(WT=="A" & MUT=="C")))
  
  return(cbind(CA,CG,CT,TA,TC,TG))
}