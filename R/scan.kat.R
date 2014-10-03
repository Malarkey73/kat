scan.kat<- function(VCF, win=5, max.threshold =500)
{
  library(RcppRoll)

  # add a rolling median window
  rmax<- VCF %>%
    group_by(chr) %>%
    mutate(intermutation = start.position-lag(start.position)) %>%
    do(., data.frame(roll_max=c(roll_max(.$intermutation, n=win), rep(NA, win-1))))

  

  # label poinst below that threshold TRUE
  kat<- mutate(VCF, below.mth = rmax$roll_max < max.threshold)

  bmth <- kat$below.mth
  burst.gr<-0
  vec=rep(NA, length(bmth))
  last=F
  
  for(i in 1:length(bmth))
  {
    if(bmth[i]==FALSE | is.na(bmth[i]))
    {
      vec[i]<-0
      last<- FALSE
    }
    
    if(bmth[i]==TRUE & last==FALSE & !is.na(bmth[i]))
    {
      burst.gr<- burst.gr +1
      last<-TRUE
      vec[i]<- burst.gr
    }
    
    if(bmth[i]==TRUE & last==TRUE & !is.na(bmth[i]))
    {
      vec[i]<- burst.gr
    }
    
  }
  return(mutate(kat, bursts=as.factor(vec)))
  
}

