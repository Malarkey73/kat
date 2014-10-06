scan.kat2<- function(VCF,  max.threshold =400)
{
  library(RcppRoll)
  
  # add 3 step lag and lead
  kat<- VCF %>%
    group_by(chr) %>%
    mutate(lag2= lag(start.position,1)-lag(start.position,2),
           lag1= start.position-lag(start.position,1),
           lead1= lead(start.position, 1)-start.position,
           lead2= lead(start.position, 2)-lead(start.position,1))
  
  # the preceding 3 gaps are below max.threshold
  kat <- mutate(kat, lagging= (lag2 < max.threshold) & (lag1 < max.threshold))
  # the next 3 gaps are below max.threshold
  kat<- mutate(kat, leading= (lead1 < max.threshold) & (lead2 < max.threshold))
  # if either leading or lagging is below mth then its within a kataegis burst
  kat <- mutate(kat, below.mth = leading | lagging)
  bmth<- (kat$below.mth)
  
  # make sure that burst are at least 3 consecutive
  bmth <- (lag(bmth) == TRUE & bmth == TRUE) | (lead(bmth) == TRUE & bmth == TRUE) 
  
  
  # label each subsequent burst with a number 1:n, the rest/baseline is 0
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
  return(as.tbl(data.frame(VCF, bursts=as.factor(vec))))
  
}

