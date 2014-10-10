scan.kat<- function(VCF, mut.gap =300)
{

  # add a rolling median window
  kat<- VCF %>%
    group_by(chr) %>%
    mutate(inburst= ifelse((start.position-lag(start.position) < mut.gap | 
                            lead(start.position)-start.position < mut.gap), 
                            TRUE, 
                            FALSE)) %>%
    mutate(inburst= ifelse((inburst==T & lead(inburst,1)==T & lead(inburst,2)==T)|
                          (inburst==T & lag(inburst,1)==T & lag(inburst,2)==T), 
                          TRUE, 
                          FALSE))  
  
  inburst<- kat$inburst
  
  # label each subsequent burst with a number 1:n, the rest/baseline is 0
  burst.gr<-0
  vec=rep(NA, length(inburst))
  last=F
  
  for(i in 1:length(inburst))
  {
    if(inburst[i]==FALSE | is.na(inburst[i]))
    {
      vec[i]<-0
      last<- FALSE
    }
    
    if(inburst[i]==TRUE & last==FALSE & !is.na(inburst[i]))
    {
      burst.gr<- burst.gr +1
      last<-TRUE
      vec[i]<- burst.gr
    }
    
    if(inburst[i]==TRUE & last==TRUE & !is.na(inburst[i]))
    {
      vec[i]<- burst.gr
    }
    
  }
  kat$bursts <- as.factor(vec)
  return(kat)
  
}

