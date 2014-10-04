scan.kat<- function(VCF, win=3, max.threshold =400, min.burst=3)
{
  library(RcppRoll)

  # add a rolling median window
  rmax<- VCF %>%
    group_by(chr) %>%
    mutate(intermutation= start.position-lag(start.position)) %>%
    do(., data.frame(lag_roll_max=c(roll_max(.$intermutation, n=win), rep(NA, win-1))))

  

  # label points below that threshol TRUE
  kat <- mutate(VCF, below.mth = rmax$lag_roll_max < max.threshold)
  
  #see http://stackoverflow.com/questions/23820491/dplyr-error-object-not-found-using-rle-in-mutate
  # the point here is that bursts of TRUE must also be at least min.burst length
  kat <- tbl_dt(kat) %>% mutate(run_len = rep( rle(below.mth)$lengths, rle(below.mth)$lengths))
  bmth<- (kat$below.mth == TRUE) & (kat$run_len > min.burst)
  
  
  # lable each subsequent burst with a number 1:n, the rest/baseline is 0
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

