scan.kat<- function(VCF, win=3, max.threshold =500)
{
  library(RcppRoll)

  # add a rolling median window
  kat<- VCF %>%
    group_by(chr) %>%
    mutate(intermutation = start.position-lag(start.position)) %>%
    do(., data.frame(roll_max=c(roll_max(.$intermutation, n=win), rep(NA, win-2))))

  

  # label poinst below that threshold TRUE
  kat<- mutate(kat, below.mth = roll_max < max.threshold)

  bmth <- kat$below.mth
  burst.gr<-0
  vec=rep(NA, length(bmth))
  last=F
  
  for(i in 1:length(bmth))
  {
    if(bmth[i]==FALSE| is.na(bmth[i]))
    {
      vec[i]<-0
      last<- FALSE
    }
    if(bmth[i]==TRUE & last==FALSE | !isna(bmth[i]))
    {
      burst.gr<- burst.gr +1
      last<-TRUE
      vec[i]<- burst.gr
    }
    if(bmth[i]==TRUE & last==TRUE | !isna(bmth[i]))
    {
      burst.gr<- burst.gr +1
      last<-TRUE
      vec[i]<- burst.gr
    }
    
  }
