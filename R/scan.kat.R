scan.kat<- function(VCF, below=100)
{
  library(RcppRoll)

  kat<- VCF %>%
    group_by(chr) %>%
    mutate(intermutation = start.position-lag(start.position)) %>%
    mutate(below.t = intermutation < below)

  #medwin = do(kat, data.frame(rollmed = c(roll_median(.$intermutation, window), rep(NA, window-1))< threshold))

  #VCF <- mutate(VCF, kat= medwin$rollmed)

  return(kat)

}
