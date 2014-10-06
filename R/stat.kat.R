stat.kat <- function(VCF, bursts="bursts", intra.burst=T)
{
  kat<- VCF %>%
    group_by(bursts, chr) %>%
    mutate(intermutation=lead(start.position)-start.position) %>%
    summarise(burst.start= min(start.position),
              burst.end=max(start.position),
              burst.length= burst.end-burst.start,
              n=n(),
              min.gap = min(intermutation, na.rm=T),
              q25.gap = quantile(intermutation, 0.25, na.rm=T),
              median.gap = median(intermutation,na.rm=T),
              mean.gap = mean(intermutation,na.rm=T),
              q75.gap = quantile(intermutation, 0.25, na.rm=T),
              max.gap = max(intermutation,na.rm=T))
           
  if(intra.burst==TRUE){
    kat <- kat %>%
      group_by(chr) %>%
      mutate(intra.burst= lead(burst.start)-burst.end) %>%
      filter(bursts != 0)
    
  }
  
  return(kat)
   
}
