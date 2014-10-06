stat.kat <- function(VCF, bursts="bursts")
{
  kat<- VCF %>%
    group_by(bursts) %>%
    
    mutate(fivenum(bursts))
}