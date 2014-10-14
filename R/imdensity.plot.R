#' Intermutation Density Plot
#' A density plot of intermutation distances with optional poissson fit. This is yet another good plot for finding chromosomes with areas of kataegis.
#' @param VCF a VCF tbl with columns containing at least "chr", "start.position", and "end.position" (a data.frame)
#' @param chr
#' @keywords rainfall, mutants, kataegesis
#' @export
#' @examples
#' imdensity.plot(VCF, chr=6)

imdensity.plot <- function(VCF, chrN=NULL, facets=NULL, type=density, adjust=0.2, binwidth=0.2, ...){
  
    type <- deparse(substitute(type))
  if(!charmatch(type, c("histogram", "density")))
    stop("unknown plot type requested")
  
  # intermutaion distance is calculated for each chromosme separately
  VCF <- VCF %>%
    filter(start.position==end.position) %>%
    arrange(chr, start.position) %>%
    group_by(chr)  %>%
    mutate(intermutation = start.position-lag(start.position)) %>%
    filter(intermutation != 0)

  library(ggplot2)
  theme_set(theme_bw())


  if(!is.null(chrN))
  {
    if(!is.character(chrN)) 
      chrN<- deparse(substitute(chrN))
    VCF <- filter(VCF, chrN %in% chr)  
  }

  m <- ggplot(VCF, aes(x = log10(intermutation)))+
    theme(strip.background=element_rect(fill="white"))


  if(!is.null(facets))
  {
    if(!is.character(facets))
      facets<- deparse(substitute(facets))
    m <- m + facet_wrap(as.formula(sprintf("~ %s", facets)))
  }
    

  if(pmatch(type, "density", nomatch=0))
    return(m+ geom_density(adjust=adjust))
  if(pmatch(type, "histogram", nomatch=0))
    return(m+ geom_histogram(binwidth = binwidth))
  
  
}
