#' Intermutation Density Plot
#'
#' A density plot of intermutation distances with optional poissson fit
#' @param VCF a VCF tbl with columns containing at least "chr", "start.position", and "end.position" (a data.frame)
#' @keywords rainfall, mutants, kataegesis
#' @export
#' @examples
#' imdensity.plot(VCF, chr=6)

imdensity.plot= function(VCF, chrN=NULL, facets=NULL, type="density"){

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
    VCF <- filter(VCF, chr %in% chrN)

  m <- ggplot(VCF, aes(x = log10(intermutation)))+
    theme(strip.background=element_rect(fill="white"))


  if(!is.null(facets))
    m <- m + facet_wrap(as.formula(sprintf("~ %s", facets)))

  if(type=="density")
    return(m+ geom_density())

  if(type=="histogram")
    return(m+ geom_histogram())
}
