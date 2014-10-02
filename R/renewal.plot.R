#' A renewal plot to detect non homogenous Poisson processes
#'
#' This idea is borrowed from the analysis of spike train dat described in the STAR package - and reapplied here to a completely different problem. The basic idea is to use the
#' rank order of the intermutation distances OJ plotted against OJ+1 (the lag). This uses the plto space better than just plotting intermutaion vs intermutaion lag.
#' Homogenous Poisson processes fill the square  space randomly, deviations from homogeneity leave gaps.
#' @param VCF
#' @keywords rainfall, mutants, kataegesis
#' @examples
#' # a minimal example
#' VCF<-read.mutations(file="1317049.tsv", chr="Chromosome", start.position="Genome.start", end.position="Genome.stop")
#' renewal.plot(VCF)

renewal.plot<-function(VCF, chr=1)
{
  kat<- VCF %>%
    group_by(chr) %>%
    mutate(intermutation=start.position-lag(start.position)) %>%
    mutate(OJ=rank(intermutation)) %>%
    mutate(OJ.l=lag(OJ))

  library(ggplot2)
  theme_set(theme_bw())
  qplot(OJ, OJ.l, data=kat)+facet_wrap(~chr, scales = "free")

}
