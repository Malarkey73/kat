#' QQ-Plots of mutation distribution vs theoretical random (Poisson) distribution
#'
#' A density plot of intermutation distances with optional poissson fit
#' @param VCF a VCF tbl with columns containing at least "chr", "start.position", and "end.position" (a data.frame)
#' @param type One of "unif" or "qexp" (character, default = "unif")
#' @keywords Poisson, mutants, kataegesis
#' @export
#' @examples
#' qqmutation(VCF, type="qexp")
#'
qqmutation<- function (VCF, type = "unif") # argument: vector of numbers
{
  theme_set(theme_bw())

  if(type=="unif")
  {
  VCF <- VCF %>%
    group_by(chr) %>%
    do(data.frame(qqnorm(.$start.position, dist=qunif, plot=F)))

  qqm <- ggplot(VCF, aes(x=x,y=y))+
    theme(strip.background=element_rect(fill="white"))+
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("Theoretical") +
    ylab("Sample") +
    facet_wrap(~chr)
  }

  if(type=="qexp")
  {
  VCF <- VCF %>%
    group_by(chr) %>%
    arrange(chr,start.position) %>%
    mutate(intermutation = start.position-lag(start.position)) %>%
    filter(intermutation > 0) %>%
    do(data.frame(qqnorm(.$intermutation, dist=qexp, plot=F)))

  qqm<- ggplot(VCF, aes(x=x,y=y))+
    theme(strip.background=element_rect(fill="white"))+
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("Theoretical ") +
    ylab("Sample") +
    scale_y_log10()+
    facet_wrap(~chr)

  }

  if(type=="rexp")
  {
    VCF2 <- VCF %>%
      group_by(chr) %>%
      arrange(chr,start.position) %>%
      mutate(intermutation = start.position-lag(start.position)) %>%
      filter(intermutation > 0) %>%
      do(data.frame(rexp= rexp(n() , 1/mean(VCF$.intermutation))))

    VCF2 <- data.frame(rexp= rexp(nrow(VCF), 1/mean(VCF$intermutation)))
    qqm<- ggplot(VCF, aes(x=x,y=y))+
      theme(strip.background=element_rect(fill="white"))+
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      xlab("Theoretical ") +
      ylab("Sample") +
      scale_y_log10()+
      facet_wrap(~chr)

  }


  return(qqm)
}
