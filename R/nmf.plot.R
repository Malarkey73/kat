#' Non-negative matrix factorisation of mutation types
#' Alexandrov et al. pioneered the extraction of mutational signatures using the NMF technique
#' @param count.df The format is a column named triplets (ATC, ATG, etc), a column named mutations(C, G, etc) and then columns of mutation counts for a set of samples
#' @param rank This is the number of mutation signatures you think maybe hiding within your data
#' @param method The NMF algorithm search algorithm. See the NMF package documentation for details of methods.
#' @param seed The NMF starting point guess. See the NMF package documentation for details of seeding.
#' @keywords NMF, non-negative matrix factorisation, mutations, kataegesis
#' @export
#' @examples
#' nmf.plot(count.df, seed=ica)


nmf.plot<- function(count.df, rank=4, method=nsNMF, seed=ica, nrun=100, ...)
{
  method=deparse(substitute(method))
  seed=deparse(substitute(seed))
  
  require(NMF)
  triplets<- count.df$triplets
  mutations<- count.df$mutations
  count.m<-count.df[,-match(c("triplets", "mutations"), colnames(count.df))]

  res <- nmf(count.m, rank = rank, method=method, seed=seed, nrun=nrun, ...)
  sig.basis <- basis(res)

  sig.df<- data.frame(vals= sig.basis[1:length(sig.basis)],
                triplets=rep(triplets,ncol(sig.basis)),
                mutations=rep(mutations, ncol(sig.basis)),
                sig = rep(1:ncol(sig.basis), each=nrow(sig.basis)))

  library(ggplot2)
  theme_set(theme_classic())
  qplot(mutations, vals, data= sig.df, geom="point")+
    facet_wrap(sig~triplets, nrow=3)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

}
