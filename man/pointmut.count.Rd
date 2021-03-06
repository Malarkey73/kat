% Generated by roxygen2 (4.0.2): do not edit by hand
\name{pointmut.count}
\alias{pointmut.count}
\title{Count of different point mutations(e.g. C>A, T>C)#'
Ignoring Indels (e.g. ATGCGC > A) this function simply counts the different types of point mutation. Note that the CtoA category  contains both CA and the reverse complement - as in most sequencing experiments we cannot distinguish which strand was actually mutated.}
\usage{
pointmut.count(VCF)
}
\arguments{
\item{VCF}{a delimited file with columns "WT" and "MUT", containing values including "A", "T", "G" or "C" (a data.frame)}
}
\description{
Count of different point mutations(e.g. C>A, T>C)#'
Ignoring Indels (e.g. ATGCGC > A) this function simply counts the different types of point mutation. Note that the CtoA category  contains both CA and the reverse complement - as in most sequencing experiments we cannot distinguish which strand was actually mutated.
}
\examples{
res <- point_mutant_count(VCF)
}
\keyword{carlo,}
\keyword{monte}
\keyword{mutations}
\keyword{mutations,}
\keyword{point}
\keyword{substitutions,}

