scan.kat<- function(VCF, win=3, median.threshold =500)
{
  library(RcppRoll)

  # add a rolling median window
  kat<- VCF %>%
    group_by(chr) %>%
    mutate(intermutation = start.position-lag(start.position)) %>%
    do(., data.frame(roll_median=c(roll_median(.$intermutation, n=win), rep(NA, win-2))))

  # remove NA values from chromosome window edges
  kat <- kat[complete.cases(kat),]

  # label poinst below that threshold TRUE
  kat<- mutate(kat, below.mt = roll_median < median.threshold)

  # make that into an RLE (run length encoding)
  rl= rle(kat$below.mt)

  v=rep(0,sum(rl[[1]]))

