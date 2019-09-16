addcaliper <- function (dmat, z, p, caliper = 0.1, penalty = 1000)
{
  #Check input
  stopifnot(is.vector(z))
  stopifnot(all((z==0)|(z==1)))
  stopifnot(is.vector(p)&is.numeric(p))
  stopifnot(length(z)==length(p))
  stopifnot(stats::sd(p)>0)
  stopifnot(is.matrix(dmat))
  stopifnot(length(z)==sum(dim(dmat)))
  stopifnot(sum(z)==(dim(dmat)[1]))
  stopifnot(is.vector(caliper)&(length(caliper)==1)&(caliper>0))
  stopifnot(is.vector(penalty)&(length(penalty)==1)&(penalty>0))


  # Add penalty for caliper violations
  sdp <- stats::sd(p)
  adif <- abs(outer(p[z == 1], p[z == 0], "-"))
  adif <- (adif - (caliper * sdp)) * (adif > (caliper * sdp))
  dmat <- dmat + adif * penalty
  dmat
}
