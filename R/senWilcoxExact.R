senWilcoxExact<-function (d, gamma = 1)
{
  # From Appendix Section 3.9 of Design of Observational Studies 2010.

  # Check input.
  stopifnot(is.vector(d)&(length(d)>1))
  stopifnot(is.vector(gamma)&(length(gamma)==1)&(gamma>=1))
  ad<-abs(d)

  if (length(ad)>length(unique(ad))){
    stop("For exact calculations, ties in |d| are not permitted.  Use senWilcox instead.")
    }

  if (sum(d==0)>0){
    d<-d[d!=0]
    warning("In computing the exact distribution, zero differences in d were removed, reducing
            the sample size.")
  }

  gconv <- function (g1, g2)
  {
    stats::convolve(g1, rev(g2), type = "o")
  }

  wilcsenexact <- function (n, gamma = 1)
  {
    p <- gamma/(1 + gamma)
    g <- c(1 - p, p)
    for (i in 2:n) {
      gi <- rep(0, i + 1)
      gi[1] <- 1 - p
      gi[i + 1] <- p
      g <- gconv(g, gi)
    }
    g
  }

  wilcsenexacttail <- function (k, n, gamma = 1)
  {
    1 - sum(wilcsenexact(n, gamma = gamma)[1:k])
  }

  a <- abs(d)
  rk <- rank(a)
  sgn <- (d > 0) * 1
  wt <- sum(sgn * rk)
  wilcsenexacttail(floor(wt), length(d), gamma = gamma)
}

