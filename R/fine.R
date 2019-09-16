fine <- function (dmat, z, f, mult = 100)
{
  #Check input
  stopifnot(is.vector(z))
  stopifnot(all((z==0)|(z==1)))
  if (is.vector(f)&!is.factor(f)){
    stopifnot(length(f)>length(unique(f)))
    f<-factor(f)
  }
  stopifnot(length(f)==length(z))
  stopifnot(is.factor(f))
  stopifnot(is.matrix(dmat))
  stopifnot(length(z)==sum(dim(dmat)))
  stopifnot(sum(z)==(dim(dmat)[1]))
  stopifnot(is.vector(mult)&(length(mult)==1)&(mult>0))

  #Expand distance matrix
  penalty <- mult * max(dmat)
  n <- dim(dmat)[1]
  m <- dim(dmat)[2]
  if (n>=m){
    stop("The distance matrix has too few columns to permit expansion for fine balance.")
  }
  t1 <- table(f[z == 1])
  t0 <- table(f[z == 0])
  L <- length(t1)
  remove <- t0 - t1
  if (min(remove) < 0) {
    stop("fine balance is infeasible.")
  }
  else {
    f<-as.integer(f)
    for (j in 1:L) {
      who <- (f[z == 0] != j) * 1
      add <- t(matrix(rep(who * penalty, remove[j]), m,
                      remove[j]))
      dmat <- rbind(dmat, add)
    }
    out <- dmat
    rownames(out)[(n + 1):m] <- (n + m + 1):(2 * m)
  }
  out
}
