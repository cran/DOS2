addalmostexact<-function(dmat,z,f,mult=10){
  #For the nominal variable f, adds a penalty for mismatch

  #Check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(f)|is.factor(f))
  stopifnot(length(z)==length(f))
  stopifnot(is.matrix(dmat))
  stopifnot(length(z)==sum(dim(dmat)))
  stopifnot(length(z)==length(f))
  stopifnot(is.vector(mult)&(mult>0))
  stopifnot(all((z==0)|(z==1)))

  #Add penalty to dist
  penalty<-mult*max(dmat)
  mismatch<-outer(f[z==1],f[z==0],"!=")
  dmat<-dmat+mismatch*penalty
  dmat
}
