cohere<-function(y,z,mpair,w=NULL,gamma=1,m=NULL,m1=NULL,m2=NULL,
                     apriori=FALSE,Scheffe=FALSE,exact=NULL){
  #Check input
  stopifnot(is.vector(gamma)&(length(gamma)==1)&(gamma >= 1))
  stopifnot((apriori==TRUE)|(apriori==FALSE))
  stopifnot((Scheffe==TRUE)|(Scheffe==FALSE))
  stopifnot(all(as.vector(!is.na(y)))) # Check for missing values
  stopifnot(all((z==0)|(z==1))) #z is 1 for treated, 0 for control
  if(!is.factor(mpair)){
    intck<-all(mpair==round(mpair))
    if(!intck){
      warning("mpair must be either be integer or a factor.")
      stopifnot(intck)
    }
  }
  mset<-mpair
  tbcheck<-table(z,mset)
  ck<-all(tbcheck[2,]==1)&all(tbcheck[1,]==1)
  if (!ck){
    warning("Every matched set must contain one treated subject and one control.
            Use the sensitivitymult package for matching with multiple controls.")
    stopifnot(ck)
  }
  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(length(z)==(dim(y)[1]))
  stopifnot(length(mset)==(dim(y)[1]))
  nvars<-dim(y)[2]
  if (nvars<2){
    warning("y must have at least two outcomes in two columns.")
    stopifnot(nvars>=2)
  }
  if (is.null(w)) w<-rep(1,nvars)
  stopifnot(length(w)==(dim(y)[2]))
  stopifnot(sum(abs(w))>0)
  if (!is.null(colnames(y))) names(w)<-colnames(y)
  if ((!is.null(m))|(!is.null(m1))|(!is.null(m2))){
    ck2<-is.null(m)|is.null(m1)|is.null(m2)
    if (ck2){
      warning("If any of (m,m1,m2) is not NULL, then all three must not be NULL.")
      stopifnot(!ck2)
    }
    stopifnot(is.vector(m)&(length(m)==1)&(m==round(m))&(m>=1))
    stopifnot(is.vector(m1)&(length(m1)==1)&(m1==round(m1))&(m1>=1))
    stopifnot(is.vector(m2)&(length(m2)==1)&(m2==round(m2))&(m2>=1))
    stopifnot((m1<=m2)&(m2<=m))
  }
  if (!is.null(exact)) stopifnot(is.logical(exact))

  if (is.null(exact)) {
    if ((dim(y)[1])<=100) exact<-TRUE
    else exact<-FALSE
  }

  mset<-as.integer(mset)
  o<-order(mset,1-z)
  y<-y[o,]
  z<-z[o]
  mset<-mset[o]
  v<-y[z==1,]-y[z==0,]
  if (!is.null(colnames(y))) colnames(v)<-colnames(y)
  rownames(v)<-mset[z==1]

  multrnks <- function (rk, m1 = 2, m2 = 2, m = 2)
  {
    n <- length(rk)
    q <- rep(0, n)

    if (exact) {
      for (l in m1:m2) {
        # This is expression (8) in Rosenbaum (2011) Biometrics 67, 1017-1027
        q <- q + (choose(rk-1,l-1)/choose(n,m))*choose(n-rk,m-l)
      }
    }
    else{
      pk <- rk/n
      q <- rep(0, n)
      for (l in m1:m2) {
        # This is expression (9) in Rosenbaum (2011) Biometrics 67, 1017-1027
        q <- q + (l*choose(m,l)*(pk^(l-1))*((1-pk)^(m-l)))
      }
    }
    q
  }

  vrk<-v
  for (j in 1:nvars){
    vark<-rank(abs(v[,j]))
    if (!is.null(m)) vrk[,j] <- multrnks(vark,m=m,m1=m1,m2=m2) * sign(v[,j])
    else vrk[,j] <- vark * sign(v[,j])
  }
  vrk<-as.matrix(vrk)
  sc<-vrk%*%w
  pr<-gamma/(1+gamma)
  TS<-sum(sc)
  eTS<-((gamma-1)/(gamma+1))*sum(abs(sc))
  vTS<-4*pr*(1-pr)*sum(sc*sc)
  deviate<-(TS-eTS)/sqrt(vTS)

  if (Scheffe){
    ScheffePVal<-1-stats::pchisq(max(0,deviate)^2,nvars)
    list(deviate=deviate,ScheffePVal=ScheffePVal,weights=w)
  }
  else if (apriori){
    aprioriPVal<-1-stats::pnorm(deviate)
    list(deviate=deviate,aprioriPVal=aprioriPVal,weights=w)
  }
  else list(deviate=deviate,weights=w)
}
