senU<-function(d,gamma=1,m=2,m1=2,m2=2,conf.int=FALSE,alpha=0.05,alternative="greater",exact=NULL){

  #Check input
  stopifnot(is.vector(d)&(length(d)>1))
  stopifnot(is.vector(gamma)&(length(gamma)==1)&(gamma>=1))
  stopifnot((alternative=="greater")|(alternative=="less")|(alternative=="twosided"))
  stopifnot(is.vector(alpha)&(length(alpha)==1)&(alpha>0)&(alpha<1))
  stopifnot(is.vector(m)&(length(m)==1)&(m==round(m))&(m>=1))
  stopifnot(is.vector(m1)&(length(m1)==1)&(m1==round(m1))&(m1>=1))
  stopifnot(is.vector(m2)&(length(m2)==1)&(m2==round(m2))&(m2>=1))
  stopifnot((m1<=m2)&(m2<=m))
  if (!is.null(exact)) stopifnot(is.logical(exact))


  if (is.null(exact)) {
    if (length(d)<=50) exact<-TRUE
    else exact<-FALSE
  }
  pr<-gamma/(1+gamma)
  if (alternative=="twosided") crit<-stats::qnorm(1-(alpha/2))
  else crit<-stats::qnorm(1-(alpha))
  int<-c(min(d),max(d))

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

  devu<-function(taus){
    ntaus<-length(taus)
    res<-rep(NA,ntaus)
    for (i in 1:ntaus){
      dt<-d-taus[i]
      adt<-abs(dt)
      rk<-multrnks(rank(adt),m1=m1,m2=m2,m=m)*(adt>0)
      sg<-1*(dt>0)
      ts<-sum(sg*rk)
      ex<-sum(pr*rk)
      va<-sum(rk*rk)*pr*(1-pr)
      res[i]<-(ts-ex)/sqrt(va)
    }
    res
  }

  devl<-function(taus){
    ntaus<-length(taus)
    res<-rep(NA,ntaus)
    for (i in 1:ntaus){
      dt<-taus[i]-d
      adt<-abs(dt)
      rk<-multrnks(rank(adt),m1=m1,m2=m2,m=m)*(adt>0)
      sg<-1*(dt>0)
      ts<-sum(sg*rk)
      ex<-sum(pr*rk)
      va<-sum(rk*rk)*pr*(1-pr)
      res[i]<-(ts-ex)/sqrt(va)
    }
    res
  }

  if (alternative=="greater") pval <- 1-stats::pnorm(devu(0))
  else if (alternative=="less") pval <- 1-stats::pnorm(devl(0))
  else {
    pvall <- 1-stats::pnorm(devl(0))
    pvalu <- 1-stats::pnorm(devu(0))
    pval <- min(1,2*min(pvall,pvalu))
    }

  estimate<-c(-Inf,Inf)
  names(estimate)<-c("low","high")
  ci<-estimate

  if ((alternative!="less")&(conf.int==TRUE)) {
    estimate[1]<-stats::uniroot(devu,int)$root
    devuCI<-function(taus){devu(taus)-crit}
    ci[1]<-stats::uniroot(devuCI,int)$root
  }

  if ((alternative!="greater")&(conf.int==TRUE)) {
    estimate[2]<-stats::uniroot(devl,int)$root
    devuCI<-function(taus){devl(taus)-crit}
    ci[2]<-stats::uniroot(devuCI,int)$root
  }

  if (conf.int==TRUE) list(pval=pval,estimate=estimate,ci=ci)
  else list(pval=pval)
}
