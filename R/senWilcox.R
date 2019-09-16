senWilcox<-function(d,gamma=1,conf.int=FALSE,alpha=0.05,alternative="greater"){

  #Check input
  stopifnot(is.vector(d)&(length(d)>1))
  stopifnot(is.vector(gamma)&(length(gamma)==1)&(gamma>=1))
  stopifnot((alternative=="greater")|(alternative=="less")|(alternative=="twosided"))
  stopifnot(is.vector(alpha)&(length(alpha)==1)&(alpha>0)&(alpha<1))

  pr<-gamma/(1+gamma)
  if (alternative=="twosided") crit<-stats::qnorm(1-(alpha/2))
  else crit<-stats::qnorm(1-(alpha))
  int<-c(min(d),max(d))

  devu<-function(taus){
    ntaus<-length(taus)
    res<-rep(NA,ntaus)
    for (i in 1:ntaus){
      dt<-d-taus[i]
      adt<-abs(dt)
      rk<-rank(adt)*(adt>0)
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
      rk<-rank(adt)*(adt>0)
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
