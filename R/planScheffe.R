planScheffe<-function(K,alpha=0.05){
  stopifnot(K>=2)
  stopifnot((alpha>0)&(alpha<1))

  jcum<-function(av,cv,K){
    #Joint probability, K>=2 variables, av and cv are scalars
    stopifnot(K>=2)
    if (cv<=0) {0}
    else{
      int<-function(x){
        stats::dnorm(x)*stats::pchisq(cv-(x*x),K-1)
      }
      stats::integrate(int,-Inf,av)
    }
  }

  diagprob<-function(cv,K){
    #Gives the cumulative probability with equal marginals
    #cv may be a vector
    stopifnot(K>=2)
    av<-stats::qnorm(stats::pchisq(cv,K))
    o<-rep(NA,length(cv))
    for (i in 1:length(cv)){
      o[i]<-as.vector(jcum(av[i],cv[i],K)$value)
    }
    names(o)<-cv
    o
  }


  up<-stats::qchisq(1-(alpha/2),K)
  lo<-stats::qchisq(1-alpha,K)
  alpha1<-1-alpha
  f<-function(cv){
    diagprob(cv,K)-alpha1
  }
  cv<-stats::uniroot(f,c(lo,up))$root
  alphacv<-1-stats::pchisq(cv,K)
  av<-stats::qnorm(1-alphacv)
  alphaav<-alphacv
  alphajnt<-1-as.vector(jcum(av,cv,K)$value)
  a<-c(alphaav,alphacv,alphajnt)
  names(a)<-c("a","c","joint")
  v<-c(av,cv)
  names(v)<-c("a","c")
  list(critical=v,alpha=a)
}
