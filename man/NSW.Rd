\name{NSW}
\alias{NSW}
\docType{data}
\title{
National Supported Work Randomized Experiment
}
\description{
National Supported Work (NSW) randomized experiment in 185
matched pairs.  Used as an example in Chapter 2 of Design
of Observational Studies.
}
\usage{data("NSW")}
\format{
  A data frame with 370 observations on the following 11 variables.
  \describe{
    \item{\code{id}}{Matched pair id, 1, 1, 2, 2, ..., 185, 185. }
    \item{\code{z}}{z=1 for treated, z=0 for control}
    \item{\code{age}}{Age in years}
    \item{\code{edu}}{Education in years}
    \item{\code{black}}{1=black, 0=other}
    \item{\code{hisp}}{1=Hispanic, 0=other}
    \item{\code{married}}{1=married, 0=other}
    \item{\code{nodegree}}{1=no High School degree, 0=other}
    \item{\code{re74}}{Earnings in 1974, a covariate}
    \item{\code{re75}}{Earnings in 1975, a covariate}
    \item{\code{re78}}{Earnings in 1978, an outcome}
  }
}
\details{
Compare with Table 2.2 of Design of Observational Studies.
}
\source{
Dehejia and Wahba (1999).
}
\references{
Couch, K.A. (1992) <doi:10.1086/298292> "New evidence on the long-term effects of employment
training programs". Journal of Labor Economics 10, 380-388.

Dehejia, R.H., Wahba, W. (1999) <doi:10.1080/01621459.1999.10473858> "Causal effects in nonexperimental
studies: Reevaluating the evaluation of training programs".
Journal of the American Statistical Association 94, 1053-1062.

LaLonde, R.J. (1986) <https://www.jstor.org/stable/1806062> "Evaluating the econometric evaluations of training
programs with experimental data". American Economic Review 76, 604-620.
}
\examples{
data(NSW)
fivepairs<-c(15,37,46,151,181)
# Table 2.2 in Design of Observational Studies (DOS)
NSW[is.element(NSW$id,fivepairs),]
# Pair differences in earnings after treatment
dif<-NSW$re78[NSW$z==1]-NSW$re78[NSW$z==0]
# Chapter 2, footnote 7 of
stats::wilcox.test(dif,conf.int=TRUE)
}
\keyword{datasets}
