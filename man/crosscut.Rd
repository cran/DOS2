\name{crosscut}
\alias{crosscut}
\concept{Crosscut test}
\concept{Observational study}
\concept{Dose-response}
\concept{Sensitivity analysis}

\title{
Crosscut Test and its Sensitivity Analysis
}
\description{
Computes the cross-cut test and its sensitivity analysis.  The cross-cut test is
a nonparametric test of dose-response correlation with good design sensitivity
when used for causal inference in observational studies.
}
\usage{
crosscut(x, y, ct = 0.25, gamma = 1, LS=FALSE)
}

\arguments{
  \item{x}{
Doses of treatment.
}
  \item{y}{
Response.
}
  \item{ct}{
The quantile that defines the cross-cut.  By default, the cross-cut
is at the outer .25 of the data, the lower 25 percent and the upper 75
percent.
}
  \item{gamma}{
Sensitivity parameter, gamma>=1.
}
  \item{LS}{
If LS=TRUE, a large sample test is performed.  If LS=FALSE,
an exact test is performed.  For LS=FALSE, the mh function
in the 'sensitivity2x2xk' package is used.  For LS=TRUE, the mhLS function
in the 'sensitivity2x2xk' package is used.
}
}
\details{
Performs a one-sided test of no association against positive association,
together with a sensitivity analysis.  The method is described in
Rosenbaum (2016), used in Rosenbaum (2017).  An adaptive cross-cut statistic is discussed in
Rosenbaum and Small (2017); it cuts at several quantiles and picks the best.  See Section 19.4 of "Design of
Observational Studies"", second edition.
}
\value{
\item{quantiles }{Quantiles that define the crosscut}
\item{table }{A 2x2 table}
\item{output }{Output from mh or mhLS when applied to table.  The
functions mh and mhLS are from the sensitivity2x2xk package.
The output includes a one-sided P-value.}
}
\references{
Rosenbaum, P. R. (2016) <doi:10.1111/biom.12373> "The crosscut statistic and its sensitivity to bias in observational studies with ordered doses of treatment". Biometrics, 72(1), 175-183.

Rosenbaum, P. R. (2017) <doi:10.1214/17-STS621> "The general structure of evidence factors in observational studies".
Statist Sci 32, 514-530.

Rosenbaum, P. R. and Small, D. S. (2017) <doi:10.1111/biom.12591> "An adaptive Mantel–Haenszel test for sensitivity analysis in observational studies". Biometrics, 73(2), 422-430.
}
\author{
Paul R. Rosenbaum
}
\note{
The 'crosscut' function makes use of 'mh' and 'mhLS' from the
'sensitivity2x2xk' package.
}

\examples{
data(periodontal)
attach(periodontal)
crosscut(cigsperday[z==1],pcteither[z==1]-pcteither[z==0],ct=.2)
crosscut(cigsperday[z==1],pcteither[z==1]-pcteither[z==0],ct=.2,gamma=1.25)
crosscut(cigsperday[z==1],pcteither[z==1]-pcteither[z==0],ct=.2,gamma=1.25,LS=TRUE)
crosscut(cigsperday[z==1],pcteither[z==1]-pcteither[z==0],ct=1/3)
detach(periodontal)

}

