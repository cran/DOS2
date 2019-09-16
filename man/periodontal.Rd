\name{periodontal}
\alias{periodontal}
\docType{data}
\title{
Smoking and Periodontal Disease
}
\description{
Data from NHANES 2011-2012 containing 441 matched pairs of a daily cigarette smoker and a never smoker, recording the extent of periodontal disease.
See Rosenbaum (2017) and Chapter 20 of "Design of Observational Studies",
second edition.
}
\usage{data("periodontal")}
\format{
  A data frame with 882 observations on the following 12 variables.
  \describe{
    \item{\code{SEQN}}{NHANES 2011-2012 sequence number}
    \item{\code{female}}{=1 for female, 0 for male}
    \item{\code{age}}{Age in years}
    \item{\code{black}}{=1 for black, 0 for other}
    \item{\code{educf}}{Education, in five categories. An ordered factor with levels \code{<9} for less than 9th grade, \code{9 to 11} for 9th to 11th grade, \code{HS/GED} for high school or GED degree,
    \code{SomeCol} for some college, \code{College} for college degree.}
    \item{\code{income}}{Ratio of family income to the poverty level, capped at 5 times
    the poverty level.}
    \item{\code{cigsperday}}{Cigarettes smoked per day for daily smokers, 0 for never smokers}
    \item{\code{either}}{Number of periodonal measurements indicative of periodontal disease.}
    \item{\code{neither}}{Number of periodonal measurements}
    \item{\code{pcteither}}{Percent indicative of periodontal disease, =100*either/neither.}
    \item{\code{z}}{Treatment indicator, 1=daily smoker, 0=never smoker}
    \item{\code{mset}}{Matched set indicator, 1 to 441.}
  }
}
\details{
Excluding wisdom teeth, 6 measurements are taken for each tooth that is present, up to 28 teeth.
Following Tomar and Asma (2000), a measurement indicates periodontal disease if either there is
a loss of attachment of at least 4mm or a pocket depth of at least 4mm.  The first individual
has 11 measurements indicative of periodontal disease, out of 106 measurements, so
pcteither is 100*11/106 = 10.38 percent.  A related data set in DOS2 with bivariate outcome is teeth.
}
\source{
Data are from the National Health and Nutrition Examination Survey 2011-2012 and were
used as an example in Rosenbaum (2017).  In the second edition of Design of Observational
Studies, these data are discussed in Chapter 20, Evidence Factors.
}
\references{

Rosenbaum, P. R. (2015) <https://obsstudies.org/two-r-packages-for-sensitivity-analysis-in-observational-studies/> "Two R packages for sensitivity analysis in observational studies". Observational Studies, 1(1), 1-17.

Rosenbaum, P. R. (2017) <doi:10.1214/17-STS621> "The general structure of evidence factors in observational studies".
Statistical Science 32, 514-530.

Tomar, S. L. and Asma, S. (2000) <doi:10.1902/jop.2000.71.5.743>  "Smoking attributable periodontitis in the US: Findings
from NHANES III".  J Periodont 71, 743-751.

"US National Health and Nutrition Examination Survey 2011-2012". www.cdc.gov/nchs/nhanes/index.htm
}
\examples{
# Figure 1 in Rosenbaum (2017)
data(periodontal)
attach(periodontal)
oldpar<-par()
m<-matrix(1:2,1,2)
layout(m,widths=c(1,2))
boxplot(pcteither[z==1]-pcteither[z==0],ylab="Smoker-Control Difference",
        main="(i)",xlab="Matched Pairs",ylim=c(-100,100))
abline(h=0,lty=2)
crosscutplot(cigsperday[z==1],pcteither[z==1]-pcteither[z==0],ylab="Smoker-Control Difference",
    xlab="Cigarettes per Day",main="(ii)",ylim=c(-100,100))
abline(h=0,lty=2)

# Sensitivity analysis in Section 2.3 of Rosenbaum (2017)
y<-pcteither[z==1]-pcteither[z==0]
x<-cigsperday[z==1]
senWilcox(y,gamma=2.76)
# The following is the same as sensitivitymw::senmw(y,gamma=2.77,method="p")
sensitivitymult::senm(pcteither,z,mset,gamma=2.77,inner=.5,trim=2)
# The following is the same as sensitivitymw::senmw(y,gamma=3.5,method="p")
sensitivitymult::senm(pcteither,z,mset,gamma=3.5,inner=.5,trim=2)
# Second evidence factor
crosscut(x,y)
crosscut(x,y,gamma=1.6)

# Note, however, that other statistics report greater insensitivity to
# bias by virtue of having larger design sensitivity:
sensitivitymult::senm(pcteither,z,mset,gamma=3.5,inner=1,trim=4)
sensitivitymult::senm(pcteither,z,mset,gamma=4.2,inner=1,trim=4)
senU(y,m1=4,m2=5,m=5,gamma=2.77)
senU(y,m1=6,m2=8,m=8,gamma=2.77)
senU(y,m1=6,m2=8,m=8,gamma=3.5)
detach(periodontal)
par(oldpar)
}
\keyword{datasets}
