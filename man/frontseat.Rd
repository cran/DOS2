\name{frontseat}
\alias{frontseat}
\concept{Causal inference}
\concept{Counterclaim}
\concept{Observational Study}
\concept{Mediation}
\docType{data}
\title{
Safety Belts in Vehicle Crashes
}
\description{
Data from the US Fatality Analysis Reporting System (FARS) in 2010-2011, as discussed in Rosenbaum (2015) and in the "Design of Observational Studies", second edition, Chapter 7.  The data concern crashes in which a driver and a right front passenger were present, following Evans (1986).  The data compare the injuries of the driver and passenger, and are particularly interesting when their safety belt use is different.  The example illustrates the analysis of a counterclaim.
}
\usage{data("frontseat")}
\format{
  A data frame with 17014 observations on the following 7 variables.
  \describe{
    \item{\code{restraint}}{Saftey belt use by (Driver,Passenger), where n=unbelted, ls=lap-shoulder belt.  A factor with levels \code{ls.ls} \code{ls.n} \code{n.ls} \code{n.n}.  Here, ls.n means the driver was belted and the passenger was not.}
    \item{\code{injury}}{Injury of (Driver,Passenger).}
    \item{\code{injurydif}}{Difference in injury scores, driver-minus-passenger, from -4 to 4.  A score of -4 means the driver was uninjured, but the passenger died.}
    \item{\code{ejection}}{Ejection from the vehicle of the (Driver,Passenger).}
    \item{\code{ejectiondif}}{1 if the driver was ejected but the passenger was not, -1 if the passenger was ejected but the driver was not, 0 if their fates were the same.}
    \item{\code{gender}}{Genders of the (Driver,Passenger).}
    \item{\code{agedif}}{Difference in ages, driver-minus-passenger.}
  }
}
\details{
This example is discussed in "Design of Observational Studies", second edition, Chapter 7.

Details are given in Rosenbaum (2015).  A crash, perhaps involving several vehicles, is recorded in FARS only if there is at least
one fatality, creating issues of ascertainment (Fisher 1934) that do not affect tests of the hypothesis of no effect,
but that do affect estimation.  Only tests of no effect are considered in this example.
}
\source{
Rosenbaum (2015)
}
\references{
Evans, L. (1986) <doi:10.1016/0001-4575(86)90007-2> "The
Effectiveness of Safety Belts in Preventing Fatalities".
Accident Analysis and Prevention, 18, 229–241.

Fisher, R.A. (1934) <doi:10.1111/j.1469-1809.1934.tb02105.x> "The Effect of Methods of Ascertainment Upon the Estimation of Frequencies". Annals of Eugenics, 6(1), 13-25.

Imai, K., Keele, L., Yamamoto, T. (2010) <doi:10.1214/10-STS321> "Identification, Inference and Sensitivity Analysis for Causal Mediation Effects". Statistical Science, 25, 51–71.

Rosenbaum, P.R. (2015) <doi:10.1080/01621459.2015.1054489> "Some Counterclaims Undermine Themselves in Observational Studies". Journal of the American Statistical Association, 110:512, 1389-1398.
}
\examples{
data(frontseat)
attach(frontseat)
use<-(!is.na(injurydif))
# Compare with Table 1 in Rosenbaum (2015), case ls.n
table(restraint[use])
use<-use&(restraint=="ls.n")
2*sensitivitymv::senmv(-injurydif[use],gamma=5,
       trim=1,lambda=.99)$pval
2*sensitivitymv::senmv(-injurydif[use],gamma=5.5,
       trim=1,lambda=.99)$pval
2*sensitivitymv::senmv(-injurydif[use],gamma=6,
       trim=1,lambda=.99,inner=.25)$pval
2*sensitivitymv::senmv(-injurydif[use],gamma=6.5,
       trim=1,lambda=.99,inner=.25)$pval

# Counterclaim analysis, one ejected individual
# Compare with Table 2 in Rosenbaum (2015), case ls.n
table(ejection,ejectiondif)
use<-use&(!is.na(ejectiondif))&(ejectiondif!=0)
sum(use)
2*sensitivitymv::senmv(-injurydif[use],gamma=9,
       trim=1,lambda=.99)$pval
2*sensitivitymv::senmv(-injurydif[use],gamma=11,
       trim=1,lambda=.99,inner=.25)$pval
detach(frontseat)
}
\keyword{datasets}
