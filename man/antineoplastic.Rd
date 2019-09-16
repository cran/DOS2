\name{antineoplastic}
\alias{antineoplastic}
\docType{data}
\title{
Biomonitoring of Workers Exposed to Antineoplastic Drugs
}
\description{
Two groups of exposed workers compared to unexposed controls, where
exposed workers prepared antineoplasitic drugs, protected either
only by gloves or by gloves and a laminar hood with vertical air flow.
The outcome is the comet assay applied to blood lymphocytes.  The comet
assay is a measure of damage to DNA.  Data from Kopjar and Garaj-Vrhovac (2001).
Illustrates the concept of evidence factors in Chapter 20 of Design of
Observational Studies, second edition.
}
\usage{data("antineoplastic")}
\format{
  A data frame with 59 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{ID number.  See Tables I and III of Kopjar and Garaj-Vrhovac (2001).}
    \item{\code{age}}{Age in years}
    \item{\code{str}}{Three age strata, cut at the thirds of age.}
    \item{\code{grp}}{Gloves = exposed nurses/doctors protected only
    by gloves, Hood = exposed nurses/doctors protected by both gloves
    and a laminar hood with vertical air flow, Control =
    students and office workers not exposed to antineoplastic
    drugs.}
    \item{\code{tailmoment}}{Tail moment of the comet assay.}
    \item{\code{taillength}}{Tail length of the comet assay mu-m.}
    \item{\code{z1}}{1 if exposed, 0 if control}
    \item{\code{z2}}{1 if Gloves, 0 if Hood, NA if control}
    \item{\code{f2}}{TRUE if in factor 2, FALSE otherwise}
  }
}
\details{
The data set is intended to illustrate evidence factors, comparing
exposed workers to controls, and workers protected only by gloves
to workers with the additional protecting of a laminar
hood with vertical air flow..
}
\source{
From Tables I and III of Kopjar and Garaj-Vrhovac (2001).
}
\references{
Kopjar, N. and Garaj-Vrhovac, V. (2001)
<doi:10.1093/mutage/16.1.71> "Application of the alkaline
comet assay in human biomonitoring for genotoxicity: a study
on Croation medical personnel handling antineoplastic drugs".
Mutagenesis 16, 71-78.
}
\examples{
data(antineoplastic)
attach(antineoplastic)
table(str)
table(grp)
oldpar<-par(mfrow=c(1,2))
boxplot(tailmoment~grp,ylab="Tail Moment",
       main="3 Groups",ylim=c(8,20))
boxplot(tailmoment~z1,names=c("Control","Exposed"),
       main="Factor 1",ylab="Tail Moment",ylim=c(8,20))
boxplot(tailmoment[f2]~z2[f2],names=c("Hood","Gloves"),
       ylab="Tail Moment",main="Factor 2",ylim=c(8,20))
y<-senstrat::hodgeslehmann(tailmoment,z1,str)
# First factor
senstrat::senstrat(y,z1,str,gamma=20)
# Second factor
senstrat::senstrat(y[f2],z2[f2],str[f2],gamma=2.75)
detach(antineoplastic)
par(oldpar)
}
\keyword{datasets}
