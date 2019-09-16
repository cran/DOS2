\name{pinto}
\alias{pinto}
\concept{Observational study}
\concept{Sensitivity analysis}
\concept{Matched pairs}
\concept{Causal inference}
\concept{Doses of treatment}
\concept{Design sensitivity}
\concept{Design of Observational Studies}
\docType{data}
\title{
Welding and DNA-Protein Crosslinks
}
\description{
This data set is from Pinto et al. (2000) and it describes 22 professional painters and 22 controls matched for age.  All are men.  The outcome is a measure of genetic damage, namely the frequency of micronuclei in 3000 oral epithelial cells scraped from the cheek, recorded as micronuclei per 1000 cells.  The data are used as an example in Chapter 5 and 18 of "Design of Observational Studies", where the data illustrate the subtle relationship between doses of treatment and sensitivity to unmeasured biases.  The dose of treatment is years of work as a professional painter, from 1.6 years to 40 years for painters, and it is highly correlated with age for painters, but it is zero for controls of all ages.
}
\usage{data("pinto")}
\format{
  A data frame with 44 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{ID number from Pinto et al. (2000).}
    \item{\code{pair}}{Pair number, 1 to 22.  The pair numbers are different from "Design of Observational Studies", where pairs were sorted by doses.}
    \item{\code{group}}{painter or control}
    \item{\code{age}}{age in years}
    \item{\code{years}}{years of work as a professional painter, 0 for controls}
    \item{\code{longEx}}{Painter in this pair had long exposure, TRUE if years>=4 years.}
    \item{\code{mn}}{Micronuclei per 1000 cells}
  }
}
\source{
The data are from Pinto et al. (2000).  The data are used as an example in Chapter 5 of
Design of Observational Studies.  Chapter 18 shows that the pattern seen in this example is expected in theory, namely focusing on high-dose pairs makes the study more insensitive to unmeasured biases, despite the loss of sample size.
}
\references{
Pinto, D., J. M. Ceballos, G. Garcia, P. Guzman, L. M. Del Razo, E. Vera, H. Gomez, A. Garcia, and M. E. Gonsebatt (2000)
<doi:10.1016/S1383-5718(00)00024-3> "Increased cytogenetic damage in outdoor painters". Mutation Research/Genetic Toxicology and Environmental Mutagenesis 467, 105-111.

Rosenbaum, P. R. (2003) <doi:10.1093/biostatistics/4.1.1> "Does a dose response relationship reduce sensitivity to hidden bias?". Biostatistics, 4, 1-10.

Rosenbaum, P. R. (2004) <doi:10.1093/biomet/91.1.153> "Design sensitivity in observational studies". Biometrika, 91, 153-164.  Does design sensitivity calculations with doses of treatment.

}
\examples{
data(pinto)
oldpar<-par(mfrow=c(1,3))
attach(pinto)
boxplot(mn~group,ylim=c(0,6),main="All",ylab="Micronuclei")
boxplot(mn[!longEx]~group[!longEx],ylim=c(0,6),main="Short Ex",ylab="Micronuclei")
boxplot(mn[longEx]~group[longEx],ylim=c(0,6),main="Long Ex",ylab="Micronuclei")

# Calculations in Table 5.5 of Design of Observational Studies (2010)
d<-mn[group=="painter"]-mn[group=="control"] # 22 pair differences
senWilcox(d,gamma=1)
senWilcox(d,gamma=2) # sensitive to gamma=2
senWilcox(d,gamma=3.3)
dLong<-d[longEx[group=="painter"]] # 12 pairs with long exposure
senWilcox(dLong,gamma=3.3) # insensitive to gamma=3.3
par(oldpar)
}
\keyword{datasets}
