\name{costa}
\alias{costa}
\concept{Observational study}
\concept{Propensity score}
\concept{Matched pairs}
\concept{Causal inference}
\concept{Multivariate matching}
\concept{Optimal matching}
\concept{Design of Observational Studies}
\docType{data}
\title{
Welding and DNA-Protein Crosslinks
}
\description{
This data set is from Costa et al. (1993) and it describes 21 welders and 26 potential controls.
All are men.  The outcome is a measure of genetic damage;
specifically, dpc is a measure of DNA-protein cross-links.  There are 3 covariates,
age, race and smoking.  This tiny example is used to illustrate the concepts of
multivariate matching in Chapter 9 of "Design of Observational Studies", second edition.  The
example is useful because its tiny size permits close inspection of the details
of multivariate matching, but its small sample size and limited number of covariates
make it highly atypical of matching in observational studies.
}
\usage{data("costa")}
\format{
  A data frame with 47 observations on the following 6 variables.
  \describe{
    \item{\code{subject}}{Within group ID number.}
    \item{\code{age}}{Age in years.}
    \item{\code{race}}{AA=African-American, C=Caucasian}
    \item{\code{smoker}}{Y=yes, N=no}
    \item{\code{welder}}{Y=yes/treated, N=no/control}
    \item{\code{dpc}}{DNA-Protein Cross-links (percent)}
  }
}
\source{
The data are from Costa et al. (1993).  The data are used as a tiny example in Chapter 9 of "Design of Observational Studies", second edition.
}
\references{
Costa, M., Zhitkovich, A. and Toniolo, P. (1993)
<https://cancerres.aacrjournals.org/content/53/3/460> "DNA-protein cross-links in welders: molecular implications". Cancer research, 53(3), 460-463.
}
\examples{
data(costa)
boxplot(costa$dpc~costa$welder,
  xlab="Control (N) or Welder (Y)",
  ylab="DNA-Protein Cross-links Percent")
}
\keyword{datasets}
