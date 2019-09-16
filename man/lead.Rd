\name{lead}
\alias{lead}
\concept{Observational study}
\concept{Sensitivity analysis}
\concept{Matched pairs}
\concept{Causal inference}
\docType{data}
\title{
Lead in Children
}
\description{
Data from Morton et al. (1982) concerning exposed children whose fathers worked in a battery plant where lead was used in the manufacture of batteries.  Exposed children were matched to controls for age and neighborhood. For exposed children, also given are the father's level of exposure to lead at work (level) and the father's hygiene upon leaving the battery plant at the end of the day.  This example is used in Chapters 7 and 9 of Rosenbaum (2017).
}
\usage{data("lead")}
\format{
  A data frame with 33 observations on the following 6 variables.
  \describe{
    \item{\code{control}}{Blood lead level for the control, micrograms of lead per decaliter of whole blood.}
    \item{\code{exposed}}{Blood lead level for the exposed/treated child, micrograms of lead per decaliter of whole blood.}
    \item{\code{level}}{Father's level of exposure to lead: a factor with levels \code{high} \code{low} \code{medium}}
    \item{\code{hyg}}{Father's hygiene before going home: a factor with levels \code{good} \code{mod} \code{poor}}
    \item{\code{both}}{A factor built from level and hyg: a factor with levels \code{high.ok} \code{high.poor} \code{low} \code{medium}}
    \item{\code{dif}}{Exposed-minus-control pair difference in blood lead levels.}
  }
}
\details{
The data were assembled from two published tables in Morton et al. (1982).  One matched pair with no control is omitted here.  Small ambiguities in assembling a complete data set from two tables were resolved by a throw of the dice; however, it is a reasonable example to illustrate statistical methods.  One table described the exposed-versus-control matched pairs, and these are as in the paper.  The second table described the exposed individuals, their level of exposure and their hygiene, and again these are as in the paper.  The two tables were linked using the blood lead level of the exposed children, and a couple of ties in these lead levels made for small ambiguities about which control responses belong with which hygienes and exposure levels for the exposed children.  See, for instance, rows 18 and 22, where both exposed children have blood lead level 34.
}
\source{
Data are from Morton et al. (1982).  They were used as an example in Rosenbaum (1991, 2002, 2011, 2017).
}
\references{
Morton, D. E., Saah, A. J., Silberg, S. L., Owens, W. L., Roberts, M. A., & Saah, M. D. (1982) <doi:10.1093/oxfordjournals.aje.a113336>  "Lead absorption in children of employees in a lead-related industry". American Journal of Epidemiology, 115(4), 549-555.

Rosenbaum, P. R. (1991) <doi:10.1214/aos/1176348141> "Some poset statistics". The Annals of Statistics, 19(2), 1091-1097.


Rosenbaum, P. R. (2002) <doi:10.1007/978-1-4757-3692-2_3> "Observational Studies" (2nd edition).  New York: Springer.  Section 4.3.

Rosenbaum, P. R. (2011) <doi:10.1111/j.1541-0420.2010.01535.x> "A new U statistic with superior design sensitivity in matched observational studies". Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2017) <https://www.hup.harvard.edu/catalog.php?isbn=9780674975576> "Observation and Experiment: An Introduction to Causal Inference".  Cambridge, MA: Harvard University Press. Chapters 7 and 9.
}
\examples{
data(lead)
# Reproduces parts of Table 2 in Rosenbaum (2011)
senU(lead$dif,gamma=5.8,m=8,m1=5,m2=8)
senU(lead$dif,gamma=5,m=5,m1=4,m2=5)

# m=2, m1=2, m2=2 is the U-statistic that closely
# resembles Wilcoxon's signed rank test.  Note
# that the results are almost the same.
senWilcox(lead$dif,gamma=5) # In Table 2
senU(lead$dif,gamma=5,m=2,m1=2,m2=2)
}
\keyword{datasets}
