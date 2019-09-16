\name{teeth}
\alias{teeth}
\docType{data}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Matched pairs}
\concept{Coherence}
\title{
Smoking and Periodontal Disease.
}
\description{
Data from NHANES 2011-2012 concerning periodonal disease in 441 matched
pairs of smokers and nonsmokers.
}
\usage{data("teeth")}
\format{
  A data frame with 882 observations on the following 4 variables.
  \describe{
    \item{\code{mset}}{Matched pair indicator: 1, 2, ..., 441.}
    \item{\code{smoker}}{Treatment indicator: 1 if current smoker, 0 if never smoker}
    \item{\code{either4up}}{Measure of periodontal disease for upper teeth; see
    Details.}
    \item{\code{either4low}}{Measure of periodontal disease for lower teeth; see
    Details.}
    \item{\code{cigsperday}}{Cigarettes smoked per day.  Zero for nonsmokers.}
  }
}
\details{
Smoking is believed to cause periodontal disease; see Tomar and Asma (2000).
Using more recent data from NHANES 2011-2012, the data describe 441 matched
pairs of a daily smoker and a never smoker.  Daily smokers smoked every day
of the last 30 days.  Never smokers smoked fewer than 100 cigarettes in their
lives, do not smoke now, and had no tobacco use in the previous five days.

Pairs are matched for education, income, age, gender and black race.

Measurements were made for up to 28 teeth, 14 upper, 14 lower, excluding
4 wisdom teeth.  Pocket depth and loss of attachment are two
complementary measures of the degree to which the gums have separated
from the teeth.  Pocket depth and loss of
attachment are measured at six locations on each tooth, providing the
tooth is present.  A measurement at a location was taken to exhibit
disease if it had either a loss of attachement >=4mm or a pocked depth
>=4mm, so each tooth contributes a score from 0 to 6.  Upper and lower
are the number of measurements exhibiting disease for upper and
lower teeth.

This example is from Rosenbaum (2016) where more information may be
found.

The data are closely related to the periodontal data set, where the outcome is univariate.  The teeth data are mentioned in Section 18.9 of "Design of Observational Studies", second edition.
}
\source{
"National Health and Nutrition Examination Survey" (NHANES), 2011-2012. cdc.gov/nchs/nhanes
}
\references{
Rosenbaum, P. R. (2016) <doi:10.1214/16-AOAS942> "Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease".  Annals of Applied Statistics, 10, 1447-1471.

Tomar, S. L. and Asma, S. (2000) <doi:10.1902/jop.2000.71.5.743> "Smoking attributable periodontitis in the United States: Findings from NHANES III". J. Periodont. 71, 743-751.

}
\examples{
data(teeth)
summary(teeth)
# See also the examples in the documentation for 'cohere'.
}
\keyword{datasets}
