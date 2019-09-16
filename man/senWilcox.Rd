\name{senWilcox}
\alias{senWilcox}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Wilcoxon}
\title{
Sensitivity Analysis for Wilcoxon's Signed-rank Statistic
}
\description{
Sensitivity analysis for Wilcoxon's signed rank statistic in observational studies.
Performs a sensitivity analysis for the P-value, the Hodges-Lehmann estimate of
an additive treatment effect, and the confidence interval for the effect.
The example is from Section 3.5 of "Design of Observational Studies", 2nd ed.,
and it also illustrates the amplification in Section 3.6 using the
'amplify' function from the 'sensitivitymult' package.
}
\usage{
senWilcox(d, gamma = 1, conf.int = FALSE, alpha = 0.05, alternative = "greater")
}
\arguments{
  \item{d}{
	A vector of treated-minus-control matched pair differences in outcomes.
}
  \item{gamma}{
	gamma >= 1 is the value of the sensitivity parameter.
}
  \item{conf.int}{
	If conf.int=TRUE, the a 1-alpha confidence interval and an interval of point estimates is returned
	in addition to the P-value testing no treatment effect.
}
  \item{alpha}{
	Coverage rate of the confidence interval.  With probability at least 1-alpha, the confidence interval
	will cover the treatment effect providing the bias in treatment assignment is at most gamma.
}
  \item{alternative}{
	If alternative = "greater" or alternative = "less", then one-sided tests and intervals are returned.
	If alternative = "twosided", then both one sided tests are done, with the smaller P-value doubled
	to yield a two-sided P-value.  If alternative = "twosided", the confidence interval is the intersection
	of two one-sided 1-alpha/2 confidence intervals.
}
}
\details{
  The senWilcox function uses a large sample Normal approximation to the distribution of Wilcoxon's signed rank statistic.  When gamma=1, it should agree with the wilcox.test() function in the stats package with exact=FALSE and correct=FALSE.

  The example reproduces the example of the large-sample approximation in Section 3.5 of "Design of Observational Studies".  Note that the confidence intervals in Table 3.3 of that book are exact, not approximate, so they are slightly different.
}
\value{
\item{pval}{
	The upper bound on the P-value testing no effect in the presence of a bias in treatment assignment of at most gamma.  If the bias in treatment assignment is at most gamma, and if there is no treatment effect, then there is at most an alpha chance that the P-value is less than alpha, this being true for all 0<alpha<1.
}
  \item{estimate}{
	If conf.int=TRUE, the interval of point estimates of an additive treatment effect in the presence of a bias in treatment assigment of at most gamma.  If gamma=1, then you are assuming ignorable treatment assignment or equivalently no unmeasured confounding, so the interval collapses to a point, and that point is the usual Hodges-Lehmann point estimate.
}
  \item{conf.int}{
	If conf.int=TRUE, the a 1-alpha confidence interval for an additive treatment effect in the presence of a bias in treatment assignment of at most gamma.  If gamma=1, then this is the usual confidence interval obtained by inverting the Wilcoxon test, and it would be appropriate in a paired randomized experiment.
}
}

\references{

Hodges Jr, J. L. and Lehmann, E. L. (1963)
<doi:10.1214/aoms/1177704172> "Estimates of location based on rank tests". The Annals of Mathematical Statistics, 598-611.

Hollander, M., Wolfe, D. and Chicken, E. (2013)
<doi:10.1002/9781119196037> "Nonparametric Statistical Methods". (3rd edition) New York: John Wiley.

Lehmann, E. L. (1975). "Nonparametrics" <https://www.springer.com/gp/book/9780387352121> San Francisco: Holden-Day.  Reprinted by Prentice-Hall and Springer.

Rosenbaum, P. R. (1987) <doi:10.1093/biomet/74.1.13> "Sensitivity analysis for certain permutation inferences in matched observational studies". Biometrika, 74(1), 13-26.

Rosenbaum, P. R. (1993) <doi:10.1080/01621459.1993.10476405> "Hodges-Lehmann point estimates of treatment effect in observational studies". Journal of the American Statistical Association, 88(424), 1250-1253.

Rosenbaum, P. R. (2002) <doi:10.1007/978-1-4757-3692-2_3>  "Observational Studies", Second edition.  New York: Springer.  Wilcoxon's test is discussed
in Section 4.3.3.

Rosenbaum, P. R. (2007) <doi:10.1111/j.1541-0420.2006.00717.x> "Sensitivity Analysis for M Estimates, Tests, and Confidence Intervals in Matched Observational Studies". Biometrics, 63(2), 456-464. R-packages 'sensitivitymult', 'sensitivitymv', 'sensitivityfull', 'sensitivitymw'

Rosenbaum, P. R. (2010) <doi:10.1198/jasa.2010.tm09570> "Design sensitivity and efficiency in observational studies". Journal of the American Statistical Association, 105(490), 692-702.

Rosenbaum, P. R. (2011) <doi:10.1111/j.1541-0420.2010.01535.x> "A new U statistic with superior design sensitivity in matched observational studies". Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2013) <doi:10.1111/j.1541-0420.2012.01821.x> "Impact of multiple matched controls on design sensitivity in observational studies". Biometrics, 69(1), 118-127. R-packages
'sensitivitymv', 'sensitivitymult' and 'sensitivityfull'

Rosenbaum, P. R. (2014) <doi:10.1080/01621459.2013.879261> "Weighted M statistics with superior design sensitivity in matched observational studies with multiple controls". Journal of the American Statistical Association, 109(507), 1145-1158.  R-package 'sensitivitymw'

Rosenbaum, P. R. (2015) <doi:10.1080/01621459.2014.960968>
"Bahadur efficiency of sensitivity analyses in observational studies". Journal of the American Statistical Association, 110(509), 205-217.
}
\author{
Paul R. Rosenbaum
}
\note{
  The Wilcoxon test is not the best test for use in sensitivity analyses -- it tends to exaggerate how sensitive a study is to bias, saying it is more sensitive than it truly is.  Other, less familiar, test statistics avoid this issue; see Rosenbaum (2010, 2011, 2013, 2014, 2015).  Learn about this using a 'shinyapp',
at rosenbap.shinyapps.io/learnsenShiny/
}

\examples{
data(werfel)
d<-werfel$serpc_p-werfel$cerpc_p

# Reproduces the approximate one-sided P-value computed
# in Section 3.5 of "Design of Observational Studies".
senWilcox(d,gamma=3)

# Amplification in Section 3.6
sensitivitymult::amplify(3,5)
sensitivitymult::amplify(3,c(5,5.8,7))

# Reproduces parts of Tables 4.3 and 4.4 in Rosenbaum (2002)
data(lead)
senWilcox(lead$dif,gamma=1,conf.int=TRUE,alternative="twosided")
senWilcox(lead$dif,gamma=2,conf.int=TRUE,alternative="twosided")

# Agrees with the usual Wilcoxon procedures when gamma=1.
senWilcox(d,gamma=1,conf.int=TRUE,alternative="twosided")
stats::wilcox.test(d,conf.int=TRUE,exact=FALSE,correct=FALSE)
}

