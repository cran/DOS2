\name{senU}
\alias{senU}
\concept{Sensitivity analysis}
\concept{Design of Observational Studies}
\concept{Observational study}
\concept{Causal inference}
\concept{Stephenson Test}
\concept{Nonresponders}
\concept{U-statistic}
\title{
Sensitivity Analysis for a New U Statistic
}
\description{
Sensitivity analysis for the new U statistic of Rosenbaum (2011).  For m=m1=m2, this is the test of Stephenson (1981).  The ranks proposed by Stephenson closely approximate the optimal ranks proposed by Conover and Salzburg (1988) for detecting a treatment that has a large effect on a small subpopulation and no effect on most of the population; see Rosenbaum (2007).  The example reproduces some results from Chapter 17 of "Design of Observational Studies", second edition, and the method is discusssed in Section 19.2.
}
\usage{
senU(d, gamma = 1, m = 2, m1 = 2, m2 = 2, conf.int = FALSE,
     alpha = 0.05, alternative = "greater", exact = NULL)
}
\arguments{
  \item{d}{
	A vector of treated-minus-control matched pair differences in outcomes.
}
  \item{gamma}{
	gamma >= 1 is the value of the sensitivity parameter.
}
  \item{m}{
	See m2.
}
  \item{m1}{
	See m2.
}
  \item{m2}{
	If (m,m1,m2) are three integers such that 1 <= m1 <= m2 <= m, then the triple (m,m1,m2) defines a U statistic.  If (m,m1,m2) = (1,1,1), then the U statistic is the sign test statistic.  If (m,m1,m2) = (2,2,2), then it is the U statistic that closely approximates Wilcoxons signed rank test.  If m=m1=m2, then the U statistic is the test of Stephenson (1981).  The general U statistic is discussed in Rosenbaum (2011).
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
\item{exact}{
	If exact is NULL, then exact is set to TRUE if length(d) <= 50, and is set to FALSE if length(d) > 50.  The ranks used by the U statistic involve combinatorial coefficiencts that grow rapidly with increasing sample size.  If exact=TRUE, these ranks are computed exactly using expression (8) in Rosenbaum (2011).  If exact=FALSE, the ranks are computed by an asymptotic approximation that does not involve large combinatorial coefficients, specifically expression (9) in Rosenbaum (2012).
}
}

\value{
\item{pval}{
	The upper bound on the P-value testing no effect in the presence of a bias in treatment assignment of at most gamma.  If the bias in treatment assignment is at most gamma, and if there is no treatment effect, then there is at most an alpah chance that the P-value is less than alpha, this being true for all 0<alpha<1.
}
  \item{estimate}{
	If conf.int=TRUE, the interval of point estimates of an additive treatment effect in the presence of a bias in treatment assigment of at most gamma.  If gamma=1, then you are assuming ignorable treatment assignment or equivalently no unmeasured confounding, so the interval collapses to a point, and that point is the usual Hodges-Lehmann point estimate.
}
  \item{conf.int}{
	If conf.int=TRUE, the a 1-alpha confidence interval for an additive treatment effect in the presence of a bias in treatment assignment of at most gamma.  If gamma=1, then this is the usual confidence interval obtained by inverting the Wilcoxon test, and it would be appropriate in a paired randomized experiment.
}
}
\details{
This method is discussed in Chapter 17 and Section 19.2 of the second edition of "Design of Observational Studies".

The general method is developed in Rosenbaum (2011).
}

\references{
Conover, W. J. and Salsburg, D. S. (1988) <doi:10.2307/2531906> "Locally most powerful tests for detecting treatment effects when only a subset of patients can be expected to" respond" to treatment". Biometrics, 189-196.

Hodges Jr, J. L. and Lehmann, E. L. (1963)
<doi:10.1214/aoms/1177704172> "Estimates of location based on rank tests". The Annals of Mathematical Statistics, 598-611.

Rosenbaum, P. R. (1993) <doi:10.1080/01621459.1993.10476405> "Hodges-Lehmann point estimates of treatment effect in observational studies". Journal of the American Statistical Association, 88(424), 1250-1253.

Rosenbaum, P. R. (2007) <doi:10.1111/j.1541-0420.2007.00783.x> "Confidence intervals for uncommon but dramatic responses to treatment". Biometrics, 63(4), 1164-1171.

Rosenbaum, P. R. (2011) <doi:10.1111/j.1541-0420.2010.01535.x> "A new U statistic with superior design sensitivity in matched observational studies". Biometrics, 67(3), 1017-1027.

Schoket, B., Phillips, D. H., Hewer, A. and Vincze, I. (1991)
<doi:10.1016/0165-1218(91)90084-Y> "32P-postlabelling detection of aromatic DNA adducts in peripheral blood lymphocytes from aluminium production plant workers". Mutation Research/Genetic Toxicology, 260(1), 89-98.

Stephenson, W. R. (1981) <doi:10.1080/01621459.1981.10477749> "A general class of one-sample nonparametric test statistics based on subsamples". Journal of the American Statistical Association, 76(376), 960-966.
}
\author{
Paul R. Rosenbaum
}
\note{
  The test of Stephenson (1981) uses ranks similar to those of Conover and Salzburg (1988) which were designed to have high power when most people are unaffected by treatment, but a small subpopulation is strongly affected; see Rosenbaum (2007).  This is the situation discussed in Chapter 17 of "Design of Observational Studies", second edition; see also Section 19.2.  Even for pair differences d that are Normal with expectation tau and constant variance, the Wilcoxon test tends to exaggerate the degree of sensitivity to unmeasured bias.  Compare the Wilcoxon test and the U statistic with (m,m1,m2) = (5,4,5) in the Normal situation.  In a randomized experiment (gamma=1), the two tests have the same Pitman efficiency.  However, as the number of pairs increases with tau=0.5, the Wilcoxon test has limiting sensitivity to bias of gamma=3.2 while (m,m1,m2) = (5,4,5) has limiting sensitivity 3.9, and (m,m1,m2) = (8,7,8) has limiting sensitivity 5.1.  See Rosenbaum (2011) for specifics.
}

\examples{
data("schoket")
d<-schoket$adductsw-schoket$adductsc

# With the defaults, m=2, m1=2, m2=2, the U-statistic is very
# similar to Wilcoxon's signed rank statistic
senWilcox(d,gamma=1)
senU(d,gamma=1)

# With m=1, m1=1, m2=1, the U-statistic is the sign test
senU(d,gamma=1,m=1,m1=1,m2=1)
prop.test(sum(d>0),length(d),p=.5,
     alternative="greater",correct=FALSE)$p.value

# With m=m1=m2, this is the test of Stephenson (1981) whose ranks are similar to
# those of Conover and Salzburg (1988); see Rosenbaum (2007).

# The calculations that follow reproduce the sensitivity analysis for the
# data of Schoket et al. (1991) in Chapter 17 of "Design of
# Observational Studies", second edition.
senU(d,gamma=1,m=2,m1=2,m2=2)
senU(d,gamma=1,m=5,m1=5,m2=5)

senU(d,gamma=1.5,m=2,m1=2,m2=2)
senU(d,gamma=1.5,m=5,m1=5,m2=5)

senU(d,gamma=1.8,m=2,m1=2,m2=2)
senU(d,gamma=1.8,m=5,m1=5,m2=5)

senU(d,gamma=2,m=2,m1=2,m2=2)
senU(d,gamma=2,m=5,m1=5,m2=5)

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

