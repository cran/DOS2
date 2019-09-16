\name{senWilcoxExact}
\alias{senWilcoxExact}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Wilcoxon}
\title{
Exact Sensitivity Analysis for Wilcoxon's Signed-rank Statistic
}
\description{
Exact sensitivity analysis for Wilcoxon's signed rank statistic in observational studies. Performs a sensitivity analysis for the one-sided P-value.  The method can be used in
small samples without ties; however, it is primarily of theoretical interest, as the large sample approximation in 'senWilcox' is fine for most samples of practical size.
}
\usage{
senWilcoxExact(d, gamma = 1)
}
\arguments{
  \item{d}{
	A vector of treated-minus-control matched pair differences in outcomes.  There must be no ties in |d| when computing the exact distribution.  If ties are present, use 'senWilcox' instead.
}
  \item{gamma}{
	gamma >= 1 is the value of the sensitivity parameter.
}
}
\details{
  The exact method is discussed in Section 3.12 of "Design of Observational Studies", second edition.  Tables 3.2 and
3.3 of Section 3.5 use these exact calculations.
}
\value{
  The upper bound on the one-sided, upper-tailed P-value testing no treatment effect in the presence of a bias in treatment assignment of at most gamma.
}
\references{

Pagano, M. and Tritchler, D. (1983) <doi:10.1080/01621459.1983.10477990> "On obtaining permutation distributions in polynomial time". Journal of the American Statistical Association, 78, 435-440.

Rosenbaum, P. R. (1987) <doi:10.1093/biomet/74.1.13> "Sensitivity analysis for certain permutation inferences in matched observational studies". Biometrika, 74(1), 13-26.

}
\author{
Paul R. Rosenbaum
}
\note{
   The 'senWilcox' function uses a large-sample approximation, adding confidence intervals and point estimates.
}

\examples{
data(werfel)
d<-werfel$serpc_p-werfel$cerpc_p

# Reproduces the exact one-sided P-value computed in Section 3.9 of
# "Design of Observational Studies".
senWilcoxExact(d,gamma=2)

# Agrees with the usual Wilcoxon procedures when gamma=1.
senWilcoxExact(d,gamma=1)
stats::wilcox.test(d,alternative="greater")

# Reproduces the one-sided confidence interval for gamma=3 in Table 3.3
# of "Design of Observational Studies".
senWilcoxExact(d-0.0935,gamma=3)
senWilcoxExact(d-0.0936,gamma=3)
}

