\name{cohere}
\alias{cohere}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Scheffe projection}
\title{
Sensitivity Analysis for a Coherent Signed Rank Statistic With Multiplicity Correction
}
\description{
For multiple outcomes in an observation study, computes
a weighted combination of signed rank statistics, one for each outcome,
and performs either a one-sided randomization test or an analysis
of sensitivity to departures from random assignment; see Rosenbaum (1997).
Each matched set contains one treated individual and one control.
The signed rank statistics can be either Wilcoxon's signed rank
statistic or the new U-statistic in Rosenbaum (2011).
The Scheffe method is described in Rosenbaum (2016).  For one outcome,
use the functions 'senWilcox' or 'senU', instead of 'cohere'.
The method is discussed in Sections 5.2.3, 18.2 and 18.3 of
"Design of Observational Studies", second edition.
}
\usage{
cohere(y,z,mpair,w=NULL,gamma=1,m=NULL,m1=NULL,m2=NULL,
                     apriori=FALSE,Scheffe=FALSE,exact=NULL)
}
\arguments{
  \item{y}{ A matrix of responses with no missing data.  Different
  columns of y are different variables.  If present, the column names
  of y are used to label output.}
  \item{z}{ Treatment indicators, z=1 for treated, z=0 for control with
  length(z)==dim(y)[1].
}
  \item{mpair}{ Matched set indicators, 1, 2, ..., sum(z) with
  length(mset)==dim(y)[1].  The vector mset may contain integers
  or may be a factor.
}
  \item{w}{ Vector of weights to be applied to the signed rank
  statistics
  for the several outcomes with length(w)==dim(y)[2].  At least
  one weight must be nonzero.  If w is NULL, then w=c(1,1,...,1)
  is used, as in Rosenbaum (1997).
}
  \item{gamma}{
gamma is the sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.  Setting
\eqn{\Gamma = 1} is equivalent to assuming randomized treatment assignment within
the matched pairs, and it performs a randomization test.
}
  \item{m}{
	See m2.
}
  \item{m1}{
	See m2.
}
  \item{m2}{
	The default sets m, m1 and m2 to NULL, and in this case Wilcoxon's signed rank statistic is used.  Otherwise, the new U-statistic in Rosenbaum (2011) is used, and the quantities (m,m1,m2) define the U-statistic, as they do for a single outcome in the senU() function.  If (m,m1,m2) are three integers such that 1 <= m1 <= m2 <= m, then the triple (m,m1,m2) defines a U statistic.  If (m,m1,m2) = (1,1,1), then the U statistic is the sign test statistic.  If (m,m1,m2) =
(2,2,2), then it is the U statistic that closely approximates Wilcoxons signed rank test.  If m=m1=m2, then the U statistic is the test of Stephenson (1981).  The general U statistic is discussed in Rosenbaum (2011).
}
  \item{apriori}{
  If Scheffe=FALSE and apriori=TRUE, then the weights w are assumed to have been chosen a priori, and a one-sided, uncorrected P-value is reported for gamma=1 or an upper bound on the one-sided, uncorrected P-value is reported for gamma>1. In
either case, this is a Normal approximation based on the central limit
theorem and equals 1-pnorm(deviate).
}
  \item{Scheffe}{
  If Scheffe=TRUE, then the weights w are assumed to have been chosen after
  looking at the data.  In this case, the P-value or P-value bound is
  corrected using Scheffe projections.  The approximate corrected P-value
  or P-value bound is 1-pchisq(max(0,deviate)^2,dim(y)[2]).  See
  Rosenbaum (2016).  The Scheffe correction permits you to look at every possible w, controlling the family-wise error rate.  In particular, a Scheffe correction entitles you to look in both tails, which you do by considering both w and -w.

  See the planScheffe() function and Rosenbaum (2019) for a combination
of an apriori and Scheffe comparisons, as discussed in section 18.3 of
"Design of Observational Studies", second edition.

  If Scheffe=FALSE and apriori=FALSE, then the deviate is returned, but no P-value is given.  This is useful with planScheffe() because the one planned comparison and the infinitely many discovered comparisons have different critical values.
}
\item{exact}{
  exact plays no role if m=NULL, m1=NULL, m2=NULL.  Otherwise, exact plays the same role in cohere() that it plays in senU().  In the new U-statistic in Rosenbaum (2011), there are both exact ranks and approximate ranks, and exact determines which will be used.  Approximate ranks are more appropriate if the sample size is large.  The exact ranks use combinatorial coefficients that become very large when the sample size is large.
	If exact is NULL, then exact is set to TRUE if length(z) <= 50, and is set to FALSE if length(z) > 50.  The ranks used by the U statistic involve combinatorial coefficiencts that grow rapidly with increasing sample size.  If exact=TRUE, these ranks are computed exactly using expression (8) in Rosenbaum (2011).  If exact=FALSE, the ranks are computed by an asymptotic approximation that does not involve large combinatorial coefficients, specifically expression (9) in Rosenbaum (2012).
}
}
\details{
With w=c(1,1,...,1) and apriori=TRUE, this is the coherent Wilcoxon signed rank
statistic in Rosenbaum (1997) and Rosenbaum (2002, Section 9.2).  Essentially, the signed rank statistics for several outcomes are added together.  It is a one-sided test, in which the several outcomes need to pointed in the same direction, so that large values of each column of y signify the expected direction of the treatment effect.  Setting a non-null value for (m,m1,m2) is the analogous statistic but with the ranks in Rosenbaum (2011) instead of Wilcoxon's ranks.

With Scheffe=TRUE, the procedure permits you to test every w that is not
(0,0,...,0), yet control the family-wise error rate.  This is essentially the method in Rosenbaum (2016), except the test-statistics are signed rank statistics rather than M-statistics.

Used in conjunction with the planScheffe() function, you can test one planned w and all possible w's controlling the family-wise error rate.  The size or alpha for the test is shared equitably between the one planned w and all the other w's, so the critical value for the one planned comparison is not greatly increased.  The one-sided critical value for the one planned w is close to the critical value for a two-sided test, although you do test in the opposite tail when you try -w in place of w, albeit with a different critical value.  See Rosenbaum (2019).
}
\value{
\item{deviate }{The upper bound on the standardized deviate
that is used to approximate P-values using the Normal or
chi-square distribution; see apriori and Scheffe in the
arguments.}
\item{aprioriPval }{If Scheffe=FALSE and apriori=TRUE, the
deviate is compared to the upper tail of the Normal
distribution to produce either a P-value for gamma=1
or an upper bound on the P-value for gamma>1.}
\item{ScheffePval }{If Scheffe=TRUE, the
deviate is compared to the upper tail of the
chi-square distribution to produce either
a P-value for gamma=1
or an upper bound on the P-value for gamma>1.}
\item{weights }{The weights possibly relabeled with
colnames of y.}
}
\references{

Rosenbaum, P. R. (1997)<doi:10.2307/2533957> "Signed rank
statistics for coherent predictions". Biometrics, 50, 368-374.

Rosenbaum, P. R. (2002) <doi:10.1007/978-1-4757-3692-2_3>
"Observational Studies" (2nd Edition).  New York: Springer.  See section 9.2.

Rosenbaum, P. R. (2011) <doi:10.1111/j.1541-0420.2010.01535.x> "A new U-statistic with superior design sensitivity in matched observational studies". Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2016) <doi:10.1214/16-AOAS942> "Using Scheffe projections for multiple outcomes in an observational study of smoking and periondontal disease".  Annals of Applied Statistics, 10, 1447-1471.

Rosenbaum, P. R. (2019) <doi:10.1093/biostatistics/kxy055> "Combining planned and discovered comparisons in observational studies".  Biostatistics, to appear.

Scheffe, H. (1953) <doi:10.1093/biomet/40.1-2.87> "A method for judging all contrasts in the analysis
of variance".  Biometrika, 40, 87-104.

}
\author{
Paul R. Rosenbaum.
}
\note{
This method is discussed in chapters 5 and 18 of the second edition of "Design of Observational Studies".

For confidence intervals for individual outcomes, use function senWilcox() or senU().
}

\examples{
# Reproduces parts of Table 5.3 from Design of Observational Studies
data(angristlavy)
attach(angristlavy)
cohere(cbind(avgmath,avgverb),1-z,pair,w=c(1,1),apriori=TRUE)
cohere(cbind(avgmath,avgverb),1-z,pair,w=c(1,0),apriori=TRUE)
cohere(cbind(avgmath,avgverb),1-z,pair,w=c(0,1),apriori=TRUE)
cohere(cbind(avgmath,avgverb),1-z,pair,w=c(1,1),gamma=1.65,apriori=TRUE)

# Uses the technique in Rosenbaum (2019)
# Rejection occurs at gamma=1.5 as 2.054>1.895
cohere(cbind(avgmath,avgverb),1-z,pair,w=c(1,1),gamma=1.5)
planScheffe(2)


detach(angristlavy)

data(teeth)
attach(teeth)
# Coherent Wilcoxon signed rank test
cohere(cbind(either4up,either4low),smoker,mset,apriori=TRUE)

# Sensitivity analysis, gamma=2
cohere(cbind(either4up,either4low),smoker,mset,
       gamma=2,apriori=TRUE)

# Upper teeth only
cohere(cbind(either4up,either4low),smoker,mset,
       w=c(1,0),gamma=2,apriori=TRUE)
# This is the same as the univariate test
y<-either4up[smoker==1]-either4up[smoker==0]
senWilcox(y,gamma=2)

# Try various weights, correcting by Scheffe's method
cohere(cbind(either4up,either4low),smoker,mset,
        w=c(1,2),gamma=2,Scheffe=TRUE)

# Replace Wilcoxon's ranks with the new U-statistic
cohere(cbind(either4up,either4low),smoker,mset,
        w=c(1,2),gamma=2,Scheffe=TRUE,m=8,m1=6,m2=8)
detach(teeth)
}

\keyword{ htest }
\keyword{ robust }
