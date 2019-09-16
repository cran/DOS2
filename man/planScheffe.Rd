\name{planScheffe}
\alias{planScheffe}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Matching}
\concept{Sheffe projection}
\concept{Multiple comparisons}
\concept{Coherence}
\title{
Combining One Planned Comparison and a Scheffe Correction For All Comparisons.
}
\description{
The function planScheffe() computes the critical values for a level
alpha test that combines one planned linear combination of
a K-dimensional multivariate Normal outcome and consideration of
all possible combinations correcting for multiple testing using
a Scheffe projection.  The function is the same as the planScheffe()
function in the sensitivitymult package, but the examples are
different.  The method is discussed in Section 18.3 of "Design of
Observational Studies", second edition, and in Rosenbaum (2019).
}
\usage{
planScheffe(K, alpha = 0.05)
}
\arguments{
  \item{K}{ An integer >=2 giving the number of outcomes to be compared.}
  \item{alpha}{ The level of the test, with 0 < alpha < 1.}
}
\details{
This method is discussed in section 18.3 of the second edition of "Design of Observational Studies".

Although the calculation uses the multivariate Normal distribution,
a typical application uses K test statistics that are
asymptotically Normal.

The method is based on Rosenbaum (2019).  The example below
reproduces some of the comparisons in that manuscript.
}
\value{
\item{critical }{ critical is a vector with two elements, a and c.
The null hypothesis is rejected at level alpha if either the
Normal deviate for the planned comparison is >= a or if the square
of the Normal deviate for any comparison is >= b.  Then the
probability of a false rejection is <= alpha.}
\item{alpha }{ alpha is a vector with three elements, a, c and
joint.  The value of joint should equal the input value of alpha
aside from numerical errors of computation: it is the probability
of a false rejection using the joint test that rejects if
either of the two critical values in critical is exceeded.
In contrast, a is the probability that the planned deviate
will be >= critical[1] when the null hypothesis is true.
Also, c is the probability that at least one comparison
will have a squared deviate >= critical[2] when the
null hypothesis is true.}
}
\references{

Rosenbaum, P. R. (2016) <doi:10.1214/16-AOAS942> "Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease".  Annals of Applied Statistics, 10, 1447-1471.

Rosenbaum, P. R. (2019) <doi:10.1093/biostatistics/kxy055> "Combining planned and discovered comparisons
in observational studies".  Biostatistics, to appear.

Scheffe, H. (1953) <doi:10.1093/biomet/40.1-2.87> "A method for judging all contrasts in the analysis
of variance".  Biometrika, 40, 87-104.

}
\author{
Paul R. Rosenbaum.
}
\note{
The method is based on Rosenbaum (2019).

The functions "cohere" may be used to calculate
the standardized deviates that are compared to the critical values
from 'planScheffe'.  The function "cohere" has options for an a priori
comparison or consideration of all possible comparisons with a
Scheffe correction.  The function "planScheffe" provides a third
option: one planned comparison plus all possible comparisons.

See also 'planScheffe' in the 'sensitivitymult' package.
}


\examples{
data(teeth)
attach(teeth)

planScheffe(2,alpha=0.05)

# Planned comparison w=c(1,1)
cohere(cbind(either4up,either4low),smoker,mset,
             w=c(1,1),gamma=2,m=8,m1=6,m2=8)

# Discovered comparison emphasizing upper teeth
cohere(cbind(either4up,either4low),smoker,mset,
        w=c(1,3),gamma=2,m=8,m1=6,m2=8)
3.465038^2 #squared deviate

# Both deviates lead to rejection, because
# 3.291909 >= 1.894915
# and 3.465038^2 = 12.00649 >= 7.077349

detach(teeth)
}

\keyword{ htest }
\keyword{ robust }
