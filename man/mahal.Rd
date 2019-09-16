\name{mahal}
\alias{mahal}
\concept{Mahalanobis distance}
\concept{Observational study}
\concept{Causal inference}
\concept{Propensity score}
\concept{Optimal matching}
\title{
Mahalanobis Distance Matrix for Optimal Matching
}
\description{
Computes a Mahalanobis distance matrix between treated individuals and potential controls; see Rubin (1980) and Rosenbaum and Rubin (1985).  The method is discussed in Section 9.3 of "Design of Observational Studies", second edition.
}
\usage{
mahal(z, X)
}
\arguments{
  \item{z}{
	z is a vector that is 1 for a treated individual and 0 for a control.
}
  \item{X}{
	A matrix of continuous or binary covariates.  The number of rows of X must equal the length of z.
}
}
\value{
  The distance matrix has one row for each treated individual (z=1) and one column for each potential control (z=0).  The row and column names of the distance matrix refer to the position in z, 1, 2, ..., length(z).
}
\note{
The method is discussed in Section 9.3 of "Design of Observational Studies", second edition.

The matching functions in the 'DOS2' package are aids to instruction or self-instruction while reading "Design of Observational Studies".  As in the book, these functions break the task of matching into small steps so they are easy to understand in depth. In practice, matching entails a fair amount of book-keeping best done by a package that automates these tasks.  Consider
R packages 'optmatch', 'rcbalance', 'DiPs', 'designmatch' or 'bigmatch'.  Section 14.10 of "Design of Observational Studies", second edition, discusses and compares several R packages for optimal matching.
}

\references{

Hansen, B. B. and Klopfer, S. O. (2006) <doi:10.1198/106186006X137047> "Optimal full matching and related designs via network flows". Journal of computational and Graphical Statistics, 15(3), 609-627. ('optmatch' package)

Hansen, B. B. (2007) <https://www.r-project.org/conferences/useR-2007/program/presentations/hansen.pdf> "Flexible, optimal matching for observational studies". R News, 7, 18-24. ('optmatch' package)

Rosenbaum, P. R. and Rubin, D. B. (1985)
<doi:10.1080/00031305.1985.10479383> "Constructing a control group using multivariate matched sampling methods that incorporate the propensity score". The American Statistician, 39, 33-38.

Rubin, D. B. (1980) <doi:10.2307/2529981> "Bias reduction using Mahalanobis-metric matching". Biometrics, 36, 293-298.
}
\author{
Paul R. Rosenbaum
}

\examples{
data(costa)
z<-1*(costa$welder=="Y")
aa<-1*(costa$race=="A")
smoker=1*(costa$smoker=="Y")
age<-costa$age
x<-cbind(age,aa,smoker)
dmat<-mahal(z,x)
# Mahalanobis distances
round(dmat[,1:6],2)
# Compare with Table 9.5 in "Design of Observational Studies", 2nd ed.
# Impose propensity score calipers
prop<-glm(z~age+aa+smoker,family=binomial)$fitted.values # propensity score
# Mahalanobis distanced penalized for violations of a propensity score caliper.
# This version is used for numerical work.
dmat<-addcaliper(dmat,z,prop,caliper=.5)
round(dmat[,1:6],2)
# Compare with Table 9.5 in "Design of Observational Studies", 2nd ed.
\dontrun{
# Find the minimum distance match within propensity score calipers.
# You must load the 'optmatch' package to try this example
optmatch::pairmatch(dmat,data=costa)
}
# Conceptual versions with infinite distances for violations of propensity caliper.
dmat[dmat>20]<-Inf
round(dmat[,1:6],2)
# Compare with Table 9.5 in "Design of Observational Studies", 2nd ed.
}

