\name{fine}
\alias{fine}
\concept{Fine balance}
\concept{Observational study}
\concept{Causal inference}
\concept{Propensity score}
\concept{Optimal matching}
\concept{Design of Observational Studies}
\title{
Expand a Distance Matrix for Matching with Fine Balance.
}
\description{
In optimal pair matching with fine balance, expand a distance matrix to become a square matrix to enforce fine balance.  The method is discussed in Chapter 11 of "Design of Observational Studies", second edition, and it is conceptually the simplest way to implement fine balance; therefore, it remains very useful for teaching and for self-study.  See details.
}
\usage{
fine(dmat, z, f, mult = 100)
}
\arguments{
  \item{dmat}{
	A distance matrix with one row for each treated individual and one column for each control.  Often, this is either the Mahalanobis distance based on covariates, 'mahal', or else a robust variant produced by 'smahal'.  The distance matrix dmat may have been penalized by 'addalmostexact' or 'addcaliper'.  An error will result unless dmat has more columns than rows.
}
  \item{z}{
	z is a vector that is 1 for a treated individual and 0 for a control.  The number of treated subjects, sum(z), must equal the number of rows of dmat, and the number of potential controls, sum(1-z), must equal the number of columns of dmat.
}
  \item{f}{
	A factor or vector to be finely balanced.  Must have length(f)=length(z).
}
  \item{mult}{
  A positive number, mult>0.  Determines the penalty used to enforce fine balance as max(dmat)*mult.  The distance matrix dmat may have been penalized by 'addalmostexact' or 'addcaliper', and in this case it makes sense to set mult=1 or mult=2, rather than the default, mult=100.  If dmat is already penalized, taking mult>1 creates a larger penalty for deviations from fine balance than the exisiting penalties.
  }
}
\value{
  An expanded, square distance matrix with "extra" treated units for use in optimal pair matching.  Any control paired with an "extra" treated unit is discarded, as are the "extra" treated units.
}

\details{
  The method is discussed in Chapter 11 of "Design of Observational Studies", second edition, and it is conceptually the simplest way to implement fine balance.  However, the expanded distance matrix can become quite large, so this simplest method is not the most efficient method in its use of computer storage.  A more compact implementation uses minimum cost flow in a network (Rosenbaum 1989, Section 3.2).  Additionally, there are several extensions of fine balance, including near-fine balance (Yang et al. 2012, in Yu's 'DiPs' package), fine balance for several covariates via integer programming (Zubizarreta 2012, 'designmatch' R-package), and refined balance (Pimentel et al. 2015, 'rcbalance' R-package).  Ruoqi Yu's 'bigmatch' R-package implements fine balance and near-fine balance in very large matching problems.
}

\note{
Fine balance is discussed in chapter 11 of "Design of Observational Studie", second edition.

The matching functions in the 'DOS2' package are aids to instruction or self-instruction while reading "Design of Observational Studies".  As in the book, these functions break the task of matching into small steps so they are easy to understand in depth. In practice, matching entails a fair amount of book-keeping best done by a package that automates these tasks.  For fine balance and similar methods, use
R packages 'rcbalance', 'DiPs', 'designmatch' or 'bigmatch'.  Section 14.10 of "Design of Observational Studies", second edition, discusses and compares several R packages for optimal matching.
}

\references{
Hansen, B. B. and Klopfer, S. O. (2006)
<doi:10.1198/106186006X137047> "Optimal full matching and related designs via network flows". Journal of Computational and Graphical Statistics, 15(3), 609-627.  The method implemented in Hansen's 'optmatch' package.

Hansen, B. B. (2007) <https://www.r-project.org/conferences/useR-2007/program/presentations/hansen.pdf> "Flexible, optimal matching for observational studies". R News, 7, 18-24. Discusses Hansen's 'optmatch' package.

Pimentel, S. D., Kelz, R. R., Silber, J. H. and Rosenbaum, P. R. (2015) <doi:10.1080/01621459.2014.997879> "Large, sparse optimal matching with refined covariate balance in an observational study of the health outcomes produced by new surgeons". Journal of the American Statistical Association, 110, 515-527.  Introduces an extension of fine balance called refined balance that is implemented in Pimentel's package
'rcbalance'.

Pimentel, S. D. (2016) "Large, Sparse Optimal Matching with R Package rcbalance" <https://obsstudies.org/large-sparse-optimal-matching-with-r-package-rcbalance/>  Observational Studies, 2, 4-23. Discusses and illustrates the use of Pimentel's 'rcbalance' package.

Rosenbaum, P. R. (1989). "Optimal matching for observational studies" <doi:10.1080/01621459.1989.10478868> Journal of the American Statistical Association, 84(408), 1024-1032.  Discusses and illustrates fine balance using minimum cost flow in a network in section 3.2.

Rosenbaum, P. R., Ross, R. N. and Silber, J. H. (2007)
<doi:10.1198/016214506000001059> "Minimum distance matched sampling with fine balance in an observational study of treatment for ovarian cancer". Journal of the American Statistical Association, 102, 75-83.  Discusses and illustrates fine balance using optimal assignment.


Yang, D., Small, D. S., Silber, J. H. and Rosenbaum, P. R. (2012)
<doi:10.1111/j.1541-0420.2011.01691.x> "Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes". Biometrics, 68, 628-636.  Extension of fine balance useful when fine balance is infeasible.  Comes as close as possible to fine balance.  Implemented as part of the 'rcbalance' and 'DiPs' packages.

Yu, Ruoqi and Rosenbaum, P.R., (2019a) <doi:10.1111/biom.13098> "Directional penalties for optimal matching in observational studies". Biometrics, to appear,   Describes the method in Yu's 'DiPs' package.

Yu, Ruoqi, Silber, J.H. and Rosenbaum, P.R. (2019b)
<https://www.imstat.org/journals-and-publications/statistical-science/statistical-science-future-papers/> "Matching methods for observational studies derived from large administrative databases". Stat Sci., to appear.  Describes the method in Yu's 'bigmatch' package.

Zubizarreta, J. R., Reinke, C. E., Kelz, R. R., Silber, J. H. and Rosenbaum, P. R. (2011) <doi:10.1198/tas.2011.11072> "Matching for several sparse nominal variables in a case-control study of readmission following surgery". The American Statistician, 65(4), 229-238.  Combines near-exact matching with fine balance for the same covariate.

Zubizarreta, J. R. (2012) <doi:10.1080/01621459.2012.703874> "Using mixed integer programming for matching in an observational study of kidney failure after surgery". Journal of the American Statistical Association, 107, 1360-1371.  Extends the concept of fine balance using integer programming.  Implemented in R in the 'designmatch' package.
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
round(dmat[,1:6],2) # Compare with Table 9.5 in "Design of Observational Studies"
# Impose propensity score calipers
prop<-glm(z~age+aa+smoker,family=binomial)$fitted.values # propensity score
# Mahalanobis distanced penalized for violations of a propensity score caliper.
# This version is used for numerical work.
dmat<-addcaliper(dmat,z,prop,caliper=.5)
round(dmat[,1:6],2) # Compare with Table 9.5 in "Design of Observational Studies"
# Because dmat already contains large penalties, we set mult=1.
dmat<-fine(dmat,z,aa,mult=1)
dmat[,1:6] # Compare with Table 11.1 in "Design of Observational Studies"
dim(dmat) # dmat has been expanded to be square by adding 5 extras, here numbered 48:52
# Any control matched to an extra is discarded.
\dontrun{
# Find the minimum distance match within propensity score calipers.
optmatch::pairmatch(dmat)
# Any control matched to an extra is discarded.  For instance, the optimal match paired
# extra row 48 with the real control in column 7 to form matched set 1.22, so that control
# is not part of the matched sample.  The harmless warning message from pairmatch
# reflects the divergence between the costa data.frame and expanded distance matrix.
}
# Conceptual versions with infinite distances for violations of propensity caliper.
dmat[dmat>20]<-Inf
round(dmat[,1:6],2) # Compare with Table 11.1 in "Design of Observational Studies"
}

