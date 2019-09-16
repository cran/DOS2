\name{addalmostexact}
\alias{addalmostexact}
\concept{Almost exact matching}
\concept{Observational study}
\concept{Causal inference}
\concept{Near exact matching}
\concept{Optimal matching}
\title{
Use a Penalty to Obtain a Near Exact Match
}
\description{
A near-exact match (also known as an almost exact match) for a nominal covariate maximizes the number
of pairs that are exactly matched for the nominal covariate.  This is done by adding a large penalty
to the covariate distance between individuals with different values of the nominal covariate.
The topic is discussed in Section 10.2 of the second edition of "Design of Observational Studies".
}
\usage{
addalmostexact(dmat, z, f, mult = 10)
}
\arguments{
  \item{dmat}{
	A distance matrix with one row for each treated individual and one column for each control.  Often, this is either the Mahalanobis distance based on covariates, 'mahal', or else a robust variant produced by 'smahal'.
}
  \item{z}{
	z is a vector that is 1 for a treated individual and 0 for a control.  length(z) must equal sum(dim(dmat)).
}
  \item{f}{
	A vector or factor giving the levels of the nominal covariate.  length(f) must equal length(z).
}
  \item{mult}{
	Let mx be the largest distance in dmat.  Then mult*mx is added to each distance in dmat for which the treated and control individuals have different values of f.  For large enough values of mult, this will ensure that a minimum distance match maximizes the number of individuals who are exactly matched for f.
}
}
\details{
  Very large values of mx may result in numerical instability or slow computations.

  Used informally, a small value of mx, say 1 or 1/10, may not maximize the number of exactly matched pairs, but it may usefully increase the number of exactly matched pairs.
}
\value{
  Returns the penalized distance matrix.
}

\note{
Section 10.2 of "Design of Observational Studies", second edition, discusses almost-exact matching, also called near-exact matching.  The topic is also discussed in section 9.2 of the first edition.

The matching functions in the 'DOS2' package are aids to instruction or self-instruction while reading "Design of Observational Studies".  As in the book, these functions break the task of matching into small steps so they are easy to understand in depth. In practice, matching entails a fair amount of book-keeping best done by a package that automates these tasks.  Consider
R packages 'optmatch', 'rcbalance', 'DiPs', 'designmatch' or 'bigmatch'.  Section 14.10 of "Design of Observational Studies", second edition, discusses and compares several R packages for optimal matching.
}
\references{

Zubizarreta, J. R., Reinke, C. E., Kelz, R. R., Silber, J. H. and Rosenbaum, P. R. (2011) <doi:10.1198/tas.2011.11072> "Matching for several sparse nominal variables in a case control study of readmission following surgery". The American Statistician, 65(4), 229-238.  Combines near-exact matching with fine balance for the same covariate.
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
# Mahalanobis distanced penalized for mismatching on smoking.
dmat<-addalmostexact(dmat, z, smoker, mult = 10)
# The first treated subject (row labeled 27) is a nonsmoker, but the
# third control (column 3) is a smoker, so there is a big penalty.
round(dmat[,1:6],2)
}

