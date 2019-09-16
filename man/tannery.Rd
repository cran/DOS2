\name{tannery}
\alias{tannery}
\docType{data}
\concept{Observational studies}
\concept{Observational study}
\concept{Evidence factor}
\concept{Causal inference}
\title{
DNA Damage Among Tannery Workers
}
\description{
Data are from Zhang et al. (2008, Table 2) who studied DNA damage among tannery workers often
exposed to trivalent chromium.  The outcome is the mean tail moment (mtm) of the comet assay,
a standard measure of DNA damage, with higher values signifying greater damage.  The
study describes 90 males in 30 blocks of 3 individuals.  There are three groups, each with
30 individuals.  The three groups had a simlar distribution of ages, and the blocks
control for smoking as closely as possible.  Group e1 consists of 30 exposed workers at the
tannery who worked in the tannery department, where the highest exposures to
trivalent chromium are expected.  Group e2 consists of 30 workers at the tannery who worked
in the finishing department, where exposure to trivalent chromium is expected to be much
lower.  Group c consists of 30 controls who did not work at the tannery.  This example is discussed in Chapter 20 of "Design of Observational Studies", 2nd ed.
}
\usage{data("tannery")}
\format{
  A data frame with 30 observations on the following 4 variables.
  \describe{
    \item{\code{block}}{Block indicator, 1 to 30.}
    \item{\code{e1mtm}}{mtm for the tannery worker from the tannery department.}
    \item{\code{e2mtm}}{mtm for the tannery worker from the finishing department.}
    \item{\code{cmtm}}{mtm for the control who did not work in the tannery}
  }
}
\details{
The comet assay is described by Collins (2004).  It is thought to measure DNA
strand breaks, producing an image that resembles the tail of a comet, a
larger, longer tail suggesting more extensive strand breaks.  This example
was discussed in Rosenbaum (2011).
}
\source{
Data from Zhang et al. (2008, Table 2).
}
\references{
Collins, A. R. (2004) <doi:10.1385/MB:26:3:249> "The comet assay for DNA damage and repair: principles,
applications, and limitations". Molecular Biotechnology 26(3), 249-261.

Rosenbaum, P. R. (2011) <doi:10.1198/jasa.2011.tm10422> "Some approximate evidence factors in observational studies".
Journal of the American Statistical Association, 106, 285-293.

Rosenbaum, P. R. (2013) <doi:10.1111/j.1541-0420.2012.01821.x> "Impact of multiple matched controls on design sensitivity in observational studies". Biometrics 69 118-127. (Introduces inner trimming.)

Rosenbaum, P. R. (2015) <https://obsstudies.org/two-r-packages-for-sensitivity-analysis-in-observational-studies/> "Two R packages for sensitivity analysis in observational studies". Observational Studies, 1(1), 1-17.

Zhang, M., Chen, Z., Chen, Q., Zou, H., Lou, J. He, J. (2008)
<doi:10.1016/j.mrgentox.2008.04.011> "Investigating DNA damage in
tannery workers occupationally exposed to trivalent chromium using comet assay".
Mutation Research/Genetic Toxicology and Environmental Mutagenesis, 654(1), 45-51.
}
\examples{
data(tannery)
boxplot(tannery[,2:4],names=c("Tannery E1","Finishing E2",
    "Control C"),ylab="Mean Tail Moment")
oldpar<-par(mfrow=c(1,2))
boxplot(tannery[,2:3],names=c("E1","E2"),ylab="Mean Tail Moment",
     main="Tannery vs. Finishing",ylim=c(0,12))
boxplot(as.vector(unlist(tannery[,2:3])),tannery[,4],
     names=c("E1+E2","C"),ylab="Mean Tail Moment",
     main="Exposed vs. Control",ylim=c(0,12))

# Stratified Wilcoxon analysis from the chapter Evidence Factors
# of Design of Observational Studies, Second Edition
# Also reproduces the F1, F2 example in Rosenbaum (2011, sec 6).
y<-tannery[,2:4]
rkc<-t(apply(y,1,rank)) # Ranks for (E1,E2,C)
sum(rkc[,3]) # Stratified rank sum for C in (E1, E2, C)
(35-60)/sqrt(20)

y<-tannery[,2:3]
rkc<-t(apply(y,1,rank)) # Ranks for (E1,E2)
sum(rkc[,2]) # Stratified rank sum for E2 in (E1, E2)

# Reorganize y for input to 'separable1v' from 'sensitivitymult'
# 'separable1v' is one-sided, looking for a large rank sum

# Factor 1
y<-tannery[,4:2]*(-1)
rkc<-t(apply(y,1,rank)) # Ranks for -y for (C,E2,E1)
sensitivitymult::separable1v(rkc,gamma=1)
# Test for C in (E1, E2, C)
(85-60)/sqrt(20)
(35-60)/sqrt(20)
sensitivitymult::separable1v(rkc,gamma=6)
# Test for C in (E1, E2, C)
p1<-sensitivitymult::separable1v(rkc,gamma=7)$pval

#Factor 2
y<-tannery[,3:2]*(-1)
rkc<-t(apply(y,1,rank)) # Ranks for -y for (E2,E1)
sensitivitymult::separable1v(rkc,gamma=1)
# Test for E2 in (E2, E1)
# Combine P-values using Fisher's method
sensitivitymv::truncatedP(c(1.134237e-08,0.001743502),trunc=1)

# Larger gammas
sensitivitymult::separable1v(rkc,gamma=1.7)
p2<-sensitivitymult::separable1v(rkc,gamma=2)$pval
# Combine P-values using Fisher's method
c(p1,p2)
sensitivitymv::truncatedP(c(p1,p2),trunc=1)

# Nearly reproduces calculations from Section 6 of Rosenbaum (2011)
# However, in Rosenbaum (2011), the second factor
# uses a pooled scale factor, whereas senm does not,
# so the result is very slightly different.
attach(tannery)
mset<-rep(block,3)
zC<-c(rep(0,60),rep(1,30))
z12<-c(rep(1,30),rep(0,30),rep(NA,30))
y<-c(e1mtm,e2mtm,cmtm)
detach(tannery)
use<-!is.na(z12)
# Factor 1
sensitivitymult::senm(y,zC,mset,gamma=1,
   alternative="less",trim=1)
sensitivitymult::senm(y,zC,mset,gamma=11.7,
   alternative="less",trim=1)
# Factor 2
sensitivitymult::senm(y[use],z12[use],mset[use],
   gamma=2,alternative="greater",trim=1)

# Combine two evidence factors
p1<-sensitivitymult::senm(y,zC,mset,gamma=12,
    alternative="less",trim=1)$pval
p2<-sensitivitymult::senm(y[use],z12[use],mset[use],gamma=3,
    alternative="greater",trim=1)$pval
c(p1,p2)
sensitivitymv::truncatedP(c(p1,p2),trunc=1)
# Combine p-values using Fisher's method

# Other psi-functions often have higher design
# sensitivity; see Rosenbaum (2013)
par(oldpar)
}
\keyword{datasets}
