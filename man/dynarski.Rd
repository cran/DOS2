\name{dynarski}
\alias{dynarski}
\concept{Observational study}
\concept{Propensity score}
\concept{Matched with multiple controls}
\concept{Causal inference}
\concept{Multivariate matching}
\concept{Optimal matching}
\concept{Caliper matching}
\concept{Design of Observational Studies}
\docType{data}
\title{
A Natural Experiment Concerning Subsidized College Education
}
\description{
The data are from Susan Dynarski (2003)'s study of the effects of a subsidy for college education provided Social Security to children whose fathers had died.  The subsidy was eliminated in 1982.  As presented here, the data are only a portion of the data used in Dynarski (2003), specifically the period from 1979-1981 when the subsidy was available.  The data set here is Table 14.1 of "Design of Observational Studies" (2nd edition), where it was used with the limited goal of illustrating matching techniques.
}
\usage{data("dynarski")}
\format{
  A data frame with 2820 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{identification number}
    \item{\code{zb}}{treatment indicator, zb=1 if 1979-1981 with father deceased,
    or zb=0 if 1979-1981 with father alive}
    \item{\code{faminc}}{family income, in units of tens of thousands of dollars}
    \item{\code{incmiss}}{indicator for missing family income}
    \item{\code{black}}{1 if black, 0 otherwise}
    \item{\code{hisp}}{1 if hispanic, 0 otherwise}
    \item{\code{afqtpct}}{test score percentile, Armed Forces Qualifications Test}
    \item{\code{edmissm}}{indicator for missing education of mother}
    \item{\code{edm}}{education of mother, 1 is <high school, 2 is high school,
    3 is some college, 4 is a BA degree or higher}
    \item{\code{female}}{1 if female, 0 if male}
  }
}
\details{
The examples reproduce steps in Chapter 14 of "Design of Observational Studies" (2nd edition).  Portions of the examples require you to load Ben Hansen's 'optmatch' package and accept its academic license; these portions of the examples do not run automatically.  Hansen's 'optmatch' package uses the Fortran code of Bertsekas and Tseng (1988).
}
\source{
Dynarski (2003).
}
\references{
Bertsekas, D. P. and Tseng, P. (1988) <doi:10.1007/BF02288322> "The relax codes for linear minimum cost network flow problems". Annals of Operations Research, 13(1), 125-190.

Dynarski, S. M. (2003) <doi:10.1257/000282803321455287> "Does aid matter? Measuring the effect of student aid on college attendance and completion". American Economic Review, 93(1), 279-288.

Hansen, B. B. (2007) <https://www.r-project.org/conferences/useR-2007/program/presentations/hansen.pdf> "Flexible, optimal matching for observational studies". R News, 7, 18-24.

Hansen, B. B. and Klopfer, S. O. (2006)
<doi:10.1198/106186006X137047> "Optimal full matching and related designs via network flows". Journal of Computational and Graphical Statistics, 15(3), 609-627.

Rosenbaum, P. R. (1989). "Optimal matching for observational studies" <doi:10.1080/01621459.1989.10478868> Journal of the American Statistical Association, 84(408), 1024-1032.

Rosenbaum, P. R. (1991) <doi:10.1111/j.2517-6161.1991.tb01848.x> A characterization of optimal designs for observational studies. Journal of the Royal Statistical Society: Series B (Methodological), 53(3), 597-610.
}
\examples{
#
data(dynarski)
# Table 14.1 of "Design of Observational Studies" (2nd edition)
head(dynarski)
# Table 14.2 of "Design of Observational Studies" (2nd edition)
zb<-dynarski$zb
zbf<-factor(zb,levels=c(1,0),labels=c("Father deceased",
    "Father not deceased"))
table(zbf)
Xb<-dynarski[,3:10]

# Estimate the propensity score, Rosenbaum (2010, Section 14.3)
p<-glm(zb~Xb$faminc+Xb$incmiss+Xb$black+Xb$hisp
    +Xb$afqtpct+Xb$edmissm+Xb$edm+Xb$female,
    family=binomial)$fitted.values
# Figure 14.1 in "Design of Observational Studies" (2nd edition)
boxplot(p~zbf,ylab="Propensity score",main="1979-1981 Cohort")

# Read about missing covariate values in section 14.4
# of "Design of Observational Studies" (2nd edition)

# Robust Mahalanobis distance matrix, treated x control
dmat<-smahal(zb,Xb)
dim(dmat)
# Table 14.3 in "Design of Observational Studies" (2nd edition)
round(dmat[1:5,1:5],2)

# Add a caliper on the propensity score using a penalty function
dmat<-addcaliper(dmat,zb,p,caliper=.2)
dim(dmat)
# Table 14.4 in "Design of Observational Studies" (2nd edition)
round(dmat[1:5,1:5],2)
\dontrun{
# YOU MUST LOAD the 'optmatch' package and accept its license to continue.
# Note that the 'optmatch' package has changed since 2010.  It now suggests
# that you indicate the data explicitly as data=dynarski.

# Creating a 1-to-10 match, as in section 14.6 of
# "Design of Observational Studies" (2nd edition)
# This may take a few minutes.
m<-fullmatch(dmat,data=dynarski,min.controls = 10,max.controls = 10,
    omit.fraction = 1379/2689)
length(m)
sum(matched(m))
1441/11 # There are 131 matched sets, 1 treated, 10 controls

# Alternative, simpler code to do the same thing
m2<-pairmatch(dmat,controls=10,data=dynarski)
# Results are the same:
sum(m[matched(m)]!=m2[matched(m2)])

# Housekeeping
im<-as.integer(m)
dynarski<-cbind(dynarski,im)
dm<-dynarski[matched(m),]
dm<-dm[order(dm$im,1-dm$zb),]

# Table 14.5 in "Design of Observational Studies" (2nd edition)
which(dm$id==10)
dm[188:198,]
which(dm$id==396)
dm[23:33,]
which(dm$id==3051)
dm[1068:1078,]
# In principle, there can be a tie, in which several different
# matched samples all minimize the total distance.  On my
# computer, this calculation reproduces Table 14.5, but were
# there a tie, 'optmatch' should return one of the tied optimal
# matches, but not any particular one.
}
}
\keyword{datasets}
