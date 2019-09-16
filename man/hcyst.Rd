\name{hcyst}
\alias{hcyst}
\docType{data}
\concept{Observational study}
\concept{Causal inference}
\concept{Multiple controls}
\concept{Sensitivity analysis}
\concept{Design sensitivity}
\concept{Crosscut test}
\title{
Homocysteine and Smoking
}
\description{
NHANES 2005-2006 data on smoking and homocysteine levels in adults.  Matched triples, one daily smoker matched to two never smokers in 548 matched sets of size 3.
}
\usage{data("hcyst")}
\format{
  A data frame with 1644 observations on the following 13 variables.
  \describe{
    \item{\code{SEQN}}{NHANES identification number}
    \item{\code{z}}{Smoking status, 1 = daily smoker, 0 = never smoker}
    \item{\code{female}}{1 = female, 0 = male}
    \item{\code{age}}{Age in years, >=20, capped at 85}
    \item{\code{black}}{1=black race, 0=other}
    \item{\code{education}}{Level of education, 3=HS or equivalent}
    \item{\code{povertyr}}{Ratio of family income to the poverty level, capped at 5 times poverty}
    \item{\code{bmi}}{BMI or body-mass-index}
    \item{\code{cigsperday30}}{Cigarettes smoked per day, 0 for never smokers}
    \item{\code{cotinine}}{Blood cotinine level, a biomarker of recent exposure to tobacco. Serum cotinine in ng/mL}
    \item{\code{homocysteine}}{Level of homocysteine. Total homocysteine in blood plasma in umol/L}
    \item{\code{p}}{An estimated propensity score}
    \item{\code{mset}}{Matched set indicator, 1, 1, 1, 2, 2, 2, ..., 548, 548, 548}
  }
}
\details{
The matching controlled for female, age, black, education, income, and BMI.  The outcome is homocysteine. The treatment is daily smoking versus never smoking.

Bazzano et al. (2003) noted higher homocysteine levels in smokers than in nonsmokers.  The NHANES data have been used several times to illustrate statistical methods; see Pimental et al (2016), Rosenbaum
(2018), Yu and Rosenbaum (2019).

This is Match 4 in Yu and Rosenbaum (2019).  The match used a propensity score estimated from these covariates, and it was finely balanced for education.  It used a rank-based Mahalanobis distance, directional penalties, and a constraint on the maximum distance.  The larger data set before matching is contained in the DiPs package as nh0506homocysteine.

The two controls for each smoker have been randomly ordered, so that, for illustration, a matched pair analysis using just the first control may be compared with an analysis of a 2-1 match.  See Rosenbaum (2013) and Rosenbaum (2017b, p222-223) for discussion of multiple controls and design sensitivity.
}
\source{
From the NHANES web page, for NHANES 2005-2006.  US National Health and Nutrition Examination Survey, 2005-2006.  From the US National Center for Health Statistics.
}
\references{
Bazzano, L. A., He, J., Muntner, P., Vupputuri, S. and Whelton, P. K. (2003)
<doi:10.7326/0003-4819-138-11-200306030-00010> "Relationship between cigarette smoking and novel risk factors for cardiovascular disease in the United States".  Annals of Internal Medicine, 138, 891-897.

Pimentel, S. D., Small, D. S. and Rosenbaum, P. R. (2016)
<doi:10.1080/01621459.2015.1076342> "Constructed second control groups and
attenuation of unmeasured biases".  Journal of the American Statistical Association, 111, 1157-1167.

Rosenbaum, P. R. (2007) <doi:10.1111/j.1541-0420.2006.00717.x> "Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies".  Biometrics, 2007, 63, 456-464.

Rosenbaum, P. R. and Silber, J. H. (2009)
<doi:10.1198/jasa.2009.tm08470> "Amplification of sensitivity analysis in observational studies".  Journal of the American Statistical Association, 104, 1398-1405.


Rosenbaum, P. R. (2013) <doi:10.1111/j.1541-0420.2012.01821.x> "Impact of multiple matched controls on design sensitivity in observational studies".  Biometrics, 2013, 69, 118-127.

Rosenbaum, P. R. (2015) <https://obsstudies.org/two-r-packages-for-sensitivity-analysis-in-observational-studies/> "Two R packages for sensitivity analysis in observational studies".  Observational Studies, v. 1.

Rosenbaum, P. R. (2016) <doi:10.1111/biom.12373> "The crosscut statistic and its sensitivity to bias in observational studies with ordered doses of treatment". Biometrics, 72(1), 175-183.

Rosenbaum, P. R. (2017a) <doi.org/10.1214/18-AOAS1153> "Sensitivity analysis for stratified comparisons in an observational study
of the effect of smoking on homocysteine levels".  Annals of Applied Statistics, 12, 2312–2334


Rosenbaum, P. R. (2017b) <https://www.hup.harvard.edu/catalog.php?isbn=9780674975576> "Observation and Experiment: An Introduction to Causal Inference".  Cambridge, MA: Harvard Univeristy Press, pp 222-223.

Yu, Ruoqi. and Rosenbaum, P. R. (2019) <doi.org/10.1111/biom.13098> "Directional penalties for optimal matching in observational studies".  Biometrics, to appear.  See Yu's 'DiPs' package.
}
\examples{
data(hcyst)
attach(hcyst)
tapply(female,z,mean)
tapply(age,z,mean)
tapply(black,z,mean)
tapply(education,z,mean)
table(z,education)
tapply(povertyr,z,mean)
tapply(bmi,z,mean)
tapply(p,z,mean)
ind<-rep(1:3,548)
hcyst<-cbind(hcyst,ind)
hcystpair<-hcyst[ind!=3,]
rm(ind)
detach(hcyst)

# Analysis of paired data, excluding second control
attach(hcystpair)
y<-log(homocysteine)[z==1]-log(homocysteine)[z==0]
x<-cotinine[z==1]-cotinine[z==0]


senWilcox(y,gamma=1)
senWilcox(y,gamma=1.53)
senU(y,m1=4,m2=5,m=5,gamma=1.53)
senU(y,m1=7,m2=8,m=8,gamma=1.53)
senU(y,m1=7,m2=8,m=8,gamma=1.7)
# Interpretation/amplification of gamma=1.53 and gamma=1.7
# See Rosenbaum and Silber (2009)
sensitivitymult::amplify(1.53,c(2,3))
sensitivitymult::amplify(1.7,c(2,3))


crosscutplot(x,y,ylab="Difference in log(homocysteine)",
   xlab="Difference in Cotinine, Smoker-minus-Control",main="Homocysteine and Smoking")
text(600,1.8,"n=41")
text(600,-1,"n=31")
text(-500,-1,"n=43")
text(-500,1.8,"n=21")

crosscut(x,y)
crosscut(x,y,gamma=1.25)


# Comparison of pairs and matched triples
# Triples increase power and design sensitivity; see Rosenbaum (2013)
# and Rosenbaum (2017b, p222-223)
library(sensitivitymult)
sensitivitymult::senm(log(homocysteine),z,mset,gamma=1.75)$pval
detach(hcystpair)
attach(hcyst)
sensitivitymult::senm(log(homocysteine),z,mset,gamma=1.75)$pval
# Inner trimming improves design sensitivity; see Rosenbaum (2013)
sensitivitymult::senm(log(homocysteine),z,mset,inner=.5,gamma=1.75)$pval

# Interpretation/amplification of gamma = 1.75
# See Rosenbaum and Silber (2009)
sensitivitymult::amplify(1.75,c(2,3,10))

# Confidence interval and point estimate
# sensitivitymult::senmCI(log(homocysteine),z,mset,inner=.5,gamma=1.5)
detach(hcyst)
}
\keyword{datasets}
