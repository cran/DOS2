\name{angristlavy}
\alias{angristlavy}
\concept{Observational study}
\concept{Sensitivity analysis}
\concept{Matched pairs}
\concept{Causal inference}
\concept{Instrumental variable}
\concept{Maimonidies Rule}
\concept{Design of Observational Studies}
\docType{data}
\title{
Class Size and Academic Performance -- Maimonidies Rule
}
\description{
This data set is from Angrist and Lavy (1999). There are 86 pairs of two Israeli schools, one with slightly more than 40 students in the fifth grade, the other with 40 or fewer in the fifth grade, together with test scores in reading and math.  This example is discussed in Chapters 1, 5 and 18 of "Design of Observational Studies".
}
\usage{data("angristlavy")}
\format{
  A data frame with 172 observations on the following 9 variables.
  \describe{
    \item{\code{scode}}{School ID}
    \item{\code{numclass}}{Number of classes in the fifth grade, 1 or 2.}
    \item{\code{cohsize}}{Total number of students in the fifth grade, near 40 for these schools.}
    \item{\code{avgmath}}{Average grade in math in the fifth grade.}
    \item{\code{avgverb}}{Average verbal grade in the fifth grade.}
    \item{\code{tipuach}}{Percent of disadvantaged students.  Used to form matched pairs.}
    \item{\code{clasz}}{Average class size in the fifth grade, equal to cohsize/numclass}
    \item{\code{z}}{1 if cohsize<=40, 0 if cohsize>40.}
    \item{\code{pair}}{pair ID, 1, 2, ..., 86}
  }
}
\details{
This example is discussed in Chapters 1, 5, and 18 of the second edition of "Design of Observational Studies".

As discussed by Angrist and Lavy (1999), Maimonidies rule requires that a class of more than 40 be divided to form two or more classes of at most 40, so there is a large discontinuity in class size at 40: at 40 students in the 5th grade, there is one class of 40, but at 41 students, there are two classes of average size 20.5.  So the enrolement of one student should cut the class size roughly in half.  Adherence to Maimonidies rule is good but not perfect.  Pairs of schools were matched for the percent of disadvantaged students (tipuach).
}
\references{
Angrist, J. D. and Lavy, V. (1999) <doi:10.1162/003355399556061> "Using Maimonides' rule to estimate the effect of class size on scholastic achievement". The Quarterly Journal of Economics, 114, 533-575.

Angrist, J. D. and Krueger, A. B. (1999)
<doi:10.1016/S1573-4463(99)03004-7> "Empirical strategies in labor economics". In Handbook of Labor Economics (Vol. 3, pp. 1277-1366). Elsevier.

}
\examples{
# Figure 1.1 in Chapter 1 of "Design of Observational Studies", 2nd edition
data(angristlavy)
attach(angristlavy)
grp<-factor(z,levels=c(1,0),labels=c("31-40","41-50"),ordered=TRUE)
oldpar<-par(mfrow=c(2,2))
boxplot(tipuach~grp,main="Disadvantaged",ylab="Percent")
boxplot(clasz~grp,main="Class Size",ylab="Students")
boxplot(avgmath~grp,main="Math",ylab="Average Score")
boxplot(avgverb~grp,main="Verbal",ylab="Average Score")
detach(angristlavy)
par(oldpar)
}
\keyword{datasets}
