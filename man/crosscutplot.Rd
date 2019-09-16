\name{crosscutplot}
\alias{crosscutplot}
\concept{Crosscut test}

\title{
Scatterplot Associated with the Cross-Cut Test
}
\description{
A scatterplot of y against x, with the points used by the
cross-cut test in black,
and the remaining points in grey.  See Figure 19.2 in
"Design of Observational Studies", second edition.
}
\usage{
crosscutplot(x, y, ct = 0.25, xlab = "", ylab = "", main = "", ylim = NULL)
}
\arguments{
  \item{x}{
Variable to be plotted on the horizontal axis.
}
  \item{y}{
Variable to be plotted on the vertical axis.
}
  \item{ct}{
The quantile that defines the crosscut.  By default, the crosscut
is at the outer .25 of the data, the lower 25 percent and the upper 75
percent.
}
  \item{xlab}{
Label for the x axis.
}
  \item{ylab}{
Label for the y axis.
}
  \item{main}{
Title of the plot.
}
  \item{ylim}{
Limits for the y axis.  See example.
}
}

\value{
A scatter plot.
}
\references{

Rosenbaum, P. R. (2016) <doi:10.1111/biom.12373> "The crosscut statistic and its sensitivity to bias in observational studies with ordered doses of treatment". Biometrics, 72(1), 175-183.

Rosenbaum, P. R. (2017) <doi:10.1214/17-STS621> "The general structure of evidence factors in observational studies".
Statistical Science 32, 514-530.

}
\author{
Paul R. Rosenbaum
}

\examples{
# Figure 1 in Rosenbaum (2017)
data(periodontal)
attach(periodontal)
m<-matrix(1:2,1,2)
layout(m,widths=c(1,2))
boxplot(pcteither[z==1]-pcteither[z==0],ylab="Smoker-Control Difference",
        main="(i)",xlab="Matched Pairs",ylim=c(-100,100))
abline(h=0,lty=2)
crosscutplot(cigsperday[z==1],pcteither[z==1]-pcteither[z==0],
    ct=.2,ylab="Smoker-Control Difference",
    xlab="Cigarettes per Day",main="(ii)",ylim=c(-100,100))
abline(h=0,lty=2)
detach(periodontal)
layout(1)
}

