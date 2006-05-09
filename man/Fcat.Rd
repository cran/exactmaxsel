\name{Fcat}
\alias{Fcat}

\title{Distribution of maximally selected statistics for multicategorical variables }
\usage{
Fcat(c, n0, n1, A, statistic)
}
\arguments{
  \item{c}{the value at which the distribution function has to be computed.}
  \item{n0}{the number of observations in class Y=0.}
  \item{n1}{the number of observations in class Y=1.}
  \item{A}{a vector of length K giving the number of observations with
  X=1,...,X=K.}
  \item{statistic}{the association measure used as criterion to select the
  best split. Currently, only \code{statistic="chi2"} (chi-square statistic)
  and \code{statistic="gini"} (the Gini-gain from machine learning) are
  implemented.}  
}
\description{
 The function \code{Fcat} computes the distribution of the maximally selected
 association criterion of interest (either the chi-square statistic or the
 Gini-gain in the current version) when Y is binary and X has unordered 
 categorical values, given \code{n0}, \code{n1} and \code{A}. 
  
}

\details{
  Suppose the response Y is
 binary (Y=0,1) and the predictor X has K unordered categorical values
 (X=1,...,K). The criterion is maximized over all the binary splittings of the
 set {1,...,K}.  For K=3, the criterion is thus maximized over the splittings
 {1,2}{3}, {1,3}{2} and {1}{2,3}.
}

\value{
  the value of the distribution function at \code{c}.
    }

\references{
 
 A.-L. Boulesteix (2006), Maximally selected chi-square statistics
 and binary splits of nominal variables, Biometrical Journal 48.
}

\author{
  Anne-Laure Boulesteix (\url{http://www.statistik.lmu.de/~socher/}) 
   
}
\seealso{\code{\link{Ford}}, \code{\link{maxsel}}.}


\examples{
# load exactmaxsel library
library(exactmaxsel)

Fcat(c=4,n0=15,n1=10,A=c(6,10,9),statistic="chi2")
Fcat(c=5,n0=15,n1=15,A=c(5,8,7,10),statistic="gini")


}
\keyword{htest}
