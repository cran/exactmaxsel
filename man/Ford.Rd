\name{Ford}
\alias{Ford}

\title{Distribution of maximally selected statistics for (at least) ordinally
scaled variables}
\usage{
Ford(c, n0, n1, A, statistic)
}
\arguments{
  \item{c}{the value at which the distribution function has to be computed.}
  \item{n0}{the number of observations in class Y=0.}
  \item{n1}{the number of observations in class Y=1.}
  \item{A}{a vector of length K giving the number of observations with
  X=1,...,X=K. In the special case of a continuous X variable taking distinct
  values in the available sample, \code{A} takes the form \code{A=rep(1,N)},
  where \code{N=n0+n1}.}
  \item{statistic}{the association measure used as criterion to select the
  best split. Currently, only \code{statistic="chi2"} (chi-square statistic)
  and \code{statistic="gini"} (the Gini-gain from machine learning) are
  implemented.}  
}
\description{
 The function \code{Ford} computes the distribution of the maximally selected
 association criterion of interest (either the chi-square statistic or the
 Gini-gain in the current version) when Y is binary and X has ordered 
 values, given \code{n0}, \code{n1} and \code{A}. Note that X must be AT LEAST
 ordinally scaled, i.e. continuous variables are also allowed as an extreme 
 special case.
  
}

\details{
  Suppose the response Y is
 binary (Y=0,1) and the predictor X has K ordered categorical values
 (X=1,...,K). The criterion is maximized over all the binary splittings of the
 set {1,...,K} that preserve the ordering.  For K=3, the criterion is thus maximized 
 over the splittings {1,2}{3} and {1}{2,3}. Note that X may also be a
 substantially continuous variable that is observed at a discrete scale and
 thus has ties.
}

\value{
  the value of the distribution function at \code{c}.
    }

\references{
 A.-L. Boulesteix (2006), Maximally selected chi-square statistics
 for ordinal variables, Biometrical Journal 48:451-462.

}

\author{
  Anne-Laure Boulesteix (\url{http://www.slcmsr.net/boulesteix}) 
}
\seealso{\code{\link{Fcat}}, \code{\link{Ford2}}, \code{\link{maxsel}}.}


\examples{
# load exactmaxsel library
library(exactmaxsel)

Ford(c=4,n0=15,n1=10,A=c(6,10,9),statistic="chi2")
Ford(c=0.02,n0=15,n1=15,A=c(5,8,7,10),statistic="gini")


}
\keyword{htest}
