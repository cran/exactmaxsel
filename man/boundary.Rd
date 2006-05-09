\name{boundary}
\alias{boundary}

\title{Computes the coordinates of the boundaries}
\usage{
boundary(x, n0, n1, c, statistic, lower=TRUE)
}
\arguments{
  \item{x}{the number of observations in the left node}
  \item{n0}{the number of observations in class Y=0}
  \item{n1}{the number of observations in class Y=1}
  \item{c}{the value of the criterion that should not be exceeded.}
  \item{statistic}{the association measure. Currently, only
  \code{statistic="chi2"} (chi-square statistic)
  and \code{statistic="gini"} (the Gini-gain from machine learning) are
  implemented.}
  \item{lower}{Should the lower boundary also be computed? }
  
}
\description{
 The function \code{boundary} is an internal function that computes the greatest (and also the smallest if
 \code{lower=TRUE})  number of observations
 from class Y=1 in the left node that lead to an association criterion <= c. 
 Remark: in general, the numbers output by \code{boundary} are not integers.
}

\details{
This function should not be called by the user in practice.
}

\value{
  a list with
   \item{upper}{the upper boundary (greatest allowed value).}
   \item{lower}{the lower boundary (smallest allowed value).}
    }

\references{
 A.-L. Boulesteix (2006), Maximally selected chi-square statistics
 for ordinal variables, Biometrical Journal 48.
 
 A.-L. Boulesteix (2006), Maximally selected chi-square statistics
 and binary splits of nominal variables, Biometrical Journal 48.
 
 C. Strobl, A.-L. Boulesteix and T. Augustin (2006), Unbiased split selection 
 for classification trees based on the Gini index, submitted.
 

  
 
}

\author{
  Anne-Laure Boulesteix (\url{http://www.statistik.lmu.de/~socher/}) 
   
}
\seealso{\code{\link{Ford}}, \code{\link{Fcat}}.}


\examples{
# load exactmaxsel library
library(exactmaxsel)

boundary(10,30,30,c=3,statistic="chi2",lower=TRUE)


}

\keyword{htest}
