\name{Ford2}
\alias{Ford2}

\title{Distribution of maximally selected statistics for (at least) ordinally
scaled variables in the two-cutpoint context}
\usage{
Ford2(c, n0, n1, A, statistic)
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
 The function \code{Ford2} computes the distribution of the maximally selected
 association criterion of interest (either the chi-square statistic or the
 Gini-gain in the current version) when Y is binary and X has ordered 
 values, given \code{n0}, \code{n1} and \code{A}, in the case of a non-monotonic 
 association represented by two cutpoints. 
  
}

\details{
  Suppose the response Y is
 binary (Y=0,1) and the predictor X has K ordered categorical values
 (X=1,...,K). The criterion is maximized over all the binary splittings of the
 set {1,...,K} that are obtained from at most two cutpoints.  For example, with K=4, 
 the criterion is maximized 
 over the splittings {1,2,3}{4}, {1,2}{3,4}, {1}{2,3,4}, {1,2,4}{3}, {1,4}{2,3} 
 and {1,3,4}{2}. 

}

\value{
  the value of the distribution function at \code{c}.
    }

\references{

 A.-L. Boulesteix and C. Strobl (2006), Maximally selected chi-square statistics
 and umbrella orderings, Computational Statistics and Data Analysis (in press).
}

\author{
  Anne-Laure Boulesteix (\url{http://www.slcmsr.net/boulesteix}) 
   
}
\seealso{\code{\link{Ford}},\code{\link{Fcat}}, \code{\link{maxsel}}.}


\examples{
# load exactmaxsel library
library(exactmaxsel)

Ford2(c=4,n0=15,n1=15,A=c(6,10,9,5),statistic="chi2")
Ford2(c=0.02,n0=15,n1=15,A=c(5,8,7,10),statistic="gini")


}
\keyword{htest}
