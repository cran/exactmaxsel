\name{maxsel.test}
\alias{maxsel.test}

\title{Test of independence based on maximally selected statistics}
\usage{
maxsel.test(x, y=NULL, type, statistic)
}
\arguments{ \item{x}{a numeric vector of length n giving the values of the variable
  X for the considered n observations. The
  classes  must be coded as 1,...,K. Alternatively, \code{x} can be a 2 x K matrix corresponding
  to a contingency table,
  where the two rows are for the values of Y (Y=0,1) and the K columns are for the values
  of X (X=1,...,K). In this case, \code{y} must be set to \code{y=NULL}.
  }
  \item{y}{a numeric vector of length n giving the class (response variable Y) of the considered observations. The
  classes  must be coded as 0 and 1. If \code{x} is a contingency table, \code{y} must be set
  to \code{y=NULL}.}
  \item{type}{the type of the considered binary splits. \code{type="ord"}
  corresponds to an ordinal X variable, \code{type="cat"} 
  corresponds to a categorical X variable with unordered categories, \code{type="ord2"} corresponds
  to an ordinal X variable with 2 cutpoints (non-monotonous association). }
  \item{statistic}{the association measure used as criterion to select the
  best split. Currently, only \code{statistic="chi2"} (chi-square statistic)
  and \code{statistic="gini"} (the Gini-gain from machine learning) are
  implemented.}  
  
}
\description{
 The function \code{maxsel.test} computes the probability that the maximally
 selected criterion is <= than the value observed
from the data, under the null-hypothesis of no association between X and 
Y, given the numbers of observations with Y=0,Y=1,X=1,...,X=K. 
The candidate binary splits over which the criterion is maximized 
depend on \code{type} (see details). If p denotes the output of the function
\code{maxsel.test}, 1-p may be seen as the p-value of an independence test.
}

\details{
For example, let us consider a variable X with the possible values \{1,2,3,4\}.
If \code{type="ord"}, the set of candidate splits consists of \{1\}\{2,3,4\}, \{1,2\}\{3,4\} and \{1\}\{2,3,4\}.
If \code{type="cat"}, the set of candidate splits consists of \{1\}\{2,3,4\}, \{1,2\}\{3,4\}, 
\{1,2,3\}\{4\}, \{1,2,4\}\{3\}, \{1,4\}\{2,3\}, \{1,3,4\}\{2\}, \{1,3\}\{2,4\}.
If \code{type="ord2"}, the set of candidate splits consists of \{1\}\{2,3,4\}, \{1,2\}\{3,4\}, 
\{1,2,3\}\{4\}, \{1,2,4\}\{3\}, \{1,4\}\{2,3\}, \{1,3,4\}\{2\}.
}

\value{
the probability that the maximally selected criterion is <= than the value observed
from the data, under the null-hypothesis of no association between \code{x} and 
\code{y}, given the numbers of observations with Y=0,Y=1,X=1,...,X=K.
    }

\references{
 A.-L. Boulesteix (2006), Maximally selected chi-square statistics
 for ordinal variables, Biometrical Journal 48:451-462.
 
 A.-L. Boulesteix (2006), Maximally selected chi-square statistics
 and binary splits of nominal variables, Biometrical Journal 48:838-848.
 
 C. Strobl, A.-L. Boulesteix and T. Augustin (2007), Unbiased split selection 
 for classification trees based on the Gini index, Computational Statistics and Data
Analysis 52:483-501.
 
 A.-L. Boulesteix and C. Strobl (2006), Maximally selected chi-square statistics and non-monotonic associations: an exact approach based on two cutpoints. Computational Statistics and Data Analysis 51:6295-6306.

}

\author{
  Anne-Laure Boulesteix (\url{http://www.ibe.med.uni-muenchen.de/organisation/mitarbeiter/020_professuren/boulesteix/}) 
   
}
\seealso{\code{\link{maxsel}}.}


\examples{
# load exactmaxsel library
library(exactmaxsel)

# First case: x and y are data vectors
# Simulate x and y
x<-sample(4,30,replace=TRUE)
y<-sample(c(0,1),30,replace=TRUE)

maxsel.test(x=x,y=y,type="ord",statistic="chi2")
maxsel.test(x=x,y=y,type="cat",statistic="gini")


# Second case: x is a contingency table, y=NULL.
x<-matrix(c(8,10,40,13,15,4),2,4,byrow=TRUE)
maxsel.test(x=x,y=NULL,type="ord",statistic="chi2")
maxsel.test(x=x,y=NULL,type="cat",statistic="gini")



}

\keyword{htest}
