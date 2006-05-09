\name{ginigain}
\alias{ginigain}

\title{Computation of the Gini gain from a 2x2 contingency table}
\usage{
ginigain(mat)
}
\arguments{
  \item{mat}{a 2x2 matrix corresponding to a two-dimensional contingency table.
  The first row and the second row correspond to Y=0 and Y=1, respectively. The
  first column and the second column correspond to the left and right nodes,
  respectively.}
  
}
\description{
 The function \code{ginigain} computes the Gini gain (also denoted as
 impurity reduction) resulting by splitting
 into the left and the right nodes whose counts are given in the contingency
 table \code{mat}. See Strobl et al. (2006) for a more precise definition.
}

\details{
Note that, in contrast to the chi-square statistic, the Gini gain does not
treat X and Y symmetrically. 
}

\value{
  the (positive) value of the Gini gain.
    }

\references{
 
 L. Breiman, J.H. Friedman, R.A. Olshen and C.J. Stone (1984), Classification
 and Regression Trees, Wadsworth, Monterey, CA.
 
 C. Strobl, A.-L. Boulesteix and T. Augustin (2006), Unbiased split selection 
 for classification trees based on the Gini index, submitted.
 
 
}

\author{
  Anne-Laure Boulesteix (\url{http://www.statistik.lmu.de/~socher/}) 
   
}
\seealso{\code{\link{maxsel.test}}.}


\examples{
# load exactmaxsel library
library(exactmaxsel)

# Define matrix
A<-matrix(c(15,20,22,12),2,2)

# Compute Gini gain
ginigain(A)


}

\keyword{htest}
