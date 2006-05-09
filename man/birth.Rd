\name{birth}
\alias{birth}

\title{Birth data set by Boulesteix (2006)}
\usage{
data(birth)
}

\description{
  A data set containing 25 qualitative and quantitative variables
  for n=501 births in 2003 and 2004. 
  
}

\details{
The missing values are of two types:
1. Some of the variables do not concern all the births. For example, the variable
\code{Episiotomy} is not relevant for cesarean births. The type \code{NA} is
used for these entries. Such observations should be removed when analysing
these variables. It does not make sense to imput these missing values.
2. Some questionnaires were returned incomplete, because the mother forgot
some questions or did not know the answer (especially for the head circumference).
}

\value{
A data frame with following variables
  \item{IndexMother}{Index of the mother (unordered)}
  \item{Sex}{Sex of the baby (1:male,2:female)} 
  \item{Weight}{Weight of the baby in g.}
  \item{Height}{Height of the baby in cm.}
  \item{Head}{Head circumference of the baby in cm.}
  \item{Month}{Month of birth (coded as 1:january,...,12:december).}
  \item{Year}{Year of birth.}
  \item{Country}{Country of birth
  (FR:France,BE:Belgium,CH:Switzerland,CA:Canada, etc).}
  \item{Term}{Duration of the pregnancy in week.}
  \item{AgeMother}{Age of the mother in year.}
  \item{Previous}{Number of previous deliveries.}
  \item{WeightBefore}{Weight of the mother before pregnancy in kg.}
  \item{HeightMother}{Height of the mother in cm.}
  \item{WeightEnd}{Weight of the mother at the end of the pregnancy.}
  \item{Twins}{Twin pregnancy? (0:no,1:yes)}
  \item{Intensive}{Number of days spent by the baby in the neonatology intensive
  care unit.}
  \item{Cesarean}{Cesarean birth? (0:no,1:yes)}
  \item{Planned}{If cesarean birth, was the cesarean planned before the
  beginning of labor? (0:no,1:yes)}
  \item{Episiotomy}{If vaginal birth, did the mother have an episiotomy? (0:no,1:yes)}
  \item{Tear}{If vaginal birth, did the mother tear? (0:no,1:yes)}
  \item{Operative}{If vaginal birth, was it a vaginal operative birth (i.e. with
  forceps,vacuum extractor,etc)?
  (0:no,1:yes)}
  \item{Induced}{If vaginal delivery or emergency cesarean, was the labor
  induced medically? (0:no,1:yes)}
  \item{Membranes}{Did the membranes rupture before the beginning of labor?
  (0:no,1:yes)}
  \item{Rest}{Was the mother prescribed bed rest for one month or more?
  (0:no,1:yes)}
  \item{Presentation}{Presentation of the baby at birth
  (1:cephalic,2:breech,3:other,e.g. transverse).}
}
  
}

\source{This data set was collected from internet users via french-speaking 
pregnancy forums in august 2004 by Anne-Laure Boulesteix.}

\references{
 A.-L. Boulesteix (2006), Maximally selected chi-square statistics for ordinal variables,
 Biometrical Journal 48:451-462.
}

\examples{
# load exactmaxsel library
library(exactmaxsel)

# load data set
data(birth)


# Display the two first observations
birth[1:2,]
}

\keyword{datasets}
