###     Computation of the maximally selected statistics from data
###
### Copyright 2006-05 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `exactmaxstat' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


maxsel<-function(x,y=NULL,type,statistic)
{
 
statvector<-c()

 
if (!is.null(y))
 {
 notNA<-which(!is.na(x)&!is.na(y))
 x<-x[notNA]
 y<-y[notNA]
 levelsx<-sort(union(x,x))
 K<-length(levelsx)
 n0<-sum(y==0)
 n1<-length(y)-n0

 
 if (type=="ord")
  {
  for (k in 1:(K-1))
   {
   c1<-sum(y==0&(x<=levelsx[k]))
   c2<-sum(y==1&(x<=levelsx[k]))
   c3<-n0-c1
   c4<-n1-c2
   if (statistic=="chi2")
    {
    statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
    }
   if (statistic=="gini")
    {
    statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
    }   
   }
  }
 
   

 if (type=="cat")
  {
  xfac<-factor(x,labels=1:K)
  x<-as.numeric(xfac)
  prop2<-as.numeric(tapply(xfac[y==1],xfac[y==1],length))/as.numeric(tapply(xfac,xfac,length))
 
  for (k in 1:(K-1))
   {
   left<-is.element(x,order(-prop2)[1:k])
   c1<-sum(y==0&left)
   c2<-sum(y==1&left)
   c3<-n0-c1
   c4<-n1-c2
   if (statistic=="chi2")
    {
    statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
    }
   if (statistic=="gini")
    {
    statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
    }   
   }
  } 
 }



if (is.null(y)&ncol(as.matrix(x))>1)
 {
 if (!is.numeric(x))
  stop("x must be given as a numeric vector") 
 if (nrow(x)!=2)
  stop("x must have 2 rows and K columns")
 
 K<-ncol(x)
 zero<-c()
 for (k in 1:K)
  {
  if (sum(x[,k])==0)
   {
   zero<-c(zero,k)
   }
  }
 if (length(zero)>0) 
  {
  x<-x[,-zero]
  }
 K<-ncol(x)
 n0<-sum(x[1,])
 n1<-sum(x[2,])

 if (K<2)
  stop("x must have at least 2 columns")
 if (K<3&(type=="cat"))
  {
  type<-"ord"
  }
 
 if (type=="ord")
  {
  for (k in 1:(K-1))
   {
   c1<-sum(x[1,1:k])
   c2<-sum(x[2,1:k])
   c3<-n0-c1
   c4<-n1-c2
   if (statistic=="chi2")
    {
    statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
    }
   if (statistic=="gini")
    {
    statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
    }   
   }
  }
  
 
  
   
  
 if (type=="cat")
  {
  prop2<-x[2,]/apply(x,FUN=sum,MARGIN=2)
  x<-x[,order(-prop2)] 
  for (k in 1:(K-1))
   {
   c1<-sum(x[1,1:k])
   c2<-sum(x[2,1:k])
   c3<-n0-c1
   c4<-n1-c2
   if (statistic=="chi2")
    {
    statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
    }
   if (statistic=="gini")
    {
    statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
    }    
   }
  }  
 }
 


if (statistic=="chi2")
 {
 statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
 }
if (statistic=="gini")
 {
 statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
 }   
    
maxsel<-max(as.numeric(statvector))
return(maxsel)

}





