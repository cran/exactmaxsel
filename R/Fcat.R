
###     Distribution of maximally selected statistics for nominal variables 
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

Fcat<-function(c,n0,n1,A,statistic)
{
N<-n0+n1
K<-length(A)
if (N!=sum(A))
 stop("error: you must have n0+n1=sum(A)")

number<-sum(as.numeric(permn(K,fun=sumloop,n0=n0,n1=n1,A=A,c=c,statistic=statistic)))
 
return(1-number/choose(n0+n1,n1))
}

###########################
loop<-function(sigma,I,k,n0,n1,A,whichboundary,c,statvector,statistic)
{
K<-length(A)
N<-n0+n1

 
upperbound<-upper(I,A,sigma)

if (k==K)
 {
 if (upperbound>=(n1-sum(I)))
  {
  number<-choose(A[K],n1-sum(I))
  }
 else
  {
  number<-0
  } 
 }
 
else 
 { 
 if (k==whichboundary)
  {
  lowerbound<-max(lower(I,A,sigma,n0,n1,whichboundary,c,statistic=statistic),floor(statvector[k]+1)-sum(I))
  upperbound<-min(upperbound,n1,A[k])
  } 
 if (k<whichboundary)
  { 
  lowerbound<-max(0,lower(I,A,sigma,n0,n1,whichboundary,c,statistic=statistic))
  upperbound<-min(upperbound,n1,A[k],floor(statvector[k])-sum(I))
  }
 if (k>whichboundary)
  {
  lowerbound<-max(0,lower(I,A,sigma,n0,n1,whichboundary,c,statistic=statistic))
  upperbound<-min(upperbound,n1,A[k])
  } 
 if (upperbound<lowerbound)
  {
  number<-0
  }
 else
  {  
 
number<-sum(sapply(as.list(lowerbound:upperbound),FUN=subloop,I=I,k=k,n0=n0,n1=n1,A=A,sigma=sigma,whichboundary=whichboundary,c=c,statvector=statvector,statistic=statistic))
  }
 }


return(number)

}

#################
subloop<-function(i,I,k,n0,n1,A,sigma,whichboundary,c,statvector,statistic)
{

number<-choose(A[k],i)*loop(I=c(I,i),k=k+1,n0=n0,n1=n1,A=A,sigma=sigma,whichboundary=whichboundary,c=c,statvector=statvector,statistic=statistic)

return(number)
}


##################

upper<-function(I,A,sigma)
{
k<-length(I)+1

if (k==1)
 return(A[1])

if (sigma[k-1]>sigma[k])
 {
 upper<-ceiling(I[k-1]*A[k]/A[k-1])-1
 }

else
 {
 upper<-floor(I[k-1]*A[k]/A[k-1])
 }

return(upper)
 
}

###############

lower<-function(I,A,sigma,n0,n1,whichboundary,c,statistic)
{
k<-length(I)+1
K<-length(A)
if (k==K)
 {
 return(0)
 }
if (k>=whichboundary)
 {
 lower<-floor(A[k]*(n1-sum(I))/sum(A[k:K]))
 }
if (k<whichboundary)
 {
 n1bound<-1+floor(boundary(x=sum(A[1:whichboundary]),n0=n0,n1=n1,c=c,statistic=statistic,lower=FALSE)$upper)
 lower<-ceiling(A[k]*(n1bound-sum(I))/sum(A[k:whichboundary]))
 } 
return(lower)
}




#########################

sumloop<-function(sigma,n0,n1,A,c,statistic)
{

K<-length(A)
statvector<-numeric(K-1)
sumloop<-numeric(K-1)
A<-A[sigma]

for (i in 1:(K-1))
 {
 statvector[i]<-boundary(sum(A[1:i]),n0=n0,n1=n1,c=c,statistic=statistic,lower=FALSE)$upper
 }

for (k in 1:(K-1))
 {
 sumloop[k]<-loop(sigma,I=c(),k=1,n0=n0,n1=n1,A=A,whichboundary=k,c=c,statvector=statvector,statistic=statistic)
 }

return(sum(sumloop))


}






