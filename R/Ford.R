###     Distribution of maximally selected statistics for ordinal variables 
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


Ford<-function(c,n0,n1,A,statistic)
{

N<-n0+n1
if (N!=sum(A))
 stop("error: you must have n0+n1=sum(A)")

K<-length(A)
paths<-path(c,n0,n1,A,statistic=statistic)
pathmin<-paths$pathmin
pathmax<-paths$pathmax
#A<-paths$A
Amin<-paths$Amin
Amax<-paths$Amax

xx<-c(Amin,Amax)
yy<-c(pathmin,pathmax)
x<-sort(xx)
y<-yy[order(xx)]

if (length(xx)==0)
 return(1)

b<-numeric(length(xx))

i<-2
b[1]<-choose(x[1],y[1])

while (i<=length(xx))
 {
 b[i]<-choose(x[i],y[i])-sum(sapply(1:(i-1),FUN=underpmax,i,x,y,b))
 i<-i+1
 }


pval<-0 
for (i in 1:length(xx))
 {
 pval<-pval+b[i]*choose(N-x[i],n1-y[i])
 }
pval<-1-pval/choose(N,n1)

pval

}



#######################
underpmax<-function(j,i,x,y,b)
{
a<-choose(x[i]-x[j],y[i]-y[j])*b[j]
a
}

##########################



path<-function(c,n0,n1,A,statistic)
{

N<-n0+n1
K<-length(A)

pathmin<-c()
pathmax<-c()
Acum<-numeric(K-1)

quot<-1/n0+1/n1

 Amin<-c()
 Amax<-c()
 for (k in 1:(K-1))
  {
  i<-sum(A[1:k])
  Acum[k]<-i
  myboundary<-boundary(x=i,n0=n0,n1=n1,c=c,lower=TRUE,statistic=statistic)
  pmk<-floor(myboundary$upper)+1
  if (pmk<=min(i,n1))
   {
   pathmaxk<-seq(pmk,min(i,n1))
   }
  else
   {
   pathmaxk<-c()
   }

  pmk<-ceiling(myboundary$lower)-1
  if (pmk>=max(0,i-n0))
   {
   pathmink<-seq(max(0,i-n0),pmk)
   }
  else
   {
   pathmink<-c()
   }

  Amin<-c(Amin,rep(Acum[k],length(pathmink)))
  Amax<-c(Amax,rep(Acum[k],length(pathmaxk)))
  pathmin<-c(pathmin,pathmink)
  pathmax<-c(pathmax,pathmaxk)
  }

list(pathmin=pathmin,pathmax=pathmax,Amin=Amin,Amax=Amax)
}

