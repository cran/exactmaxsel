###    Test of independence based on maximally selected statistics
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



maxsel.test<-function(x,y=NULL,type,statistic)
{

maxselcrit<-maxsel(x,y,type=type,statistic=statistic)


if (!is.null(y)&ncol(as.matrix(x))==1)
 {
 notNA<-which(!is.na(x)&!is.na(y))
 x<-x[notNA]
 y<-y[notNA]
 levelsx<-sort(union(x,x))
 K<-length(levelsx)
 n0<-sum(y==0)
 n1<-length(y)-n0
 
 A<-numeric(K)
 for (k in 1:K)
  {
  A[k]<-sum(x==levelsx[k])
  }
 }
else 
 {
 n0<-sum(x[1,])
 n1<-sum(x[2,])
 A<-apply(x,FUN=sum,MARGIN=2)
 } 

if (type=="ord")
 {
 p<-Ford(c=maxselcrit,n0=n0,n1=n1,A=A,statistic=statistic)
 return(p)
 }
 

if (type=="cat")
 {
 p<-Fcat(c=maxselcrit,n0=n0,n1=n1,A=A,statistic=statistic)
 return(p)
 }

}

