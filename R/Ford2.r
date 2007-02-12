###     Distribution of maximally selected statistics for ordinal variables with one cutpoints
###
### Copyright 2006-09 Anne-Laure Boulesteix
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

Ford2<-function(c,n0,n1,A,statistic)
{
N<-n0+n1
if (N!=sum(A))
 stop("error: you must have n0+n1=sum(A)")
K<-length(A)

if (K<=3)
 stop("when there are 3 categories, use Fcat")

b<-numeric(K-1)
b[1]<-choose(N,n1)*(1-Ford(c,n0=n0,n1=n1,A=A[c(1,K:2)],statistic=statistic))

uppervector<-numeric(N)
lowervector<-numeric(N)
for (i in 1:N)
 {
 boundaryi<-boundary(i,n0=n0,n1=n1,c=c,statistic=statistic,lower=TRUE)
 uppervector[i]<-boundaryi$upper
 lowervector[i]<-boundaryi$lower
 }

for (k in 2:(K-1))
 {
 Ak<-numeric(K)
 Ak[1:k]<-A[1:k]
 Ak[(k+1):K]<-A[K:(k+1)]

 bb<-numeric(K-k)
 # bb[i] is the number of paths crossing at k-1+i
 for (i in 1:(K-k))
  {
  bb[i]<-loop_ord2(A1=c(),k=k,n0=n0,n1=n1,A=Ak,whichboundary=k-1+i,c=c,statistic=statistic,uppervector=uppervector,lowervector=lowervector)
  }
 b[k]<-sum(bb)
 }

return(1-sum(b)/choose(N,n1))
}

##################



loop_ord2<-function(A1,k,n0,n1,A,whichboundary,c,statistic,uppervector,lowervector)
{

iter<-length(A1)+1
K<-length(A)
n_p<-sum(A[as.numeric(iter>1):(iter-1)])
n1_p<-sum(A1[as.numeric(iter>1):(iter-1)])


forbidden<-c()

candidate<-max(0,n1-n1_p-sum(A[-(1:iter)])):min(A[iter],n1-n1_p)

if (iter>k)
 {
 for (i in 2:k)
  {
  if (iter>(k+1))
   {
   a<-sum(A1[c(1:(i-1),(k+1):(iter-1))])
   }
  else
   {
   a<-sum(A1[c(1:(i-1))])
   }
  u<-floor(uppervector[sum(A[c(1:(i-1),(k+1):iter)])])+1-a
  l<-ceiling(lowervector[sum(A[c(1:(i-1),(k+1):iter)])])-1-a

  if (l>=0)
   {
   forbidden<-union(forbidden,0:l)
   }
  if (u<=A[iter])
   {
   forbidden<-union(forbidden,u:A[iter])
   }
  }
 }


if (iter==K)
 {
 if (is.element(n1-n1_p,forbidden))
  {
  number<-0
  }
 else
  {
  number<-choose(A[iter],n1-n1_p)
  }
 return(number)
 }


if (iter<whichboundary)
 {
 u<-floor(uppervector[sum(A[1:iter])])+1-n1_p
 l<-ceiling(lowervector[sum(A[1:iter])])-1-n1_p
 if (l>=0)
  {
  forbidden<-union(forbidden,0:l)
  }
 if (u<=A[iter])
  {
  forbidden<-union(forbidden,u:A[iter])
  }
 }

if (iter==whichboundary)
 {
 forbidden<-union(forbidden,(ceiling(lowervector[sum(A[1:iter])])-n1_p):(floor(uppervector[sum(A[1:iter])])-n1_p))
 }



authorized<-setdiff(candidate,forbidden)
if (length(authorized)==0)
  {
  number<-0
  }
 else
  {
  number<-sum(sapply(authorized,FUN=subloop_ord2,iter=iter,A1=A1,k=k,n0=n0,n1=n1,A=A,whichboundary=whichboundary,c=c,statistic=statistic,uppervector=uppervector,lowervector=lowervector))
  }
return(number)
}


#########################

subloop_ord2<-function(i,iter,A1,k,n0,n1,A,whichboundary,c,statistic,uppervector,lowervector)
{
number<-choose(A[iter],i)*loop_ord2(A1=c(A1,i),k=k,n0=n0,n1=n1,A=A,whichboundary=whichboundary,c=c,statistic=statistic,uppervector=uppervector,lowervector=lowervector)

return(number)
}

