
###     Computes the coordinates of the boundaries 
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

boundary<-function(x,n0,n1,c,statistic, lower=TRUE)
{
N<-n0+n1
xN<-x/N
boundary<-list()
if (statistic=="chi2")
 {
 if (lower==FALSE)
  {
  boundary$upper<-n1*x/N+n0*n1*sqrt(c)/N*sqrt(xN*(1-xN)*(1/n0+1/n1))
  return(boundary)
  }
 else
  {
  a<-n1*x/N
  b<-n0*n1*sqrt(c)/N*sqrt(xN*(1-xN)*(1/n0+1/n1))
  boundary$upper<-a+b
  boundary$lower<-a-b
  return(boundary)
  } 
 }

if (statistic=="gini")
 {
 if (lower==FALSE)
  {
  boundary$upper<-n1*x/N+x*(N-x)*sqrt(8*c/(x*(N-x)))/4
  return(boundary)
  }
 else
  {
  a<-n1*x/N
  b<-x*(N-x)*sqrt(8*c/(x*(N-x)))/4
  boundary$upper<-a+b
  boundary$lower<-a-b
  return(boundary)
  } 
 } 
 
return(boundary)
}
