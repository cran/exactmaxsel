###     Computation of the Gini gain from a 2x2 contingency table
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


ginigain<-function(mat)
{

N<-sum(mat)
N1<-sum(mat[1,])
N2<-sum(mat[2,])
Nl<-sum(mat[,1])
Nr<-sum(mat[,2])
gain<-2/N*(N2*N1/N-mat[2,2]*mat[1,2]/Nr-mat[2,1]*mat[1,1]/Nl)
return(gain)

}
