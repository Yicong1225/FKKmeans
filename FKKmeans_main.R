###################
### auxiliary functions
###################
library(vegan)
library(GUniFrac)
library(cluster)
library(dirmult)
library(MASS)

gen_solve <-  function(V) return(ginv(V))
sqrt_matrix<-function (V) {
  V.eig<-eigen(V)
  V.eigvalues<-V.eig$values
  Vsqrt.eigvalues<-V.eigvalues  
  for(i in 1:length(V.eigvalues)){
    if (V.eigvalues[i]<=0) {Vsqrt.eigvalues[i]=0
    } else {Vsqrt.eigvalues[i]=sqrt(V.eigvalues[i])}
  }
  V.sol=V.eig$vectors %*% diag(Vsqrt.eigvalues) %*% solve(V.eig$vectors)
  return(V.sol)
}

nor_matrix<-function (n) {return(diag(n)-matrix(rep(1,n),n,1)%*%t(matrix(rep(1,n),n,1))/n)}

gen_solve<-function(V) return(ginv(V))
tr = function(M){sum(diag(M))}
Trace = function(M){sum(diag(M))}

#' select the tuning parameter lambda
#'
#' @param y the outcome of interest
#' @param K kernel matrix
#' @param subID individuals' ID  
#' @param Lambda candidates of lambda
#'
#' @return A list of the selection of optimal lambda based on BIC and the corresponding 
#'
#' @examples
KMR.lamda.new = function(y, K, subID, Lambda){
  ###################################################
  ## tuning function of lambda parameter in KMR
  ## Input Lambda now is a vector of possible lambda values
  ## Need to figure out possible candidates
  m=length(Lambda)
  BICs=rep(NA,m)
  out=list()
  for(i in 1:m){
    KMR.tem <- KMR(y, K, subID, Lambda[i])
    BICs[i]=KMR.tem$BIC
    out[[i]]=KMR.tem$alphahat
  }
  index=order(BICs)[1]
  optlambda=Lambda[index]
  alphamat=out[[index]]
  return(list(optlambda = optlambda, BICs=BICs,  alphamat=alphamat))
}

mysolve<-function (V) {
  V.eig<-eigen(V)
  V.eigvalues<-V.eig$values
  Vi.eigvalues<-V.eigvalues
  for(i in 1:length(V.eigvalues)){
    if (V.eigvalues[i]==0) {Vi.eigvalues[i]=0
    } else {Vi.eigvalues[i]=1/V.eigvalues[i]}
  }
  V.sol=V.eig$vectors %*% diag(Vi.eigvalues) %*% solve(V.eig$vectors)
  return(V.sol)
}

get_K_tilde <- function(K, n_i_div) {
  N <- nrow(K); n <- length(n_i_div)
  K_tilde <- matrix(0, nrow=N, ncol=n*N)
  a <- 1; b <- 1 # initialize index for K_tilde
  a_K <- 1 # initialize index for K
  for (i in 1:n){
    c <- a + n_i_div[i] - 1; d <- b + N - 1 # current index for K_tilde
    b_K <- a_K + n_i_div[i] - 1 # current index for K
    
    K_tilde[a:c,b:d] <- K[a_K:b_K,]
    
    a_K <- b_K + 1 # update index for K
    a <- c + 1; b <- d + 1 # update index for K_tilde
  }
  return(K_tilde)
}

get_D <- function(K,n) {
  N <- nrow(K)
  D <- matrix(0, nrow=n*N, ncol=n*N)
  for (i in 1:n) {
    a <- N*(i-1) + 1; b <- N*i
    D[a:b,a:b] <- K
  }
  return(D)
}

mydsit2euclid<-function(D,sigma){
  nn = NROW(D)
  for (l in 1:nn)
    for (j in 1:nn)
    { if (l!=j){
      D[l,j]=sqrt(D[l,j]^2+2*abs(sigma))
    }
    }
  return (D)
}

scale2 = function(x) as.numeric(scale(x))  # This needs only once
#gauss kernel
sqare.diff.mat = function(x){  # Each column is a sample
  Kmat = matrix(NA, nrow=ncol(x), ncol=ncol(x))
  for (ii in 1:ncol(x)){Kmat[ii,] = diag(t(x-x[,ii])%*%(x-x[,ii]))}
  return(Kmat)
}
kernel.Gauss = function(Kmat, rho){
  return(exp(-Kmat/rho))
}

norm_minmax <- function(x){
  (x- min(x)) /(max(x)-min(x))
}

