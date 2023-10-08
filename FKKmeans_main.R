
# source("Function-copy.R")
# source("function_x_highdim.R")
# library(vegan)
# library(GUniFrac)
# library(cluster)
# library(dirmult)
# library(MASS)
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
# redundant case
KMR = function(y, K, subID, lambda){
  N=length(y)
  n=length(unique(subID))
  out=matrix(NA, nrow=n,ncol=N)  ## output of alpha-coef matrix
  res=NULL
  trhat=0
  for(i in 1:n){
    index=which(subID==unique(subID)[i])
    Ki=K[index,]
    yi=y[index]
    # get alphahat by individual's block 
    alpha=mysolve(t(Ki)%*%Ki + lambda*K)%*%t(Ki)%*%yi
    #alpha=mysolve(Ki%*%t(Ki) + lambda*K)%*%Ki%*%yi
    out[i,]=alpha
    res=c(res,yi-Ki%*%alpha)
    hati=Ki%*%mysolve(t(Ki)%*%Ki + lambda*K)%*%t(Ki) ## Hat-matrix
    #hati=Ki%*%mysolve(Ki%*%t(Ki) + lambda*K)%*%Ki
    trhat=trhat+Trace(hati)
  }
  RSS=sum(res^2)
  BIC=log(RSS)+trhat*log(N)/N
  return(list(alphahat = out, BIC=BIC))
}

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
  ## ?? Need to figure out possible candidates
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

############################################################################
#   SIMULATION under 3 Scenarios in the paper
#
#  You have to choose one scenario parameter 、
#  (default:  SC = 1)
# 
# 
############################################################################

SC <- 1 ## Scenario Num (1:normal, 2:X causal, 3:T causal)
s <- 3; noise.sd <- 0.2 ## Noise level
m <- 10 	## No. of individual (load with common file)
k <- 3 ## No. of repeated measures within a cluster/subject
n <- m * k	## Total sample size
subID=gl(m,k)
p <- 856 ## always use 856 and select col randomly inside rep loop
pt <- 1
rep <- 1 # replication times

cluster_number <- c(rep(0,n/4),rep(1,n/4),rep(2,n/4),rep(3,n/4))

# make preparations for generating X 
if (SC == 1) { # Scenario 1
  theta <- rep(0, k)
  Sigma <- matrix(1, nrow = k, ncol = k)
  
  for (ii in 1:k) {
    for (jj in 1:k) {
      Sigma[ii, jj] <- 0.5^abs(ii - jj)
    }
  }
  message("Normal")
} else if (SC == 2) { # Scenario 2
  data(throat.tree)
  data(throat.otu.tab)
  nClus = 20
  depth = 10000
  tree = throat.tree  # tree = midpoint(tree)
  tree.dist = cophenetic(tree)
  obj <- pam(tree.dist, nClus)
  clustering <- obj$clustering
  otu.ids <- tree$tip.label
  load("DirMultOutput.RData")
  p.est = dd$pi						## estimation of OTU proportions from a real data
  p.est1 = p.est[order(p.est,decreasing=T)[1:p]]  	## take the first p of 856 OTUs as our OTU-design matrix
  p.est1 = p.est1/sum(p.est1)
  names(p.est1) <- names(dd$pi[order(p.est,decreasing=T)[1:p]])
  theta <- dd$theta
  gplus <- (1 - theta) / theta
  g.est <- p.est1 * gplus
  comm <- matrix(NA, n , length(g.est)) 
  rownames(comm) <- 1:nrow(comm)
  colnames(comm) <- names(g.est)
  
  message("X as causal")
} else if (SC == 3) { # Scenario 2
  message("T as causal")
} else {
  stop("Scenario Num not supported")
}


ri <- matrix(0,rep,2)# ARI for L2-norm/RKHS approach
ru <- matrix(0,rep,1)# ARI for K-means on Y approach
r_admm <- matrix(0,rep,1) # ARI for ADMM approach
ri_MLD <- matrix(0,rep,1) # ARI for MLD approach
time_collection <- matrix(0,rep,4) # collect the computation time
num_group_collection <- matrix(0,rep,4) # collect the estimated number of groups

k_out <- matrix(0,2*rep,m)
m_out <- matrix(0,1*rep,m)
admm_out <- matrix(0,1*rep,m)
admm_collection <- matrix(0,1*rep,2)

aa<-permutations(4,4,seq(1,4))
myclusters<-matrix(rep(0,24*m),ncol=m,nrow=24)
for(r in 1:24) {
  myclusters[r,] <- rep(aa[r,], times = 1, length.out = NA, each = m/4)
}

i <- 1
while (i <= rep) {
  try({
    # generate X 
    if (SC == 1) { 
      X.matrix <- MASS::mvrnorm(m, theta, Sigma)
      Xt <- matrix(t(X.matrix), nrow = n, ncol = 1)
      Xt <- cbind(Xt,rep(0,n)) # n by 1 matrix couldn't be handled by the following normalization method
      ## normalize
      ## normalize comm by columns & every subjects
      for(r in 1:m){
        s1 <- k*(r-1)+1
        s2 <- k*r
        Xt[s1:s2,] <- apply(Xt[s1:s2,],2,norm_minmax)
      }
      Xt <- Xt[,1]
    } else if (SC == 2) {
      ## Data generation scheme 2: A proportion X-design matrix
      for (h in 1:m) {
        comm.p = rdirichlet(1, g.est)[1, ]
        comm[(k*(h-1)+1),]=comm.p
        for(j in 2:k){
          perb=rep(1,p)+rnorm(p,0,0.3)
          perb[perb<0]=0
          comm.p1=comm.p*perb/sum(comm.p*perb)
          comm[(k*(h-1)+j),] =comm.p1
          comm.p=comm.p1
        }
      }
      
      ## normalize
      ## normalize comm by columns & every subjects
      for(r in 1:m){
        s1 <- k*(r-1)+1
        s2 <- k*r
        comm[s1:s2,] <- apply(comm[s1:s2,],2,norm_minmax)
      }
      
      random_col <- sample(1:ncol(comm),1)
      Xt=comm[,random_col]
      
    } else if (SC == 3) {
      time.x <- rep(seq(0,1,length.out=k),m)
      Xt <- time.x
    }

    X=Xt
    ## generate Y
    Y = case_when(cluster_number == 0 ~ as.numeric(-1.5*X),
                  cluster_number == 1 ~ as.numeric(-1.5*X+1.5),
                  cluster_number == 2 ~ as.numeric(cos(2*pi*(X))),
                  cluster_number == 3 ~ as.numeric(1-2*exp(-6*X)))+rnorm(n,0,noise.sd)
    y=Y
    
    
    ## calculate cluster
    public.start <- Sys.time()
    
    ## Gaussian kernel
    dMat = sqare.diff.mat(t(Xt)) ## Remark. P=1, adjustment needed here
    rho=median(dMat)
    K=kernel.Gauss(dMat,rho)
    
    min=0
    max=1
    Lambda=seq(min,max,by=(max-min)/10)
  
    alpha_0_mat <-KMR.lamda.new(y,K,subID,lambda_1)$alphamat
    lambda_1 = KMR.lamda.new(y,K,subID,lambda_1)$optlambda 
    alpha_0 <- matrix(t(alpha_0_mat), ncol=1) ## convert matrix form alpha_0 to vector
    
    ## L2-norm, RKHS, K-means on Y
    ## each sample fit separately
    alpha.l2<-matrix(rep(0,1*n*m),ncol=n,nrow=m)
    
    ## l2 norm k-means
    for (r in 1:m){
      c1<-k*(r-1)+1
      c2<-k*r
      alpha.l2[r,]<-gen_solve(t(K[c1:c2,])%*%K[c1:c2,] + lambda_1*K)%*%t(K[c1:c2,])%*%y[c1:c2]
    }
    
    public.end <- Sys.time()
    
    l2.start <- Sys.time()
    gap_stat.l2 <- clusGap(alpha.l2, FUN = kmeans, nstart = 20, K.max = 10, B = 50)
    best_k.l2 <- maxSE(gap_stat.l2$Tab[, "gap"], gap_stat.l2$Tab[, "SE.sim"], method = "globalSEmax")
    #centers = best_k.l2 
    
    kmeans1 <- kmeans(alpha.l2,centers=best_k.l2,nstart=50,iter.max = 200)
    l2.end <- Sys.time()
    
    ## rkhs
    rkhs.start <- Sys.time()
    
    alpha.rkhs<-matrix(rep(0,1*n*m),ncol=n,nrow=m)
    alpha.rkhs=t(sqrt_matrix(K)%*%t(alpha.l2))
    
    gap_stat.rkhs <- clusGap(alpha.rkhs, FUN = kmeans, nstart = 20, K.max = 10, B = 50)
    best_k.rkhs <- maxSE(gap_stat.rkhs$Tab[, "gap"], gap_stat.rkhs$Tab[, "SE.sim"], method = "globalSEmax")
    #centers = best_k.rkhs
    
    kmeans2 <- kmeans(alpha.rkhs,centers=best_k.rkhs,nstart = 50,iter.max = 200)
    rkhs.end <- Sys.time()
    
    clu1 <- kmeans1[["cluster"]]
    clu2 <- kmeans2[["cluster"]]

    k_out[2*i-1,] <- clu1 #cluster outcome for l2 norm
    k_out[2*i,] <- clu2 #cluster outcome for RKHS norm
  })
}

# ## get ARI
# for (i in 1:rep) {
#   b1 <- rep(0,24);b2 <- b1;b3 <- b1; b4 <- b1
#   for (j in 1:24) b1[j]=adjustedRandIndex(k_out[2*i-1,],myclusters[j,])
#   ri[i,1]=max(b1)
#   for (l in 1:24) b2[l]=adjustedRandIndex(k_out[2*i,],myclusters[l,])
#   ri[i,2]=max(b2)
#   for (r in 1:24) b4[j]=adjustedRandIndex(m_out[i,],myclusters[r,])
#   ri_MLD[i,1]=max(b4)
# }


