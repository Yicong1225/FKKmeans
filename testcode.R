############################################################################
#   a toy example for Scenario 1
############################################################################
source("FKKmeans_main.R")
s <- 3
noise.sd <- 0.2 ## Noise level
m <- 10 	## No. of individual (load with common file)
k <- 3 ## No. of repeated measures within a cluster/subject
n <- m * k	## Total sample size
subID=gl(m,k)

cluster_number <- c(rep(0,n/4),rep(1,n/4),rep(2,n/4),rep(3,n/4))
aa<-permutations(4,4,seq(1,4))
myclusters<-matrix(rep(0,24*m),ncol=m,nrow=24)
for(r in 1:24) {
  myclusters[r,] <- rep(aa[r,], times = 1, length.out = NA, each = m/4)
}

theta <- rep(0, k)
Sigma <- matrix(1, nrow = k, ncol = k)

for (ii in 1:k) {
  for (jj in 1:k) {
    Sigma[ii, jj] <- 0.5^abs(ii - jj)
  }
}

X.matrix <- MASS::mvrnorm(m, theta, Sigma)
Xt <- matrix(t(X.matrix), nrow = n, ncol = 1)
Xt <- cbind(Xt,rep(0,n)) # n by 1 matrix couldn't be handled by the following normalization method

## normalize
for(r in 1:m){
  s1 <- k*(r-1)+1
  s2 <- k*r
  Xt[s1:s2,] <- apply(Xt[s1:s2,],2,norm_minmax)
}
Xt <- Xt[,1]
X=Xt
## generate Y
Y = case_when(cluster_number == 0 ~ as.numeric(-1.5*X),
              cluster_number == 1 ~ as.numeric(-1.5*X+1.5),
              cluster_number == 2 ~ as.numeric(cos(2*pi*(X))),
              cluster_number == 3 ~ as.numeric(1-2*exp(-6*X)))+rnorm(n,0,noise.sd)
y=Y

## Gaussian kernel
dMat = sqare.diff.mat(t(Xt)) ## Remark. P=1, adjustment needed here
rho=median(dMat)
K=kernel.Gauss(dMat,rho)

# tune lambda
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

# Kmeans method (run if necessary)
# gap_stat.l2 <- clusGap(alpha.l2, FUN = kmeans, nstart = 20, K.max = 10, B = 50)
# best_k.l2 <- maxSE(gap_stat.l2$Tab[, "gap"], gap_stat.l2$Tab[, "SE.sim"], method = "globalSEmax")
# centers = best_k.l2 
# kmeans1 <- kmeans(alpha.l2,centers=best_k.l2,nstart=50,iter.max = 200)

alpha.rkhs<-matrix(rep(0,1*n*m),ncol=n,nrow=m)
alpha.rkhs=t(sqrt_matrix(K)%*%t(alpha.l2))

gap_stat.rkhs <- clusGap(alpha.rkhs, FUN = kmeans, nstart = 20, K.max = 10, B = 50)
best_k.rkhs <- maxSE(gap_stat.rkhs$Tab[, "gap"], gap_stat.rkhs$Tab[, "SE.sim"], method = "globalSEmax")
kmeans2 <- kmeans(alpha.rkhs,centers=best_k.rkhs,nstart = 50,iter.max = 200)
clu2 <- kmeans2[["cluster"]]

## get ARI
b1 <- rep(0,24);b2 <- b1;b3 <- b1; b4 <- b1
for (j in 1:24) b1[j]=adjustedRandIndex(clu2,myclusters[j,])
ri=max(b1)


