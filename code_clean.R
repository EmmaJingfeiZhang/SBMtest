library(mgcv)
library(RMTstat)
library(rARPACK)
library(reliaR)
library(maxLik)
library(blockmodels)


Generate.A <- function(P) {
  # generates random adjacency matrix
  n <- dim(P)[1]
  upper.tri.ind <- upper.tri(P)
  p.upper <- P[upper.tri.ind]
  A.upper <- rbinom(n*(n-1)/2, 1, p.upper)
  A <- matrix(0, ncol = n, nrow = n)
  A[upper.tri.ind] <- A.upper
  A <- A + t(A)
  return(A)
}


ClustVec2Mat <- function(clust.vec, K) {
  # convert a membership vector to a membership matrix, the number of different values in
  # clust.vec must be at most K
  clust.mat <- matrix(0, ncol = K, nrow = length(clust.vec))
  for (i in 1:K) {
    clust.mat[clust.vec==i,i] <- 1
  }
  return(clust.mat)
}


SpecClust <- function(A, K, nstart = 100) {
  # spectral clustering
  n <- nrow(A)
  if (K == 1) {
    return(rep(1, n))
  }
  U <- slanczos(A, K)$vectors
  return(kmeans(U, K, nstart = nstart)$cluster)
}



EmpB<- function(A,c) {
  # empirical connectivity matrix
  k <- max(c)
  C<-matrix(0,length(c),k)
  for(k in 1:k){C[c==k,k]<-1}
  size<-apply(C,2,sum)
  msize<-size%o%size
  diag(msize)<-size*(size-1)
  if (k == 1) {
    n <- nrow(A)
    B <- sum(A)/(n*(n-1))
  } else {
    B <- t(C)%*%A%*%C/msize
  }
  return(B)
}


llf <- function(param) {
  mu <- param[1 ]
  sigma <- param[2]
  llValue <- dgumbel(x, mu, sigma, log=TRUE)
  return(sum(llValue))
}


GoFTest<- function(A, c, K0) {
#function for calcuating the test statistics
  n <- nrow(A)
  clust.size <- rep(0, K0)
  B.hat <- EmpB(A,c)
  P.hat<-matrix(0,n,n)
  clust.size<-table(c)
  if(K0==1) P.hat=matrix(B.hat,n,n)
  if(K0>1) P.hat=B.hat[c,c]
  tilde.A <- (A - P.hat+diag(diag(P.hat))) / sqrt(P.hat*(1-P.hat))
  L.test<-matrix(0,n,K0)
  for(j in 1:K0){
      L.test[,j]<-abs(apply(tilde.A[,c==j],1,sum))/sqrt(clust.size[j]-1)
  }
  L.test[is.na(L.test)]<-0
  max.L.test=max(L.test)
  gumbel.test<-max.L.test^2-2*log(2*K0*n)+log(log(2*K0*n))
  gumbel.test
}


GoFTestBoot <- function(A, c, K.test=3, n.b = 100) {
#function for calcuating the bootstrap corrected test statistics
  N <- nrow(A)
  B.hat <- EmpB(A,c)
  P.hat <- B.hat[c,c]
  gumbel.test<-GoFTest(A,c,K.test)

  mu.gum<- -2.531024
  sd.gum <- 2
  test.vec <- rep(0, n.b)
  for (i in 1:n.b) {
    A.i<-Generate.A(P.hat) ## adjacency matrix
    diag(A.i)<-0
    test.vec[i] <- GoFTest(A.i,c,K.test)
  }
  llf <- function(param) {
    mu <- param[1 ]
    sigma <- param[2]
    llValue <- dgumbel(x, mu, sigma, log=TRUE)
    return(sum(llValue))
  }
  x<- test.vec
  A <- matrix(c(0, 1), 1, 2); B <- 0
  ml <- maxLik(llf, start = c(mu=-mean(test.vec), sigma=sd(test.vec)),constraints=list(ineqA=A, ineqB=B))
  mu.test.hat <- ml[[2]][1]
  sig.test.hat <- ml[[2]][2]
  test.stat.boot <- mu.gum + sd.gum*(gumbel.test - mu.test.hat)/sig.test.hat
  return(list(test.stat = gumbel.test, test.stat.boot = test.stat.boot))
}



######################################################
#################### SBM Examples  ###################
######################################################

## generate a network with three communities of sizes (400,400,500)
K=3
npc<-c(400,400,500)
r<-4;rho <- 0.1
srange<-rep(npc,K)
c<-unlist(sapply(1:K, function(x) rep(x,srange[x]))) 
P = matrix(1, K, K) + diag(rep(r, K))
P = rho*P
N = length(c)
A = matrix(0, ncol = N, nrow = N)
C<-matrix(0,N,K); for(k in 1:K){C[c==k,k]<-1}
A<-Generate.A(C%*%P%*%t(C))


## cluster under K.test
K.test=3
chat <- SpecClust(A, K.test)

## p-value of Tn ##
gumbel.test1<-GoFTest(A, chat, K.test)
pgumbel(gumbel.test1, -2.531024, 2, lower.tail = F)


## p-value of Tn_boot ##
gumbel.test2<-GoFTestBoot(A, chat, K.test, n.b=100)[[2]]
pgumbel(gumbel.test2, -2.531024, 2, lower.tail = F)



## calculate augmented test statistics
## the augmented test statistic is recommended for the planted partition model
Bhat<-EmpB(A,chat)
npc_aug<-round(min(table(chat))/2)
A12<-matrix(rbinom(N*npc_aug, 1, min(Bhat)/2),npc_aug,N)
A22.upper <- rbinom(npc_aug*(npc_aug-1)/2, 1, max(Bhat))
A22 <- matrix(0, ncol = npc_aug, nrow = npc_aug)
A22[upper.tri(A22)] <- A22.upper
A22 <- A22 + t(A22)
A_aug<-cbind(rbind(A22,t(A12)),rbind(A12,A))


## cluster under K.test+1
chat2 <-SpecClust(A_aug, K.test+1)

## p-value of Tn+ ##
gumbel.test3<-GoFTest(A_aug, chat2, K.test+1)
pgumbel(gumbel.test3, -2.531024, 2, lower.tail = F)

## p-value of Tn_boot+ ##
gumbel.test4<-GoFTestBoot(A_aug, chat2, K.test+1, n.b=100)[[2]]
pgumbel(gumbel.test4, -2.531024, 2, lower.tail = F)

