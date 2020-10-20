setwd("/home/diana/Workspace/cnww/cancer_networks/")
S <- read.csv("data/gene_depedency_rank_corr.csv", row.names = 1)
S <- data.matrix(S)
make.affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])
  
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A  
}

A <- make.affinity(S, 3)  # use 3 neighboors (includes self)
A[1,]
D <- diag(apply(A, 1, sum))
D[1,]
U <- D - A
round(U[1:12,1:12],1)
evL <- eigen(U, symmetric=TRUE)
signif(evL$values,2)

library(ggplot2)
plot(1:10, rev(evL$values)[1:10], pch=20, 
     xlab = "Eigenvalues of Graph Laplacian", ylab = "")
