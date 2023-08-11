library('rARPACK')
library('Matrix')
library(roxygen2)
library(quadprog)

source("r/vertex_hunting_functions.R")
source("r/simplex_dist.R")



score <- function(D, K, scatterplot=FALSE, K0=NULL, m=NULL, N=NULL, threshold=FALSE,
                  Mquantile=0.05, VHMethod = 'SP', normalize="none",
                  alpha=0.5, max_K=150, returnW=FALSE){
  #' This function computes the estimates for the A and W matrix based on the algorithm proposed in Ke and Wang's work: https://arxiv.org/pdf/1704.07016.pdf
  #' 
  #'
  #' @param D the document term matrix. Can be a dense or sparse matrix. The dimensions of this matrix should be p x n, with
  #'          p the number of words in the dictionary, and n the number of documents. D is assumed to have been appropriately 
  #'          preprocessed beforehand.
  #' @param K the number of topics. If null, it is computed from the data using the SVD decomposition of D.
  #' @param normalize Normalization method used on the topic document matrix. Choices include "none" (no normalization), "norm" (row-wise mean normalization), 
  #'                  and "norm_score_N" (mean/N).
  #' 
  #' @param N the document length. Only necessary for the prenormalization option "norm_score_N" and
  #' @param threshold boolean. Should the words in the dictionary be thresholded or not (corresponds to a sparsity assumption on A)
  #' @param VHMethod vertex hunting method. Choice between SVS (sketched vertex search), SP (successive projections) and SVS-SP. 
  #'         Note that the paper by Tracy Ke recommends using SVS, as it tends to be more robust to noise and outliers. SVS requires additional parameters K0 and m.
  #' @param m parameter for the SVS vertex hunting procedure. Corresponds to the number of clusters used in the denoising step of SVS (clusters are taken as
  #'          a proxy for the r computed in the previous step, and the averaging is supposed to make them less noisy).
  #'          Ke recommends using m >>> K, but smaller than p. The value m = 10 * K is used in the paper.
  #' @param K0 parameter for the SVS vertex hunting procedure. Corresponds to the number of clusters closest to the overall 
  #'           cluster grand mean (amongst the m) that are kept. The value K0= ceiling(1.5 * K) is recommended in the paper and used as default here.
  #' @param Mquantile quantile of M to cap the lowest values of M to. In real data analysis, it is sometimes beneficial to use a positive value of quantile tau, 
  #'                  to avoid over-weighting those extremelylow-frequency words in the pre-SVD normalization
  #'
  #' @return The square of the input value.
  #'
  #' @examples
  #' square(4)
  #' 
  #'
  #' @export
  #'
  #' 
  
  p <- dim(D)[1]
  n <- dim(D)[2]
  print(c(n, p, K, N))
  
  M <- rowMeans(D)   #### average frequency at which each word appears in each document
  M_trunk <- sapply(M,function(x){max(quantile(x, Mquantile))})
  if(normalize == "norm_score_N"){
      M2 <- rowSums(D/t(matrix(rep(N,p),n,p)))
  }
  
  if (is.null(K0)){
    K0 <- ceiling(1.5 * K) 
  }
  if (is.null(m)){
    m <- ceiling(10 * K) 
  }
  
  ### Step 0: Pre-SVD thresholding and selection of words (optional) 
  if (threshold){
    threshold_J = alpha * sqrt(log(max(p,n))/(N *n))
    print(sprintf("Threshold for alpha = %f  is %f ", alpha, threshold_J))
    setJ = which(M > threshold_J)
    print(sprintf("Nb of elected words = %i  ( %f percent) ", length(setJ), 100 * length(setJ)/p))
    if (length(setJ) < 0.1 * length(M)){
      setJ = sort(M, decreasing=TRUE, index.return=TRUE)$ix[1:ceiling( 0.1 * length(M))]
    }
    newD = D[setJ,]
    M = M[setJ]
    if(normalize == "norm_score_N"){
      M2 = M2[setJ]
    }
    M_trunk = M_trunk[setJ]
    print(paste0(p-length(setJ), " words were thresholded (", (p-length(setJ))/p * 100, "%)"))
    print(paste0(length(setJ), " words remain"))
    new_p <- length(setJ)
    
  }else{
    new_p = p
    newD=D
    setJ = 1:length(M)
  }
  
  
  newD <- switch(normalize, 
                 "norm" = diag(sqrt(M_trunk^(-1))) %*% newD,
                 "norm_score_N" = diag(sqrt(M2^(-1)))%*% newD,
                 "huy" = newD %*% t(newD) - n/N * diag(M),
                 "none" = newD)
  if (K >= min(dim(newD))){
    obj = svd(newD, min(K, min(dim(newD))))
  }else{
    obj = svds(newD, K)
  }
  if (max_K >= min(dim(newD))){
    obj_full = svd(newD, min(max_K, min(dim(newD))))
  }else{
    obj_full = svds(newD, max_K)
  }
  eigenvalues = obj_full$d
  Xi <- obj$u
  
  
  #Step 1: SVD
  Xi[,1] <- abs(Xi[,1]) ### to get rid of some small irregularities, but should be positive and bounded away from 0 (this is guaranteed by Perron's theorem)
  R <- apply(as.matrix(Xi[,2:K]),2,function(x) x/Xi[,1])
  
  #Step 2: Post-SVD normalization
  print("Start VH")
  if (VHMethod == 'SVS'){
    vertices_est_obj <- vertices_est(R, K0, m)
    V <- vertices_est_obj$V
    theta <- vertices_est_obj$theta
  } else if (VHMethod == 'SP'){
    vertices_est_obj <- successiveProj(R, K)
    V <- vertices_est_obj$V
    theta <- NULL
  } else if (VHMethod == 'SVS-SP'){
    vertices_est_obj <- vertices_est_SP(R, m)
    V <- vertices_est_obj$V
    theta <- NULL
  }else if (VHMethod == 'AA'){
    vertices_est_obj <- ArchetypeA(r_to_py(R),as.integer(K))
    V <-vertices_est_obj$V
    theta<-NULL

  }
  
  if (scatterplot){
    par(mar=c(1,1,1,1))
    plot(R[,1],R[,2])
    points(V[,1],V[,2],col=2,lwd=5)
  }
  
  
  #Step 3: Topic matrix estimation
  print("Start Step 3")
  if(rankMatrix(cbind(V,rep(1,K)))[1]<K){
    Pi <- cbind(R, rep(1,new_p)) %*% MASS::ginv(cbind(V,rep(1,K)))
  }else{
    Pi <- cbind(R, rep(1,new_p)) %*% solve(cbind(V,rep(1,K)))
  }

  Pi <- pmax(Pi,matrix(0,dim(Pi)[1],dim(Pi)[2])) ### sets negative entries to 0 
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 3b 
  if (normalize %in% c("norm", "norm_score_N" )){
    A_hat <- switch(normalize, 
                    "norm" = diag(sqrt(M_trunk)),
                    "norm_score_N" = diag(sqrt(M2))) %*% (Xi[,1]*Pi)
    
  }else{
    A_hat <- Xi[,1]*Pi
  }

  #Step  3c: normalize each column to have a unit l1 norm
  temp <- colSums(A_hat)
  A_hat <- t(apply(A_hat,1,function(x) x/temp))
  
  A_hat_final = matrix(0, p, K )
  if(threshold){
      A_hat_final[setJ, ] = A_hat
      if(returnW){
        print("Start estimation of W")
        W_hat <- compute_W_from_AD(A_hat, D[setJ,])
      }else{
        W_hat <- NULL
      }
    }else{
      A_hat_final = A_hat
      if(returnW){
        W_hat <- compute_W_from_AD(A_hat, D)
      }else{
        W_hat <- NULL
      }
    }
  
  
  

  #Step 4
  
  return(list(A_hat=A_hat_final, R=R,V=V, Pi=Pi, theta=theta, W_hat = W_hat,
              eigenvalues = eigenvalues, thresholded = 1- length(setJ)/ dim(D)[1]))
}



compute_W_from_AD <- function(A_hat, D){
  #' This function computes the estimates for the  W matrix when given an estimate of A based on a weighted least-squares algorithm, 
  #' where the weights come from the normalizing factors in the pre-SVD normalization and aim to tackle severe frequency heterogeneity
  #'
  #' @param D the document term matrix. Can be a dense or sparse matrix. The dimensions of this matrix should be p x n, with
  #'          p the number of words in the dictionary, and n the number of documents. D is assumed to have been appropriately 
  #'          preprocessed beforehand.
  #' @param A_hat the estimated A_hat matrix
  #' 
  
  K <-dim(A_hat)[2]
  n <- dim(D)[2]
  
  W_hat <- matrix(0, K, n)
  M <- rbind(diag(K-1), rep(-1,K-1))
  bM <- diag(K)[,K]
  Dmat <- 2*t(A_hat %*% M) %*% (A_hat %*% M)
  Amat <- t(M)
  bvec <- -bM
  
  AM <- A_hat %*% M
  AbM <- A_hat %*% bM
  for (i in 1:n){
    dvec <- 2*t(D[,i]-AbM) %*% AM
    # Dmat <- matrix(nearPD(Dmat)$mat, nrow(Dmat), ncol(Dmat))
    # Dmat <- nearPD(Dmat)
    qp_sol <- solve.QP(Dmat, dvec, Amat, bvec)$solution
    W_hat[,i] <- c(qp_sol, 1-sum(qp_sol))
  }
  W_hat <- pmax(W_hat,0)
  
  return(W_hat)
}







