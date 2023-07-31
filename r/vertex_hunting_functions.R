replaceWithLeastPositive <- function(vec){
  vec[vec<=0] = min(vec[vec>0])
  return(vec)
}

nearPD <- function(mat){
  mat <- (mat+t(mat))/2
  eigenObj <- eigen(mat)
  values <- eigenObj$values
  vectors <- eigenObj$vectors
  values <- replaceWithLeastPositive(values)
  return(vectors%*%diag(values)%*%t(vectors))
}

successiveProj <- function(R, K){
  # succesive projection on rows of R
  n <- dim(R)[1]
  
  Y <- cbind(rep(1,n),R)
  indexSet = c()
  it = 1
  valid_indices <- 1:dim(Y)[1]
  while (length(indexSet) < K & dim(Y)[1]>0){
    print(paste0("here it is ", it))
    print(sprintf("Length set is %d and dim(Y)[1]=%d", length(indexSet), dim(Y)[1] ))
    l2Norms <- apply(Y,1,function(x) sqrt(sum(x^2)))
    #### check if they are in the same document
    index <- which(l2Norms == max(l2Norms))
    
    if (length(index) >1){
      print("yoho")
      r = rankMatrix(R[valid_indices[index], ])[1]
      if (r < length(index)){
        ### only select 1
        chosen = sample(index,  r)
        u <- Y[chosen,] / sqrt(sum(Y[chosen,]^2))
        indexSet <- c(indexSet, valid_indices[chosen])
        Y <- Y[-setdiff(index, c(chosen)), ]
        valid_indices <- setdiff(valid_indices, valid_indices[setdiff(index, c(chosen))])
        Y <- t(apply(Y, 1, function(x) x-sum(x%*%(u))%*%t(u)))
      }else{
        u <- Y[index,] / sqrt(sum(Y[index,]^2))
        indexSet <- c(indexSet, valid_indices[index])
        Y <- t(apply(Y, 1, function(x) x-sum(x%*%(u))%*%t(u)))
      }
    }else{
      indexSet <- c(indexSet, valid_indices[index])
      u <- Y[index,] / sqrt(sum(Y[index,]^2))
      Y <- t(apply(Y,1,function(x) x-sum(x*u)*u))
    }
    it = it + 1

    
  }
  if(length(indexSet)>K){
       indexSet = sample(indexSet, K)
  }
  
  return(list(V=R[indexSet,], indexSet=indexSet))
}

vertices_est_SP <- function(R,m){
  library(quadprog)
  K <- dim(R)[2] + 1
  
  obj <- kmeans(R,m,iter.max=K*100,nstart = K*10)
  theta <- as.matrix(obj$centers)
  return(successiveProj(theta, K))
}


vertices_est <- function(R, K0, m){
  library(quadprog)
  K <- dim(R)[2] + 1
  #Step 2a
  print(paste0("min is ", min(m, nrow(R))))
  print(c(m, dim(R)))
  obj <- kmeans(R, min(m, nrow(R)), iter.max=K*100, nstart = K*10)
  
  theta <- as.matrix(obj$centers)
  theta_original <- theta
  plot(R[,1],R[,2])
  points(theta[,1], theta[,2], col=2,lwd=4)
  
  #Step 2b'
  inner <- theta%*%t(theta)
  distance <- diag(inner)%*%t(rep(1,length(diag(inner)))) + rep(1,length(diag(inner)))%*%t(diag(inner)) - 2*inner
  top2 <- which(distance==max(distance),arr.ind=TRUE)[1,]
  theta0 <- as.matrix(theta[top2,])
  theta <- as.matrix(theta[-top2,])
  
  if (K0 > 2){
    for (k0 in 3:K0){
      inner <- theta%*%t(theta)
      distance <- rep(1,k0-1)%*%t(diag(inner))-2*theta0%*%t(theta)
      ave_dist <- colMeans(distance)
      index <- which(ave_dist==max(ave_dist))[1]
      theta0 <- rbind(theta0, theta[index,])
      theta <- as.matrix(theta[-index,])
    }
    theta <- theta0
  }
  
  #Step 2b
  if (K0<=K){
    #plot(theta[,1],theta[,2])
    #points(theta[comb[,min_index],1],theta[comb[,min_index],2],col=2,pch=2)
    
    return(list(V=theta, theta=theta_original))
  }else{
      comb <- combn(1:K0, K)
      max_values <- rep(0, dim(comb)[2])
      for (i in 1:dim(comb)[2]){
        for (j in 1:K0){
          max_values[i] <- max(simplex_dist(as.matrix(theta[j,]), as.matrix(theta[comb[,i],])), max_values[i])
        }
      }
      min_index <- which(max_values == min(max_values))
      new_theta = theta[comb[,min_index[1]],]
    
    #plot(theta[,1],theta[,2])
    #points(theta[comb[,min_index],1],theta[comb[,min_index],2],col=2,pch=2)
    
    return(list(V=new_theta, theta=theta_original))
  }
 
}