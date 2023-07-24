select_K <- function(eigenvals, p,n, N, method="huy"){
  if (method == "huy"){
    thresh = 8 * log(min(p, n)) * sqrt( n * log(min(p, n))/ N)
    return(K)
  }
  if(method == "olga"){
    thresh = 4 *  sqrt( n * log(p + n)/ N)
  }
  if(method == "tracy"){
    if (n > max(N *p^2, p^3, N^2 * p^5)){
      beta_n = 1 + min(p/N, p^2/N^(3/2))
    }else{
      beta_n = 1 + p^2/N^(3/2)
    }
    g_n = log(log(n))
    thresh = n/N  +  beta_n * sqrt( n * p * log(n)/ N) * gn
  }else{
    print("Not implemented yet")
    return(NULL)
  }
  K = sort(eigenvals[which(eigenvals > thresh)])[1]
}