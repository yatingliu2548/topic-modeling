
select_K <- function(eigenvals, p,n, N, method="huy"){
  if (method == "huy"){
    thresh = 8 * log(min(p, n)) * sqrt( n * log(min(p, n))/ N)
  }else{
    if(method == "olga"){
      thresh = 4 * sqrt( n * log(p + n)/ N)
    }else{
      if(method == "tracy"){
        if (n > max(N *p^2, p^3, N^2 * p^5)){
          beta_n = 1 + min(p/N, p^2/N^(3/2))
        }else{
          beta_n = 1 + p^2/N^(3/2)
        }
        g_n = log(log(n))
        thresh = sqrt(n/N  +  beta_n * sqrt( n * p * log(n)/ N) * g_n)
        #print(c(n/N, beta_n, g_n, sqrt( n * p * log(n)/ N), beta_n * sqrt( n * p * log(n)/ N) * g_n ))
      }else{
        print("Not implemented yet")
        return(NULL)
      }
    }
  }
  print(thresh)
  if(length(which(eigenvals> thresh))==0){
    return(c(list(Khat=1, thresh= thresh)))
  }else{
    K = sort(which(eigenvals> thresh), decreasing = TRUE)[1]
    return(c(list(Khat=K, thresh= thresh)))
  }
}
