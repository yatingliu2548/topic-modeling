#### Generate Sigma 
import numpy as np

def Cal_lbd(Ns, X_freq, D_X, P, K, anchor_vec):
  p = len(D_X)
  fixed_rate = K * np.sqrt(K * np.log(P))

  if len(Ns) == 1:
    return fixed_rate * np.sqrt(p / Ns[0] / min(D_X[anchor_vec]))
  else:
    X_Ni = X_freq[anchor_vec,] / Ns
    sum_weighted = np.sum(X_Ni, 1)
    return fixed_rate * np.sqrt(min(p ** 2 * sum_weighted / D_X[anchor_vec] ** 2))



def GetSigma(X, Ns):
    [n, p] = [len(Ns), X.shape[0]]
    ratio = Ns / (Ns - 1)
    X_Ns = X * ratio
    diag_X_Ns = np.mean(X / (Ns - 1), 1)
    return X_Ns.dot(np.transpose(X_Ns)) / n - p * np.diag(diag_X_Ns)



def Obtain_eta_Q(N, X, D_X, Sigma):
  [p, n] = X.shape 
  M = max(n, p, max(N))
  mu = D_X / n
  m = np.amax(X, axis = 1)
  if len(N) == 1:
      inv_N = float(1) / N[0]
      r1 = np.sqrt(np.log(M) / n * inv_N)
      r2 = r1 ** 2
      r3 = np.sqrt(np.log(M) ** 4 / n) * inv_N ** 1.5

      T1 = (np.sqrt(Sigma * m[:, None]) + np.sqrt(Sigma * m)) * r1 * 3 * np.sqrt(6)
      T2 = (np.repeat(m[:,None], p, 1) + np.repeat(m[None, :], p, 0)) * r2 * 2
      T3 = np.sqrt(np.repeat(mu[:,None], p, 1) + np.repeat(mu[None, :], p, 0)) * r3 * 31
      
      eta_mat = T1 + T2 + T3
      Sigma_norm = Sigma / mu[:, None] / mu
      Q = eta_mat / mu[:, None] / mu + (Sigma_norm / np.sqrt(mu[:,None])\
                      + Sigma_norm / np.sqrt(mu)) * r1 * 2
  else:
      X_Ni = X / N[None,:]
      inv_mean_N = np.mean(N ** (-1))
      Sigma_N = (X_Ni).dot(X.transpose()) / n
      mu_Ni3 = np.mean(X_Ni / N[None,:]**2, 1)
      mu_Ni = np.mean(X_Ni, 1)
      
      r1 = np.sqrt(np.log(M) / n)
      r2 = r1 ** 2 * inv_mean_N
      r3 = np.sqrt(np.log(M) ** 4 / n)
      
      T1 = (np.sqrt(Sigma_N * m[:, None]) + np.sqrt(Sigma_N * m)) * r1 * 3 * np.sqrt(6)
      T2 = (np.repeat(m[:,None], p, 1) + np.repeat(m[None, :], p, 0)) * r2 * 2
      T3 = np.sqrt(np.repeat(mu_Ni3[:,None],p,1) + np.repeat(mu_Ni3[None,:],p,0)) * r3 * 31
      
      eta_mat = T1 + T2 + T3
      Sigma_norm = Sigma / mu[:, None] / mu
      mu_norm = np.sqrt(mu_Ni) / mu
      Q = eta_mat / mu[:, None] / mu + (Sigma_norm * mu_norm[:,None]\
                      + Sigma_norm * mu_norm) * r1 * 2
  return [eta_mat, Q]
