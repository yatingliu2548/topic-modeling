#####   Cross validation  ######

import numpy as np
import Sp_Anchor
import Sp_Utils

def CV_C1(X, Ns, C1_grid):
  # Cross validation for choosing C1. For each C1 from the given grids, 
  # split the data into two parts. Estimate I, AI and C from the first
  # dataset and calculate AICAI^T. Choose C1 which minimizes the criterion 
  # AICAI' - Sigma^(2) from the second dataset.
  #
  # Args: 
  #   X: p by n frequency data matrix.
  #   Ns: n document lengths 
  #   C1_grid: vector of numerical constants.
  # 
  # Returns: 
  #   the selected optimal C1

  [p, n] = X.shape
  sampInd = np.random.choice(range(n), n / 2, replace = False)
  sampInd_comp = [x for x in range(n) if x not in sampInd]
  
  X1, X2 = X[:,sampInd], X[:,sampInd_comp]
  D_X1, D_X2 = np.sum(X1, axis = 1), np.sum(X2, axis = 1)
  n1, n2 = X1.shape[1], X2.shape[1]
  Sigma1 = Sp_Utils.GetSigma(X1, Ns[sampInd])
  Sigma2 = Sp_Utils.GetSigma(X2, Ns[sampInd_comp])
  R1 = n1 ** 2 * Sigma1 / D_X1[:, None] / D_X1
  R2 = n2 ** 2 * Sigma2 / D_X2[:, None] / D_X2
  
  M, arg_M = np.amax(R1, axis = 1), np.argmax(R1, axis = 1)
  Q = Sp_Utils.Obtain_eta_Q(Ns[sampInd], X1, D_X1, Sigma1)[1]
  
  loss = np.zeros(len(C1_grid))
  for i in range(len(C1_grid)):
    [fit, anchor_vec] = CalFittedSigma(R1, C1_grid[i], M, arg_M, Q)
    subR2 = R2[anchor_vec,:][:, anchor_vec]
    loss[i] = np.mean(abs(subR2 - fit))

  return np.median(C1_grid[loss == min(loss)])



def CalFittedSigma(Sigma, C1, M, arg_M, Q):
  # Calculate the fitted value of A_ICA_I^T for given Sigma and C1.
  # 
  # Args: 
  #   Sigma: p by p matrix. 
  #   C1: given parameter. 
  #   Ms: the calculated maximal values of Sigma per row.
  #
  # Returns: 
  #     anchor_vec: vector of the indices of estimated pure variables.
  #     fitted: fitted value of A_ICA_I'

  p = Sigma.shape[0]
  anchor_group = Sp_Anchor.SelectAnchor(Sigma, C1, M, arg_M, Q)
  anchor_vec = [anchor for group in anchor_group for anchor in group]
  K = len(anchor_group)
  A = np.zeros((p, K))

  for i in range(K):
    A[anchor_group[i], i] = 1   

  AI = A[anchor_vec, ]
  P = AI.dot(np.linalg.solve(AI.transpose().dot(AI), AI.transpose()))
  fitted = P.dot(Sigma[anchor_vec,][:,anchor_vec]).dot(P.transpose()) 

  return [fitted, anchor_vec]
