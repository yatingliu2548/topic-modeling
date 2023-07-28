#################################################################################
#####                                                                       #####
#####       Code for Step 2: regression to estimate the non-anchor rows     #####
#####                                                                       #####
#################################################################################
import cvxpy as cvx
import numpy as np

def Est_A(D_X, R, anchor_group, anchor_vec, C0, lbd, thresh):
  p = len(D_X)
  K = len(anchor_group)
  non_anchor_vec = [x for x in range(p) if x not in anchor_vec]

  B = np.zeros((p, K))

  for k in range(K):
    anchor_k = anchor_group[k]
    B[anchor_k, k] = 1

  BI = B[anchor_vec, ]
  PI = np.linalg.solve(BI.transpose().dot(BI), BI.transpose())
  C_tilde = PI.dot(R[anchor_vec,][:,anchor_vec]).dot(PI.transpose())
  lbd_step = C0 * lbd

  counter = 0
  pd_flag = False
  while not pd_flag:
    try:
      np.linalg.cholesky(C_tilde)
      pd_flag = True
    except np.linalg.LinAlgError:
      counter += 1
      C_tilde += lbd_step * np.identity(K)

  Theta_tilde = PI.dot(R[anchor_vec,])

  for j in non_anchor_vec:
    if D_X[j] >= thresh:
      B[j, ] = Est_AJ_row(C_tilde, Theta_tilde[:,j])

  A_tilde = B * D_X[:, None] / float(p)
  return A_tilde / np.sum(A_tilde, 0)


def Est_AJ_row(C_tilde, theta_tilde):
  K = C_tilde.shape[0]
  one_vec = np.ones(K)

  # Construct convex problem
  beta = cvx.Variable(K)
  objective = cvx.Minimize(cvx.quad_form(beta, C_tilde) - 2 * beta.T @ theta_tilde)
  constraints = [beta >= 0, one_vec.T @ beta == 1]
  prob = cvx.Problem(objective, constraints)
  prob.solve(solver = cvx.OSQP)

  return beta.value
