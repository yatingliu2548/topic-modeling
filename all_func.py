import numpy as np
from numpy.linalg import lstsq, norm
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from scipy.spatial.distance import euclidean
from sklearn.preprocessing import normalize
import random
import logging
logging.disable(logging.INFO)

import time

import pylab
pylab.rcParams['figure.figsize'] = (10.0, 8.0)

## Data simul
def gen_data(V=1000, K=10, M=1000, Nm=100, eta=0.1, alpha=0.1, M_test=500, anchor=True, equi=False,sparseW=False,sparseA=False):#number of nonzero S_w=int(1/2)
    #beta is A
    s=3
    beta_t = np.random.dirichlet(np.ones(V)*eta, K)
    if equi:
        if K!=V:
            print('DONT DO IT')
        beta_t = np.eye(K)
    theta_t = np.random.dirichlet(np.ones(K)*alpha, M)
    theta_test = np.random.dirichlet(np.ones(K)*alpha, M_test)
    anchor_set = []
    #if anchor:
    #    for i in range(K):
    #        cand = np.argmax(beta_t[i,:])
    #        anchor_set.append(cand)
    #        beta_t[np.arange(K)!=i,cand] = 0
    #    beta_t = np.apply_along_axis(lambda x: x/x.sum(), 1, beta_t)
    cardset=random.sample(range(1000), k=int(K*(K-1)/2))
    if anchor:
        r=0
        for i in range(K):
            for j in [x for x in range(i,K) if x != i]:
           
                if i != j:
                    x_ij = np.random.uniform(low=0, high=1/K)  # Generate random value 0 <= x_ij < 1/K
            column = np.random.uniform(low=0, high=1/K) * np.eye(K)[:, j] + x_ij * np.eye(K)[:, i]
            beta_t[:, cardset[r]] = column
            r=r+1
        beta_t = np.apply_along_axis(lambda x: x/x.sum(), 1, beta_t)
    if sparseA:
        
        for ii in range(V):
            beta_copy=beta_t[:,ii].copy()
            beta_copy.sort()
            card_list = [np.argwhere(beta_t[:,ii] == temp) for temp in beta_copy[-s:]]
            beta_t[~np.isin(np.arange(K),card_list),ii] = 0
        beta_t = np.apply_along_axis(lambda x: x/x.sum(), 1, beta_t)
    if sparseW:
        for j in range(M):
            theta_copy=theta_t[j,:].copy()
            theta_copy.sort()
            card_listw = [np.argwhere(theta_t[j,:] == temp) for temp in theta_copy[-s:]]
            theta_t[j,~np.isin(np.arange(K),card_listw)] = 0
        theta_t = np.apply_along_axis(lambda x: x/x.sum(), 1, theta_t)
    simplex = np.dot(theta_t, beta_t)
    wdf = np.apply_along_axis(lambda x: np.random.multinomial(Nm, x), 1, simplex).astype('float')
#    sum(wdf.sum(axis=0)==0)
    full_set = wdf.sum(axis=0)>0
    new_V = [i for i in range(V) if full_set[i]]
    anchor_set = [new_V.index(i) for i in anchor_set]
    wdf = wdf[:,full_set]
    beta_t = normalize(beta_t[:,full_set], 'l1')
    simplex_test = np.dot(theta_test, beta_t)
    wdf_test = np.apply_along_axis(lambda x: np.random.multinomial(Nm, x), 1, simplex_test).astype('float')
    return wdf, beta_t, theta_t, simplex[:,full_set], anchor_set, theta_test, wdf_test
##Generate sample method 2
def gen_data_sparse(V=1000, K=10, M=1000, Nm=100, s=20, eta_1=0.3,eta_2=3, alpha=0.1, M_test=500, anchor=True):#number of nonzero S_w=int(1/2)
    #beta is A

    beta_t = np.zeros((V, K))
    
    beta_t[0:s,:]=np.abs(np.random.normal(0,100,size=(s,K)))#np.random.dirichlet(np.ones(s)*eta_1, K).T
    beta_t[s:,:]=np.random.dirichlet(np.ones(V-s)*eta_2, K).T*0.00001
    
    beta_t=beta_t.T
    beta_t = np.apply_along_axis(lambda x: x/x.sum(), 1, beta_t)
    
    theta_t = np.random.dirichlet(np.ones(K)*alpha, M)
    theta_test = np.random.dirichlet(np.ones(K)*alpha, M_test)
    anchor_set = []
    #if anchor:
    #    for i in range(K):
    #        cand = np.argmax(beta_t[i,:])
    #        anchor_set.append(cand)
    #        beta_t[np.arange(K)!=i,cand] = 0
    #    beta_t = np.apply_along_axis(lambda x: x/x.sum(), 1, beta_t)
    cardset=random.sample(range(1000), k=int(K*(K-1)/2))
    if anchor:
        r=0
        for i in range(K):
            for j in [x for x in range(i,K) if x != i]:
           
                if i != j:
                    x_ij = np.random.uniform(low=0, high=1/K)  # Generate random value 0 <= x_ij < 1/K
            column = np.random.uniform(low=0, high=1/K) * np.eye(K)[:, j] + x_ij * np.eye(K)[:, i]
            beta_t[:, cardset[r]] = column
            r=r+1
        beta_t = np.apply_along_axis(lambda x: x/x.sum(), 1, beta_t)
    simplex = np.dot(theta_t, beta_t)
    
    wdf = np.apply_along_axis(lambda x: np.random.multinomial(Nm, x), 1, simplex).astype('float')
#    sum(wdf.sum(axis=0)==0)
    full_set = wdf.sum(axis=0)>=0
    new_V = [i for i in range(V) if full_set[i]]
    anchor_set = [new_V.index(i) for i in anchor_set]
    wdf = wdf[:,full_set]
    beta_t = normalize(beta_t[:,full_set], 'l1')
    simplex_test = np.dot(theta_test, beta_t)
    wdf_test = np.apply_along_axis(lambda x: np.random.multinomial(Nm, x), 1, simplex_test).astype('float')
    return wdf, beta_t, theta_t, simplex[:,full_set], anchor_set, theta_test, wdf_test

## GDM algorithms
def get_beta(cent, centers, m):
    betas = np.array([cent + m[x]*(centers[x,:] - cent) for x in range(centers.shape[0])])
    betas[betas<0] = 0
    betas = normalize(betas, 'l1')
    return betas

def gdm(wdfn, K, ncores=-1):
    glob_cent = np.mean(wdfn, axis=0)
    kmeans = KMeans(n_clusters=K, n_jobs=ncores, max_iter=1000).fit(wdfn)
    centers = kmeans.cluster_centers_
    labels = kmeans.labels_
    m = []
    for k in range(K):
        k_dist = euclidean(glob_cent, centers[k])
        Rk = max(np.apply_along_axis(lambda x: euclidean(glob_cent, x), 1, wdfn[labels==k,:]))
        m.append(Rk/k_dist)
    
    beta_means = get_beta(glob_cent, centers, m)
    
    return beta_means

## Geometric Theta
def proj_on_s(beta, doc, K, ind_remain=[], first=True, distance=False):
    if first:
        ind_remain = np.arange(K)
    s_0 = beta[0,:]
    if beta.shape[0]==1:
        if distance:
            return norm(doc-s_0)
        else:
            theta = np.zeros(K)
            theta[ind_remain] = 1.
            return theta
    beta_0 = beta[1:,:]
    alpha = lstsq((beta_0-s_0).T, doc-s_0)[0]
    if np.all(alpha>=0) and alpha.sum()<=1:
        if distance:
            p_prime = (alpha*(beta_0-s_0).T).sum(axis=1)
            return norm(doc-s_0-p_prime)
        else:
            theta = np.zeros(K)
            theta[ind_remain] = np.append(1-alpha.sum(), alpha)
            return theta
    elif np.any(alpha<0):
        ind_remain = np.append(ind_remain[0], ind_remain[1:][alpha>0])
        return proj_on_s(np.vstack([s_0, beta_0[alpha>0,:]]), doc, K, ind_remain, False, distance)
    else:
        return proj_on_s(beta_0, doc, K, ind_remain[1:], False, distance)

## Evaluation
def min_match(beta, beta_t):
    b_to_t = np.apply_along_axis(lambda x: np.sqrt(((beta_t-x)**2).sum(axis=1)), 1, beta)
    return max([max(np.min(b_to_t, axis=0)), max(np.min(b_to_t, axis=1))])

def perplexity(docs, beta, theta='geom', scale=True):
  if type(theta)==str:
      theta = np.apply_along_axis(lambda x: proj_on_s(beta, x, beta.shape[0]), 1, normalize(docs, 'l1'))  
      scale = True
  est = np.dot(theta, beta)
  if scale:
      est = np.log(normalize(np.apply_along_axis(lambda x: x + x[x>0].min(), 1, est), 'l1'))
      mtx = docs * est
  else:
      est = np.log(est)
      mtx = docs * est
      mtx[np.isnan(mtx)] = 0.
  return np.exp(-mtx.sum()/docs.sum())

def get_stat(t_s, beta_est, beta_t, data_test, theta_test, hdp_gibbs = False):
    mm = min_match(beta_est, beta_t)
    if type(hdp_gibbs)!=bool:
        pp = perplexity(data_test, hdp_gibbs, theta=theta_test, scale=False)
        if pp == np.inf:
            pp = perplexity(data_test, hdp_gibbs, theta=theta_test, scale=True)
    else:
        pp = perplexity(data_test, beta_est, theta=theta_test, scale=False)
        if pp == np.inf:
            pp = perplexity(data_test, beta_est, theta=theta_test, scale=True)
    t = time.time() - t_s
    return [mm, pp, t]
