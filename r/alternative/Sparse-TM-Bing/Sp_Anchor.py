### 	Anchor word recovering 
import numpy as np

def SelectAnchor(Sigma, C1, rowMax, argRowMax, Q):
	group = []
	p = Sigma.shape[0]
	for i in range(p):
		lbd_rowi = C1 * (Q[i, argRowMax[i]] + Q[i, ])
		Si = np.where(rowMax[i] - Sigma[i, ] <= lbd_rowi)[0].tolist()
		if len(Si) == 0:
			next
		anchor_flag = TestAnchor(Sigma, i, Si, rowMax, argRowMax, C1, Q)
		if anchor_flag:
			group = Merge(group, Si)
	return group



def TestAnchor(Sigma, rowInd, Si, rowMax, argRowMax, C1, Q):
  for i in range(len(Si)): 
    j = Si[i]
    lbd = C1 * (Q[rowInd, j] + Q[j, argRowMax[j]])
    if abs(Sigma[rowInd, j] - rowMax[j]) > lbd:
      return False
  return True



def Merge(group, Si):
  # merge the new group with the previous ones which have common nodes
  newgroup = list(group)
  if (len(newgroup) != 0):
    for i in range(len(newgroup)):
      common_nodes = list(set(newgroup[i]).intersection(Si))
      if len(common_nodes) != 0:
      	newgroup[i] = common_nodes
      	return newgroup
  newgroup.append(Si)
  return newgroup






