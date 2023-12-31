from __future__ import division

import numpy as np
import sys
import os
import errno
from numpy.random import RandomState
import random_projection as rp
import gram_schmidt_stable as gs

def findAnchors(Q, K, candidates, seed=1234,
                checkpoint_prefix=1, new_dim=16):
    # Random number generator for generating dimension reduction
    prng_W = RandomState(seed)

    # row normalize Q
    row_sums = Q.sum(1)
    for i in range(len(Q[:, 0])):
        Q[i, :] = Q[i, :]/float(row_sums[i])

    # Reduced dimension random projection method for recovering anchor words
    Q_red = rp.Random_Projection(Q.T, new_dim, prng_W)
    Q_red = Q_red.T
    (anchors, anchor_indices) = gs.Projection_Find(Q_red, K, candidates)

    # restore the original Q
    for i in range(len(Q[:, 0])):
        Q[i, :] = Q[i, :]*float(row_sums[i])

    return anchor_indices
