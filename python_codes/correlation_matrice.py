#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:46:04 2025
@author: axel
"""

import matplotlib.pyplot as plt
import cooler
import numpy as np
import hicstuff
from scipy import sparse
import pandas as pd

%matplotlib

# loading of data 
# cool_file1="/media/axel/RSG5/disk/copy_diverse_yeast/bordelet_paper/SRR12284705_G1/contacts/matrices/rescued_map.cool"
cool_file1 = "/media/axel/EVO/yeast/XVIfMmycot1/out_hicberg_AT882/contacts/pairs/tmp/valid_idx_pcrfree.pairs.5000.cool"

c1 = cooler.Cooler(cool_file1)
cooler.balance_cooler(c1,store=True, mad_max=1000)   # Normalisation

# chr1="chrVI"
# chr2="chrVI"

chr1='chrXVI-Mmmyco-translocation1'
chr2='chrXVI-Mmmyco-translocation1'

m1 = c1.matrix(balance=True, sparse=False).fetch(chr1, chr2)

def correlation_mat(m) :
    """ Compute the correlation matrice with detrending
    """
    ms = sparse.csr_matrix(m)
    N=ms 
    N = N.tocoo()
    # Detrend by the distance law
    dist_bins, dist_vals = hicstuff.distance_law_from_mat(N, log_bins=False)
    N.data /= dist_vals[abs(N.row - N.col)]
    N = N.tocsr()  
    ND=N.todense()  
    
    N = hicstuff.corrcoef_sparse(N)
    
    # N=N.todense()
    # N[np.isnan(N)] = 0.0   
    
    # [eigen_vals, pr_comp] = hicstuff.eig(N.todense())
    return N , ND

# computation of correlation mat
m1[np.isnan(m1)] = 0
N,ND=correlation_mat(m1)

# plots:
plt.figure(1)
expo=0.15
plt.imshow(m1**expo, cmap="afmhot_r")
plt.colorbar()
plt.title(cool_file1)

plt.figure(2)
value_bar = 0.8
plt.imshow(ND,vmin=0, vmax=2.0)
plt.colorbar()
plt.title(cool_file1)

plt.figure(3)
value_bar = 0.8
plt.imshow(N.todense(), cmap="seismic",vmin=-value_bar, vmax=value_bar)
plt.colorbar()
plt.title(cool_file1)

