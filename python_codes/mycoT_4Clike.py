#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 09:40:56 2025
@author: axel
"""

import cooltools
# Import python package for working with cooler files and tools for analysis
import cooler

import re
import matplotlib.gridspec as gridspec
import serpentine as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy
import hicstuff
from scipy import sparse

# from hicberg

cool_file1 = "/media/axel/EVO/yeast/XVIfMmycot1/out_hicberg_AT881/contacts/pairs/tmp/valid_unrescued_idx_pcrfree.pairs.2000.cool"
cool_file2 = "/media/axel/EVO/yeast/XVIfMmycot1/out_hicberg_AT881/contacts/pairs/tmp/valid_idx_pcrfree.pairs.2000.cool"

# Computation

clr1 = cooler.Cooler(cool_file1)
clr2 = cooler.Cooler(cool_file2)

clr1.chromnames

clr1.info
clr2.info

cooler.balance_cooler(clr1, store=True, mad_max=1000) # Normalisation 
cooler.balance_cooler(clr2, store=True, mad_max=1000) # Normalisation 

chr_chosen1 = 'chrXVI-Mmmyco-translocation1'
chr_chosen2 = 'chrXVI-Mmmyco-translocation1'

m1 = clr1.matrix(balance=False).fetch(chr_chosen1, chr_chosen2)
m2 = clr1.matrix(balance=False).fetch(chr_chosen1, chr_chosen2)
 
#  matrices
m11 = clr1.matrix(balance=True).fetch(chr_chosen1, chr_chosen2)
m22 = clr2.matrix(balance=True).fetch(chr_chosen1, chr_chosen2)

expo=0.1    # for hicberg data 
# expo=0.1   # fro processed files 
 
plt.figure(1)
plt.imshow(m11**expo, cmap="afmhot_r")
plt.colorbar()
plt.title(cool_file1)
plt.ylabel(chr_chosen1)
plt.xlabel(chr_chosen2)


plt.figure(2)
plt.imshow(m22**expo, cmap="afmhot_r")
plt.colorbar()
plt.title(cool_file2)
plt.ylabel(chr_chosen1)
plt.xlabel(chr_chosen2)


# # 4C like plot from telomeric region 

plt.figure(3)
telo_signal = m22[1077,:]+m22[1076,:]+m22[1075,:]+m22[1074,:]+m22[1073,:]
plt.plot(telo_signal)

# telo_signal1 = m11[431,:]+m11[430,:]
# telo_signal2 = m22[431,:]+m22[430,:]


# plt.plot(telo_signal1)
# plt.plot(telo_signal2)

plt.xlabel("Genomic coordinates")
plt.ylabel("4C like signal from telomeric region")
plt.semilogy()

