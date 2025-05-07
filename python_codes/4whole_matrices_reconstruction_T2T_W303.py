#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 14:49:46 2022
@author: axel
"""
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import cooler
import numpy as np
import pyBigWig
import pandas as pd 
import random 
import re   # regular expressions
from matplotlib.ticker import FuncFormatter

# from 2 bw files Chip-seq: 
bw = pyBigWig.open("/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR15041199/contacts/pairs/unrescued_Chip_Scc1_over_input.bw")
bw2 = pyBigWig.open("/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR15041199/contacts/pairs/rescued_Chip_Scc1_over_input.bw")

bw = pyBigWig.open("")
bw2 = pyBigWig.open("")

bw2.chroms()

# on one chr 
chr1='chrIV'
chr2="chrIV"

i=np.array(bw.intervals(chr1) )
i2=np.array(bw2.intervals(chr1) )

x1=i[:,0]
x2=i2[:,0]
v1=i[:,2]
v2=i2[:,2]

x1=x1[~np.isnan(v1)]
x2=x2[~np.isnan(v2)]

v1=v1[~np.isnan(v1)]
v2=v2[~np.isnan(v2)]

# #-------------------------------------------------
# cool_file1="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR12284705/contacts/matrices/unrescued_map.cool"
# cool_file2="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR13736508/contacts/matrices/unrescued_map.cool"
# cool_file3="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR12284705/contacts/matrices/rescued_map.cool"
# cool_file4="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR13736508/contacts/matrices/rescued_map.cool"

# with filtering of uncuts: 
cool_file1="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR12284705/contacts/pairs/tmp/valid_unrescued_idx_pcrfree.pairs.2000.cool"
cool_file2="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR13736508/contacts/pairs/tmp/valid_unrescued_idx_pcrfree.pairs.2000.cool"
cool_file3="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR12284705/contacts/pairs/tmp/valid_rescued_idx_pcrfree.pairs.2000.cool"
cool_file4="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR13736508/contacts/pairs/tmp/valid_rescued_idx_pcrfree.pairs.2000.cool"


# # # with lib of Freds: 
# cool_file1="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR12284705/contacts/pairs/tmp/valid_unrescued_idx_pcrfree.pairs.2000.cool"
# cool_file2="/media/axel/EVO/yeast/T2T_W303/out_hicberg_NB9/contacts/matrices/unrescued_map.cool"
# cool_file3="/media/axel/EVO/yeast/T2T_W303/out_hicberg_SRR12284705/contacts/pairs/tmp/valid_rescued_idx_pcrfree.pairs.2000.cool"
# cool_file4="/media/axel/EVO/yeast/T2T_W303/out_hicberg_NB9/contacts/matrices/rescued_map.cool"

name_bank1="SRR12284705_G1"
name_bank2="SRR13736508_cdc20"
name_bank3="SRR12284705_G1"
name_bank4="SRR13736508_cdc20"

# name_bank1="G1 (alpha-factor)\n before reconstruction"
# name_bank2="G2/M (nocodazole) \n before reconstruction"
# name_bank3="G1 (alpha-factor)\n after reconstruction"
# name_bank4="G2/M (nocodazole) \n after reconstruction"


c1 = cooler.Cooler(cool_file1)
c2 = cooler.Cooler(cool_file2)
c3 = cooler.Cooler(cool_file3)
c4 = cooler.Cooler(cool_file4)
#
# cooler.balance_cooler(c1,store=True, mad_max=1000)   # Normalisation
# cooler.balance_cooler(c2,store=True, mad_max=1000)   # Normalisation
# cooler.balance_cooler(c3,store=True, mad_max=1000)   # Normalisation
# cooler.balance_cooler(c4,store=True, mad_max=1000)   # Normalisation

# c is a cool file 
#m1 = c1.matrix(balance=False, sparse=False)
#m2 = c2.matrix(balance=False, sparse=False)
#
#m1=m1[:]   #  whole matrice
#m2=m2[:]  


m1 = c1.matrix(balance=True).fetch(chr1, chr2)
m2 = c2.matrix(balance=True).fetch(chr1, chr2)

m3 = c3.matrix(balance=True).fetch(chr1, chr2)
m4 = c4.matrix(balance=True).fetch(chr1, chr2)

#  -------------------------------------------------  plot  


#
expo=0.12
va=0
vb=np.nanmax(m4**expo)

current_cmap = mpl.cm.get_cmap('afmhot_r')
current_cmap.set_bad(color='white')

facteur = 2
def formatter_x(value, _):
# 'value' est la valeur réelle sur l’axe.
# On renvoie la valeur multipliée par 'facteur', formatée à la décimale près.
    return f"{value * facteur:.1f}"


# Plot with coverages 
gs = gridspec.GridSpec(3,2,height_ratios=[10,10,1],width_ratios=[1,1])
BIN = 2000
bin_matrice= 2000


ax1 = plt.subplot(gs[0])
ax1.imshow(m1**expo,interpolation='none',vmin= va, vmax= vb,cmap=current_cmap)
ax1.set_title(name_bank1+" "+chr1)
ax1.set_yticks([])
ax1.set_xticks([])

ax1 = plt.subplot(gs[1], sharex=ax1,sharey=ax1)
ax1.imshow(m3**expo,interpolation='none',vmin= va, vmax= vb,cmap=current_cmap)
ax1.set_title(name_bank3+" "+chr1)
# ax1.set_yticks([])
# ax1.set_xticks([])

ax2 = plt.subplot(gs[2], sharex=ax1,sharey=ax1)
ax2.imshow(m2**expo,interpolation='none',vmin= va, vmax= vb,cmap=current_cmap)
ax2.set_title(name_bank2+" "+chr1)
ax2.set_yticks([])
ax2.set_xticks([])

ax2 = plt.subplot(gs[3], sharex=ax1,sharey=ax1)
ax2.imshow(m4**expo,interpolation='none',vmin= va, vmax= vb,cmap=current_cmap)
ax2.set_title(name_bank4+" "+chr1)
ax2.set_yticks([])
ax2.set_xticks([])

# ChIP-seq data
ax12 = plt.subplot(gs[4], sharex=ax1)
#ax12.plot(i[:,0]/bin_matrice, i[:,2], color="royalblue",label="unrescued")
ax12.fill_between(x1/bin_matrice,v1, facecolor='royalblue', alpha=1., label="unrecued")
ax12.set_ylim([-5, 10.0])
ax12.set_title("Scc1 occupancy with reconstruction")
ax12.set_xlabel("Genomic coordinates")
ax12.set_ylabel("ChIP/input")

# ax12.xaxis.set_major_formatter(FuncFormatter(formatter_x))

ax22 = plt.subplot(gs[5], sharex=ax12,sharey=ax12)
#ax22.plot(i2[:,0]/bin_matrice, i2[:,2], color="limegreen",label="rescued")
ax22.fill_between(x2/bin_matrice,v2, facecolor='royalblue', alpha=1., label="recued")
ax22.set_ylim([-7, 10.0])
ax22.set_title("Scc1 occupancy with reconstruction")
ax22.set_xlabel("Genomic coordinates")
ax22.set_ylabel("ChIP/input")



# plot of genes and other elements:
genes = pd.read_csv('/media/axel/EVO/yeast/T2T_W303/all_objects_genes_W303.asm01.HP0.tidy.gff3.txt', 
                    sep="\t", header= None)

genes_chr = genes[genes[0]==chr1]
genes_chr=np.array(genes_chr)
for j in range(len(genes_chr)) :
#    if genes_chr[j][1]=="DYN1":
    if genes_chr[j][3]-genes_chr[j][2]>100 or genes_chr[j][3]-genes_chr[j][2]<-100:
        y_rand=random.uniform(0, 1)
        p1=int(genes_chr[j][2])
        delta_p=float( (genes_chr[j][3]-genes_chr[j][2]))
        
        if delta_p>0:
            col_arrow="royalblue"
        else :
            col_arrow="red"
        if bool(re.search(r'TY', genes_chr[j][1])):
#        if re.findall(genes_chr[j][1], r'TY\d') :      
            col_arrow="orange"
            print(genes_chr[j][1])
            
        ax12.arrow(p1/bin_matrice, -5-y_rand, delta_p/bin_matrice, 0,
                  width=1,length_includes_head=True,head_length=0.4,
                  color=col_arrow)     
        ax12.arrow(p1/bin_matrice, -5-y_rand, delta_p/bin_matrice, 0,
                  width=1,length_includes_head=True,head_length=0.4,
                  color=col_arrow) 
        ax22.arrow(p1/bin_matrice, -5-y_rand, delta_p/bin_matrice, 0,
                  width=1,length_includes_head=True,head_length=0.4,
                  color=col_arrow)     
        ax22.arrow(p1/bin_matrice, -5-y_rand, delta_p/bin_matrice, 0,
                  width=1,length_includes_head=True,head_length=0.4,
                  color=col_arrow) 
        # ax22.text(p1/bin_matrice+(delta_p/bin_matrice)/2.0, y_rand, genes_chr[j][1], rotation=45)
    

plt.tight_layout()

## plot sum of columns
#signal_sum=[sum(x) for x in zip(*m1)]
#plt.plot(signal_sum)
#
#
## ratio plot
#np.max(m2-m1)
#m=m2-m1
#v = np.reshape(m, m.shape[0]*m.shape[0])
#np.median(v)
#np.max(v)
#
#plt.hist(v, np.max(v))
#
#plt.imshow(m2-m1,interpolation='none',vmin= 0, vmax= 100,cmap="Blues")
#plt.colorbar()
#
#
