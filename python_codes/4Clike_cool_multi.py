#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 10:34:57 2022
@author: axel
with a cool file, generates pdf files for each chromosome 
with the 4C like signal for several area of interest (multi, ex:TY1)
"""

import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import os
from scipy.signal import find_peaks
from scipy import interpolate
import numpy as np
import cooler

# part to fill 
cool_file="/media/axel/RSG51/diverse_yeast_data_copy/test_hicberg/AC1_pairs_files/all_group.pairs.header.200.cool"
name_bank="AC1_with_Ty1"

cool_file="/home/axel/Bureau/hicberg_project/test_20240105/cellular_cycle/sc_45m/matrices/rescued_map.cool"
name_bank="sc_45m"

cool_file="/home/axel/Bureau/hicberg_project/test_20240105/sacCer_no-array/matrices/rescued_map.cool"
name_bank="sacCer_no-array"

cool_file="/home/axel/Bureau/hicberg_project/test_20240105/LM171/matrices/rescued_map.cool"
name_bank="LM171"

cool_file="/media/axel/RSG51/diverse_yeast_data_copy/TY1_potential_integration/out_FG0129/tmp/valid_idx.pairs.2000.cool"
name_bank="FG0129"

cool_file="/media/axel/RSG51/diverse_yeast_data_copy/myco_2023/FG0095_myco_redone/contacts/matrices/rescued_map.cool"
name_bank="FG0095_rescued_redone"

cool_file="/media/axel/RSG51/diverse_yeast_data_copy/dChr_2024/dChr_F8/contacts/matrices/rescued_map.cool"
name_bank="dChr_F8_with_TY1"

cool_file="/home/axel/Bureau/hicberg_project/test_20240105/sacCer_no-array/matrices/rescued_map.cool"
name_bank="sacCer_no-array_with-TY1"

cool_file="/media/axel/RSG51/diverse_yeast_data_copy/SynEc/LM54_SynEc/contacts/matrices/unrescued_map.cool"
name_bank="LM54_SynEc_with_TY1_unrescued"

cool_file= "/media/axel/RSG51/diverse_yeast_data_copy/heat_shock_2024/FG128/matrices/rescued_map.cool"
name_bank="FG128_after_HS_with_rDNA"


cool_file= "/media/axel/RSG51/diverse_yeast_data_copy/3Cseq/SRR1649436_hicberg02/contacts/matrices/rescued_map.cool"
name_bank="SRR1649436_with_TY1"


cool_file= "/media/axel/RSG51/diverse_yeast_data_copy/AC1/AC1_hicberg02/contacts/matrices/rescued_map.cool"
name_bank="AC1_hicberg02_with_TY1"

cool_file= "/home/axel/Bureau/hicberg_project/AC1_update_ps/rescued_map.cool"
name_bank="AC1_update_ps"

cool_file= "/media/axel/RSG51/diverse_yeast_data_copy/TY1_shot_gun_signal/out_FG0129_mq1/tmp/valid_idx.pairs.2000.cool"
name_bank="FG0129_mq1_shot_gun"

cool_file= "/home/axel/Bureau/hicberg_project/test_20240105/AC1_random/matrices/rescued_map.cool"
name_bank="AC1_random"

# input group of genomic positions to sum from wich we want the 4Clike
ty1 = pd.read_table("/home/axel/Bureau/YEAST/TY/Ty1.txt",header=None)

# here positions of TY1 or all TY  
ty1 = pd.read_csv('/home/axel/Bureau/YEAST/TY/positions_TY.txt2', 
                    sep="\t", header= None)

ty1 = pd.read_csv('/home/axel/Bureau/YEAST/TY/seq_telomeres.txt2', 
                    sep="\t", header= None)

ty1 = pd.read_csv('/media/axel/RSG51/diverse_yeast_data_copy/dChr_2024/pHM240.txt', 
                    sep="\t", header= None)

ty1 = pd.read_csv('/media/axel/RSG51/diverse_yeast_data_copy/heat_shock_2024/heat_shock_plots_rDNA/rDNA_coordinates.txt', 
                    sep="\t", header= None)

ty1 = pd.read_csv('/media/axel/RSG51/plasmodium/Pf3D7_MIT_v3.txt', 
                    sep="\t", header= None)
#  
c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=1000)   # Normalisation 
d=c.info
total_reads = d['sum']
# SC288 assembly:
BIN=2000
list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16","chrM"]
list_chr = c.chromnames
list_chr1 = c.chromnames

total_number_reads=c.info['sum']
reads_plasmide=0
for chr3 in list_chr1:
    ty1_chr = ty1[ty1[0] ==chr3]
    for ind in ty1_chr.index:
        chr_chosen=ty1_chr[0][ind]
        pos1=ty1_chr[1][ind]
        pos2=ty1_chr[2][ind]
        chr2 = (chr_chosen,pos1,pos2) # region of focus 
        
        m1 = c.matrix(balance=False).fetch(chr2, chr3)
        if chr_chosen==chr3:
            reads_plasmide += m1.sum()/2+np.trace(m1)/2
        else:
            reads_plasmide += m1.sum()

perc_plasmid=reads_plasmide/total_number_reads
print("Proportion of focus area")
print(perc_plasmid)

#
limit_y =   2.1
limit_x =  int(1532000/BIN)
expo=0.12

def nan_sum(row):
    total = np.nanmean(row)*len(row[~np.isnan(row)])
    return total


if not os.path.exists(name_bank+"_files") :
    os.makedirs(name_bank+"_files")

list_all_contact = []
list_all_chip = []

list_all_contact_lg = []
list_all_chip_lg = []

list_all_chr = []
list_all_pos = []

df_total = pd.DataFrame()


for chr1 in list_chr1 :
    ty1_chr = ty1[ty1[0] ==chr1]
    # 1) previous method
    matscn = c.matrix(balance=True).fetch(chr1, chr1)
    coverage1 = np.apply_along_axis(nan_sum, 1, matscn) # Hi-C coverage intra only
    
    k=0
    for ind in ty1.index:
        print(ind)
        chr_chosen=ty1[0][ind]
        pos1=ty1[1][ind]
        pos2=ty1[2][ind]
        chr22 = (chr_chosen,pos1,pos2) # region of focus
        mat = c.matrix(balance=False).fetch(chr1, chr22)  # already normalised 4C signal
        if k==0:
            coverage_plasmid1 = np.apply_along_axis(nan_sum, 1, mat)
        else: 
            coverage_plasmid1 = coverage_plasmid1 + np.apply_along_axis(nan_sum, 1, mat)
        k+=1
        
    # to work with the raw number of reads
    reads_chr1=0
    coverage=0
    for chr3 in list_chr1:
        m1 = c.matrix(balance=False).fetch(chr1, chr3)
        reads_chr1 += m1.sum()
        coverage_temp=np.apply_along_axis(nan_sum, 1, m1)
        coverage_temp[np.isnan(coverage1)] = np.nan
        coverage += coverage_temp              # Hi-C general coverage

    # 2) for each member of the group of interest
    j=0
    for ind in ty1.index:
        print(ind)
        chr_chosen=ty1[0][ind]
        pos1=ty1[1][ind]
        pos2=ty1[2][ind]
        chr22 = (chr_chosen,pos1,pos2) # region of focus  
        
        m2 = c.matrix(balance=True).fetch(chr22, chr22)
        reads_chr2= m2.sum()
        contact_of_plasmid_i= np.apply_along_axis(nan_sum, 1, m2)
#        if j==0:
#            contact_of_plasmid = contact_of_plasmid_i
#        else :
#            contact_of_plasmid += contact_of_plasmid_i
        
        m12 = c.matrix(balance=False).fetch(chr1, chr22)
        reads_chr12= m12.sum()
        coverage_plasmid_i = np.apply_along_axis(nan_sum, 1, m12)
        if j==0:
            coverage_plasmid=coverage_plasmid_i 
        else: 
            coverage_plasmid=coverage_plasmid+coverage_plasmid_i  
        j+=1  
        
    coverage_plasmid[np.isnan(coverage_plasmid1)] = np.nan

    coverage_plasmid=coverage_plasmid/coverage # normalisation 
    coverage_plasmid=coverage_plasmid/perc_plasmid # taking into accournt plasmid proportion

    # computation of correlation:
    maxi= matscn.shape[0]*BIN
    def fill_nan(A):
        '''
        interpolate to fill nan values
        '''
        inds = np.arange(A.shape[0])
        good = np.where(np.isfinite(A))
        f = interpolate.interp1d(inds[good], A[good],bounds_error=False)
        B = np.where(np.isfinite(A),A,f(inds))
        return B

    if BIN == 2000:
        proportion_nan = np.count_nonzero(np.isnan(coverage_plasmid))/ len(coverage_plasmid)
    if proportion_nan  < 0.9:
        coverage_plasmid_new = fill_nan(coverage_plasmid)
    else :
        coverage_plasmid_new = np.zeros(len(coverage_plasmid))
    hot_spots = find_peaks(coverage_plasmid_new, 
                           height = 0.0006,
                           width = 2.,
                           distance=2)  # for 2000 bp resolution    

    if BIN == 200:
        coverage_plasmid[np.isnan(coverage_plasmid)] = 0.0
        coverage_plasmid_new = coverage_plasmid
        hot_spots = find_peaks(coverage_plasmid_new,
                               height = 0.0011,
                               width = 5.,
                               distance=70)  # for 200 bp resolution

    col_chr = [chr1] * len(hot_spots[0])
    col_pos = list(hot_spots[0])
    col_pos = [int(e*BIN) for e in col_pos]
    df = pd.DataFrame({'chrom1': col_chr, 'start1': col_pos})
    df_total = pd.concat([df_total, df])

    # Multiplot:
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8, 12.8)
    fig.tight_layout(pad=5.0)
    gs = gridspec.GridSpec(4, 1,height_ratios=[7,1,1,1])

    ax1 = plt.subplot(gs[0])
    ax1.imshow(matscn**expo, interpolation="none", cmap="afmhot_r",
           vmin = 0.0)
    locs, labels = plt.xticks()            # Get locations and labels
    locs2 = [int(i*BIN) for i in locs]
    plt.xticks(locs,locs2)
    fig.tight_layout(pad=2.)
    plt.title("\n\n\n\n"+chr1+" "+str(total_reads)+" "
              +str(reads_chr1)+" "+str(reads_chr2)+" "+str(reads_chr12)+"  "
              +name_bank)

    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.plot(coverage)
    plt.title("Proportion of the signal in INTRA (in %)")

#    coverage_plasmid = coverage_plasmid /(float(sizes_dist[chr2]))
    median_plasmid = np.nanmedian(coverage_plasmid)
    print(median_plasmid)
    ax3 = plt.subplot(gs[2], sharex=ax1)
    ax3.plot(coverage_plasmid, color="royalblue" )
    plt.yscale('log')
#    plt.ylim(10**-1,30)
    
#    plt.xticks([])
    # plot with rectangles of TY1: 
    pos1=np.array(ty1_chr[1])
    pos2=np.array(ty1_chr[2])
    for i in range(len(pos1)):
        rectangle = plt.Rectangle( (min(int(pos1[i]),int(pos2[i]))/BIN, 0.0),
                                  abs(int(pos2[i])-int(pos1[i]))/BIN, 
                                  limit_y/4, 
                                  fc="orange",
                                  label="TY1")
        plt.gca().add_patch(rectangle)

    plt.axhline(y=median_plasmid,ls='--',color="Black")
#    plt.xticks([])
#    plt.ylim(0, limit_y)
#    plt.xlim(0, limit_x)
#    plt.plot(hot_spots[0],hot_spots[0]/hot_spots[0]*limit_y,'v', color="orange")
    h=hot_spots[1]
    left_ips = h['left_ips']
    right_ips = h['right_ips']
#    plt.plot(left_ips,left_ips*0.0,'|')
#    plt.plot(right_ips,right_ips*0.0,'|')
    plt.title(chr2[0]+" Median="+str(round(median_plasmid*10**4,2))+"x 10^-4")

    ax4 = plt.subplot(gs[3], sharex=ax1)
    ax4.plot(coverage_plasmid_new, color="grey")
    plt.ylim(10**-1,30)
#    plt.axhline(y=median_chrM,ls='--',color="Black")
#    plt.ylim(0, limit_y)
#    ax4.plot(b_ip, values_chip, color="royalblue")
    plt.xlabel("Position along the chromosome (bins "+ str(BIN) +"bp)")
#    plt.title(str(chrM)+" Median="+str(round(median_chrM*10**8,2))+"x 10^8")
    plt.title("Plasmid contact with 1D interpolation")

    list_all_contact=np.concatenate((list_all_contact, coverage_plasmid), axis=0)
    list_all_chr=np.concatenate((list_all_chr, [chr1] * len(coverage_plasmid)), axis=0)
    list_all_pos=np.concatenate((list_all_pos, range(len(coverage_plasmid))), axis=0)
    
    # only for long genes:
    plt.savefig(name_bank+"_files"+"/MAT_SCN_"+chr1+"_"+
            name_bank+"_"+str(BIN/1000)+"kb"+".pdf")

 
#
#
## writting of all hot spots of contact in one file:
#df_total.to_csv("area_in_"+name_bank+".txt", sep='\t',
#                index=False)
#
#df_all = pd.DataFrame(
#    {'chr': list_all_chr,
#     'bin': list_all_pos,
#     'contact': list_all_contact
#    })
#
#df_all.to_csv("contact_signal2_"+chr2[0]+"_"+name_bank+".txt", sep='\t',
#                index=False,na_rep='NA')

