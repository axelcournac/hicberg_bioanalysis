#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 14:49:46 2022
@author: axel
mappability: adding of the signal of the mappability
"""
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import cooler
import numpy as np
import pandas as pd

# cool_file1= "/home/axel/Bureau/Apollo_project/test_yeast_20230217/S_cerevisiae_pattern_exp_bench_map.cool"
# cool_file2= "/home/axel/Bureau/Apollo_project/YEAST_20230217/S_cerevisiae_pattern_exp_bench_map_rescued.cool"

# cool_file1= "/home/axel/Bureau/Apollo_project/test_yeast_20230203/CH_128_filtered_full_map.cool"
# cool_file2= "/home/axel/Bureau/Apollo_project/test_yeast_20230203/CH_128_filtered_full_map_rescued.cool"

# cool_file1="/home/axel/Bureau/AC1/matrices/unrescued_map.cool"
# cool_file2="/home/axel/Bureau/AC1/matrices/rescued_map.cool"

# cool_file1="/home/axel/Bureau/hicberg_bacteries/OS_46_quad_map.cool"
# cool_file2="/home/axel/Bureau/hicberg_bacteries/OS_46_quad_map_rescued.cool"

# cool_file1="/home/axel/Bureau/hicberg_bacteries/unrescued_map.cool"
# cool_file2="/home/axel/Bureau/hicberg_bacteries/rescued_map.cool"

# cool_file1="/home/axel/Bureau/matrices_epichloe_festucase/unrescued_map.cool"
# cool_file2="/home/axel/Bureau/matrices_epichloe_festucase/rescued_map.cool"

cool_file1="/media/axel/RSG5/disk/copy_bureau/hicberg_project/test_20240105/AC1/matrices/unrescued_map.cool"
cool_file2="/media/axel/RSG5/disk/copy_bureau/hicberg_project/test_20240105/AC1/matrices/rescued_map.cool"

c1 = cooler.Cooler(cool_file1)
c2 = cooler.Cooler(cool_file2)

cooler.balance_cooler(c1, ignore_diags=False,mad_max=1000, store=True)   # Normalisation
cooler.balance_cooler(c2, ignore_diags=False,mad_max=1000, store=True)

# c is a cool file 
# m1 = c1.matrix(balance=True, sparse=False)
# m2 = c2.matrix(balance=True, sparse=False)

# m12=m1[:]   #  whole cool 
# m2=m2[:]   

c1.info
c2.info
BIN=2000
expo=0.12
# expo=0.5   #  for the inter maps 

# mappability
#bedGraph = pyBedGraph.BedGraph('/home/axel/Bureau/Apollo_project/sacCer3_with_plasmid_2micron_output/SC288_with_micron.genmap.chrom.sizes', '/home/axel/Bureau/Apollo_project/sacCer3_with_plasmid_2micron_output/SC288_with_micron.genmap.bedgraph')
#bed_chr=bedGraph.load_chrom_data("chr1")

#file_bed="/home/axel/Bureau/Apollo_project/sacCer3_with_plasmid_2micron_output_W303_out/W303.genmap.bedgraph"
file_bed="/media/axel/RSG5/disk/copy_bureau/hicberg_project/sacCer3_with_plasmid_2micron_output/SC288_with_micron.genmap.bedgraph"
df=pd.read_table(file_bed,header=None, delimiter="\t")


list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16","chrM"]


list_chr =  ["chr1"]

# list_chr =  c1.chromnames
factor_bin = 2  #  number of kb for the bin 

# chr_chosen2="chr5"
i_fig=0
for chr1 in list_chr :
    chr_chosen=chr1
    i_fig=i_fig+1
    df2=df.loc[(df[0] == chr_chosen)]
    x_values = [df2[1]/BIN, df2[2]/BIN]
    y_values = [df2[3], df2[3]]
    #plt.plot(x_values, y_values, linewidth=3.)
    
    fig =plt.figure(i_fig)
    fig.set_size_inches(12, 8)
    gs = gridspec.GridSpec(2,2,height_ratios=[10,1],width_ratios=[1,1])
    BIN = 2000

    # Plot for one chr
    m12_a=c1.matrix(balance=False).fetch(chr_chosen, chr_chosen)
    m12_b=c2.matrix(balance=False).fetch(chr_chosen, chr_chosen)
    
    m12_1=c1.matrix(balance=True).fetch(chr_chosen, chr_chosen)
    m12_2=c2.matrix(balance=True).fetch(chr_chosen, chr_chosen)
    
    # # plot for whole genome 
    # m12_a=c1.matrix(balance=False)[:]
    # m12_b=c2.matrix(balance=False)[:]
    
    # m12_1=c1.matrix(balance=True)[:]
    # m12_2=c2.matrix(balance=True)[:]
    
    v1=0
    v2=np.nanmax(m12_1**expo)
    
    ax1 = plt.subplot(gs[0])
    ax1.imshow(m12_1**expo,interpolation='none',vmin= v1, vmax= v2,cmap="afmhot_r")
    ax1.set_title("Contacts map "+chr_chosen)
    plt.xticks(ticks=plt.xticks()[0][1:], labels=factor_bin * np.array(plt.xticks()[0][1:], dtype=np.float64))
    plt.yticks(ticks=plt.xticks()[0][1:], labels=factor_bin * np.array(plt.xticks()[0][1:], dtype=np.float64))
    
    ax1.set_xlim(-1,np.shape(m12_1)[0])
    ax1.set_ylim(np.shape(m12_1)[0],-1)
    
    ax3 = plt.subplot(gs[2], sharex=ax1)
    #ax3.set_ylim([min(gc_signal(BIN)/100), max(gc_signal(BIN)/100)]);
    ax3.plot(x_values, y_values, linewidth=3.)
    ax3.set_xlabel("Position along the genome (in kb)")
    ax3.set_ylabel("Mappability")
    plt.xticks(ticks=plt.xticks()[0][1:], labels=factor_bin * np.array(plt.xticks()[0][1:], dtype=np.float64))
    
    #ax3.set_title("red=rRNA operons,HT=yellow")
    
    #for i in range(0,len(r)) :
    #    rectangle = plt.Rectangle(   (r[i,0]/BIN,min(gc_signal(BIN)/100)),(r[i,1]-r[i,0])/BIN,0.2 , fc='r'  )
    #    plt.gca().add_patch(rectangle)
    #
    #for i in range(0,len(p)) :
    #    rectangle = plt.Rectangle(   (p[i,0]/BIN,min(gc_signal(BIN)/100)),(p[i,1]-p[i,0])/BIN,0.2 , fc='yellow'  )
    #    plt.gca().add_patch(rectangle)
    
    ax2 = plt.subplot(gs[1], sharex=ax1,sharey=ax1)
    ax2.imshow(m12_2**expo,interpolation='none',vmin= v1, vmax= v2,cmap="afmhot_r")
    #ax2.imshow( (MAT_INT1+MAT_INT2)**expo,interpolation='none');
    #ax2.imshow( (matscn1)**expo,interpolation='none');
    ax2.set_title("reconstruction with Hicberg "+chr_chosen)
    plt.xticks(ticks=plt.xticks()[0][1:], labels=factor_bin * np.array(plt.xticks()[0][1:], dtype=np.float64))
    plt.yticks(ticks=plt.xticks()[0][1:], labels=factor_bin * np.array(plt.xticks()[0][1:], dtype=np.float64))
    
    ax2.set_xlim(-1,np.shape(m12_1)[0])
    ax2.set_ylim(np.shape(m12_1)[0],-1)
    
    ax4 = plt.subplot(gs[3], sharex=ax1);
    ax4.plot(  (m12_a).sum(axis=0) ,linewidth=2.0,label="coverage before", color="blue")
    #ax4.plot(  (MAT_INT1+MAT_INT2).sum(axis=0) ,linewidth=2.0,label="coverage after");
    ax4.plot(  (m12_b).sum(axis=0) ,linewidth=2.0,label="coverage after", color="green")
    #ax4.axhline(y= np.mean( m1.sum(axis=0)),color="red",linewidth=3.0)
    #ax4.legend(loc=4);
    ax4.set_title("Hi-C coverage: blue=before, green=after")
    ax4.set_ylabel("Coverage")
    ax4.set_xlabel("Position along the genome (in kb)")
 
    # ax4.plot(292000/BIN,0.0,"|")
    # ax4.plot(295000/BIN,0.0,"|")
    
    # ax4.plot(11000/BIN,0.0,"|")
    # ax4.plot(14000/BIN,0.0,"|")
    
    # ax4.plot(200442/BIN,0.0,"|")
    # ax4.plot(200969/BIN,0.0,"|")
    
    # ax4.plot(198671/BIN,0.0,"|")
    # ax4.plot(201177/BIN,0.0,"|")
    # ax4.plot(212535/BIN,0.0,"|")
    # ax4.plot(212720/BIN,0.0,"|")   

    plt.tight_layout(pad=3.0)
    plt.savefig("CM_"+chr1+".svg", dpi=600, format="svg")



# ratio plot
np.max(m12_2-m12_1)
m=m12_2-m12_1
v = np.reshape(m, m.shape[0]*m.shape[0])
np.median(v)
np.max(v)

plt.hist(v, np.max(v))

plt.imshow(m12_2-m12_1,interpolation='none',vmin= 0, vmax= 100,cmap="Blues")
plt.colorbar()


