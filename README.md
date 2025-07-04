# Codes developped for the bioanalysis with Hicberg

Codes and scripts to reproduce the analyses and plots. 

### Dependencies

Scripts and codes can be run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=3.10)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `SRA tools` / [SRA tools](https://github.com/ncbi/sra-tools)
* `Bowtie2` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `hicstuff` / [hicstuff](https://github.com/koszullab/hicstuff)
* `cooler` / [cooler](https://github.com/open2c/cooler)
* `tinyMapper` / [tinyMapper](https://github.com/js2264/tinyMapper)
* `deepTools` / [deepTools](https://deeptools.readthedocs.io/en/develop/)
* `pyBigWig`   /   [pyBigWig](https://github.com/deeptools/pyBigWig)


### Raw data extraction and alignment
#### Data extraction
Data can be dowloaded on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.

A SRA executable called fastq-dump from SRA toolkit can be used to download and split both mates of a NGS library: 
 
```bash
./fastq-dump SRR639031 --split-3 -O /home/data/
```

#### Alignment of the Hi-C libraries with hicberg
To align the reads and generate the contact files in cool format, we used hicstuff pipeline: 
```bash
hicberg pipeline -f -e DpnII -e HinfI -c plasmid_p2-micron -c chrMT   --cpus 8 -o T2T_S288C/  -n out_hicberg_AC1_T2T/ T2T_BY4742.asm01.HP0_2micron.fa   AC1.end1.fastq.gz   AC1.end2.fastq.gz
```

### Processing of genomic data like Mnase-seq, ChIP-seq after hicberg
We used bedtools and bedGraphToBigWig: 
```bash

awk '{if($2==$5) {if(sqrt(($6-$3)^2)<1000) {if($3<$6) {printf "%s\t%s\t%s\t%d\n" , $2,$3,$6,"1";} if($3>$6) {printf "%s\t%s\t%s\t%d\n" , $2,$6,$3,"1";}}}}' group1.pairs > group1.pairs2

bedtools coverage -a /home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt3 -b group1.pairs2 -d > group1.pairs2.bp

awk '{if($3!=$4) print $0}' group1.pairs2.bp > group1.pairs2.bp2  # remove of last bp of each chr

awk '{print $1"\t"$4"\t"$4+1"\t"$5}'  group1.pairs2.bp2 > group1.pairs2.bp.bedgraph  # conversion in bedgraph format 

sort -k1,1 -k2,2n group1.pairs2.bp.bedgraph > group1.pairs2.bp.sorted.bedgraph      # sort 

/home/axel/Bureau/tools/./bedGraphToBigWig group1.pairs2.bp.sorted.bedgraph /home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt  group1.pairs_sorted.bw   # conversion in bw  
```


#### Computation of contact signal with ty1 of 2 micron molecular parasite without and with hicberg
```bash
python3 /home/axel/Bureau/z_python_scripts_copy/plasmid_HSC_1D_agglo.py /home/axel/Bureau/hicberg_project/test_20240105/AC1/matrices/rescued_map.cool  plasmid_p2-micron AC_rescued_with_ty1_rescued /home/axel/Bureau/YEAST/ty1/ty1_scerevisiae.txt5
```

Negatif control was computed using chrM (mitochondria molecule)
```bash
python3 /home/axel/Bureau/z_python_scripts_copy/plasmid_HSC_1D_agglo.py /home/axel/Bureau/hicberg_project/test_20240105/AC1/matrices/rescued_map.cool chrM AC_rescued_chrM_with_ty1_rescued /home/axel/Bureau/YEAST/ty1/ty1_scerevisiae.tx
```

#### Computation of cohesin signals around TY elements

```bash
bw_file1="unrescued_Chip_over_input.bw"
bw_file2="rescued_Chip_over_input.bw"

computeMatrix scale-regions -p 8 --averageTypeBins sum \
-S  $bw_file1 $bw_file2 \
-R /home/axel/Bureau/YEAST/TY/TY_elements.renamed.bed  --beforeRegionStartLength 5000 --regionBodyLength 5000 --afterRegionStartLength 5000 --outFileName signal_scc2.gz

name="ChIP_of_Scc1_LM96_TY"
name2="ChIP_of_Scc1"

plotProfile -m signal_scc2.gz -out Scc1_peaks_on_$name.pdf     --perGroup --plotTitle "ChIP-seq Scc1" --colors royalblue salmon --samplesLabel "unrescued" "rescued"   --regionsLabel "on TY in S288C"  --startLabel "TY"  --endLabel " "  --yAxisLabel $name2  --yMin 0 0 --yMax 100 100

plotHeatmap  -m signal_scc2.gz  -out map_Scc1_peaks_on_$name.pdf   --colorMap 'Blues'  --perGroup --plotTitle "ChIP-seq Scc1"  --samplesLabel "unrescued" "rescued"  --regionsLabel "on TY in S288C"  --startLabel "TY"  --endLabel " "  --yAxisLabel $name2   --yMin 0 0  --yMax 100 100   --zMin  0 --zMax 100 
```





