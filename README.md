# Nanom6A

**Quantitative profiling of N6-methyladenosine at single-base resolution using Nanopore direct RNA sequencing**

The fast5 raw file of nanopore direct rna sequence  in this method is below:

rep1:

Under uploading to ncbi.

rep2:

https://sra-download.ncbi.nlm.nih.gov/traces/sra68/SRZ/012822/SRR12822922/Ptr-WT-SDX-20200827.tar

Indexs

**Download scripts and run test files**

**Details about how this pipeline works**

**Train your own model**
## Download  pre-compiled bianry or source code or docker image

**(1)**
The scripts was pre-compiled into binary. You can download from google drive.
https://drive.google.com/drive/folders/1Dodt6uJC7lBihSNgT3Mexzpl_uqBagu0?usp=sharing

**(2)**
 If the binary was not work, you can also compile by yourself. The Dependence was show blow.

`git clone https://github.com/gaobaibai/nanom6A.git`

**Source code dependence**

###### extract_raw_and_feature_fast.py
soft or module | version
---|---
python                               |2.7.15
h5py                               |2.9.0
statsmodels                        |0.9.0
numpy                              |1.16.6
tqdm                               |4.32.1
######  predict_sites.py and plot_m6A.py


soft or module | version
---|---
bedtools | v2.29.2
samtools | 1.3.1
minimap2 | 2.17-r941
python                               |3.7.3
joblib                        |0.14.1
xgboost                       |0.90
pysam                         |0.16.0.1
tqdm                          |4.39.0
pycairo                       |1.19.1
pyBigWig                      |0.3.17
**(3)**
 The easiest way is running with docker.
 
**Under construction.**

## Run test files

binary python executable file
```
unzip binary_2020_10_15.tar.gz
tar -xvzf c.tar.gz 
export PATH=`pwd`/binary_2020_10_15.tar.gz:$PATH
```
source code and example file
```
git clone https://github.com/gaobaibai/nanom6A.git
cd nanom6A
cd model && gunzip *gz && cd ../
cd data && gunzip *gz && cd ../
#binary
sh run.sh
#source code
run_source_code.sh
```



## Details about how this pipeline works
##### 1. Preprocess

###### Nanopore Basecalling using guppy (version 3.6.1)

`guppy_basecaller -i $f5 -s guppy --num_callers 40 --recursive  --fast5_out --config rna_r9.4.1_70bps_hac.cfg `

###### single big fast5 to small size fast5 file

`single_to_multi_fast5 -i guppy -s single -t 40 --recursive -n 8000`

###### resquiggle raw signals

`tombo resquiggle --overwrite --basecall-group Basecall_1D_001 single referance.transcript.fa --processes 40 --fit-global-scale  --include-event-stdev`

###### list all fast5 file

`find single -name "*.fast5" >all_f5.list `

##### 2.Extarct normalized raw signals and predict

###### extract signals

`extract_raw_and_feature_fast --cpu=8  --fl=all_f5.list -o extrat.feature --clip=10` 

###### predict m6A site

`predict_sites -i extrat.feature -o result_final -r $refbed6 -g $genome`
                        
The main output is the ratio.x.tsv and genome_abandance.x.bed in the output dir.

The header of ratio.x.tsv.
gene\|chrom | coordinate\|mod number\|total number\|mod ratio 
---|---
ACTB|chr7	5566813\|162\|639.0\|0.2535211267605634	


The header of genome_abandance.x.bed.

chrom | coordinate|gene| read id|read pos|kmer
---|---|---|---|---|---
chr7|5567320|ACTB	|e88129423ae1.fast5	|1257	|AAACA
##### 3. Visualization of m6A sites in single-base and single DRS read resolution

`plot_m6A -c nanom6A/conf_plot.txt`

### Train your own model

`python train.py`
##### All suggestions are welcome to lfgu@fafu.edu.cn or yubanggaofafu@gmail.com 
