**Nanom6A**

**Quantitative profiling of N6-methyladenosine at single-base resolution using Nanopore direct RNA sequencing**

**Install quick-start**

In order to make it easy to install Nanom6A, we provided three different methods for users.

**(1). Installing a pre-compiled binary release (the first method to use nanom6A)**

To use the binary package, simply download the pre-compiled Linux binary from following link:
 https://drive.google.com/drive/folders/1Dodt6uJC7lBihSNgT3Mexzpl_uqBagu0?usp=sharing


Users can untar binary_2020_10_29.tar.gz, and make sure the binaries in your PATH environment variable.   
Testing the pre-compiled binary installation:


```
tar -xvzf binary_2020_10_29.tar.gz
cd binary_2020_10_29
sh run_binary.sh
```


We tested pre-compiled binary release in **ubuntu** and **centos**.

For **ubuntu 16.04**, user may need to install libgomp1 
and libxcb1:

```
apt install libgomp1
apt install libxcb1
```

For **ubuntu 20.10**, user may need to install libncurses5.


```
apt install libncurses5
```


For **centos 8.2**, user may need to install ncurses.

```
yum install ncurses*
```


**(2). Testing Nanom6A from source (the second method to use nanom6A)**


**In order to test nanopore from source code**, you can install the dependence through conda. 

**Firstly, user can install miniconda or conda**

```
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh
chmod 777 Miniconda3-py37_4.8.3-Linux-x86_64.sh
./Miniconda3-py37_4.8.3-Linux-x86_64.sh
# After installation, please close the shell terminal and open a new one.
```


**Install conda environment**


```
tar -xvzf binary_2020_10_29.tar.gz
cd binary_2020_10_29
conda env create -f step1.yml #install conda environment for extract_raw_and_feature_fast.py
conda env create -f step2.yml #install conda environment for predict_sites.py and nanoplot.py
```



Following list was the detailed dependence:

**Source code dependence of extract_raw_and_feature_fast.py**

soft or module | version
---|---
python                               |2.7.15
h5py                               |2.9.0
statsmodels                        |0.10.0
numpy                              |1.16.6
tqdm                               |4.32.1

**Source code dependence of predict_sites.py and nanoplot.py**

soft or module | version
---|---
bedtools | v2.29.2
samtools | 1.3.1
minimap2 | 2.17-r941
python                               |3.7.3
joblib                        |0.16.0
xgboost                       |0.80
pysam                         |0.16.0.1
tqdm                          |4.39.0
pycairo                       |1.19.1
scikit-learn              |0.22

Testing the installation (Please make sure the dependence was installed).

```
tar -xvzf binary_2020_10_29.tar.gz
cd binary_2020_10_29
sh run_source_code.sh
```


**(3). Testing docker (the third method to use nanom6A)**

This is the simplest way to use nanom6A and save users from occasionally frustrating process

```
sudo docker pull gaoyubang/nanom6a:v0
```


Testing the Docker:


```
tar -xvzf binary_2020_10_29.tar.gz
cd binary_2020_10_29
sudo docker run -it -v `pwd`:/data gaoyubang/nanom6a:v0 /bin/bash
cd /data/
sh run_docker.sh
```



**FAQ**

1. The screen/nohup log file might show following output:


```
cat: result_final/AGACA.mod: No such file or directory
cat: result_final/AGACT.mod: No such file or directory
cat: result_final/GAACA.mod: No such file or directory
cat: result_final/GAACT.mod: No such file or directory
```


The tested ACTB gene  in the package did not contains m6A modification of above four kmers based on current Nanopore DRS reads.

2. The screen/nohup log file might show following output:
```
Fontconfig error: Cannot load default config file
```

The incompatibility from cairo library caused this problem,which donot obstruct m6A identification. User can fix this issue using following method: https://github.com/alacritty/alacritty/issues/2675 

**Manual for Nanom6A**

**1. Preproccess**

**Basecalling using guppy (version 3.6.1)**


```
guppy_basecaller -i $f5 -s guppy --num_callers 40 --recursive --fast5_out --config rna_r9.4.1_70bps_hac.cfg
```

**merged single big fast5 into small size fast5 file**

```
single_to_multi_fast5 -i guppy -s single -t 40 --recursive -n 8000
```



**resquiggle raw signals**



```
tombo resquiggle --overwrite --basecall-group Basecall_1D_001 single referance.transcript.fa --processes 40 --fit-global-scale --include-event-stdev
```

**list all fast5 file**

```
find single -name "*.fast5" >files.txt
```

2. **identification of m6A sites based on DRS reads**

**extracting signals**


```
extract_raw_and_feature_fast --cpu=20 --fl=files.txt -o result --clip=10
```

**predicting m6A site**

```
predict_sites --cpu 20 -i result -o result_final -r data/anno.bed -g data/anno.fa
```


The main output is the ratio.x.tsv and genome_abandance.x.bed in the output dir.
The header of ratio.x.tsv.

gene\|chrom | coordinate\|mod number\|total number\|mod ratio 
---|---
ACTB\|chr7|	5566813\|162\|639.0\|0.2535211267605634	

The header of genome_abandance.x.bed.

chrom | coordinate|gene| read id|read pos|kmer
---|---|---|---|---|---
chr7|5567320|ACTB	|e88129423ae1.fast5	|1257	|AAACA

**3. Visualization of m6A sites**


```
nanoplot --input result_final -o plot_nano_plot
```


**4. Train your own model**


```
python train.py
```



**The fast5 raw file of nanopore direct RNA sequence in our unpublished studies:**

repeat1:

https://sra-download.ncbi.nlm.nih.gov/traces/sra36/SRZ/012881/SRR12881185/poplar_guppy_recall.tar.gz

repeat2:

https://sra-download.ncbi.nlm.nih.gov/traces/sra68/SRZ/012822/SRR12822922/Ptr-WT-SDX-20200827.tar

All suggestions are welcome to lfgu@fafu.edu.cn or yubanggaofafu@gmail.com

