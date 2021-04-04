**Nanom6A**

**Quantitative profiling of N6-methyladenosine at single-base resolution using Nanopore direct RNA sequencing**

**Install quick-start**

In order to make it easy to install Nanom6A, we provided three different methods for users.

**(1). Installing a pre-compiled binary release (the first method to use nanom6A)**

To use the binary package, simply download the pre-compiled Linux binary from following link:
 https://drive.google.com/drive/folders/1Dodt6uJC7lBihSNgT3Mexzpl_uqBagu0?usp=sharing


Users can untar nanom6A_2021_3_18.tar.gz, and make sure the binaries in your PATH environment variable.   
Testing the pre-compiled binary installation:



```
tar -xvzf nanom6A_2021_3_18.tar.gz
cd nanom6A_2021_3_18
sh run_binary.sh
```

User may still need sam2tsv in your $PATH (after 2021_3_18 version), you can install it through conda.



```
conda install -c hcc jvarkit-sam2tsv
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

user may need to install libxrender and libxext.

```
apt install -y libxrender-dev
apt install -y libxext-dev
```

```
tar -xvzf nanom6A_2021_3_18.tar.gz
cd nanom6A_2021_3_18
conda env create -f conda.yml #install conda environment
```



Following list was the detailed dependence:



**Source code dependence**

soft or module | version
---|---
bedtools | v2.29.2
samtools | 1.3.1
minimap2 | 2.17-r941
python                               |3.7.3
h5py                               |2.9.0
statsmodels                        |0.10.0
joblib                        |0.16.0
xgboost                       |0.80
pysam                         |0.16.0.1
tqdm                          |4.39.0
pycairo                       |1.19.1
scikit-learn              |0.22

Testing the installation (Please make sure the dependence was installed).

```
tar -xvzf nanom6A_2021_3_18.tar.gz
cd nanom6A_2021_3_18
sh run_source_code.sh
```


**(3). Testing docker (the third method to use nanom6A)**

This is the simplest way to use nanom6A and save users from occasionally frustrating process

```
sudo docker pull gaoyubang/nanom6a:v1
```


Testing the Docker:


```
tar -xvzf nanom6A_2021_3_18.tar.gz
cd nanom6A_2021_3_18
sudo docker run -it -v `pwd`:/data gaoyubang/nanom6a:v1 /bin/bash
cd /data/
sh run_docker.sh
```



**FAQ**

1. The screen/nohup log file might show following output:
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

**Convert merged single big fast5 into small size fast5 file**

When the fast5 file was stored in multi_read formats, this step is required (mostly seuqenced with SQK-RNA002 kit).
You can check the size of one fast5 file, if it's about ~300MB, this step is required. If it's about ~100KB, You have to skip this step to avoid error.

```
multi_to_single_fast5 -i guppy -s single -t 40 --recursive 
```



**resquiggle raw signals**

The tombo resquiggle referance.transcript.fa should not be genome file, it should be the referance gene fasta file.

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

the -g option should be the genome file.

```
predict_sites --cpu 20 -i result -o result_final -r data/anno.bed -g data/anno.fa
```
(1) the -r parameter is file of referance corrd of transcripts.
It should contains 6 colums

chrom | start| end| name| .| strand
---|---|---|---|---|---
chr7|	5566779	|5570232|	ACTB|	.|	-

(2) please check your genome file index, make sure you index with samtools index ref.fa and picard CreateSequenceDictionary R=ref.fa O=ref.dict 




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



**The fast5 raw file of nanopore direct RNA sequence in our published studies:**

repeat1:

https://sra-download.ncbi.nlm.nih.gov/traces/sra36/SRZ/012881/SRR12881185/poplar_guppy_recall.tar.gz

repeat2:

https://sra-download.ncbi.nlm.nih.gov/traces/sra68/SRZ/012822/SRR12822922/Ptr-WT-SDX-20200827.tar


**Acknowledgement**

**2020.11.6 16:56 Fuzhou**

LiuFuyuxaing helped me with testing the code and improvement of the Manual!

**2021.3.4 22:00 Fuzhou**

Fixed bugs due to difference between tombo aligned sequence and minimap2 aligned sequence.

Fixed bugs due to difference between Sam2tsv and samtools depth.

Thank you to Hang Qin from Institute of Plant Physiology and Ecology for bringing this to our attention!

**2021.3.5 12:00 Fuzhou**

Update binary and docker version due to bugs finding at 2021.3.4!

**2021.3.11 22:00 Fuzhou**

Fixed bugs lead to 1-based m6A sites in negative strand , Thank you to Yan Xin from The University of Hong Kong for bringing this to our attention!

**2021.3.18 19:00 Fuzhou**

Fixed bugs due to samtools depth default 8000 maximum coverage!


All suggestions are welcome to lfgu@fafu.edu.cn or yubanggaofafu@gmail.com


