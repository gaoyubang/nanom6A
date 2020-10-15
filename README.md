# Nanom6A

**Quantitative profiling of N6-methyladenosine at single-base resolution using Nanopore direct RNA sequencing**
Two step detail of nanom6A
**Download scripts and run test files**
**Details about how this pipeline works**

## Download scripts and run test files
`git clone https://github.com/gaobaibai/nanom6A.git`
The scripts was pre compile into binary. If the binary was not work, you can also compile by yourself
`wget http://forestry.fafu.edu.cn/tool/nanom6A/*gz`
## Details about how this pipeline works
###1. Preprocess
######Nanopore Basecalling using guppy (version 3.6.1)
`guppy_basecaller -i $f5 -s guppy --num_callers 40 --recursive  --fast5_out --config rna_r9.4.1_70bps_hac.cfg `
######single big fast5 to small size fast5 file
`single_to_multi_fast5 -i guppy -s single -t 40 --recursive -n 8000`
######resquiggle raw signals
`tombo resquiggle --overwrite --basecall-group Basecall_1D_001 single referance.transcript.fa --processes 40 --fit-global-scale  --include-event-stdev`
######list all fast5 file
`find single -name "*.fast5" >all_f5.list `
###2.Extarct normalized raw signals and predict
######extract signals
`extract_raw_and_feature_fast --cpu=8  --fl=all_f5.list -o extrat.feature --clip=10` 
######predict m6A site
`predict_sites -i extrat.feature -o result_final -r $refbed6 -g $genome`
                        
###3. Visualization of m6A sites in single-base and single DRS read resolution
`plot_m6A -c nanom6A/conf_plot.txt`

#####All suggestions are welcome to lfgu@fafu.edu.cn or yubanggaofafu@gmail.com 
