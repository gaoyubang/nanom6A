export PATH=/media/gyb/big_enough1/H/nanom6A/human/nanom6A_pre_compile/test_files/bin:$PATH
data=data/f5
output=result

cd model && gunzip *gz && cd ../
cd data && gunzip *gz && cd ../

find $data  -name "*.fast5"  >files.txt
extract_raw_and_feature_fast --cpu=20  --fl=files.txt -o result --clip=10
predict_sites --cpu 20  -i result -o result_final -r data/anno.bed -g data/anno.fa
nanoplot --input result_final  -o plot_nano_plot 
