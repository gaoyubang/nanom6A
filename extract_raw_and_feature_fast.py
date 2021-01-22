from __future__ import absolute_import
import argparse
import os,sys,re,h5py
import multiprocessing
from statsmodels import robust
import numpy as np
#import joblib
from tqdm import tqdm
##################################
def get_label_raw(fast5_fn, basecall_group, basecall_subgroup,reverse = False):
	##Open file
	try:
		fast5_data = h5py.File(fast5_fn, 'r')
	except IOError:
		raise IOError('Error opening file. Likely a corrupted file.')

	# Get raw data
	try:
		raw_dat = list(fast5_data['/Raw/Reads/'].values())[0]
		# raw_attrs = raw_dat.attrs
		raw_dat = raw_dat['Signal'][()]
		# ~ .value
	except:
		raise RuntimeError(
			'Raw data is not stored in Raw/Reads/Read_[read#] so ' +
			'new segments cannot be identified.')

	# Read corrected data
	try:
		corr_data = fast5_data['/Analyses/'+basecall_group +'/' + basecall_subgroup + '/Events']
		corr_attrs = dict(list(corr_data.attrs.items()))
		corr_data = corr_data[()]
		# ~ .value
	except:
		raise RuntimeError(('Corrected data not found.'))

	#fast5_info = fast5_data['UniqueGlobalKey/channel_id'].attrs
	# sampling_rate = fast5_info['sampling_rate'].astype('int_')

	# Reading extra information
	corr_start_rel_to_raw = corr_attrs['read_start_rel_to_raw']  #
	if len(raw_dat) > 99999999:
		raise ValueError(fast5_fn + ": max signal length exceed 99999999")
	if any(len(vals) <= 1 for vals in (corr_data, raw_dat)):
		raise NotImplementedError(('One or no segments or signal present in read.'))
	event_starts = corr_data['start'] + corr_start_rel_to_raw
	event_lengths = corr_data['length']
	event_bases = corr_data['base']
	fast5_data.close()
	# ~ label_data = np.array(
		# ~ list(zip(event_starts, event_lengths, event_bases)),
		# ~ dtype=[('start', '<u4'), ('length', '<u4'), ('base', 'S1')])
	return (raw_dat, event_bases, event_starts, event_lengths)
	#######################################
def search_RRACH(signal,start,length,base,fn_string):
	uniq_arr=np.unique(signal)
	signal = (signal - np.median(uniq_arr)) / np.float(robust.mad(uniq_arr))
	raw_signal = signal.tolist()
	kmer_fillter="[AG][AG]AC[ACT]"
	line=""
	total_seq="".join([x.decode() for x in base])
	clipnum=int(FLAGS.clip)
	# ~ print(length)
	# ~ print(base)
	# ~ print(len(length))
	for indx in range(len(length)):
		if 2+clipnum<=indx<=len(length)-3-clipnum:
			base0,base1,base2,base3,base4=[base[indx+x].decode() for x in [-2,-1,0,1,2]]
			kmer_now_t="%s%s%s%s%s"%(base0,base1,base2,base3,base4)
			# ~ print(kmer_now_t)
			# ~ print(indx,kmer_now_t)
			list_have=[x.start() for x in re.finditer(kmer_fillter,kmer_now_t)]
			if len(list_have)==0:
				continue
			raw_signal_every=[raw_signal[start[indx+x]:start[indx+x]+length[indx+x]] for x in [-2,-1,0,1,2]]
			mean=[np.mean(x) for x in raw_signal_every]
			std=[np.std(x) for x in raw_signal_every]
			md_intense = [np.median(x) for x in raw_signal_every]
			length2=[length[indx+x] for x in [-2,-1,0,1,2]]
			#############
			line+="%s|%s|%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(str(fn_string).split("/")[-1],indx,"N",base2,kmer_now_t,"|".join([str(x) for x in mean]),"|".join([str(x) for x in  std]),"|".join([str(x) for x in  md_intense]),"|".join([str(x) for x in length2]),kmer_now_t)
	# ~ print(line)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
	return line
def extract_file(input_file):
	try:
		(raw_data, raw_label, raw_start, raw_length) = get_label_raw(input_file, FLAGS.basecall_group,FLAGS.basecall_subgroup)
	except Exception as e:
		# ~ print(str(e))
		return False, (None, None)
	raw_data = raw_data[::-1]
	# ~ print(input_file,raw_start,raw_length,raw_label)
	total_seq="".join([x.decode() for x in raw_label])
	ids=input_file.split("/")[-1]
	total_seq=">%s\n%s\n"%(ids,total_seq)
	line=search_RRACH(raw_data,raw_start,raw_length,raw_label,input_file)
	return line,total_seq
    # ~ print(line)
def subcon(fls):
	# ~ if Multile_processing=="True":
	if True:
		results=[]
		##########
		pool = multiprocessing.Pool(processes = int(FLAGS.cpu))
		##########	
		for fl in fls:
			# ~ print(fl)
			result=pool.apply_async(extract_file,(fl,))
			results.append(result)
		pool.close()
		############################
		pbar=tqdm(total=len(fls),position=0, leave=True)
		nums=[]
		for result in results:
			num,seq=result.get()
			if num and seq:
				nums.append([num,seq])
			pbar.update(1)
		#############################
		pool.join()
		#############################
	else:
		nums=[]
		for fl in fls:
			# ~ print(fl)
			num,seq=extract_file(fl)
			if num and seq:
				nums.append([num,seq])
	output=open(FLAGS.output+".feature.tsv","w")
	output.write("".join([str(x[0]) for x in nums]))
	output.close()
	##########################
	# ~ lines="".join([str(x[0]) for x in nums]).rstrip().split("\n")
	# ~ X,Y=[],[]
	# ~ for line in lines:
		# ~ ele=line.rstrip().split()
		# ~ GXB01170_2018.fast5|233|23,3,1,00,2,3	0	0|0|0|0|0	0.47062142444662086|0.8176029853529686|0.7531814474848483|-0.0835983106934529|-0.5416802793696001	0.14888963355157137|0.1977707177676313|0.11147011292496822|0.18032028688587404|0.15884083915957087	0.5053374754088856|0.9240774715516673|0.7871817035819118|-0.06103521168167164|-0.5495652071815835	15|6|6|138|64
		# ~ insert=[]
		# ~ for item in ele[3],ele[4],ele[5],ele[6]:
			# ~ for itemsub in item.split("|"):
				# ~ insert.append(float(itemsub))
		# ~ X.append(insert)
		# ~ Y.append(ele[0])
	# ~ joblib.dump(X,"%s.joblib.X"%(FLAGS.output))
	# ~ joblib.dump(Y,"%s.joblib.Y"%(FLAGS.output))
	#############################
	output=open(FLAGS.output+".feature.fa","w")
	output.write("".join([str(x[1]) for x in nums]))
	output.close()
#########################################################################################################################
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Extract fast5 files.')
	# ~ parser.add_argument('-i', '--input', required = True,
						# ~ help="Directory that store the fast5 files, multiple directories separted by commas.")
	parser.add_argument('-o', '--output', required = True, help="Output file")
	parser.add_argument('--basecall_group',default = "RawGenomeCorrected_000",
                        help='The attribute group to extract the training data from. e.g. RawGenomeCorrected_000')
	# ~ parser.add_argument('-f', '--tffile', default="train.tfrecords",
						# ~ help="tfrecord file")
	parser.add_argument('--basecall_subgroup', default='BaseCalled_template',help='Basecall subgroup Nanoraw resquiggle into. Default is BaseCalled_template')
	parser.add_argument('--cpu', default=1,help='cpu number usage')
	parser.add_argument('--clip', default=10,help='reads first and last N base signal drop out')
	parser.add_argument('--fl',required = True,help='files comtained fast5 path, one line, one fast5 file')
	args = parser.parse_args(sys.argv[1:])
	# ~ run(args)
	global FLAGS
	FLAGS = args
	total_fl=[]
	for i in open(FLAGS.fl,"r"):
		total_fl.append(i.rstrip())
	subcon(total_fl)
##################################################################################
