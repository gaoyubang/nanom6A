#!/usr/bin/env python
import argparse
import joblib,sys,os,re,gzip
from collections import defaultdict
from xgboost.sklearn import XGBClassifier
from pysam import FastaFile
from tqdm import tqdm
import subprocess
###########################################
def predict(tpp,fl1,dir1):
# ~ ######################
	# ~ ins_match={}
	# ~ fl="signal_match.tsv"
	# ~ for i in open(fl,"r"):
		# ~ ele=i.rstrip().split()
		# ~ ins_match[ele[0]]=float(ele[-2])
	###########################################################
	# ~ FLAGS.abspath="/media/gyb/big_enough1/H/nanom6A/epinano_rep1"
	model = joblib.load(filename="%s/%s.m"%(FLAGS.model.rstrip("/"),tpp))
	X=[]
	Y=[]
	for i in open(fl1,"r"):
		ele=i.rstrip().split()
		# ~ GXB01170_2018.fast5|233|23,3,1,00,2,3	0	0|0|0|0|0	0.47062142444662086|0.8176029853529686|0.7531814474848483|-0.0835983106934529|-0.5416802793696001	0.14888963355157137|0.1977707177676313|0.11147011292496822|0.18032028688587404|0.15884083915957087	0.5053374754088856|0.9240774715516673|0.7871817035819118|-0.06103521168167164|-0.5495652071815835	15|6|6|138|64
		# ~ if ele[-1]!=tpp:
			# ~ continue
		# ~ ids=ele[0].split("|")[0]
		# ~ cv=ins_match[ids]
		insert=[]
		for item in ele[3],ele[4],ele[5],ele[6]:
			# ~ item=item.decode("utf-8")
			for itemsub in item.split("|"):
				insert.append(float(itemsub))
		#ele[0]="|".join(ele[0].split("|")[:-1])
		# ~ insert.append(cv)
		X.append(insert)
		Y.append(ele[0])
	if X and Y:
		pass
	else:
		return False
	##########################################################
	output1=open("%s/%s.mod"%(dir1,tpp),"w")
	output2=open("%s/%s.unmod"%(dir1,tpp),"w")
	model.nthread=-1
	results=model.predict_proba(X)
	for i in range(len(results)):
		punmod=results[i][0]
		pmod=results[i][1]
		if pmod>punmod:
			output1.write("%s\t%s\t%s\n"%(punmod,pmod,Y[i]))
		else:
			output2.write("%s\t%s\t%s\n"%(punmod,pmod,Y[i]))
	output1.close()
	output2.close()
	return True
#################################################################
def sampare(i,storepos,store):
	# ~ for i in open(sys.argv[1],"r"):
	if i.startswith("@"):
		return False
	dd=i.rstrip().split()
	if dd[2]=="*":
		return False
# ~ GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_545_ch_160_strand.fast5        16      Chr05   13495549        60      274M71N205M     *       0       0       AATTCCAATCTCTCAAGTATGTTACAATATGTATAATCAATCATACATGCTTTCAGTGCTGAAGAGTTCCCCTCAACGACGACAAAGACAACAGCAAAACGAACAAGAGAGTTACAGAGACACTGTCAAAGAGAAACAGCGTCGATTCCATTATTTACAAAAATCCAACTTTATTTTGCCTTTGCAGCAGATGCAGTGAAAGCAGAAAACAGTCTACAAATCTCTTTTTATTCCGGGAGAAATTAAGACAACAAACACGAACTTAAACGCAAGGCATTTCCTTTGGGGTTATTGGGATCGACACCAGGCCATTCACGGCAGAGACAGTGTTGCCCCTCAGGGATGCCCCTATAAACGGAAGGAACATGCCTGCCACAACCACCCCAGGTGGTCTTGCCGCACGTTGAGCATCTCTCCCTGTAACACATCTTTCGGCTCTGAGAGTACTGTTGATTGTTTCTTTTTATGTTTTCTTCTCT *       NM:i:0  ms:i:479        AS:i:447        nn:i:0  ts:A:+  tp:A:P  cm:i:158        s1:i:467        s2:i:0  de:f:0  rl:i:0
	ids=dd[0]
	###############
	# check modes
	if not ids in storepos:
		return False
	#####################
	(flag,pos,cigar,seq)=(int(dd[1]),int(dd[3]),dd[5],dd[9])
	if (flag&16)==0:
		strand='+'
	elif (flag&16)==16:
		strand='-'
	want=["D","I","S","H"]
	mark=0
	for subitem in want:
		if subitem in cigar:
			mark=1
			# ~ sys.stderr.write("bad cigar %s\t%s\n"%(i,cigar))
			continue
	if mark==1:
		# ~ sys.stderr.write("bad cigar %s\t%s\n"%(cigar,ids))
		return False
	######################
	reflect={}
	startpos=pos-1
	start=-1
	ele=[int(x) for y,x in enumerate(re.split('M|N',cigar)) if int(y)!=(len(re.split('M|N',cigar))-1)]
	if strand=="-":
		startpos=pos+sum(ele)
		ele=ele[::-1]
	if "M" in cigar and "N" not in cigar:
		if strand=="+":
			for subpos in range(ele[0]):
				start+=1
				startpos+=1
				reflect[start]=startpos
		elif strand=="-":
			for subpos in range(ele[0]-1,-1,-1):
				start+=1
				startpos-=1
				reflect[start]=startpos
		########
	elif "M" in cigar and "N" in cigar:
		for j in range(0,len(ele),2):
			if strand=="+":
				for exon in range(ele[j]):
					start+=1
					startpos+=1
					reflect[start]=startpos
				if j!=len(ele)-1:
					for intron in range(ele[j+1]):
						startpos+=1
			elif strand=="-":
				########################################################
				for exon in range(ele[j]-1,-1,-1):
					start+=1
					startpos-=1
					reflect[start]=startpos
				##########
				if j!=len(ele)-1:
					for intron in range(ele[j+1]-1,-1,-1):
						startpos-=1
				##########################
	#######
	#check option
	#######s
	if len(reflect)!=len(store[ids]):
		# ~ print(cigar,pos,startpos,reflect,len(store[ids]))
		# ~ linerror="%s\t%s\t%"%(cigar,pos,startpos,reflect,len(store[ids]))
		# ~ sys.stderr.write(reflect)
		# ~ sys.stderr.write(ids)
		# ~ sys.exit("not equel")
		return False
	################
	# ~ want_base=ids.split("|")[-1].split(",")
	want_base=sorted(storepos[ids].keys())
	line="%s\t%s\t%s\t%s\t%s\t"%(dd[2],ids,strand,pos,startpos)
	for subwantpos in want_base:
		subwantpos=int(subwantpos)
		mer=store[ids][subwantpos-2:subwantpos+3]
		line+="%s|%s|%s\t"%(subwantpos,reflect[subwantpos],mer)
	return line
	# ~ print(line)
def tsvparese(basefl):
	pass
def read2genome2(basefl):
	cmd="cat {0}/AAACA.mod  {0}/AAACC.mod  {0}/AAACT.mod  {0}/AGACA.mod  {0}/AGACC.mod  {0}/AGACT.mod  {0}/GAACA.mod  {0}/GAACC.mod  {0}/GAACT.mod  {0}/GGACA.mod  {0}/GGACC.mod  {0}/GGACT.mod >{0}/total_mod.tsv".format(basefl)
	os.system(cmd)
	cmd="cat {0}/*.*mod >{0}/total_prediction.tsv".format(basefl)
	os.system(cmd)
	fl="{0}/total_prediction.tsv".format(basefl)
	storepos=readprediction(fl)
	fl=FLAGS.input+".feature.fa"
	store=readfasta(fl)
	fl="{0}/extract.reference.bed12".format(basefl)
	# ~ head result_final/extract.reference.bed12
	# ~ NM_001197125.1	30	1440	1e0208b3-8061-451f-9415-73df9654f9da.fast5	0	+	30	1440	255,0,0	1	1410	0
	readgene={}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		readgene[ele[3]]=ele[0]
	fl="%s/extract.sort.bam.tsv.gz"%(basefl)
	##########################################
	##########################################
	# ~ chr04	W_003002_20180416_FAH83697_MN23410_sequencing_run_20180415_FAH83697_mRNA_WT_Col0_2918_23801_read_57_ch_290_strand.fast5	-	10201753	10201753	236|10203065|GAACA	275|10202917|AGACC	991|10202201|GGACA	1003|10202189|AGACC	1373|10201819|AAACA
	# ~ c1,c2,c3=0,0,0
	total_m6A_reads=defaultdict(dict)
	pbar=tqdm(total=len(store.keys()),position=0, leave=True)
	pre1,pre2="",""
	for i in gzip.open(fl,"r"):
		i=i.decode("utf-8").rstrip()
		# ~ NR_002323.2	0	chr22	7541	A	I	31375381	a	M
		if i.startswith("#"):
			pre1="#"
			continue
		ele=i.rstrip().split()
		if ele[3]=="." or ele[6]==".":
			continue
		ids,chro,idspos,gpos=ele[0],ele[2],int(ele[3]),ele[6]
		if ids!=pre1:
			pbar.update(1)
			pre1=ids
		if ele[1]=="0":
			strand="+"
		elif ele[1]=="16":
			strand="-"
			lens=len(store[ids])
			idspos=lens-idspos-1
		if ids in storepos and idspos in storepos[ids] and ids in readgene:
			kmer=store[ids][idspos-2:idspos+3]
			line="%s|%s|%s"%(idspos,gpos,kmer)
			total_m6A_reads["%s\t%s\t%s\t%s\tNA\t"%(chro,ids,strand,readgene[ids])][line]=1
	output=open("%s/sam_parse2.txt"%(basefl),"w")
	for item in total_m6A_reads:
		sorts=sorted(total_m6A_reads[item].keys(), key=lambda d:(int(d.split("|")[0])))
		output.write("%s\t%s\n"%(item,"\t".join(sorts)))
	output.close()
	# ~ print("sucess parse ids ratio number %s, ratio %s"%(c1,c1/(c1+c2+0.0)))
	########################################################################################
def read2genome1(basefl):
	cmd="cat {0}/AAACA.mod  {0}/AAACC.mod  {0}/AAACT.mod  {0}/AGACA.mod  {0}/AGACC.mod  {0}/AGACT.mod  {0}/GAACA.mod  {0}/GAACC.mod  {0}/GAACT.mod  {0}/GGACA.mod  {0}/GGACC.mod  {0}/GGACT.mod >{0}/total_mod.tsv".format(basefl)
	os.system(cmd)
	fl="{0}/total_mod.tsv".format(basefl)
	storepos=readprediction(fl)
	fl=FLAGS.input+".feature.fa"
	store=readfasta(fl)
	fl="%s/extract.sam"%(basefl)
	##########################################
	c1,c2=0,0
	output=open("%s/sam_parse.txt"%(basefl),"w")
	for i in open(fl,"r"):
		i=i.rstrip()
		status=sampare(i,storepos,store)
		if status:
			output.write(status+"\n")
			c1+=1
		else:
			c2+=1
	output.close()
	print("sucess parse ids ratio number %s, ratio %s"%(c1,c1/(c1+c2+0.0)))
	########################################################################################
def readprediction(fl):
	storepos=defaultdict(dict)
	for i in open(fl,"r"):
		# ~ c9a0d84d-42f6-456c-895f-e957d1172623|6|AAGGGAAAGACTCCAGAGGAAATTAGGAAGACCTTTAACATCAAGAATGACTTTACACCTGAGGAGGAGGAGGAAGTTCGCCGTGAGAACCAGTGGGCATTTGAATGAAGTGCGTCTGATGGTTTCATGGAAGGAATGTTGTTCTAATGCCAAATGAATGCTGTGGGTTATCTTAGCGTAGACAAGACTATGTTTCTATGACTTTATTGTGAACCTGTGAGCACATTGACTGTAAATAATACTTGTATTCTGGGGAGGGGATTGGTAGTAGTTTCCTGCAATCAATCCTCTGCTTGTGGGCAAATGTTATTTGTTGCAGACTTGCAGTGATCCTTATCTGTTGTATCTGTTTTCCCTCTGTGTTCCTGCCAAGTTTGTTTCTTGGACATAATCATCAAGTCTTGGTGTCTCTT	1.0	0.06883740425109862	0.9311625957489014
		# ~ 0.09221864	0.90778136	GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_17200_ch_182_strand.fast5|177|2,1,3,2,0,1,3,2,0,3,1,1,0,0,0,3,1,1,3,2,0,3,2,0,1,1,1,3,1,3,3,2,3,0,1,1,0,2,0,2,0,3,3,2,1,1,1,0,3,0,3,2,3,0,1,0,0,2,0,1,3,2,0,1,1,2,0,2,1,1,0,0,0,3,0,3,2,0,0,2,1,1,0,1,3,2,1,3,1,2,1,0,2,3,3,2,2,0,1,0,1,0,2,0,0,2,3,0,3,2,1,1,0,3,2,2,2,1,3,2,0,2,2,1,3,2,0,2,2,3,3,0,3,2,2,0,0,0,1,0,0,3,2,3,0,3,0,0,0,3,2,2,1,3,3,2,1,2,3,3,3,0,3,3,2,3,3,0,3,2,3,2,3,3,2,0,0,0,1,0,3,2,2,3,1,3,2,3,3,3,0,1,3,1,3,3,3,3,2,2,2,2,3,3,2,2,3,3,3,3,2,3,2,0,2,2,2,3,3,3,2,0,0,3,3,3,1,0,3,0,0,2,0,0,3,2,0,0,3,2,0,3,0,3,3,3,1,2,3,2,1,0,2,1,3,1,1,0,0,0,1,3,0,3,2,0,3,3,3,2,2,2,2,2,3,3,2,0,0,3,2,2,0,0,0,3,0
		ele=i.rstrip().split()
		unmod,mod,mark=ele
		ids,spos,seq=mark.split("|")
		# ~ ids=namechange[ids]
		# ~ print(ids,int(spos))
		storepos[ids][int(spos)]=1
	return storepos
def readfasta(fl):
	store={}
	lines=open(fl,"r").readlines()
	for index,i in enumerate(lines):
		if i.startswith(">"):
			ids=i.rstrip().lstrip(">")
			read=lines[index+1].rstrip()
			store[ids]=read
	return store
###################
def pare_sam_site(fl,geneids):
	store=defaultdict(dict)
# ~ Chr06	GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_61758_ch_424_strand.fast5	-	19873404	19873404	230|19874613|GGACA	532|19873617|GGACA	646|19873503|AAACT	714|19873435|GGACC	
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		ref,ids,strand=ele[0],ele[1],ele[2]
		if ids in geneids:
			genename=geneids[ids]
		else:
			genename="NA"
		for item in ele[5:]:
			spos,gpos,gbase=item.split("|")
			store[ids][spos]=gpos,ref,gbase,genename
	return store
def pare_annotation2(fl):
	geneids=defaultdict(dict)
	#NM_001197123.2	52	1650	13d4649d-79d3-4593-9cfc-14fac5bee959.fast5	0	+	52	1650	255,0,0	2	175,1020	0,578
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		# ~ for item in ele[-1].split(","):
			# ~ item=item.split(";")[0]
		geneids[ele[3]]=ele[0]
	return geneids
###############################################################
def pare_annotation(fl):
	geneids=defaultdict(dict)
	# ~ chr1	564442	564813	ENSG00000225972.1	+	278b2b79-27eb-47bb-94a6-0cf34c53cd47;0,68435139-630b-4d3d-8308-fa2eada3987e;0,
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		for item in ele[-1].split(","):
			item=item.split(";")[0]
			geneids[item]=ele[3]
	return geneids
###############################################################
def site2corrd(basefl):
	# ~ fl="%s/extract.bed6.gene"%(basefl)
	# ~ geneids=pare_annotation(fl)
	fl="%s/extract.reference.bed12"%(basefl)
	geneids=pare_annotation2(fl)
	fl="%s/sam_parse2.txt"%(basefl)
	sites=pare_sam_site(fl,geneids)
	fl="{0}/total_mod.tsv".format(basefl)
	limit=float(FLAGS.proba)
	line1,line2={},{}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		unmod,mod,mark=ele
		ids,spos,seq=mark.split("|")
		mod=float(mod)
		####################
		if mod>limit and ids in sites and spos in sites[ids]:
			# ~ if spos not in sites[ids]:
				# ~ print(ids,spos)
				# ~ continue
			gpos,ref,gbase,genename=sites[ids][spos]
			line1["%s\t%s\t%s"%(ref,gpos,gpos)]=1
			line2["{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(ref,gpos,genename,ids,spos,gbase)]=1
	########################
	output1=open("%s/genome.%s.bed"%(basefl,limit),"w")
	output1.write("\n".join(line1.keys()))
	output1.close()
	#################
	output2=open("%s/genome_abandance.%s.bed"%(basefl,limit),"w")
	output2.write("\n".join(line2.keys()))
	output2.close()
	#################
def paresread_sites(fl):
	read=defaultdict(dict)
	small_memory={}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		#######
		#################
		# ~ read["%s:%s"%(ele[0],ele[1])]=ele[2],ele[-1]
		# ~ Chr02	2218542	POTRI.002G034400.2.V3.0	GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_498_ch_413_strand.fast5	502	AAACA
		# ~ if not ele[2].startswith("POTRI"):
		# ~ continue
		name=ele[2]
		# ~ .split(".")[1]
		read[name+"|"+ele[0]][i.rstrip()]=1
		small_memory["%s%s"%(ele[0],ele[1])]=1
	return read,small_memory
def parse_depth(fl,small_memory):
	readfeature={}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		# ~ Chr01	8408	1
		# ~ Chr10   20468025        CAAGG|+|CAAGG   1       c698ff1b-921e-40a0-bfd8-c1c76c05fb06
		# ~ Chr06	8086941	TGACA	3	GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_3121_ch_398_strand.fast5|GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_22161_ch_467_strand.fast5|GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_22479_ch_393_strand.fast5
		# ~ if "%s%s"%(ele[0],int(ele[1])-1) in small_memory:
			# ~ readfeature["%s|%s"%(ele[0],int(ele[1])-1)]=int(ele[2])+0.0
		if "%s%s"%(ele[0],int(ele[1])) in small_memory:
			readfeature["%s|%s"%(ele[0],int(ele[1]))]=int(ele[2])+0.0
	return readfeature
def parse_depth2(fl,small_memory):
	readfeature={}
	#      9 NM_001040668.1
	for i in open(fl,"r"):
		ele=i.strip().split()
		# ~ Chr01	8408	1
		# ~ Chr10   20468025        CAAGG|+|CAAGG   1       c698ff1b-921e-40a0-bfd8-c1c76c05fb06
		# ~ Chr06	8086941	TGACA	3	GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_3121_ch_398_strand.fast5|GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_22161_ch_467_strand.fast5|GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_22479_ch_393_strand.fast5
		# ~ if "%s%s"%(ele[0],int(ele[1])-1) in small_memory:
			# ~ readfeature["%s|%s"%(ele[0],int(ele[1])-1)]=int(ele[2])+0.0
			
		# ~ if "%s%s"%(ele[0],int(ele[1])) in small_memory:
		readfeature[ele[1]]=int(ele[0])
	return readfeature
def parse_depth3(fl):
	readfeature=defaultdict(dict)
	#chr19	ef3b0b8e-623a-4e62-8d31-19e275acdf65.fast5	-	IRF3	NA		113|50168934|AGACA	165|50168098|GGACC	174|50168089|GAACC	217|50168046|GGACC	250|50168013|GAACA	382|50166722|AGACC	451|50166653|GGACC	463|50166641|GGACC	499|50166605|GAACT	514|50166506|GGACT
	#      9 NM_001040668.1
	for i in open(fl,"r"):
		ele=i.strip().split()
		# ~ Chr01	8408	1
		# ~ Chr10   20468025        CAAGG|+|CAAGG   1       c698ff1b-921e-40a0-bfd8-c1c76c05fb06
		# ~ Chr06	8086941	TGACA	3	GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_3121_ch_398_strand.fast5|GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_22161_ch_467_strand.fast5|GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_22479_ch_393_strand.fast5
		# ~ if "%s%s"%(ele[0],int(ele[1])-1) in small_memory:
			# ~ readfeature["%s|%s"%(ele[0],int(ele[1])-1)]=int(ele[2])+0.0
			
		# ~ if "%s%s"%(ele[0],int(ele[1])) in small_memory:
		for item in ele[5:]:
			gpos=item.split("|")[1]
			readfeature["%s|%s|%s"%(ele[3],ele[0],gpos)][ele[1]]=1
	return readfeature
def readprediction2(fl):
	storepos={}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		#######
		#################
		# ~ read["%s:%s"%(ele[0],ele[1])]=ele[2],ele[-1]
		# ~ Chr02	2218542	POTRI.002G034400.2.V3.0	GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_498_ch_413_strand.fast5	502	AAACA
		# ~ if not ele[2].startswith("POTRI"):
		# ~ continue
		storepos["%s|%s"%(ele[0],ele[1])]=1
	return storepos
def establish_ratio(i,read,readfeature):
	genename,chro=i.split("|")
	poss_summary=defaultdict(dict)
	numlimit=int(FLAGS.support)
	poss=[]
	######
	# ~ Chr10	21971841	Potri.010G244500.1	GXB01149_20180715_FAH87828_GA10000_sequencing_run_20180715_NPL0183_I1_33361_read_252_ch_88_strand.fast5	162	CTACA
	for x in read[i]:
		ele=x.split("\t")
		poss.append(int(ele[1]))
		poss_summary[int(ele[1])][ele[3]]=1
	# ~ poss=set([int(x.split("\t")[1]) for x in read[i]])
	poss=sorted(set(poss))
	#######################
	add=[]
	for subpos in poss:
		# ~ modids=[x.split("\t")[3] for x in read[i] if int(x.split("\t")[1])==subpos]
		nummod=len(poss_summary[subpos].keys())
		if "%s|%s|%s"%(genename,chro,subpos) in readfeature:
			numtotal=len(readfeature["%s|%s|%s"%(genename,chro,subpos)].keys())
			if nummod<numlimit:
				continue
			if numtotal!=0:
				# ~ print("checks","%s|%s"%(chro,subpos),subpos,nummod,numtotal)
				fre=nummod/numtotal
				add.append("%s|%s|%s|%s"%(subpos,nummod,int(numtotal),fre))
	if add:
		return("%s\t%s"%(i,"\t".join(add)))
	else:
		return False
def ratio(basefl):
	limit=float(FLAGS.proba)
	fl="%s/genome_abandance.%s.bed"%(basefl,limit)
	read,small_memory=paresread_sites(fl)
	# ~ fl="%s/extract.depth"%(basefl)
	# ~ readfeature=parse_depth(fl,small_memory)
	#      9 NM_001040668.1
	# ~ readfeature=parse_depth2(fl,small_memory)
	fl="%s/sam_parse2.txt"%(basefl)
	readfeature=parse_depth3(fl)
	###########
	# ~ for pos in ["Chr10|21619580","Chr10|21620214"]:
		# ~ print("bamdepth",pos,readfeature[pos])
	###########
	output=open("%s/ratio.%s.tsv"%(basefl,limit),"w")
	for i in read:
		stats=establish_ratio(i,read,readfeature)
		if stats:
			output.write(stats+"\n")
	output.close()
	#################################################################################################
def replace_gene(basefl):
	fl=FLAGS.isoform
	# ~ head gene2transcripts.txt 
	# ~ DDX11L1	NR_046018.2
	rs_gene={}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		for item in ele[1:]:
			rs_gene[item]=ele[0]
	fl="%s/extract.reference.isoform.bed12"%(basefl)
	output=open("%s/extract.reference.bed12"%(basefl),"w")
	#NM_001197125.1	30	1440	1e0208b3-8061-451f-9415-73df9654f9da.fast5	0	+	30	1440	255,0,0	1	1410	0
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		if ele[0] in rs_gene:
			ele[0]=rs_gene[ele[0]]
			output.write("\t".join(ele)+"\n")
	output.close()
	####################
def run_main():
	cmd="mkdir %s"%(FLAGS.output)
	os.system(cmd)
	basefl=FLAGS.output.rstrip("/")
	fl1=FLAGS.input+".feature.tsv"
	# ~ fl1=FLAGS.input+".feature.tsv"
	########################################
	print("1.start predict")
	########################################
	#split file
	########################################
	for tpp in ["AAACA","AAACT","AGACC","GAACA","GAACT","GGACC","AAACC","AGACA","AGACT","GAACC","GGACA","GGACT"]:
		print("extract %s"%(tpp))
		cmd="grep %s %s >%s/%s.tsv"%(tpp,fl1,basefl,tpp)
		os.system(cmd)
		print("predicting %s"%(tpp))
		predict(tpp,"%s/%s.tsv"%(basefl,tpp),basefl)
	# ~ ####################################
	print("2.start mapping")
	fl1=FLAGS.input+".feature.fa"
	cmd="%ssamtools faidx %s"%(FLAGS.abspathexe,fl1)
	os.system(cmd)
	cmd="%sminimap2 --secondary=no -ax splice -uf -k14 -t %s %s  %s|%ssamtools view -@ %s -bS - |%ssamtools sort -@ %s - >%s/extract.sort.bam"%(FLAGS.abspathexe,FLAGS.cpu,FLAGS.genome,fl1,FLAGS.abspathexe,FLAGS.cpu,FLAGS.abspathexe,FLAGS.cpu,basefl)
	os.system(cmd)
	cmd="%ssamtools index %s/extract.sort.bam"%(FLAGS.abspathexe,basefl)
	os.system(cmd)
	cmd='sam2tsv -r {1} {0}/extract.sort.bam|gzip -c >{0}/extract.sort.bam.tsv.gz'.format(basefl,FLAGS.genome)
	os.system(cmd)
	cmd="%sbedtools bamtobed -bed12 -split -i %s/extract.sort.bam >%s/extract.bed12"%(FLAGS.abspathexe,basefl,basefl)
	os.system(cmd)
	print("gene annotation")
	cmd="%sminimap2 --secondary=no -ax splice -uf -k14 -t %s %s  %s|%ssamtools view -@ %s -bS - |%ssamtools sort -@ %s - >%s/extract.reference.sort.bam"%(FLAGS.abspathexe,FLAGS.cpu,FLAGS.referance,fl1,FLAGS.abspathexe,FLAGS.cpu,FLAGS.abspathexe,FLAGS.cpu,basefl)
	os.system(cmd)
	cmd="%sbedtools bamtobed -bed12 -split -i %s/extract.reference.sort.bam >%s/extract.reference.isoform.bed12"%(FLAGS.abspathexe,basefl,basefl)
	os.system(cmd)
	replace_gene(basefl)
	print("parse bam")
	# ~ ################################################################################
	print("3.m6A site to genome sites")
	# ~ read2genome1(basefl)
	read2genome2(basefl)
	site2corrd(basefl)
	ratio(basefl)
	##################################
	# ~ method2(basefl)
	####################################################################################
def dependence_check():
	#genome file
	fa=FLAGS.genome
	if fa.endswith("fa") or fa.endswith("fasta"):
		pass
	else:
		sys.exit("please check your genome file, make shure it's end with fa or fasta!\n")
	sys.stderr.write("genome file ok!\n")
	#############
	ref=FLAGS.referance
	if ref.endswith("fa") or ref.endswith("fasta"):
		pass
	else:
		sys.exit("please check your referance transcripts sequence file, make shure it's end with fa or fasta!\n")
	sys.stderr.write("referance transcripts sequence file ok!\n")
	#############
	#gene to trans
	fl=FLAGS.isoform
	trans={}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		for item in ele[1:]:
			trans[item]=ele[0]
	#############
	ref=FLAGS.referance
	svs_id={}
	for i in open(ref,"r"):
		ele=i.rstrip().split()
		if i.startswith(">"):
			id1=ele[0].lstrip(">")
			if id1 not in trans:
				sys.exit("please check %s transcripts not find in %s file!\n"%(id1,fl))
	#############
	#genome file index
	faidx1=".".join(fa.split(".")[:-1])+".dict"
	faidx2=fa+".fai"
	if os.path.isfile(faidx1) and os.path.isfile(faidx2):
		pass
	else:
		sys.exit("please check your genome file index, make shure you index with samtools index and picard CreateSequenceDictionary R=ref.fa O=ref.dict !\n")
	sys.stderr.write("genome file index ok!\n")
	##############
	for com in ["samtools","bedtools","minimap2","sam2tsv"]:
		cmd="which %s"%(com)
		pids=subprocess.getstatusoutput(cmd)
		if pids[1]:
			sys.stderr.write("%s ok!\n"%(com))
		else:
			sys.exit("please check %s in your $PATH !\n"%(com))
	#############
	sys.stderr.write("check finsh!\n")
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Predict to genome sites.')
	parser.add_argument('-i', '--input', required = True,help="features_extract")
	parser.add_argument('-o', '--output', required = True, help="Output file")
	parser.add_argument('-g', '--genome', required = True, help="genome file for mapping")
	parser.add_argument('-r', '--referance', required = True, help="referance transcripts sequence file")
	parser.add_argument('-b', '--isoform', required = True, help="gene to referance transcripts information")
	parser.add_argument('--cpu', default=8,help='cpu number usage,default=8')
	parser.add_argument('--support', default=20,help='one m6A site supported read number,default=20')
	parser.add_argument('--proba', default=0.5,help='m6A site predict probability,default=0.5')
	parser.add_argument('--model',required = True, help='model dir')
	args = parser.parse_args(sys.argv[1:])
	global FLAGS
	FLAGS = args
	FLAGS.abspathexe=""
	# ~ folder_path, file_name = os.path.split(os.path.abspath(__file__))
	# ~ FLAGS.abspath=folder_path
	# ~ FLAGS.abspathexe=os.path.dirname(os.path.realpath(sys.executable))
	dependence_check()
	run_main()
###########################################################
