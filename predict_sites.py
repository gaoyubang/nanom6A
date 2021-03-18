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
	fl="{0}/total_mod.tsv".format(basefl)
	storepos=readprediction(fl)
	fl=FLAGS.input+".feature.fa"
	store=readfasta(fl)
	fl="%s/extract.sort.bam.tsv.gz"%(basefl)
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
		if ids in storepos and idspos in storepos[ids]:
			kmer=store[ids][idspos-2:idspos+3]
			line="%s|%s|%s"%(idspos,gpos,kmer)
			total_m6A_reads["%s\t%s\t%s\tNA\tNA\t"%(chro,ids,strand)][line]=1
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
	fl="%s/extract.bed6.gene"%(basefl)
	geneids=pare_annotation(fl)
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
		if "%s|%s"%(chro,subpos) in readfeature:
			numtotal=readfeature["%s|%s"%(chro,subpos)]
			
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
	fl="%s/extract.depth"%(basefl)
	readfeature=parse_depth(fl,small_memory)
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
def method2(basefl):
	fa=FLAGS.input+".feature.fa"
	loader=FastaFile(fa)
	fl1=FLAGS.input+".feature.tsv"
	output=open("%s/20bp.fa"%(basefl),"w")
	for i in open(fl1,"r"):
		ele=i.rstrip().split()
		ids,pos=ele[0].split("|")[:-1]
		pos=int(pos)
		try:
			seq=loader.fetch(ids,pos-30,pos+30)
			output.write(">%s|%s\n%s\n"%(ids,pos,seq))
		except:
			print("ids %s %s,error"%(ids,pos))
	output.close()
	align_hisat2()
	# ~ write_ratio()
	# ~ oneline()
def oneline():
	basefl=FLAGS.output.rstrip("/")
	fl="%s/20bp.fillter.ratio"%(basefl)
	output=open("%s/20bp.fillter.oneline.ratio"%(basefl),"w")
	alls=defaultdict(dict)
	for i in open(fl,"r"):
		# ~ TUG1	chr22|31372152	1	16.0	0.0625
		name,item,n1,n2,rs=i.rstrip().split()
		# ~ RPLP0|chr2	38708960|101|689.0|0.14658925979680695	38709006|30|709.0|0.04231311706629055	38709090|46|734.0|0.06267029972752043	38709165|33|747.0|0.04417670682730924	38709204|65|751.0|0.08655126498002663	38709223|59|751.0|0.07856191744340879	38709273|74|759.0|0.09749670619235837	38709291|257|759.0|0.3386034255599473	38709372|87|767.0|0.11342894393741851	38709395|63|772.0|0.08160621761658031	38709479|182|795.0|0.2289308176100629	38709486|117|798.0|0.14661654135338345	38709519|219|802.0|0.2730673316708229	38709620|107|811.0|0.1319358816276202	38709987|167|786.0|0.21246819338422393
		chro,pos=item.split("|")
		alls["%s|%s"%(name,chro)]["%s|%s|%s|%s"%(pos,n1,n2,rs)]=1
	for item in alls:
		output.write("%s\t%s\n"%(item,"\t".join(alls[item].keys())))
	output.close()
def write_ratio():
	basefl=FLAGS.output.rstrip("/")
	fl="%s/20bp.fillter.bed12"%(basefl)
	alls=defaultdict(dict)
	read_site={}
	# ~ cc6m_2709_T7_ecorv	2703	2723	MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_103979_ch_204_strand.fast5|355	60	+	2703	2723	255,0,0	1	20	0
	# ~ cc6m_2459_T7_ecorv	875	895	MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_23855_ch_465_strand.fast5|109	60	+	875	895	255,0,0	1	20	0
	############################
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		# ~ MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_2305_ch_326_strand.fast5_785.mrna2	cc6m_2459_T7_ecorv	+	813	833	833	833	1	813,	833,	0	MinION_2_20181108_FAK35406_MN29046_sequencing_run_RNA081120181_29618_read_2305_ch_326_strand.fast5_785.path2	none	none	-1,
		chro,s,e=ele[0],int(ele[1]),int(ele[2])
		ids,pos=ele[3].split("|")
		if e-s==30:
			pos=int((s+e)/2+1)
			alls["%s|%s"%(chro,pos)][ids]=1
			read_site[ele[3]]="%s|%s"%(chro,pos)
	##############################
	fl="%s/extract.bed6.gene"%(basefl)
	geneids=pare_annotation(fl)
	##############################
	alls_mod=defaultdict(dict)
	output=open("%s/20bp.fillter.gpos"%(basefl),"w")
	fl="%s/total_mod.tsv"%(basefl)
	gpos2gene={}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		ids,pos=ele[2].split("|")[:-1]
		if "%s|%s"%(ids,pos) in read_site:
			gpos=read_site["%s|%s"%(ids,pos)]
			alls_mod[gpos][ids]=1
			output.write("%s\t%s\t%s\n"%(i.rstrip(),gpos,geneids[ids]))
			gpos2gene[gpos]=geneids[ids]
	output.close()
	################################
	output=open("%s/20bp.fillter.ratio"%(basefl),"w")
	for gpos in alls_mod:
		n1=len(alls_mod[gpos].keys())
		n2=len(alls[gpos].keys())+0.0
		output.write("%s\t%s\t%s\t%s\t%s\n"%(gpos2gene[gpos],gpos,n1,n2,n1/n2))
		# ~ output.write("%s\t%s\t%s\n"%(gpos,"|".join(alls_mod[gpos].keys()),"|".join(alls[gpos].keys())))
	output.close()
	################################
def align_hisat2():
	basefl=FLAGS.output.rstrip("/")
	# ~ cmd="mkdir -p %s/db"%(basefl)
	# ~ os.system(cmd)
	# ~ cmd="%s/hisat_pre/hisat2-build -p %s %s %s/db/fa_idx"%(FLAGS.abspath,FLAGS.cpu,FLAGS.genome,basefl)
	# ~ os.system(cmd)
	cmd="%s/hisat_pre/hisat2 -x %s/db/hg19/genome  -f -U %s/20bp.fa -S %s/20bp.sam -p %s"%(FLAGS.abspath,basefl,basefl,basefl,FLAGS.cpu)
	os.system(cmd)
	cmd="%s/samtools_pre view -@ %s -bS %s/20bp.sam  >%s/20bp.bam"%(FLAGS.abspath,FLAGS.cpu,basefl,basefl)
	os.system(cmd)
	cmd="%s/samtools_pre view -@ %s -F 256 -b %s/20bp.bam  >%s/20bp.fillter.bam"%(FLAGS.abspath,FLAGS.cpu,basefl,basefl)
	os.system(cmd)
	cmd="%s/bedtools_pre bamtobed -i %s/20bp.fillter.bam -bed12 -split -cigar >%s/20bp.fillter.bed12"%(FLAGS.abspath,basefl,basefl)
	os.system(cmd)
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
	cmd="%sminimap2 --secondary=no -ax splice -uf -k14 -t %s %s  %s|%ssamtools view -@ %s -bS - |%ssamtools sort -@ %s - >%s/extract.sort.bam"%(FLAGS.abspathexe,FLAGS.cpu,FLAGS.genome,fl1,FLAGS.abspathexe,FLAGS.cpu,FLAGS.abspathexe,FLAGS.cpu,basefl)
	os.system(cmd)
	cmd="%ssamtools index %s/extract.sort.bam"%(FLAGS.abspathexe,basefl)
	os.system(cmd)
	cmd="%ssamtools view %s/extract.sort.bam >%s/extract.sam"%(FLAGS.abspathexe,basefl,basefl)
	os.system(cmd)
	cmd="%ssamtools depth -d 100000000 %s/extract.sort.bam >%s/extract.depth"%(FLAGS.abspathexe,basefl,basefl)
	os.system(cmd)
	cmd="%ssamtools faidx %s"%(FLAGS.abspathexe,fl1)
	os.system(cmd)
	print("gene annotation")
	cmd="%sbedtools bamtobed -bed12 -split -i %s/extract.sort.bam >%s/extract.bed12"%(FLAGS.abspathexe,basefl,basefl)
	os.system(cmd)
	cmd="cut -f 1,2,3,4,5,6 %s/extract.bed12 >%s/extract.bed6"%(basefl,basefl)
	os.system(cmd)
	cmd="%sbedtools intersect  -a %s -b %s/extract.bed6 -wo|bedtools groupby -g 1,2,3,4,6 -c 10 -o collapse >%s/extract.bed6.gene"%(FLAGS.abspathexe,FLAGS.referance,basefl,basefl)
	os.system(cmd)
	cmd='sam2tsv -r {1} {0}/extract.sort.bam|gzip -c >{0}/extract.sort.bam.tsv.gz'.format(basefl,FLAGS.genome)
	os.system(cmd)
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
	parser.add_argument('-r', '--referance', required = True, help="referance corrd of transcripts")
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
