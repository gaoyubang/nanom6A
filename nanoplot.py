#!/usr/bin/env python
from __future__ import division
import argparse
from  collections import defaultdict
import cairo,random
import os,sys,collections,re,multiprocessing
# ~ ,commands,shutil
from argparse import ArgumentParser
from math import pi as p
#import pyBigWig
from tqdm import tqdm
def hex2rgb(hexstring, digits=2):
    """Converts a hexstring color to a rgb tuple.

    Example: #ff0000 -> (1.0, 0.0, 0.0)

    digits is an integer number telling how many characters should be
    interpreted for each component in the hexstring.
    """
    if isinstance(hexstring, (tuple, list)):
        return hexstring

    top = float(int(digits * 'f', 16))
    r = int(hexstring[1:digits+1], 16)
    g = int(hexstring[digits+1:digits*2+1], 16)
    b = int(hexstring[digits*2+1:digits*3+1], 16)
    return r / top, g / top, b / top
def parsepos(fl):
	returnfl=[]
	for i in open(fl,"r"):
		# ~ m6A_Seq merip.pos.txt.format #1E90FF
		if i.startswith("#"):
			continue
		name,fl1,colr=i.rstrip().split()
		dic={}
		for i in open(fl1,"r"):
			ele=i.rstrip().split()
			dic[ele[0]]=ele[1]
		returnfl.append([name,dic,colr])
	return returnfl
def wig_count(bw,chro,s,e,libn):
	count=dict()
	# ~ bw=pyBigWig.open(lib)
	oks=bw.intervals(chro,s,e)
	for start,end,high in oks:
		if high==0:
			continue
		for pos in range(start,end+1):
			count[pos]=high*10000000/libn
			#print pos,high
	return count
def parsebw(fl):
	# ~ SRR8244749	data/SRR8244749.bw	#8B6914
	# ~ IP	/home/gyb/Desktop/G/other_frame2/m6A_h2o/output_guppy/merip_Seq_2020_9_30_v2/18.rmdup.bam.bw	9494429	#8B6914
	# ~ name,colour,lib=libs
	returnlist=[]
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		if i.startswith("#"):
			continue
		bw=pyBigWig.open(ele[1])
		returnlist.append([ele[0],ele[3],bw,int(ele[2])])
	return returnlist
def parebed(bed):
	drsbed={}
	for i in open(bed,"r"):
		ele=i.rstrip().split()
		# ~ chr1	12000	14406	DESKTOP_0G0ETDR_20181109_FAK25608_MN30022_sequencing_run_Kris_HEK_20181109_55198_read_52220_ch_182_strand.fast5	10	+	12000	14406	255,0,0	3	227,109,1186	0,612,1220
		drsbed[ele[3]]=i.rstrip()
	return drsbed
def parsebedrr(bed):
	# ~ /home/gyb/Desktop/previos_work/Part3_Circular/input/bed/Ptrichocarpa_210_v3.0.transcript.bed12
	drsbed=defaultdict(dict)
	for i in open(bed,"r"):
		ele=i.rstrip().split()
		# ~ Chr01	8390	12209	Potri.001G000400.3	60	-	8390	12209	255,0,0	4	276,73,96,174	0,2509,2691,3645
		name=".".join(ele[3].split(".")[:-1])
		drsbed[name][i.rstrip()]=1
	return drsbed
def parebedgene(bedgene):
	ds=defaultdict(dict)
	for i in open(bedgene,"r"):
		ele=i.rstrip().split()
		# ~ chr1	761586	762902	LINC00115	-	DESKTOP_0G0ETDR_20181109_FAK25608_MN30022_sequencing_run_Kris_HEK_20181109_55198_read_5047_ch_441_strand.fast5,DESKTOP_0G0ETDR_20181109_FAK25608_MN30022_sequencing_run_Kris_HEK_20181109_55198_read_12129_ch_39_strand.fast5,DESKTOP_0G0ETDR_20181109_FAK25608_MN30022_sequencing_run_Kris_HEK_20181109_55198_read_12134_ch_39_strand.fast5
		for item in ele[-1].split(","):
			ds[ele[3]][item]=1
	return ds
def parsemod(modbed,fl):
	oks_region={}
	for i in open(fl,"r"):
		ele=i.rstrip().split()
		name,chro=ele[0].split("|")
		# ~ Potri.005G192100|Chr05	20941130|24|39|0.6153846153846154	20941586|11|53|0.20754716981132076
		for idx in ele[1:]:
			pos,p1,p2=[int(float(x)) for x in idx.split("|")[:3]]
			if p1>=20:
				oks_region["%s%s"%(chro,pos)]=1
	ds=defaultdict(dict)
	modgs={}
	for i in open(modbed,"r"):
		ele=i.rstrip().split()
		if "%s%s"%(ele[0],ele[1]) in oks_region:
			# ~ chr1	32508228	KHDRBS1	DESKTOP_0G0ETDR_20181109_FAK25608_MN30022_sequencing_run_Kris_HEK_20181109_55198_read_17951_ch_126_strand.fast5	101	AAACA
			ds[ele[3]][ele[1]]=1
			modgs[ele[2]]=1
	return ds,modgs
def example(cr,h,scale,line,arrayline,bar_high,genename):
	cr.move_to(50,bar_high)
	cr.rectangle(50,bar_high-h,10,h)
	cr.set_font_size(8)
	cr.set_source_rgba(0,0,1,0.5)
	cr.fill()
	cr.move_to(50+10+2,bar_high)
	cr.set_source_rgb(0,0,0)
	cr.show_text("nanopore exon")
	###############################################
	cr.move_to(50,bar_high+10)
	cr.line_to(50+10,bar_high+10)
	(r,g,b)=hex2rgb("#000000")
	cr.set_source_rgb(r,g,b)
	cr.set_line_width(line)
	cr.stroke()
	cr.move_to(50+10+2,bar_high+10)
	cr.show_text("Intron")
	cr.fill()
	######################
	cr.move_to(50,bar_high+10+10+10)
	cr.line_to(50+50*scale,bar_high+10+10+10)
	cr.set_source_rgb(0,0,0)
	cr.set_line_width(line)
	cr.stroke()
	cr.move_to(50+50*scale+2,bar_high+10+10+10)
	cr.show_text("50 bp")
	cr.fill()
	#########################################
	cr.move_to(250-30,20)
	cr.show_text("%s"%(genename))
	cr.fill()
	###############################
	##############################################
	step=-30
	for item in ["m6a"]:
		step+=30
		cr.set_line_width(line)
		cr.move_to(50+step,bar_high+50)
		cr.line_to(50+step,bar_high+50-10)
		#cr.rectangle(50,bar_high-h,h,h)
		cr.set_font_size(8)
		r,g,b=hex2rgb("#FF0000")
		cr.set_source_rgba(r,g,b,0.5)
		cr.stroke()
		cr.move_to(50+step+5,bar_high+50)
		cr.set_source_rgb(0,0,0)
		dicts={"m6a":"m6A","pa":"PA","IP":"IP","Input":"Input"}
		cr.show_text(dicts[item])
	#############################################
	bar_high+=20
def transcripts_nanopore(cr,h,scale,line,arrayline,sort_dict,minpos,maxpos,genename):
	cr.set_font_size(4)
	global hh
	################################
	# ~ pa[ids]=pos
	#######################################
	# ~ print(len(sort_dict))
	# ~ cnn=0
	for j in sort_dict:
		# ~ ###########################################
		# ~ cnn+=1
		# ~ cr.move_to(25,hh)
		# ~ cr.show_text("%s"%(cnn))
		# ~ cr.fill()
		############################################
		info=j.split()
		strand=info[5]
		###
		hh+=h*1.5
		############
	# ~ chr1	12000	14406	DESKTOP_0G0ETDR_20181109_FAK25608_MN30022_sequencing_run_Kris_HEK_20181109_55198_read_52220_ch_182_strand.fast5	10	+	12000	14406	255,0,0	3	227,109,1186	0,612,1220
		######
		# ~ ids=info[5].split(";")[1]
		ids=info[3]
		# ~ [0].lower()
		#hh=100	
		############
		total_polya_len=""
		####
		if ids in modbeddict:
			r,g,b=hex2rgb("#FF0000")
			cr.set_source_rgb(r,g,b)
			cr.set_line_width(line/4.0)
			# ~ scaffold_263;1	scaffold_263	11451	12056	-	cufflinks;ed5cb864-b892-4110-8790-3c0cdeb7129e	11451:12056
			total_modification_number=len(modbeddict[ids].keys())
			for m6apos in modbeddict[ids]:
				# ~ print(m6apos)
				m6apos=int(m6apos)
				cr.move_to((m6apos-minpos)*scale+50,hh)
				cr.line_to((m6apos-minpos)*scale+50,hh+h)
				cr.stroke()
		####################################
		###########
		lens=[int(x) for x in info[-2].split(",")]
		starts=[int(x) for x in info[-1].split(",")]
		idss,idee=int(info[1]),int(info[2])
		#########################
		for st in range(len(lens)):
			(s,e)=idss+starts[st],idss+starts[st]+lens[st]
			s=int(s)
			e=int(e)
			#############################################
			if strand == '+' or strand=="-":
				cr.set_source_rgba(0,0,1,0.5)
				cr.rectangle(50+(s-minpos)*scale,hh,(e-s)*scale,h)
				cr.fill()
				if  st!=(len(lens)-1):
					(sn,en)=idss+starts[st+1],idss+starts[st+1]+lens[st+1]
					sn=int(sn)
					en=int(en)
					cr.set_source_rgb(0,0,0)
					cr.set_line_width(line)
					cr.move_to((e-minpos)*scale+50,hh+1*h/2)
					cr.line_to((sn-minpos)*scale+50,hh+1*h/2)
					cr.stroke()
			if strand=="+" and st==(len(lens)-1):
				cr.set_source_rgb(1,0,0)
				cr.set_line_width(arrayline)
				cr.move_to(50+(e-minpos)*scale,hh+1*h/2)
				cr.line_to(50+(e-minpos)*scale+8,hh+1*h/2)
				cr.stroke()
				cr.move_to(50+(e-minpos)*scale+8,hh+1*h/2)
				cr.rel_line_to(-2,-0.5)
				cr.rel_line_to(0,1)
				cr.close_path()
				cr.fill_preserve ()
				cr.stroke()
			if strand == '-' and st==0:
				cr.set_source_rgb(1,0,0)
				cr.set_line_width(arrayline)
				cr.move_to(450-(maxpos-s)*scale,hh+1*h/2)
				cr.line_to(450-(maxpos-s)*scale-8,hh+1*h/2)
				cr.stroke()
				cr.move_to(450-(maxpos-s)*scale-8,hh+1*h/2)
				cr.rel_line_to(2,-0.5)
				cr.rel_line_to(0,1)
				cr.close_path()
				cr.fill_preserve ()
				cr.stroke()
def transcripts_nanoporerr(cr,h,scale,line,arrayline,sort_dict,minpos,maxpos,genename):
	cr.set_font_size(4)
	global hh
	################################
	# ~ pa[ids]=pos
	#######################################
	# ~ print(len(sort_dict))
	# ~ cnn=0
	for j in sort_dict:
		# ~ ###########################################
		# ~ cnn+=1
		# ~ cr.move_to(25,hh)
		# ~ cr.show_text("%s"%(cnn))
		# ~ cr.fill()
		############################################
		info=j.split()
		strand=info[5]
		###
		hh+=h*1.5
		############
	# ~ chr1	12000	14406	DESKTOP_0G0ETDR_20181109_FAK25608_MN30022_sequencing_run_Kris_HEK_20181109_55198_read_52220_ch_182_strand.fast5	10	+	12000	14406	255,0,0	3	227,109,1186	0,612,1220
		######
		# ~ ids=info[5].split(";")[1]
		ids=info[3]
		# ~ [0].lower()
		#hh=100	
		############
		total_polya_len=""
		####
		if ids in modbeddict:
			r,g,b=hex2rgb("#FF0000")
			cr.set_source_rgb(r,g,b)
			cr.set_line_width(line/4.0)
			# ~ scaffold_263;1	scaffold_263	11451	12056	-	cufflinks;ed5cb864-b892-4110-8790-3c0cdeb7129e	11451:12056
			total_modification_number=len(modbeddict[ids].keys())
			for m6apos in modbeddict[ids]:
				# ~ print(m6apos)
				m6apos=int(m6apos)
				cr.move_to((m6apos-minpos)*scale+50,hh)
				cr.line_to((m6apos-minpos)*scale+50,hh+h)
				cr.stroke()
		###############################################################################
		###########
		lens=[int(x) for x in info[-2].split(",")]
		starts=[int(x) for x in info[-1].split(",")]
		idss,idee=int(info[1]),int(info[2])
		#########################
		for st in range(len(lens)):
			(s,e)=idss+starts[st],idss+starts[st]+lens[st]
			s=int(s)
			e=int(e)
			#############################################
			if strand == '+' or strand=="-":
				# ~ cr.set_source_rgba(0,0,1,0.5)
				r,g,b=hex2rgb("#1E90FF")
				cr.set_source_rgba(r,g,b,0.5)
				cr.rectangle(50+(s-minpos)*scale,hh,(e-s)*scale,h)
				cr.fill()
				if  st!=(len(lens)-1):
					(sn,en)=idss+starts[st+1],idss+starts[st+1]+lens[st+1]
					sn=int(sn)
					en=int(en)
					cr.set_source_rgb(0,0,0)
					cr.set_line_width(line)
					cr.move_to((e-minpos)*scale+50,hh+1*h/2)
					cr.line_to((sn-minpos)*scale+50,hh+1*h/2)
					cr.stroke()
			if strand=="+" and st==(len(lens)-1):
				cr.set_source_rgb(1,0,0)
				cr.set_line_width(arrayline)
				cr.move_to(50+(e-minpos)*scale,hh+1*h/2)
				cr.line_to(50+(e-minpos)*scale+8,hh+1*h/2)
				cr.stroke()
				cr.move_to(50+(e-minpos)*scale+8,hh+1*h/2)
				cr.rel_line_to(-2,-0.5)
				cr.rel_line_to(0,1)
				cr.close_path()
				cr.fill_preserve ()
				cr.stroke()
			if strand == '-' and st==0:
				cr.set_source_rgb(1,0,0)
				cr.set_line_width(arrayline)
				cr.move_to(450-(maxpos-s)*scale,hh+1*h/2)
				cr.line_to(450-(maxpos-s)*scale-8,hh+1*h/2)
				cr.stroke()
				cr.move_to(450-(maxpos-s)*scale-8,hh+1*h/2)
				cr.rel_line_to(2,-0.5)
				cr.rel_line_to(0,1)
				cr.close_path()
				cr.fill_preserve ()
				cr.stroke()
def wig(cr,h,scale,line,arrayline,minpos,maxpos,chro,wigs):
	global hh
	maxheight=[]
	#######################
	#cal maxpos
	for libs in wigs:
		hh+=50
		name,colour,lib,libn=libs
		# ~ print(lib,chro,minpos,maxpos)
		wiggg=wig_count(lib,chro,minpos,maxpos,libn)
		if not wiggg:
			maxacnum=0
		else:
			maxacnum=max(wiggg.values())
		maxheight.append(maxacnum)
	if max(maxheight)!=0:
		maxacnum=max(maxheight)
		scaleh=32/maxacnum
	else:
		return 0
	########################
	hh+=50
	for libs in wigs:
		# ~ hh+=50
		name,colour,lib,countn=libs
		# ~ print(lib,chro,minpos,maxpos)
		wiggg=wig_count(lib,chro,minpos,maxpos,countn)
		# ~ if not wiggg:
			# ~ maxacnum=0
			# ~ scaleh=0
		# ~ else:
			# ~ maxacnum=max(wiggg.values())
			# ~ scaleh=32/maxacnum
		cr.set_source_rgb(0,0,0)
		cr.move_to(8,hh)
		cr.set_font_size(5.5)
		cr.save()
		cr.rotate(p*270/180)
		cr.show_text(name)
		cr.fill()
		cr.restore()
		cr.move_to(50-2,hh)
		cr.line_to(50-2,hh-maxacnum*scaleh)
		cr.stroke()
		####################
		aesmin=maxacnum/4
		if True:
			for stick in range(4+1):
			#print stick,aesmin
			#print hh-aesmin*scaleh*(4-stick)
				cr.move_to(50-2,hh-aesmin*scaleh*stick)
				cr.line_to(50-2-2,hh-aesmin*scaleh*stick)
				cr.stroke()
				cr.move_to(50-2-2-4*len(str("%.1f"%(aesmin*stick))),hh-aesmin*scaleh*stick+0.4*h)
				cr.set_font_size(5)
				cr.show_text("%.1f"%(aesmin*stick))
		##########################
		for gp in wiggg:
			(r,g,b)=hex2rgb(colour)
			cr.set_source_rgba(r,g,b,0.05)
			if gp<minpos or maxpos<gp:
				continue
			cr.rectangle(50+(gp-minpos)*scale,hh-wiggg[gp]*scaleh,1*scale,wiggg[gp]*scaleh)
			cr.fill()
def plot_single(gene,mapdict,beddict):
	basefl=FLAGS.output.rstrip("/")
	# ~ if gene not in modgs:
		# ~ return False
	# ~ chr1	12000	14406	DESKTOP_0G0ETDR_20181109_FAK25608_MN30022_sequencing_run_Kris_HEK_20181109_55198_read_52220_ch_182_strand.fast5	10	+	12000	14406	255,0,0	3	227,109,1186	0,612,1220
	minpos  = min(int(beddict[x].split()[1]) for x in mapdict[gene])
	maxpos  = max(int(beddict[x].split()[2]) for x in mapdict[gene])
	# ~ print(mapdict[gene])
	# ~ print(beddict)
	# ~ print(mapdict[gene])
	# ~ print(mapdict[gene].keys()[0],beddict[mapdict[gene].keys()[0]])
	strand = beddict[list(mapdict[gene].keys())[0]].split()[5]
	chro=beddict[list(mapdict[gene].keys())[0]].split()[0]
	# ~ sort_dict= sorted([beddict[x] for x in mapdict[gene]], key=lambda d:(int(d.split()[2])-int(d.split()[1])), reverse = True)
	if strand=="+":
		sort_dict= sorted([beddict[x] for x in mapdict[gene]], key=lambda d:(int(d.split()[2])), reverse = True)
	else:
		sort_dict= sorted([beddict[x] for x in mapdict[gene]], key=lambda d:(int(d.split()[1])), reverse = True)
	# ~ h=70/(len(sort_dict)+2)
	h=5
	hight_of=70+8+1.5*h*(len(sort_dict)+5)+100
	################################################
	surface = cairo.PDFSurface("%s/%s.pdf"%(basefl,gene),500,hight_of)
	# ~ sys.stderr.write("	%s has processed\n"%(gene))
	cr = cairo.Context(surface)
	cr.select_font_face('Helvetica Neue LT Pro', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
	#################################
	line=1*h/3
	arrayline=1*line/2
	scale=400/(maxpos-minpos+1)
	bar_high=35
	scale=400/(maxpos-minpos+1)
	posnum=(maxpos-minpos+1)//10
	example(cr,h,scale,line,arrayline,bar_high,gene)
	global hh
	hh=70+10+20
	#######################
	#############
	transcripts_nanopore(cr,h,scale,line,arrayline,sort_dict,minpos,maxpos,gene)
	#################################
	#################################
	cr.show_page()
	surface.finish()
	# ~ return True
	#######################
def run_main():
	cmd="mkdir %s"%(FLAGS.output)
	os.system(cmd)
	basefl=FLAGS.output.rstrip("/")
	inputfl=FLAGS.input.rstrip("/")
	########
	#parse bed
	bed=inputfl+"/extract.bed12"
	beddict=parebed(bed)
	########
	#parse bed gene
	bedgene=inputfl+"/extract.bed6.gene"
	bedgenedict=parebedgene(bedgene)
	########
	global modbeddict,modgs,bwlist
	modbed=inputfl+"/genome_abandance.%s.bed"%(FLAGS.proba)
	modbedrs=inputfl+"/ratio.%s.tsv"%(FLAGS.proba)
	modbeddict,modgs=parsemod(modbed,modbedrs)
	###############
	###############
	#############################
	nums=[x for x in bedgenedict if x in modgs]
	pbar=tqdm(total=len(nums),position=0, leave=True)
	for i in nums:
		# ~ if i in ["Potri.002G113300","Potri.008G203100"]:
		# ~ if i in ["Potri.009G061500","Potri.003G181400","Potri.007G031700","Potri.001G068100","Potri.002G107300","Potri.005G108900"]:
		plot_single(i,bedgenedict,beddict)
		pbar.update(1)
	##########################################
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Plot modfication sites on DRS reads')
	parser.add_argument('--input', required = True,help="prediction step output")
	parser.add_argument('--proba', default=0.5,help='m6A site predict probability')
	parser.add_argument('-o', '--output', required = True, help="Output file")
	args = parser.parse_args(sys.argv[1:])
	global FLAGS
	FLAGS = args
	folder_path, file_name = os.path.split(os.path.abspath(__file__))
	FLAGS.abspath=folder_path
	run_main()
###############################################

