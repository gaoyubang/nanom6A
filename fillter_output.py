import sys,os,re
fl=sys.argv[1]
# ~ fl=sys.argv[1].rstrip("/")
# ~ fl1="%s/ratio.0.5.tsv"%(fl)
save_dict={}
# ~ result_45CM/ratio.0.5.tsv 
# ~ AT4G38840|chr04	18124995|419|671|0.624441132637854	18125012|292|794|0.3677581863979849	18125050|380|835|0.4550898203592814	18125056|453|845|0.5360946745562131	18125080|273|851|0.3207990599294947	18125301|283|815|0.347239263803681	18125342|22|799|0.02753441802252816	18125414|554|783|0.7075351213282248	18125480|185|725|0.25517241379310346
for i in open(fl,"r"):
	ele=i.rstrip().split()
	name,chro=ele[0].split("|")
	for pos in ele[1:]:
		pos=pos.split("|")[0]
		save_dict["%s|%s|%s"%(chro,pos,name)]=1
# ~ fl2="%s/genome_abandance.0.5.bed"%(fl)
fl=sys.argv[2]
output=open(fl+".filler.bed","w")
for i in open(fl,"r"):
	ele=i.rstrip().split()
	# ~ chr04	18125080	AT4G38840	00c33dc7-0df5-44ff-a5e2-3ab513f58b63.fast5	424	AAACA
	if "%s|%s|%s"%(ele[0],ele[1],ele[2]) in save_dict:
		output.write(i)
output.close()
####################
