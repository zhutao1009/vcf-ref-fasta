#-*-coding:utf-8-*-
import re
import argparse
import linecache
from pyfasta import Fasta

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description = 'conver a vcf file to fasta file')
parser.add_argument('-v','--vcf',help="The vcf file",type=str,required=True)
parser.add_argument('-r','--ref',help="The reference squence file",type=str,required=True)
parser.add_argument('-k','--keep',help="keep the intermediate output files",type=str)
args = parser.parse_args()

def usage():
	print ("************************************************")
	print ("python vcf-ref-fa.py -r /path/to/refrence -m /path/to/intermediate/file -k [T|F] ")
	print ("************************************************")

ambi={
"AC":"M",
"AG":"R",
"AT":"W",
"CG":"S",
"CT":"Y",
"GT":"K",
"N":"N",
"CA":"M",
"TA":"W",
"GA":"R",
"GC":"S",
"TC":"Y",
"TG":"K"
} #The IUPAC Ambiguity Codes

def vcf2csv(vcf,csv):#Convert a vcf file to a table file and return samples ID
	vcf=open(vcf,"r")
	result=open(csv,"w")
	pattern=re.compile(r'./.')
	linedic={}
	sample=[]
	count=0
	for line in vcf:
		if line.startswith("#CHROM"):
			sample=line.strip("\n").split("\t")[9:]
			count+=1
			result.write("#chrom"+"\t"+"pos"+"\t"+"\t".join(sample)+"\n")
		elif line.startswith("##"):
			pass
		else:
			l=line.strip("\n").split("\t")
			genetype=pattern.findall(line)
			chrom=l[0]
			pos=l[1]
			convered=[]
			ref=l[3]
			alt=l[4]
			if len(alt)<2:
				count+=1
				conver={"0/0":ref,"1/1":alt,"0/1":ambi[ref+alt],"./.":"N"}
				for s in genetype:
					convered.append(conver[s])
				result.write(chrom+"\t"+pos+"\t"+"\t".join(convered)+"\n")
			else:
				pass
			if chrom not in linedic:
				linedic[chrom]=[0,1]
				linedic[chrom][0]=count
				linedic[chrom][1]=count
			else:
				linedic[chrom][-1]=count
	result.close()
	return (sample,linedic)

def fastawrit(list,listname,file):#Write sequence to a file
	f=open(file,"a+")
	newlist="".join(list)
	a=0
	f.write(">"+str(listname)+"\n")
	l=len(newlist)
	while a < l:
		cache=newlist[a:(a+80)]
		a+=80
		f.write(cache+'\n')
	f.close()

def replace_base(reference,linedic,samplename,table):#Read table file and replace the base pairing in reference squence
	for chromID in linedic:
		refseq=list(str(reference[chromID]))
		linerange=range(linedic[chromID][0],linedic[chromID][-1]+1)
		for sample in samplename:
			seq=open(sample+".fa","a+")
			csv=open(table,"r")
			for index in linerange:
				line=linecache.getline(table,index)
				head=line.strip("\n").split("\t")
				chrom=head[0]
				pos=int(head[1])
				getindex=samplename.index(sample)
				bp=head[getindex+2]
				refseq[pos-1]=bp
			fastawrit(refseq,chromID,sample+".fa")
		seq.close()

if __name__ == '__main__':
	usage()
	ref= Fasta(args.ref)
	vcf=args.vcf
	csv=args.vcf+".out"
	sample_and_line=vcf2csv(vcf,csv)
	sample=sample_and_line[0]
	lin_num=sample_and_line[-1]
	replace_base(ref,lin_num,sample,csv)
	if args.keep == "F":
		os.remove(csv)