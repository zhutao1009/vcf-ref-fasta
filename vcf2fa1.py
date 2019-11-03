#-*-coding:utf-8-*-
import re
import argparse
import numpy as np
import os, sys

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
	sample=[]
	for line in vcf:
			if line.startswith("#CHROM"):
				sample=line.strip("\n").split("\t")[9:]
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
					conver={"0/0":ref,"1/1":alt,"0/1":ambi[ref+alt],"./.":"N"}
					for s in genetype:
						convered.append(conver[s])
					result.write(chrom+"\t"+pos+"\t"+"\t".join(convered)+"\n")
				else:
					pass
	result.close()
	return sample

def readfasta(input):#Read reference sequence
	with open(input,'r') as f:
		fasta={}
		for line in f:
			line = line.strip("\n")
			if line[0] == '>':
				h = line[1:].split(" ")
				header=h[0]
			else:
				sequence = np.array(list(line))
				if header in fasta:
					fasta[header]=np.append(fasta[header],sequence)
				else:
					fasta[header]=sequence		
	return fasta

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


def replace_base(refrence,alist,table):#Read table file and replace the base pairing in reference squence
	for sample in alist:
		seq=open(sample+".fa","a+")
		csv=open(table,"r")
		for line in csv:
			if line.startswith("#chrom"):
				pass
			else:
				head=line.strip("\n").split("\t")
				chrom=head[0]
				pos=int(head[1])
				getindex=alist.index(sample)
				bp=head[getindex+2]
				refrence[chrom][pos-1]=bp
		for scaffold in refrence:
			fastawrit(refrence[scaffold],scaffold,sample+".fa")
		seq.close()
		csv.close()

if __name__ == '__main__':
	usage()
	refseq=readfasta(args.ref)
	vcf=args.vcf
	csv=args.vcf+".txt"
	samplename=vcf2csv(vcf,csv)
	replace_base(refseq,samplename,csv)
	if args.keep == "F":
		os.remove(csv)
	end_time = time()
	run_time = end_time-begin_time