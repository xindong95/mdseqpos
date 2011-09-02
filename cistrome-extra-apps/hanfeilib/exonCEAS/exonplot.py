
# ------------------------------------
# Python modules
# ------------------------------------
import sys,operator,warnings
from math import *
from array import *
import copy
import logging
from logging import info
import re
from corelib import _mean,_median,_std,array_extractor,find_nearest

def readDump ( fhd,num,leftbins,rightbins,flag ):
	data = {}
	for l in fhd:

		(chrom,start,end,name,scores) = l.strip().split()
		key = chrom+":"+start+".."+end
		s = scores.replace("nan","0").strip(",").split(",")
		if len(s) == num:

			pool = map(float,s)
			if flag == "in":
				data[key] = pool[(len(pool)-1)/2-leftbins:(len(pool)-1)/2+rightbins+1]
			elif flag == "out":
				data[key] = pool[(len(pool)-1)/2-rightbins:(len(pool)-1)/2+leftbins+1]
	return data

def statistic ( data ):
	values = data.values()
	s = [None]*len(values[0])
	n = 0
	for i in xrange(len(values[0])):
		t = [x[i] for x in values]
		if len(t)==1:
			s[i]=(t[0],t[0],0,1,t[0],t[0])
		else:
			s[i] = (_mean(t),_median(t),_std(t),len(t),min(t),max(t))
	return s
		
def get_exons_byindex(GeneTable,idx,chr,conditions):
	result=[]
	
	if idx==[[]]:
		return False
	for ix in idx:
		ES=array_extractor(GeneTable[chr]['exonStarts'],ix)
		EE=array_extractor(GeneTable[chr]['exonEnds'],ix)
		ED=array_extractor(GeneTable[chr]['strand'],ix)
		for i in zip(ES,EE,ED):
		#i includes exons in one gene
			if conditions.exon_utr==False:
				if len(i[0])>2:
					(ESs,EEs)=(i[0][1:-1],i[1][1:-1])
				else:
					continue
			else:
				(ESs,EEs)=(i[0],i[1])             
			EDs=[i[2]]*len(ESs)
			if i[2]=='+':
				result.extend(zip(ESs,EEs,EDs))
			else:
				result.extend(zip(EEs,ESs,EDs))
	return [filter for filter in result if abs(filter[0]-filter[1])>=conditions.ex_espan]


def enumerate_genes_region(genes_regions):
	"""Enumerate genes"""

	schema={}
	for i,x in enumerate(genes_regions): 
#		print i,x
		if schema.has_key(x):
			schema[x].append(i)
		else:
			schema[x] = [i]
#	print "scccc"
#	print schema
	return schema

def read_wiggle(wiggle_path):
		
		wiggle_file=open(wiggle_path,"r")
		wiggle_content=wiggle_file.read()
		wiggle_chroms_re=re.compile(r'chrom=(\S+)\s')
		
		wiggle_chroms=wiggle_chroms_re.findall(wiggle_content)
		
		wiggle_splited=wiggle_content.split("chrom=")[1:]
		wiggle_profile_re=re.compile(r'(\d+)\t([\-]*\d+\.\d*)')
		
		wiggle_dictionary={}
		
		for (one_section,one_chrom) in zip(wiggle_splited,wiggle_chroms):
#			print "shokou"
			
			profile_list=wiggle_profile_re.findall(one_section)
#			print profile_list[0:10]
			transposed_section=zip(*((int(i),float(j)) for (i,j) in profile_list))
			wiggle_dictionary[one_chrom]=transposed_section
		
		return wiggle_dictionary
			
	
def  get_gene_indicies_by_wiggleprofile(genetable,wiggle_dict,chrom,cutoff):
	#===========================================================================
	# return gene indexes of peak-enriched genes and none peak-enriched genes
	#===========================================================================
	gene_start=list(genetable[chrom]['cdsStart'])
	gene_end=list(genetable[chrom]['cdsEnd'])
	genes_regions=zip(gene_start,gene_end)
	schema=enumerate_genes_region(genes_regions)
	highpeak_gene_idx=[]
	lowpeak_gene_idx=[]
	highpeak_gene_count,lowpeak_gene_count=(0,0)
	sorted_keys=schema.keys()
	sorted_keys.sort()
	start=0
	for (gs,ge) in  sorted_keys:
		if gs==ge:
			continue
		
		wiggle_s=find_nearest(wiggle_dict[chrom][0],gs,start)
		wiggle_e=find_nearest(wiggle_dict[chrom][0],ge,start)
		start=wiggle_s
		try:
			gene_average_profile=_mean(wiggle_dict[chrom][1][wiggle_s:wiggle_e+1])
		except:
			print "ERROR"
		if gene_average_profile>1.5:
			highpeak_gene_idx.append(schema[gs,ge])
			highpeak_gene_count+=1
		else:
			lowpeak_gene_idx.append(schema[gs,ge])
			lowpeak_gene_count+=1
			
	return (highpeak_gene_idx,lowpeak_gene_idx,highpeak_gene_count,lowpeak_gene_count)

def get_gene_indicies_by_bedpeaks(genetable,bed_list,chrom,cutoff=0):
	print chrom
	gene_start=list(genetable[chrom]['cdsStart'])
	gene_end=list(genetable[chrom]['cdsEnd'])	
	genes_regions=zip(gene_start,gene_end)
	schema=enumerate_genes_region(genes_regions)
	has_peak_gene_idx=[]
	no_peak_gene_idx=[]
	has_peak_gene_count,no_peak_gene_count=(0,0)
		
	sorted_keys=schema.keys()
	sorted_keys.sort()

	for (gs,ge) in  sorted_keys:
		
		HasPeak=False
		if gs==ge:

			continue
		for peak in bed_list:
			ps=int(peak[1])
			pe=int(peak[2])
			if  peak[0]!=chrom:
				continue
			elif (ps>gs and ps<ge) or( pe>gs and pe<ge):
				HasPeak=True
				break
			else:
				continue
		
		if HasPeak==True:
			has_peak_gene_idx.append(schema[gs,ge])
			has_peak_gene_count+=1
		else:
			no_peak_gene_idx.append(schema[gs,ge])
			no_peak_gene_count+=1
		if (has_peak_gene_count+no_peak_gene_count)%500 ==0:
			info("##processed genes has peaks/all processed<%d/%d>"\
			     %(has_peak_gene_count,has_peak_gene_count+no_peak_gene_count))
	return (has_peak_gene_idx,no_peak_gene_idx,has_peak_gene_count,no_peak_gene_count)

def get_gene_indicies_by_groups(genes,subsets):
    """Return the indicies of genes belonging to subsets"""
    
    # get the table
    schema=enumerate_genes_region(genes)
    
    ixs=[]
    missing_genes=[]
    for sub in subsets:
        ix=[]
        missing_gene=[]
        for s in sub:
            try:
                ix.extend(schema[s])
            except KeyError,e:
#                warnings.warn("%s does not exist in the gene annotation table" %e)
                missing_gene.append(s)
                pass
        ix = list(set(ix))

        ixs.append(ix)
        missing_genes.append(missing_gene)
    
    return ixs,missing_genes

def paired_bed_make(paired_list,file_a,file_b,chr):
	
	
	lines1=''
	lines2=''
	i=1
	for j in paired_list:
	    lines1+=chr+'\t'+`j[0]`+'\t'+`j[0]+1`+'\t''exon_peak_%s_%d'%(chr,i)+'\t0\t'+j[2]+'\n'
	    lines2+=chr+'\t'+`j[1]`+'\t'+`j[1]+1`+'\t''exon_peak_%s_%d'%(chr,i)+'\t0\t'+j[2]+'\n'
	    i+=1
	file_a.writelines(lines1)
	file_b.writelines(lines2)

def plot ( names, ie, oe, leftbins, rightbins, resolution ,expname):
	l = len(ie)
	m = int(l)/2
	x5prime = [i*resolution for i in xrange(-1*leftbins,rightbins+1) ]
	x3prime = [i*resolution for i in xrange(-1*rightbins,leftbins+1) ]

	txt =  "library(gplots)\n"
	txt += "pdf(\"%s_CI.pdf\")\n" % (expname)
	txt += "par(mfrow=c(1,2))\n"
	txt += "cr <- colorRampPalette(col=c(\"#C8524D\", \"#BDD791\", \"#447CBE\", \"#775A9C\"), bias=1)\n"


	txt += "linecols <- cr(%d)\n" % len(names)
	
	txt +="iey=list();iel=list();iee=list();oey=list();oel=list();oee=list()\n"
	
	i=1
	for a_name in names:
	
		txt += "iey[[%d]] <- c(%s)\n" %(i, ",".join( [str(x[0]) for x in ie[a_name]] ))
		txt += "iel[[%d]] <- c(%s)\n" %(i, ",".join( [str(x[3]) for x in ie[a_name]] ))
		txt += "iee[[%d]] <- qt(0.975,iel[[%d]]-1)*c(%s)/sqrt(iel[[%d]])\n" % (i,i,",".join( [str(x[2]) for x in ie[a_name]] ),i)
	
		txt += "oey[[%d]] <- c(%s)\n" % (i,",".join( [str(x[0]) for x in oe[a_name]] ))
		txt += "oel[[%d]] <- c(%s)\n" % (i,",".join( [str(x[3]) for x in oe[a_name]] ))
		txt += "oee[[%d]] <- qt(0.975,oel[[%d]]-1)*c(%s)/sqrt(oel[[%d]])\n" % (i,i,",".join( [str(x[2]) for x in oe[a_name]] ),i)
		
		i+=1
	
	
	txt += "x  <- c(%s)\n" % ",".join( [str(x) for x in x5prime])
	txt += "x2  <- c(%s)\n" % ",".join( [str(x) for x in x3prime])
	
	txt += "yt=c()\n"
	txt += "yb=c()\n"
	i=1
	for a_name in names:
		txt +="yt[%d]=max(iey[[%d]]+iee[[%d]],oey[[%d]]+iee[[%d]])\n"%(i,i,i,i,i)
		txt +="yb[%d]=min(iey[[%d]]-iee[[%d]],oey[[%d]]-iee[[%d]])\n"%(i,i,i,i,i)
		i+=1
	
	txt+="yt=max(yt)\n"
	txt+="yb=min(yb)\n"
	txt+="yu=yt+(yt-yb)*0.24\n"
	txt+="yd=yb-(yt-yb)*0.15\n"
	i=1
	for a_name in names:
		txt += "par(new=T,mfg=c(1,1),font=2)\n"
		txt += "plotCI(x,iey[[%d]],uiw=iee[[%d]],type=\"p\",pch=NA,cex=0,xlab=\"Around 5' boundary (bp)\",ylab=\"Profile\",col=linecols[%d],ylim=c(yd,yu),lty=2,lwd=0.75,gap=0.1,sfrac=0.001,axes=F)\n" % (i,i,i)
		txt += "lines(x,iey[[%d]],lty=1,lwd=2.43,col=linecols[%d])\n"%(i,i)
		i+=1
	txt += "axis(1)\n"
	txt += "axis(2)\n"
	txt += "abline(v=0,col='black',lty=2)\n"
	
	txt += "names=c(%s)\n"% ",".join( [repr(str(x)) for x in names])
	txt += "legend(\"topleft\",names,col=linecols,lty=c(1,1),pch=c(20,20),lwd=2.43,bg=gray(0.9),text.col=linecols,text.width=strwidth(names,font=2),box.lwd=0,box.col=gray(0.9))\n"
	
	i=1
	for a_name in names:
		txt += "par(new=T,mfg=c(1,2))\n"
		txt += "plotCI(x2,oey[[%d]],oee[[%d]],type=\"p\",pch=NA,cex=0,xlab=\"Around 3' boundary (bp)\",ylab=\"\",col=linecols[%d],ylim=c(yd,yu),lty=2,lwd=0.75,gap=0.1,sfrac=0.001,axes=F)\n" % (i,i,i)
		txt += "lines(x2,oey[[%d]],lty=1,lwd=2.43,col=linecols[%d])\n"	%(i,i)
		i+=1
	txt += "axis(1)\n"
	txt += "abline(v=0,col='black',lty=2)\n"	
	txt += "title(\"%s - %s\",outer=T,line=-3)\n" % ("Average Exon Boundary Profile",expname)
	
	return txt

def make_Rscript_with_CI(gn_names,exp_name,lbin, rbin,res,span):

	num=span/res*2+1
	
	into_exon_stat={}
	outof_exon_stat={}
	for a_group in gn_names:
		into_exon_fhd=open(exp_name+'_'+a_group+'_es_dump.txt','r')
		outof_exon_fhd=open(exp_name+'_'+a_group+'_ee_dump.txt','r')
		print "Start cal dump"
		
		try:
			into_exon = readDump(into_exon_fhd,num,lbin,rbin,"in")
		except:
			print "error reading into_exon"
	
		try:
			outof_exon = readDump(outof_exon_fhd,num,lbin,rbin,"out")
			
		except:
			print "error reading outof_exon"
	
		print into_exon
		into_exon_stat[a_group] = statistic(into_exon)
		outof_exon_stat[a_group] = statistic(outof_exon)
	text = plot(gn_names,into_exon_stat,outof_exon_stat,lbin,rbin,res,exp_name)
	

	return text

if __name__ == '__main__':
	main()


