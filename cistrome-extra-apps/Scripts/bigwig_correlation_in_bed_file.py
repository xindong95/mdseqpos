#!/usr/bin/env python
# Time-stamp: <2011-07-22 13:45:10 Jian Ma>

"""Description: Draw correlation plot for many wiggle files for a given bed file.

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
#import re
import logging
from optparse import OptionParser
#import urllib2
#import tempfile
#import gzip
import subprocess
#from CistromeAP.taolib.CoreLib.Parser import WiggleIO, BedIO
from CistromeAP.taolib.CoreLib.BasicStat.Func import * 
from jianlib.BwReader import BwIO
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

# ------------------------------------
# Misc functions
# ------------------------------------

error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info

def BedInput(fn=''):
    """Read a bed file, return a list"""
    
    f=open(fn,'r')
    standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
    bedlist = []
    for line in f:
        if line.startswith('track') or line.startswith('#') or line.startswith('browser') or not line.strip():
            continue
        l=line.strip().split()
        
        try:
            l[0]=standard_chroms[l[0]]
        except KeyError:
            pass
        bedlist.append(l)
    
    f.close()
    return bedlist

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog [options] <-r rfile> <-b bed file> <-w bigwig file>(>=2)"
    description = """Draw correlation plot for many bigwig files at regions by a bed file.

Method: It will calculate a value for each region defined in a bed
file based on each bigwig files. The method can be chosen from -m
option.
    """
    
    optparser = OptionParser(version="%prog 0.2",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    #optparser.add_option("-d","--db",type="str",dest="dbname",help="UCSC db name for the assembly. Default: ce4",default="ce4")
    optparser.add_option("-z","--imgsize",dest="imgsize",type="int",
                         help="image size. default: 10, minimal: 10",default=10)
    optparser.add_option("-f","--format",dest="imgformat",type="string",
                         help="image format. PDF or PNG",default='PDF')
    #optparser.add_option("-m","--method",dest="method",type="string",default="median",
    #                     help="method to process the paired two sets of data in the sampling step. Choices are 'median', 'mean', and 'sample' (just take one point out of a data set). Default: median")
    optparser.add_option("-r","--rfile",dest="rfile",
                         help="R output file. If not set, do not save R file.")
    optparser.add_option("-b","--bed",dest="bed",type="string",
                         help="the bed file you want to include in the calculation.")
    optparser.add_option("-w","--bw",dest="wig",type="string",action="append",
                         help="the bigwig file you want to include in the calculation. This option should be used for at least twice.")
    optparser.add_option("-l","--wig-label",dest="wiglabel",type="string",action="append",
                         help="the bigwig file labels in the figure. No space is allowed. This option should be used same times as -w option, and please input them in the same order as -w option. default: will use the wiggle file filename as labels.")
    optparser.add_option("--min-score",dest="minscore",type="float",default=0,
                         help="minimum score included in calculation. Points w/ score lower than this will be discarded.")
    optparser.add_option("--max-score",dest="maxscore",type="float",default=10000,
                         help="maximum score included in calculation. Points w/ score larger than this will be discarded.")    
    optparser.add_option("-H","--heatmap",dest="heatmap",action="store_true",default=False,
                         help="If True, a heatmap image will be generated instead of paired scatterplot image.")
    (options,args) = optparser.parse_args()

    imgfmt = options.imgformat.upper()
    if imgfmt != 'PDF' and imgfmt != 'PNG':
        print "unrecognized format: %s" % imgfmt
        sys.exit(1)
    
    medfunc = mean
    #method = options.method.lower()
    #if method == 'median':
    #    medfunc = median
    #elif method == 'mean':
    #    medfunc = mean
    #elif method  == 'sample':
    #    medfunc = lambda u: u[-1]
    #else:
    #    print "unrecognized method: %s" % (method)
    #    sys.exit(1)

    # must provide >=2 wiggle files
    if not options.wig or len(options.wig) < 2 or not options.rfile or not options.bed:
        optparser.print_help()
        sys.exit(1)
    # check the files
    if not os.path.isfile(options.bed):
        error("%s is not valid!" % options.bed)
        sys.exit(1)
    for f in options.wig:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)

    # wig labels
    if options.wiglabel and len(options.wiglabel) == len(options.wig):
        wiglabel = options.wiglabel
    else:        # or use the filename
        wiglabel = map(lambda x:os.path.basename(x),options.wig)
        
    wigfilenum = len(options.wig)
    info("number of wiggle files: %d" % wigfilenum)

    #get chromosome length from optins.wig[0]:
    p=BwIO(options.wig[0])
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
        
    # get the common chromosome list:
    chrset = set([t['key'] for t in p.chromosomeTree['nodes']])
    for bw in options.wig[1:]:
        p=BwIO(bw)
        chrset = chrset.intersection(set([t['key'] for t in p.chromosomeTree['nodes']]))
        
    info("common chromosomes are %s..." % ",".join(chrset))

    # open the R script file handler
    rfhd = open(options.rfile,"w")
    rfhd.write('''require("RColorBrewer") ## from CRAN\n''')

    info("read bed file %s" % os.path.basename(options.bed))        
    bedregion = BedInput(options.bed)
    index_all = range(len(bedregion))
    index_keep = index_all[:] #filter None and keep these regions
    profiles = []
    # for each wig file, sample...
    for i in range(len(options.wig)):
        bw = BigWigFile(open(options.wig[i], 'rb'))
        
        profile = []
        for k in index_all:
            summary = bw.summarize(bedregion[k][0], int(bedregion[k][1]), int(bedregion[k][2]), 1)
            if not summary:
                try:
                    index_keep.remove(k)
                except ValueError:
                    pass
                profile.append(None)
            else:
                profile.append((summary.sum_data / summary.valid_count)[0])
        profiles.append(profile)
    
    for i in range(len(profiles)):
        str_profile = ','.join([str(profiles[i][t]).replace('nan', 'NA') for t in index_keep])
        info("write values to r file")
        rfhd.write('p%d <- c(%s)\n' %(i, str_profile))
        
    rfhd.write("c <- cbind(%s)\n" %','.join(['p%d'%t for t in range(wigfilenum)]))

    rfhd.write("c <- c[ c[,1]<=%f & c[,1]>=%f " % (options.maxscore,options.minscore))
    for i in range(wigfilenum-1):
        rfhd.write("& c[,%d]<=%f & c[,%d]>=%f " % (i+2,options.maxscore,i+2,options.minscore))
    rfhd.write(",]\n")

    if imgfmt == 'PDF':
        rfhd.write("pdf(\"%s.pdf\",width=%d,height=%d)\n" % (options.rfile,options.imgsize,options.imgsize))
    elif imgfmt == 'PNG':
        rfhd.write("png(\"%s.png\",units=\"in\",res=150,width=%d,height=%d)\n" % (options.rfile,options.imgsize,options.imgsize))
    
    if options.heatmap:
        rfhd.write('library(gplots)\n')
        rfhd.write('''
m <- cor(c, method="pearson", use="pairwise.complete.obs")
''')
        labels = ",".join(map(lambda x:"\""+x+"\"",wiglabel))
        rfhd.write("rownames(m) <- c(%s)\n" % labels)
        rfhd.write("colnames(m) <- c(%s)\n" % labels)         
        rfhd.write('# draw the heatmap using gplots heatmap.2\n') 
#        rfhd.write('bitmap("%s.bmp",width=%d,height=%d)\n' % (options.rfile,options.imgsize,options.imgsize))
        rfhd.write('mn <- -1\n')
        rfhd.write('mx <- 1\n')
        rfhd.write('n <- 98\n')
        rfhd.write('bias <- 1\n')
        rfhd.write('mc <- matrix(as.character(round(m, 2)), ncol=dim(m)[2])\n')
        rfhd.write('breaks <- seq(mn, mx, (mx-mn)/(n))\n')
        rfhd.write('cr <- colorRampPalette(colors = c("#2927FF","#FFFFFF","#DF5C5C"), bias=bias)\n')
        rfhd.write('heatmap.2(m, col = cr(n), breaks=breaks, trace="none", cellnote=mc, notecol="black", notecex=1.8, keysize=0.5, density.info="histogram", margins=c(27.0,27.0), cexRow=2.20, cexCol=2.20, revC=T, symm=T)\n')
    else:
        rfhd.write('''
panel.plot <- function( x,y, ... )
{
  par(new=TRUE)
  m <- cbind(x,y)
  plot(m,col=densCols(m),pch=20)
  lines(lowess(m[!is.na(m[,1])&!is.na(m[,2]),]),col="red")  
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,use="complete.obs")
  txt <- format(round(r,2),width=5,nsmall=2)
  #txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * abs(r))
  text(0.5, 0.5, txt, cex = cex.cor)
}
''')
#        rfhd.write("bitmap(\"%s.bmp\",width=%d,height=%d)\n" % (options.rfile,options.imgsize,options.imgsize))
        labels = ",".join(map(lambda x:"\""+x+"\"",wiglabel))
        rfhd.write('''
pairs(c, lower.panel=panel.plot, upper.panel=panel.cor, labels=c(%s))
''' % (labels))
    rfhd.write("dev.off()\n")
    rfhd.close()

    # try to call R
    try:
        subprocess.call(['Rscript',options.rfile])
    except:
        info("Please check %s" % options.rfile)
    else:
        info("Please check %s" % (options.rfile+'.bmp'))
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
