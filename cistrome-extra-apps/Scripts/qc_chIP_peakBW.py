#!/usr/bin/env python
# Time-stamp: <2010-05-27 14:31:30 Tao Liu>

"""Module Description

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
import re
import logging
import math
from optparse import OptionParser
#from CistromeAP.taolib.CoreLib.Parser import XLSIO,WiggleIO,BedIO
from CistromeAP.taolib.CoreLib.BasicStat.Func import * 
from CistromeAP.jianlib.BwReader import BwIO
import pybedtools

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

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog [options]"
    description = "QC for two replicates of ChIP-chip experiment. (big wig)"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-p","--peak1",dest="peak1",type="string",
                         help="peak file in bed or xls format for #1 replicate")
    optparser.add_option("-q","--peak2",dest="peak2",type="string",
                         help="peak file in bed of xls format for #2 replicate")
    optparser.add_option("-x","--wig1",dest="wig1",type="string",
                         help="bigwig file for #1 replicate")
    optparser.add_option("-y","--wig2",dest="wig2",type="string",
                         help="bigwig file for #2 replicate")
    optparser.add_option("-r","--rfile",dest="rfile",
                         help="R output file")
    #optparser.add_option("-f","--format",dest="format",type="string",
    #                     help="bed, ma2c, mat or macs, default: ma2c",default="bed")
    #optparser.add_option("-s","--step",dest="step",type="int",
    #                     help="number of steps to calculate cor based on some score ranges, default: 5",default=5)
    optparser.add_option("--res",dest="res",type="int",
                          help="resolution to calculate scores. default: 160 (bp)",default=160)
    #optparser.add_option("-m","--method",dest="method",type="string",default="sample",
    #                     help="method to process the paired two sets of data in the sampling step. Choices are 'median', 'mean','sum', and 'sample' (just take one point out of a data set). Default: sample")

    (options,args) = optparser.parse_args()
    if not options.peak1 or not options.peak2 or not options.wig1 or not options.wig2 or not options.rfile:
        optparser.print_help()
        sys.exit(1)

    rfhd = open(options.rfile,"w")

    info("#1 read peaks from bed file")
    peak1 = pybedtools.BedTool(options.peak1)
    peak_merge = peak1.cat(options.peak2)
    info("#1 read score from bigwig file")
    bw1 = BigWigFile(open(options.wig1,'rb'))
    bw2 = BigWigFile(open(options.wig2,'rb'))
    p1 = []
    p2 = []
    for peak in peak_merge:
        summary = bw1.summarize(peak[0], int(peak[1]), int(peak[2]), (int(peak[2])-int(peak[1]))/options.res+1 )
        if not summary:
            continue
        dat1 = summary.sum_data / summary.valid_count
        summary = bw2.summarize(peak[0], int(peak[1]), int(peak[2]), (int(peak[2])-int(peak[1]))/options.res+1 )
        if not summary:
            continue
        dat2 = summary.sum_data / summary.valid_count
        for i in range(len(dat1)):
            if math.isnan(dat1[i]) or math.isnan(dat2[i]):
                continue
            else:
                p1.append(dat1[i])
                p2.append(dat2[i])

    p1name = os.path.basename(options.wig1.rsplit(".bw",2)[0])
    p2name = os.path.basename(options.wig2.rsplit(".bw",2)[0])
    rfhd.write('''library("geneplotter")  ## from BioConductor
    require("RColorBrewer") ## from CRAN
    ''')
    rfhd.write("p1 <- c(%s)\n" % (",".join(map(str,p1))) )
    rfhd.write("p2 <- c(%s)\n" % (",".join(map(str,p2))) )
    rfhd.write("pdf(\"%s.pdf\",width=10,height=10)\n" % options.rfile)

    info("#2 centering pscore1")
    pscore1 = centering(p1)
    info("#2 centering pscore2")    
    pscore2 = centering(p2)
    # overall correlation coefficient
    cor = sum(map(lambda x:x[0]*x[1],zip(pscore1,pscore2)))/(len(pscore1)-1)/std(pscore1)/std(pscore2)

    rfhd.write("smoothScatter(p1,p2,main=\"cor=%.2f\",xlab=\"%s\",ylab=\"%s\")\n" % (cor,p1name,p2name))
    #rfhd.write("plot(p1,p2,cex=0.5,main=\"cor=%.2f\",xlab=\"%s\",ylab=\"%s\")\n" % (cor,p1name,p2name))
    rfhd.write("dev.off()\n")
    rfhd.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
