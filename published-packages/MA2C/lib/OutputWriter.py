# Time-stamp: <2012-01-10 16:04:49 Tao Liu>

"""Module Description: Write output files like peak information and
wiggle files.

"""

# ------------------------------------
# python modules
# ------------------------------------
import gzip
import os
import Consts
from math import log

# ------------------------------------
# Misc functions
# ------------------------------------

log10_times_neg_10 = lambda x:-10*log(x,10)

def print_FDR_table (FDRtable, tagFileInfo):
    """format:

    (fdr, scoreValue, pValueCut, posLen, negLen)
    """
    fdr_filename = tagFileInfo["output"]["FDRtable"]
    f = open (fdr_filename,"w")
    f.write("FDR\tMA2Cscore\tPValue\tPosNum-NegNum\n")
    for (fdr, scoreValue, pValueCut, posLen, negLen) in FDRtable:
        f.write("%.2f\t%.6f\t%.4e\t%d\n" % (fdr,scoreValue,pValueCut,posLen-negLen))

    f.close()

    # write R script file
    raw_name = fdr_filename.rstrip("txt")
    r_filename = raw_name+"r"
    pdf_filename = os.path.basename(raw_name)+"pdf"
    f = open(r_filename,"w")
    f.write("# R script generated by MA2C. Usage: $ R --vanilla < %s\n" % (os.path.basename(r_filename)))
    f.write("d <- read.table(\"%s\",header=T)\n" % (os.path.basename(fdr_filename)))
    f.write("pdf(\"%s\",width=8,height=6)\n" % (pdf_filename))
    f.write("plot(d$FDR,d$PosNum.NegNum,xlab=\"FDR (%)\",ylab=\"#Positive peaks - # of Negative peaks\",pch=20)\n")
    f.write("dev.off()")

def print_peak_to_xls ( ppeaks, npeaks, tagFileInfo ):

    f = open (tagFileInfo["output"]["PeakXLS"],"w")
    _print_peak_to_xls(ppeaks, f)
    f.close()
    
    #f = open (tagFileInfo["output"]["NegPeakXLS"],"w")
    #_print_peak_to_xls(npeaks, f)
    #f.close()


def _print_peak_to_xls ( peaks, fhd):
    """format:

    (startPosition, stopPosition, summitPosition-startPosition, lenOfPeak, lenOfProbes, summitScore, pValue, fdr)
    """
    fhd.write("Chr\tStart\tEnd\tSummit\tLength\tName\tMA2Cscore\t-10*log10(PValue)\tFDR\n")
    chroms = peaks.keys()
    chroms.sort()
    i = 0
    for chrom in chroms:
        for (startPosition, stopPosition, summitPosition, lenOfPeak, lenOfProbes, summitScore, pValue, fdr) in peaks[chrom]:
            i+=1
            name = "MA2C_peak_%d" % i
            pValue = log10_times_neg_10(pValue)
            fhd.write("%s\t%d\t%d\t%d\t%d\t%s\t%.6f\t%.2f\t%.2f\n" % (chrom,startPosition, stopPosition, summitPosition, lenOfPeak, name,summitScore, pValue, fdr))

def print_peak_to_bed ( peaks, tagFileInfo ):
    """format:

    (startPosition, stopPosition, summitPosition-startPosition, lenOfPeak, lenOfProbes, summitScore, pValue, fdr)
    """

    f = open (tagFileInfo["output"]["PeakBED"],"w")
    chroms = peaks.keys()
    chroms.sort()
    i = 0
    for chrom in chroms:
        for (startPosition, stopPosition, summitPosition, lenOfPeak, lenOfProbes, summitScore, pValue, fdr) in peaks[chrom]:
            i+=1
            name = "MA2C_peak_%d" % i
            f.write("%s\t%d\t%d\t%s\t%.4f\n" % (chrom,startPosition, stopPosition, name, summitScore))
    f.close()

def save_score_to_wig ( peak_input, tagFileInfo ):
    """Save to wiggle file using variableStep defination.

    """
    f = gzip.open (tagFileInfo["output"]["wiggle"],"w")
    # the header
    f.write("track type=wiggle_0 name=\"%s MA2Cscore\" description=\"MA2Cscore for %s\"\n" % (tagFileInfo["rawName"],tagFileInfo["rawName"])) # data type line

    # probe length from tpmap
    span = 1
    fhd = open(tagFileInfo["output"]["tpmap"])
    for l in fhd:
        if not l.startswith("#") and not l.startswith("PROBE"):
            span = len(l.split()[0])
            break

    # for each chrom
    chroms = peak_input[0].keys()
    chroms.sort()
    for chrom in chroms:
        f.write("variableStep chrom=%s span=%d\n" % (chrom,span))
        # for each probe

        for i in range(len(peak_input[0][chrom])):
            p = peak_input[0][chrom][i]-span/2
            s = peak_input[1][chrom][i]
            if s != Consts.NEG_INF:
                f.write("%d\t%f\n" % (p,s))
    f.close()
