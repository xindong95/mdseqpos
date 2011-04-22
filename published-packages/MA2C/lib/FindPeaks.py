from MA2C import Functions
from MA2C import CallPeaksCalc

def find(peakInput, negPeakInput, fdrTable, peakDetectionMethod, threshold, maxGap, nullMean, nullStdErr):
    
    cutoff = 0.0
    
    if peakDetectionMethod == "fdr":
        cutoff = _makeCutOffWithFdr(fdrTable, threshold)
    elif peakDetectionMethod == "pvalue":
        cutoff = nullMean + nullStdErr * Functions.qnorm(threshold, True)
    elif peakDetectionMethod == "ma2c":
        cutoff = threshold
    else:
        raise Exception("unknown method!")
        pass # peak detection method is checked while reading the tag file.

    pbed = CallPeaksCalc.calculate(peakInput, cutoff, maxGap)
    nbed = CallPeaksCalc.calculate(negPeakInput, cutoff, maxGap)
    
    return (pbed, nbed)


def _makeCutOffWithFdr(fdrTable, threshold):
    
    if len(fdrTable) == 0:
        raise Exception('FdrTable is empty!')
    
    result = fdrTable[0][1] # the first item's scoreValue
    
    for record in fdrTable:
        fdr = record[0]
        
        if fdr <= threshold:
            result = record[1]
        else:
            break
    
    return result
