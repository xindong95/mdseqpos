from MA2C import Functions
from MA2C import CallPeaksCalc
from MA2C import Consts
from MA2C import ScorePeakConvert

def calculate(peakInput, negPeakInput, nullMean, nullStdErr, maxGap):
    
    result = []
    cut = Functions.qnorm(Consts.INITIAL_PVALUE_IN_FDR, True) * nullStdErr + nullMean

    # call peaks to get negative score list
    nbed = CallPeaksCalc.calculate(negPeakInput, cut, maxGap)

    scoreList = _makeScoreList(nbed)
    
    for scoreValue in scoreList:
        #print "score: ",scoreValue
        pValueCut = Functions.pnorm((scoreValue - nullMean) / nullStdErr, True)
        
        pbeds = CallPeaksCalc.calculate(peakInput, scoreValue, maxGap)
        nbeds = CallPeaksCalc.calculate(negPeakInput, scoreValue, maxGap)
        posLen = 0
        for chrom in pbeds.keys():
            posLen = posLen + len(pbeds[chrom])
        
        negLen = 0
        for chrom in nbeds.keys():
            negLen = negLen + len(nbeds[chrom])
        
        fdr = 100.0 * negLen / max(1, posLen)
        
        fdr = min(fdr, Consts.UPPER_LIMIT_FDR)
        if (fdr > Consts.QUIT_FDR and posLen > Consts.QUIT_P_PEAK_NUM):
            break
        #print "FDR: ", fdr
        result.append((fdr, scoreValue, pValueCut, posLen, negLen))
    #result.sort(cmp = lambda x, y: cmp(x[0], y[0]))
    
    return result


def _makeScoreList(nbed):
    
    result = []
    
    for seqId in nbed:
        for record in nbed[seqId]:
            (startIndex, stopIndex, scoreValue) = record
            result.append(scoreValue)
    
    result.sort(reverse = True)
    
    return result
