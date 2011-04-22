import math
from MA2C import Consts

def _computeMean(ip, input):

    nElements = len(ip)
        
    ipmean = sum(ip) / float(nElements)
    inputmean = sum(input) / float(nElements)
    variance = 0.0
    
    for i in xrange(nElements):
        tmp = (ip[i] - input[i] - ipmean + inputmean)
        variance += (tmp*tmp)
        
    variance /= nElements
 
    return (ipmean, inputmean, math.sqrt(variance))

def compute(ipGcBin, inputGcBin):
    
    ipMean = {}
    inputMean = {}
    std = {}
    
    gcKeys = ipGcBin.keys()
    
    for gc in gcKeys:
        ip = ipGcBin[gc]
        input = inputGcBin[gc]
        
        if len(ip) < Consts.NORMALIZE_MINBIN:
            (ipMean[gc], inputMean[gc], std[gc]) = (None, None, None)
        else:
            (ipMean[gc], inputMean[gc], std[gc]) = _computeMean(ip, input)
    
    return (ipMean, inputMean, std) 