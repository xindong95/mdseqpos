import math
from MA2C import Consts

def _absmedian(d):
    t = map(abs, d)
    t.sort()
    return t[int((len(t) - 1)/ 2)]

def _computeMean(ip, input, C):
    error = 1.0
    
    ipmean_centered = []
    inputmean_centered = []
    d = []
    count = len(ip)
    sumw = float(count)
    
    ipmean = sum(ip) / sumw
    inputmean = sum(input) / sumw
    w = [1.0] * count
    
    ipmean_centered = map(lambda x: x - ipmean, ip)
    inputmean_centered = map(lambda x: x - inputmean, input)
    sig11 = sum(map(lambda x: x * x, ipmean_centered)) / sumw
    sig22 = sum(map(lambda x: x * x, inputmean_centered)) / sumw
    sig12 = sum(map(lambda x, y: x * y, ipmean_centered, inputmean_centered)) / sumw
    
    det_inv = 1.0 / (sig11 * sig22 - sig12 * sig12)
    
    while (error > Consts.ROBUST_CONVERGENCE):
        ipmean_centered = map(lambda x: x - ipmean, ip)
        inputmean_centered = map(lambda x: x - inputmean, input)
        
        d = map(
            lambda x, y: det_inv * (sig22 * x * x - 2 * sig12 * x * y + sig11 * y * y),
            ipmean_centered, inputmean_centered)
        
        mad = _absmedian(d)
        tmp2 = C * mad
        sumw = 0.0
        
        for i in xrange(count):
            tmp = d[i] / tmp2
            if tmp > 1:
                w[i] = 0.0
            else:
                tmp3 = 1 - tmp * tmp
                w[i] = tmp3 * tmp3
                sumw += w[i]
                
        old_ipmean = ipmean
        old_inputmean = inputmean
        
        ipmean = sum(map(lambda x, y: x * y, w, ip)) / sumw
        inputmean = sum(map(lambda x, y: x * y, w, input)) / sumw
        
        sig11 += sum(map(lambda w, x: w * x * x, w, ipmean_centered))
        sig22 += sum(map(lambda w, x: w * x * x, w, inputmean_centered))
        sig12 += sum(map(lambda w, x, y: w * x * y, w, ipmean_centered, inputmean_centered))
        
        sig11 /= sumw
        sig22 /= sumw
        sig12 /= sumw
        
        det_inv = 1.0 / (sig11 * sig22 - sig12 * sig12)
        
        error = max(abs((ipmean - old_ipmean) / (old_ipmean + Consts.ROBUST_CONVERGENCE / 10.0)),
            abs((inputmean - old_inputmean) / (old_inputmean + Consts.ROBUST_CONVERGENCE / 10.0)))
        
    return (ipmean, inputmean, math.sqrt(sig11 + sig22 - 2 * sig12))


def compute(ipGcBin, inputGcBin, C):
    
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
            (ipMean[gc], inputMean[gc], std[gc]) = _computeMean(ip, input, C)
    
    return (ipMean, inputMean, std)