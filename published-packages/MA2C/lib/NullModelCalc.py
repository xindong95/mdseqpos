import math

def _calcStdErr(nullList):
    mean = sum(nullList) / len(nullList)
    
    sd = sum(map(lambda x: (x - mean) * (x - mean), nullList))
    
    return math.sqrt(sd / (len(nullList) - 1))


def calculate(score, bandwidth):
    
    # build nullList
    nullList = []
    
    for seqId in score:
        seqData = score[seqId]
        
        lastPos = None
        
        for record in seqData:
            (position, lenOfProbeSeq, normalRatio, numOfNearby, scoreValue) = record
            
            if scoreValue == None:
                continue
            
            if lastPos == None:
                nullList.append(scoreValue)
                lastPos = position
                continue
            
            if (position - lastPos) > 2 * bandwidth:
                nullList.append(scoreValue)
                lastPos = position
    
    # find nullMean
    nullListLen = len(nullList)
    
    if nullListLen < 2:
        raise Exception("The elements in NullList is not enough! NullList need 2 elements at least.")
    
    nullList.sort()
    meanIndex1 = (nullListLen - 1) / 2
    meanIndex2 = nullListLen / 2
    nullMean = (nullList[meanIndex1] + nullList[meanIndex2]) / 2.0
    
    # make right part of nullList
    for i in xrange(nullListLen / 2):
        nullList[nullListLen - i - 1] = 2 * nullMean - nullList[i]
    
    # calculate nullStdErr
    nullStdErr = _calcStdErr(nullList)
    
    return (nullMean, nullStdErr)
