def calculate(score, bandwidth, minProbes):
    for seqId in score:
        seqData = score[seqId]
        _dealSeqData(seqData,bandwidth, minProbes)
        
    return score 

def _dealSeqData(seqData,bandwidth, minProbes):
    for i in xrange(len(seqData)):
        _dealOneRecord(seqData, i, bandwidth, minProbes)
        
def _dealOneRecord(seqData, i, bandwidth, minProbes):
    
    nearByNr = []
    nearByCount = 0
    
    position = seqData[i][0]
    lenOfProbeSeq = seqData[i][1]
    normalRatios = seqData[i][2]
    
    # find left side
    for testIndex in xrange(i - 1, 0, - 1): # -1->0
        testPosition = seqData[testIndex][0]
        testLenOfProbeSeq = seqData[testIndex][1]
        testNormalRatios = seqData[testIndex][2]
    
        distance = abs(position + lenOfProbeSeq / 2 - testPosition - testLenOfProbeSeq / 2)
        
        if distance >= bandwidth:
            break
        
        nearByNr.extend(testNormalRatios)
        nearByCount = nearByCount + 1
        
    # add self
    nearByNr.extend(normalRatios)
    nearByCount = nearByCount + 1
    
    # find right side
    for testIndex in xrange(i + 1, len(seqData)):
        (testPosition, testLenOfProbeSeq, testNormalRatios) = seqData[testIndex]
    
        distance = abs(position + lenOfProbeSeq / 2 - testPosition - testLenOfProbeSeq / 2)
        
        if distance >= bandwidth:
            break
        
        nearByNr.extend(testNormalRatios)
        nearByCount = nearByCount + 1

    # calculate score
    if nearByCount < minProbes:
        scoreValue = None
    else:
        scoreValue = _calcMatValue(nearByNr)
        
    seqData[i].append(nearByCount)
    seqData[i].append(scoreValue)

def _calcMatValue(nearByNr):
    k=len(nearByNr)
    if (k==1):
        return nearByNr[0]
    nearByNr.sort()    
    if (k%2==0):
        return (nearByNr[k/2]+nearByNr[k/2-1])/2    
    else:
        return nearByNr[int(len(nearByNr)/2)]
