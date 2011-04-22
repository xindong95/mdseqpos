import math
from MA2C import Consts

def pool(ipMean, inputMean, std):
    gcKeys = ipMean.keys()
    gcKeys.sort()
    
    lowest = gcKeys[0]
    highest = gcKeys[-1]
    middle = (lowest + highest) / 2.0
    
    # lowest to middle
    for gcIndex in xrange(len(gcKeys)):
        gc = gcKeys[gcIndex]
        if gc > middle:
            break
        
        if ipMean[gc] != None:
            continue
        
        poolSourceIndex = None
        for poolIndex in xrange(gcIndex + 1, len(gcKeys)):
            poolValue = ipMean[gcKeys[poolIndex]]
            if poolValue != None:
                poolSourceIndex = poolIndex
                break
            
        if poolSourceIndex == None:
            raise Exception("All value is None!")
            
        # found
        sourceGc = gcKeys[poolSourceIndex]
        ipMean[gc] = ipMean[sourceGc]
        inputMean[gc] = inputMean[sourceGc]
        std[gc] = std[sourceGc]
        
    # highest to middle
    for gcIndex in xrange(len(gcKeys) - 1, - 1, - 1):
        gc = gcKeys[gcIndex]
        if gc < middle:
            break
        
        if ipMean[gc] != None:
            continue
        
        poolSourceIndex = None
        for poolIndex in xrange(gcIndex - 1, - 1, - 1):
            poolValue = ipMean[gcKeys[poolIndex]]
            if poolValue != None:
                poolSourceIndex = poolIndex
                break
            
        if poolSourceIndex == None:
            raise Exception("All value is None!")
            
        # found
        sourceGc = gcKeys[poolSourceIndex]
        ipMean[gc] = ipMean[sourceGc]
        inputMean[gc] = inputMean[sourceGc]
        std[gc] = std[sourceGc]
