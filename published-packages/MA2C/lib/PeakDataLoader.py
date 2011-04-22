import csv
from MA2C import CsvFilter

def _matchNormalizedRecord(normalizedReaders, x, y):
    normalRatios = []
    
    skip = False
    
    for readerRecord in normalizedReaders:
        reader = readerRecord["reader"]
        current = readerRecord["current"]
        
        if current == None:
            current = reader.next()
            
            readerRecord["current"] = current
            
        nx = int(current["X"])
        ny = int(current["Y"])
        
        if (nx, ny) != (x, y):
            skip = True
            continue            # mismatch, so 'current' is kept
        else:                   # a match has been found
            normalRatio = float(current["NormalizedLog2Ratio"])
            normalRatios.append(normalRatio)
            current = None      # the match has been added, so empty 'current'
            readerRecord["current"] = current
    
    if skip:
        return None
    else:
        return normalRatios
    
def _appendRecord(store, seqId, record):
    
    if not store.has_key(seqId):
        store[seqId] = []
        
    store[seqId].append(record)
    
    
def read(tagFileInfo):
    
    tpmapFileName = tagFileInfo["output"]["tpmap"]
    normalizedFileNames = map(lambda x:tagFileInfo["output"]["chips"][x]["normalized"],
                              tagFileInfo["output"]["chips"].keys())
                              
    tpmapReader = csv.DictReader(CsvFilter.Filter(tpmapFileName), delimiter = '\t')
    
    normalizedReaders = []
    
    for normalizedFileName in normalizedFileNames:
        normalizedReaders.append(
            {
                "reader": csv.DictReader(CsvFilter.Filter(normalizedFileName), delimiter = '\t'),
                "current": None
                }
            )
    
    result = {}

    for tpmapRow in tpmapReader:
        x = int(tpmapRow["X"])
        y = int(tpmapRow["Y"])
        seqId = tpmapRow["SEQ_ID"]
        position = int(tpmapRow["POSITION"])
        lenOfProbeSeq = len(tpmapRow["PROBE_SEQUENCE"])
        
        try:
            normalRatios = _matchNormalizedRecord(normalizedReaders, x, y)
        except StopIteration:
            break
        
        if normalRatios != None:
            _appendRecord(result, seqId, [position, lenOfProbeSeq, tuple(normalRatios)])
        
    return result
