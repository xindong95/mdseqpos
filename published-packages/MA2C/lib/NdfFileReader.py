import csv
from MA2C import CsvFilter

def compareTp(x, y):
    result = cmp(x[1], y[1])
    if result != 0:
        return result
    
    return cmp(x[2], y[2])

def read(fileName, posFileContent = {}):
    reader = csv.DictReader(CsvFilter.Filter(fileName), delimiter = '\t')
    
    result = {}
    result["tpmap"] = []
    
    maxX = 0
    maxY = 0
    
    for row in reader:
        probeId = row.setdefault("PROBE_ID", None)
        if probeId == None:
            raise IOError("PROBE_ID column does not exist in %s!" % fileName)
        
        seqId = row.setdefault("SEQ_ID", None)
        if seqId == None:
            raise IOError("SEQ_ID column does not exist in %s!" % fileName)
        
        probeSequence = row.setdefault("PROBE_SEQUENCE", None)
        if probeSequence == None:
            raise IOError("PROBE_SEQUENCE column does not exist in %s!" % fileName)
        
        x = row.setdefault("X", None)
        if x == None:
            raise IOError("X column does not exist in %s!" % fileName)
        
        y = row.setdefault("Y", None)
        if y == None:
            raise IOError("Y column does not exist in %s!" % fileName)
        
        position = row.setdefault("POSITION", None)
        if position == None:
            raise IOError("POSITION column does not exist in %s!" % fileName)
        
        probe_class = row.setdefault("PROBE_CLASS", None)
        if probe_class != None and probe_class and probe_class.find("experimental") == -1:
            continue

        container = row.setdefault("CONTAINER",None)
        if container != None and container and container.find("BLOCK") == -1:
            continue

        x = int(x)
        y = int(y)
        position = int(position)
        
        maxX = max(maxX, x)
        
        maxY = max(maxY, y)
            
        if posFileContent.has_key(seqId):
            seqId = posFileContent[seqId]
        
        result["tpmap"].append((probeSequence, seqId, position, x, y))
    
    result["tpmap"].sort(cmp = compareTp)
    result["maxX"] = maxX
    result["maxY"] = maxY
    
    return result
        
