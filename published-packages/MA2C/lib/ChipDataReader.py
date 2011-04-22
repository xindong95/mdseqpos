# Time-stamp: <2008-06-03 17:04:32 Tao Liu>

"""Module Description: Read raw data

"""

# ------------------------------------
# python modules
# ------------------------------------
import csv
from MA2C import CsvFilter

# ------------------------------------
# Misc functions
# ------------------------------------


def read(dataFileNames):
    result = {}
    result["ip"] = {}
    result["input"] = {}
    
    # ip
    fileName = dataFileNames["ipFile"]
    
    reader = csv.DictReader(CsvFilter.Filter(fileName), delimiter = '\t')
    store = result["ip"]
    
    for row in reader:
        
        if not row.has_key("X"):
            raise IOError("X column does not exist in %s!" % fileName)
        
        if not row.has_key("Y"):
            raise IOError("Y column does not exist in %s!" % fileName)
        
        if not row.has_key("PM"):
            raise IOError("PM column does not exist in %s!" % fileName)
        
        x = int(row["X"])
        y = int(row["Y"])
        pm = row["PM"]
        
        key = (x, y)
        
        store[key] = pm
        
    del reader
        
    # input
    fileName = dataFileNames["inputFile"]
    
    reader = csv.DictReader(CsvFilter.Filter(fileName), delimiter = '\t')
    store = result["input"]
    
    for row in reader:
        
        if not row.has_key("X"):
            raise IOError("X column does not exist in %s!" % fileName)
        
        if not row.has_key("Y"):
            raise IOError("Y column does not exist in %s!" % fileName)
        
        if not row.has_key("PM"):
            raise IOError("PM column does not exist in %s!" % fileName)
        
        x = int(row["X"])
        y = int(row["Y"])
        pm = row["PM"]
        
        key = (x, y)
        
        store[key] = pm
        
    del reader
    
    return result
        
